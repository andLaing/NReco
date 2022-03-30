if abspath(PROGRAM_FILE) == @__FILE__
    using Pkg
    Pkg.activate(normpath(joinpath(@__DIR__, "..")))
end

using ATools
using NReco

using ArgParse
using Configurations
using DataFrames
using Glob
using HDF5

import Base.@kwdef
@kwdef struct SiPMEvt <: ATools.OutputDataset
    nevt ::Int64   = one(Int64)
    nvtx ::Int64   = zero(Int64)
    nsipm::Int64   = zero(Int64)
    tmin ::Float32 = zero(Float32)
    tmax ::Float32 = zero(Float32)
    k2sil::Float32 = zero(Float32)
    k3sil::Float32 = zero(Float32)
    k4sil::Float32 = zero(Float32)
end

using Clustering
using Distances
using Statistics
function lor_kmeansTest(hitdf::DataFrame, k::Int64)
    Mhits = transpose(Matrix(hitdf[!, [:x, :y, :z]]))
	kr    = kmeans(Mhits, k, weights=hitdf.q)
	#ka    = assignments(kr)

    distances = pairwise(SqEuclidean(), Mhits)
	mean(silhouettes(kr, distances))
end

function count_vertices(vertices::GroupedDataFrame, vols::Vector{Int64}, evts::Vector{Int64})
    vtx_count = 0
    for evt in evts
        try
            vrtx = vertices[(evt,)] 
        catch
            continue
        end
        vrtx = filter(df -> in(vols).(df.volume_id), vrtx)
        vtx_count += length(unique(vrtx.track_id))
    end
    vtx_count
end

function sipm_rate(file_path::String, outfile::String, dconf::NReco.DetConf,
                   source_rate::Float32, evt_window::Float32)

    simd_events = 0
    evt_sipm    = SiPMEvt[]
    time_bins   = NReco.calculate_timebins(source_rate, evt_window)
    for fn in glob("*.h5", file_path)
        pdf = read_abc(fn)

        simd_events    += nrow(pdf.primaries)
        mixed_waveforms = time_bins(pdf.primaries[:, [:event_id]], pdf.waveform)
        grp_vrtx = groupby(filter(df -> df.parent_id .== 0, pdf.vertices), :event_id)
        volumes  = findall((pdf.volume_names .==     "LXe") .|
						   (pdf.volume_names .== "Steel_1")   ) .- 1
        for wvf in groupby(mixed_waveforms, :time_bin)
            kvals = Float32[0.0, 0.0, 0.0]
            nvtx  = 0
            reco_hits  = NReco.select_sensors(wvf, dconf.ecut, dconf.pde, dconf.sigma_tof)
            if !isnothing(reco_hits) && !isempty(reco_hits)
                nevt     = length(unique(wvf.event_id))
                nsens    = length(keys(reco_hits))
                min_maxt = combine(reco_hits, :mtime => (x -> [extrema(x)]) => [:minv, :maxv])
                tmin     = minimum(min_maxt.minv) 
                tmax     = maximum(min_maxt.maxv)
                xyzqt    = NReco.sensor_positions(reco_hits, pdf.sensor_xyz)
                if nrow(xyzqt) > 4
                    evts = unique(wvf.event_id)
                    nvtx = count_vertices(grp_vrtx, volumes, evts)
                    #_, _, dotdf1, dotdf2 = NReco.lor_maxq(xyzqt)
                    #println("nevent $nevt, dot product split df1 sensors = $(nrow(dotdf1)), df2 sensors = $(nrow(dotdf2))")
                    for k in 2:4
                        kvals[k-1] = lor_kmeansTest(xyzqt, k)
                        #println("Test k = $k, mean silhouette = ", lor_kmeansTest(xyzqt, k))
                    end
                end
                push!(evt_sipm, SiPMEvt(nevt  = nevt    , nsipm = nsens   ,
                                        tmin  = tmin    , tmax  = tmax    ,
                                        k2sil = kvals[1], k3sil = kvals[2],
                                        k4sil = kvals[3], nvtx  = nvtx    ))
            end
        end
    end

    h5open(outfile, "w") do h5out
        grp    = create_group(h5out, "SiPMRate")
        attributes(grp)["total_generations"] = simd_events

        dtype  = ATools.generate_hdf5_datatype(SiPMEvt)
        dspace = dataspace(evt_sipm)
        dset   = create_dataset(grp, "counts", dtype, dspace)
        write_dataset(dset, dtype, evt_sipm)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--dir", "-d"
            help     = "directory with nema simulations"
            arg_type = String
        "--ofile", "-x"
            help     = "output file"
            arg_type = String
            default  = "sipm_counts.h5"
        "--config", "-c"
			help     = "detector configuration"
			default  = "default"
			arg_type = String
        "--rate", "-r"
            help     = "Source rate in Bq"
            default  = 2000.0f0
            arg_type = Float32
        "--window", "-w"
            help     = "Event coincidence window in ns"
            default  = 1000.0f0
            arg_type = Float32
    end
    parsed_args = parse_args(s)

    dconf = parsed_args["config"]
    if dconf != "default"
        dconf = from_toml(NReco.DetConf, dconf)
    else
        dconf = NReco.DetConf()
    end
    sipm_rate(parsed_args["dir"], parsed_args["ofile"], dconf,
              parsed_args["rate"], parsed_args["window"])
end