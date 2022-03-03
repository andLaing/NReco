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
    nsipm::Int64 = zero(Int64)
    tmin ::Float32 = zero(Float32)
    tmax ::Float32 = zero(Float32)
end

function sipm_rate(file_path::String, outfile::String, dconf::NReco.DetConf)

    simd_events = 0
    evt_sipm    = SiPMEvt[]
    for fn in glob("*.h5", file_path)
        pdf = read_abc(fn)

        simd_events += nrow(pdf.primaries)
        for wvf in groupby(pdf.waveform, :event_id)
            reco_hits  = NReco.select_sensors(wvf, dconf.ecut, dconf.pde, dconf.sigma_tof)
            if !isnothing(reco_hits) && !isempty(reco_hits)
                nsens    = length(keys(reco_hits))
                min_maxt = combine(reco_hits, :mtime => minimum => :minv, :mtime => maximum => :maxv)
                tmin     = minimum(min_maxt.minv) 
                tmax     = maximum(min_maxt.maxv)
                push!(evt_sipm, SiPMEvt(nsipm = nsens, tmin = tmin, tmax = tmax))
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
    end
    parsed_args = parse_args(s)

    dconf = parsed_args["config"]
    if dconf != "default"
        dconf = from_toml(NReco.DetConf, dconf)
    else
        dconf = NReco.DetConf()
    end
    sipm_rate(parsed_args["dir"], parsed_args["ofile"], dconf)
end