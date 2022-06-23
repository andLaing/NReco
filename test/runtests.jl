using NReco
using ATools

using DataFrames
using Distributions
using HDF5
using Logging
using Statistics

using Test

# Lower function verbosity
logger = global_logger(SimpleLogger(stdout, Logging.Warn))

fname = "testdata/n3-window-1m-LXe-20mm-1-20.h5"
pdf   = NReco.read_abc(fname)
dconf = NReco.DetConf(0.3f0, 0.05f0, 2.0f0, 100.0f0, 5000.0f0, 7, 3)
sxyz  = pdf.sensor_xyz
wfm   = pdf.waveform
vdf   = pdf.vertices




@testset "util" begin
    @test NReco.select_by_index(sxyz,
    "sensor_id", 0) == NReco.select_by_column_value(sxyz, "sensor_id", 0)

    wevt = NReco.select_event(wfm, 4995000)
    @test mean(wevt.time) ≈ 42.39652f0
    @test NReco.select_by_column_value(wevt, "sensor_id", 1783).sensor_id[1] == 1783
    @test mean(NReco.select_by_column_value_lt(wevt, "time", 5.0).time) < 5.0
    @test mean(NReco.select_by_column_value_gt(wevt, "time", 5.0).time) > 5.0
    @test mean(NReco.select_by_column_value_interval(wevt, "time", 5.0, 10.0).time) >5.0
    @test mean(NReco.select_by_column_value_interval(wevt, "time", 5.0, 10.0).time) <10.0

    _, imx, _ = NReco.find_max_xy(wevt, "sensor_id", "time")
    m, i      = findmax(wevt.time)
    @test wevt.sensor_id[i] == imx
end

@testset "images" begin
    img_file = "testdata/imgtest.raw"
    img_pix, img_size, img_data = NReco.read_img(img_file)
    @test all(img_pix .== [51, 51, 101])
    @test all(isapprox.(img_size, [50.0f0, 50.0f0, 100.0f0]))
    @test all(img_pix .== size(img_data))

    fom_file = "testdata/fom_smagicFV_vac_steel"
    fom_df   = NReco.read_foms(fom_file)
    @test nrow(fom_df) == 6
    @test ncol(fom_df) == 4

    allfom_file = "testdata/allFoms_smagicFV_vac_steel"
    allfom_df   = NReco.read_allfoms(allfom_file)
    @test nrow(allfom_df) == 30
    @test ncol(allfom_df) == 19
end

@testset "recof" begin
    test_rad = 360.0f0
    test_xyz = Float32[test_rad * cos(pi / 4), test_rad * sin(pi / 4), 30.0f0]
    test_hit = NReco.Hit(test_xyz..., 2.0f0)
    rad_corr = NReco.radial_correction(test_hit, test_rad)
    @test all(isapprox.(rad_corr, test_xyz))

    wvf1 = first(groupby(wfm, :event_id))
    wvf1 = combine(groupby(wvf1, :sensor_id), nrow => :q)
    wvf1 = leftjoin(wvf1, sxyz, on=:sensor_id)
    disallowmissing!(wvf1, [:x, :y, :z])
    @test NReco.phistd(wvf1) > 0.0
    @test NReco.xyzstd(wvf1, "x") > 0.0

    phis_nadj = NReco.transverse_angle(wvf1, false)
    @test all(abs.(phis_nadj) .<= pi)
    phis_adj  = NReco.transverse_angle(wvf1)
    @test count(phis_adj .> pi) > 0

    @test all(NReco.primary_in_lxe(vdf).parent_id .== 0)

    df          = DataFrame(:q1=>Float32[100.0, 120.5, 132.6],
                            :q2=>Float32[123.6, 122.9,  99.9])
    filtered_df = NReco.filter_energies(df, 100.0f0, 150.0f0)
    @test nrow(filtered_df) == 1
    @test all(in(Array(filtered_df[1, :])).([120.5f0, 122.9f0]))

    df[!, :zstd1] = Float32[2.5, 3.2, 4.9]
    df[!, :zstd2] = Float32[7.9, 4.5, 2.3]
    NReco.calculate_interaction_radius!(df, x -> 2*x, "zstd")
    @test all(in(propertynames(df)).([:r1x, :r2x]))
    @test all(isapprox.(2 * df[!, :zstd1], df[!, :r1x]))
    @test all(isapprox.(2 * df[!, :zstd1], df[!, :r1x]))
end

@testset "lor_space" begin
    ## First a simple one.
    p1 = [ 10.0f0, 325.0f0, 0.0f0]
    p2 = [-10.0f0, 325.0f0, 0.0f0]
    r, phi, z, theta = NReco.lor_from_interaction(p1, p2)
    @test r     == 325.0
    @test isapprox(phi, pi / 2)
    @test z     == 0.0f0
    @test isapprox(theta, 0.0f0)
    ## More complicated
    p1 = [  0.0f0,  350.0f0,  10.0f0]
    p2 = [-10.0f0, -325.0f0, -20.0f0]
    r, phi, z, theta = NReco.lor_from_interaction(p1, p2)
    @test r > 0
    @test 0.0 <= phi < pi
    @test z == (p1[3] + p2[3]) / 2

    ## Now with a selected value. !Functions should be in a file somewhere from here
    function first_vertex(vrts::SubDataFrame)::Bool
        any(vrts.volume_id .== 0 .&& vrts.track_id .== 1) &&
        any(vrts.volume_id .== 0 .&& vrts.track_id .== 2)
    end
    first_vrt  = first(filter(first_vertex, groupby(vdf, :event_id)))
    first_vrt  = filter(row -> row.volume_id == 0 && (row.track_id == 1 || row.track_id == 2), first_vrt)
    ## Up to here needs to be put somewhere central.
    first_prim = pdf.primaries[pdf.primaries.event_id .== unique(first_vrt.event_id), :]
    prim_point = Array(first_prim[1, [:x, :y, :z]])
    prim_dir   = Array(first_prim[1, [:vx, :vy, :vz]])

    ## From primaries: point + unit vector
    track1 = Array(first_vrt[first_vrt.track_id .== 1, [:x, :y, :z]])
    track2 = Array(first_vrt[first_vrt.track_id .== 2, [:x, :y, :z]])
    lor_from_vrt  = NReco.lor_from_interaction(vec(track1), vec(track2))
    lor_from_prim = NReco.lor_from_primary(prim_point, prim_dir)

    # Check on z much less strict as definitions not the same.
    @test all(isapprox.(lor_from_prim[1:end .!= 3], lor_from_vrt[1:end .!= 3]))
    # atol of 3 mm is sufficient here but difference could be more
    # since true interaction points inside the xenon at arbitrary
    # depth up to the outer edge.
    @test isapprox(lor_from_prim[3], lor_from_vrt[3], atol=3)
end

@testset "recohits" begin
    grp_wvf = groupby(wfm, :event_id)
    wvf1    = first(grp_wvf)

    slct_sens = NReco.select_sensors(wvf1, dconf.ecut, dconf.pde, dconf.sigma_tof)
    @test typeof(slct_sens) <: GroupedDataFrame
    @test all(in([:sensor_id]).(groupcols(slct_sens)))
    expected_columns = ["event_id", "sensor_id", "time", "dt", "mtime"]
    @test all(in(expected_columns).(names(slct_sens)))
    @test all(combine(slct_sens, nrow).nrow .> dconf.ecut)

    xyzqt = NReco.sensor_positions(slct_sens, sxyz)
    @test typeof(xyzqt) <: DataFrame
    expected_columns = [:sensor_id, :tmin, :trmin, :q, :x, :y, :z]
    @test all(in(expected_columns).(propertynames(xyzqt)))
    @test length(unique(xyzqt.sensor_id)) == length(xyzqt.sensor_id)

    mean_time = NReco.average_first_hits(slct_sens, xyzqt.sensor_id[1:20], dconf)
    min_max = extrema(filter(row -> in(xyzqt.sensor_id[1:20])(row.sensor_id), xyzqt).trmin)
    @test (mean_time >= min_max[1]) && (mean_time < min_max[2])

    inbound = ATools.range_bound(dconf.qmin, dconf.qmax, ATools.OpenBound)
    hemis1  = NReco.split_hemispheres(wvf1, sxyz, dconf,
                                      inbound, NReco.lor_maxq)
    ## The first event is rejected due to low charge.
    @test isnothing(hemis1)
    hemis2 = NReco.split_hemispheres(grp_wvf[(4995005,)], sxyz, dconf,
                                     inbound, NReco.lor_maxq)
    @test !isnothing(hemis2)
    @test length(hemis2) == 6
    @test typeof(hemis2[1]) <: NReco.Hit && typeof(hemis2[2]) <: NReco.Hit
    @test typeof(hemis2[3]) <: Float32   && typeof(hemis2[4]) <: Float32
    @test typeof(hemis2[5]) <: DataFrame && typeof(hemis2[6]) <: DataFrame
end

@testset "nemareco" begin
    exp_keys = [:event_id, :phot1, :phot2, :nsipm1, :nsipm2, :q1, :q2,
	            :E1, :E2, :r1,  :r2, :r1x, :r2x,
                :phistd1, :zstd1, :widz1, :widphi1, :corrzphi1,
                :phistd2, :zstd2, :widz2, :widphi2, :corrzphi2,
			    :xs, :ys, :zs, :ux, :uy, :uz, :xt1, :yt1, :zt1,
                :ti1, :t1, :xt2, :yt2, :zt2, :ti2, :t2, :x1, :y1, :z1,
                :x2, :y2, :z2, :xr1, :yr1, :zr1, :tr1,
                :xr2, :yr2, :zr2, :tr2, :xb1, :yb1, :zb1, :ta1,
			    :xb2, :yb2, :zb2, :ta2]
    _, result = NReco.nemareco([fname], dconf)
    result_fields = fieldnames(typeof(result[1]))
    @test length(result_fields) == length(exp_keys)
    @test all(in(exp_keys).(result_fields))
    filter_func = ismissing, isnothing, isnan
    corrzphi1   = filter(c -> !any(f -> f(c), filter_func), getfield.(result, :corrzphi1))
    corrzphi2   = filter(c -> !any(f -> f(c), filter_func), getfield.(result, :corrzphi2))
    @test all((corrzphi1 .<= 1.0) .& (corrzphi1 .>= -1.0))
    @test all((corrzphi2 .<= 1.0) .& (corrzphi2 .>= -1.0))
end

@testset "classify_events" begin
    evt_keys   = [:total, :empty, :single, :prompt, :single_prompt, :good_prompt]
    evt_counts = NReco.event_classifier([fname], dconf)
    evt_fields = fieldnames(typeof(evt_counts))
    @test length(evt_fields) == length(evt_keys)
    @test all(in(evt_keys).(evt_fields))

    @test evt_counts.total         == 20
    @test evt_counts.empty         ==  7
    @test evt_counts.single        == 11
    @test evt_counts.prompt        ==  2
    @test evt_counts.single_prompt == 11
    @test evt_counts.good_prompt   ==  2
end


## Tests for executables
include("../scripts/makenema.jl")
@testset "makenema" begin
    config = Dict(
        "dir"     => "testdata/"  ,
        "pattern" => fname[10:end],
        "odir"    => tempdir()    ,
        "ofile"   => "evtdf.h5"   ,
        "filei"   => 1            ,
        "filel"   => 1            ,
        "loralgo" => "lor_kmeans" ,
        "detconf" => "default"
    )
    makenema(config)

    outfile = joinpath(config["odir"], config["ofile"])
    @test isfile(outfile)
    h5open(outfile) do h5test
        @test haskey(h5test                   , "configuration"   )
        @test haskey(h5test["configuration"]  , "RunConfiguration")
        @test haskey(h5test                   , "selected_events" )
        @test haskey(h5test["selected_events"], "EventParameters" )
    end
end


include("../scripts/makeEventTable.jl")
@testset "makeEventTable" begin
    config = Dict(
        "dir"     => "testdata/"  ,
        "pattern" => fname[10:end],
        "odir"    => tempdir()    ,
        "ofile"   => "evtdf.h5"   ,
        "filei"   => 1            ,
        "filel"   => 1            ,
        "config"  => "default"
    )
    makezoo(config)

    outfile = joinpath(config["odir"], config["ofile"])
    @test isfile(outfile)
    h5open(outfile) do h5test
        @test haskey(h5test               , "EventCounts")
        @test haskey(h5test["EventCounts"], "counts"     )
    end
end


include("../scripts/calibrate_lors.jl")
@testset "calibrate_lors" begin
    config = Dict(
        "input-file"     => "n3-window-1m-LXe-20mm-1-20_reduced.h5",
        "output-file"    => "n3-window-1m-LXe-20mm-1-20_reduced"   ,
        "time-parameter" => "ta"                                   ,
        "trueout"        => false
    )
    cal_conf = NReco.CalConfig(
        input_dir = "testdata/",
        conf_dir  = ""         ,
        plot_dir  = tempdir()  ,
        qmin      =  600.0f0   ,
        qmax      = 2200.0f0   ,
        cal_func  = NReco.CalFunction(
                        cal_grp = "n3-20mm",
                        cal_std = "cstd"
                    )
    )
    path_in, path_out = define_paths(cal_conf)
    calibrate_lors(config, cal_conf, path_in, path_out)

    out_file = joinpath(path_out, config["output-file"] * "_mlor.h5")
    @test isfile(out_file)
    expected = DataFrame(:dt => Float32[   1.0083368, 0.20166937],
                         :x1 => Float32[ 357.57245, 103.550095],
                         :y1 => Float32[  33.63522 , -343.89935],
                         :z1 => Float32[ 166.94383, -219.65816],
                         :x2 => Float32[-230.92096, 143.3677],
                         :y2 => Float32[-272.25385, -329.29483],
                         :z2 => Float32[-197.30713, -231.62079],
                         :q1 => Float32[1985.0, 607.0],
                         :q2 => Float32[814.0, 919.0],
                         :E1 => Float32[511.0, 511.0],
                         :E2 => Float32[232.63318, 105.336716])
    h5open(out_file) do h5test
        @test haskey(h5test, "reco_info")
        @test haskey(h5test["reco_info"], "lors")
        lors = readh5_todf(h5test, "reco_info", "lors")
        @test nrow(        lors ) ==  2
        @test length(names(lors)) == 11
        @test all(all.(eachrow(isapprox.(lors, expected))))
    end
end
