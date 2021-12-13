using Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "..")))

using ATools
using NReco

using ArgParse
using Configurations
using DataFrames
using Glob
using HDF5
using Unitful

import Unitful: ns, mm


function define_paths(config::NReco.CalConfig)
    in_path  = joinpath(config.input_dir, config.conf_dir)
    out_path = joinpath(config.plot_dir , config.conf_dir)
    if !isdir(out_path)
        mkpath(out_path)
    end
    return in_path, out_path
end


function calibrate_lors()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--conf", "-c"
            help     = "Calibration configuration"
			arg_type = String
        "--trueout", "-t"
            help     = "Output true position LORs?"
            arg_type = Bool
            default  = false
        "--time-parameter", "-m"
            help     = "Time parameter used: ta or tr"
            arg_type = String
            default  = "ta"
        "--output-file", "-o"
            help     = "Name of output file"
            arg_type = String
            default  = "default"
    end
    confArgs     = parse_args(s)

    conf         = from_toml(NReco.CalConfig, confArgs["conf"])


    (path_in ,
     path_out)   = define_paths(conf)

     if confArgs["output-file"] == "default"
         outfile = joinpath(path_out, conf.conf_dir[1:end-1])
     else
         outfile = joinpath(path_out, confArgs["output-file"])
     end

    (nsim, rmin,
     rmax, ndf ) = read_evtpar(glob("evt*.h5", path_in))

    cal_func     = NReco.calibration_function(conf.cal_func, rmin, rmax)

    ndfq         = NReco.filter_energies(ndf, conf.qmin, conf.qmax)

    NReco.calculate_interaction_radius!(ndfq, cal_func, conf.cal_func.cal_std, rmax)

    # We need/want units here. More generalised use of units should be implemented.
    ## Get the true interaction radius (Why not saved?)
    transform!(ndfq, [:xt1, :yt1] => ByRow(rxy) => :rt1, [:xt2, :yt2] => ByRow(rxy) => :rt2)
    units_ndfq = ATools.set_units(ndfq)

    # Calculate the x, y, z of the interaction.
    xint1, yint1, zint1 = ATools.radial_correction(units_ndfq.xr1 / mm, units_ndfq.yr1 / mm,
                                                   units_ndfq.zr1 / mm, units_ndfq.r1x / mm)
    xint2, yint2, zint2 = ATools.radial_correction(units_ndfq.xr2 / mm, units_ndfq.yr2 / mm,
                                                   units_ndfq.zr2 / mm, units_ndfq.r2x / mm)
    # Calculate the interaction time.
    time1 = Symbol(confArgs["time-parameter"] * "1")
    time2 = Symbol(confArgs["time-parameter"] * "2")
    t1    = ATools.interaction_time(units_ndfq, :r1x, time1, rmax, conf.cal_func.nLXe)
    t2    = ATools.interaction_time(units_ndfq, :r2x, time2, rmax, conf.cal_func.nLXe)

    # Calculate LORs.
    dt   = uconvert.(ns, t2 - t1) ./ ns
    mLor = ATools.MlemLor.(dt, xint1, yint1, zint1, xint2, yint2, zint2,
        units_ndfq[!, :q1], units_ndfq[!, :q2], units_ndfq[!, :E1], units_ndfq[!, :E2])
    # Will want a dataset/table with metadata too, to be decided.
    mlor_filename = outfile * "_mlor.h5"
    ATools.write_lors_hdf5(mlor_filename, mLor)
    if confArgs["trueout"]
        flight1   = ATools.time_of_flight(units_ndfq, [:xs, :ys, :zs], [:xt1, :yt1, :zt1])
        flight2   = ATools.time_of_flight(units_ndfq, [:xs, :ys, :zs], [:xt2, :yt2, :zt2])
        dt        = uconvert.(ns, flight2 - flight1) ./ ns
        true_mLor = ATools.MlemLor.(dt, ndfq.xt1, ndfq.yt1, ndfq.zt1,
            ndfq.xt2, ndfq.yt2, ndfq.zt2,
            units_ndfq[!, :q1], units_ndfq[!, :q2], units_ndfq[!, :E1], units_ndfq[!, :E2])
        mlor_filename = outfile * "_Truemlor.h5"
        ATools.write_lors_hdf5(mlor_filename, true_mLor, "true_info")
    end
end
calibrate_lors()
