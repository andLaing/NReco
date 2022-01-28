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
using Logging
using Printf

logger = SimpleLogger(stdout, Logging.Warn)
old_logger = global_logger(logger)


function makenema(args)

	lorf      = args["loralgo"]
	detconf   = args["detconf"]
	dr        = args["dir"]
	f_pattern = args["pattern"]
	outd      = args["odir"]
	outf      = args["ofile"]
	file_i    = args["filei"]
	file_l    = args["filel"]

	lor_algo = NReco.lor_kmeans
	if lorf == "lor_qmax"
		lor_algo = NReco.lor_maxq
	end

	if isdir(outd) == false
		mkdir(outd)
	end
	output = joinpath(outd, outf)

	## Patch here so we can read actually select in a directory which
	## is already partially processed. Needs to be improved in general.
<<<<<<< HEAD
	files  = sort(glob(f_pattern, dr), by=x->parse(Int64, split(x, "-")[end][1:end-3]))
=======
	#files  = sort(glob("*.h5", dr), by=x->parse(Int64, split(x, "-")[end][1:end-3]))
	files = glob("*.h5", dr)
>>>>>>> 89b8515 (Remove sort)

	if detconf != "default"
		dconf = from_toml(NReco.DetConf, detconf)
	else
		dconf = NReco.DetConf()
	end

	println("makenema configuration")
	println("detector configuration", dconf)
	println(" lor_algo  = $lorf")

	println("number of files in data dir = $(length(files))")
	println("reading = $(file_l - file_i + 1) files")
	println("output file  = $output")

	h5open(output, "w") do h5out
		pars_grp    = create_group(h5out, "selected_events")
		pars_dtype  = ATools.generate_hdf5_datatype(ATools.EventParameters)
		conf_grp    = create_group(h5out, "configuration")
		conf_dtype  = ATools.generate_hdf5_datatype(ATools.SimConfiguration)

		nevt, evts  = NReco.nemareco(files, dconf, file_i, file_l, lor_algo)

		rmin, rmax  = extrema(getfield.(evts, :r1))
		# Maybe combine the DetConf and SimConfiguration types?
		config_pars = ATools.SimConfiguration(SiPMPDE=dconf.pde,
						ElecSTD=dconf.sigma_tof, SiPMThr=dconf.ecut,
						SumQmin=dconf.qmin, SumQmax=dconf.qmax,
						Ntof=dconf.ntof, NEvent=nevt, Rmin=rmin,
						Rmax=rmax)#TODO: Calibration option implementation

		## Save to file.
		pars_dspace = dataspace(evts)
		pars_dset   = create_dataset(pars_grp, "EventParameters", pars_dtype, pars_dspace)
		write_dataset(pars_dset, pars_dtype, evts)
		conf_dspace = dataspace([config_pars])
		conf_dset   = create_dataset(conf_grp, "RunConfiguration", conf_dtype, conf_dspace)
		write_dataset(conf_dset, conf_dtype, [config_pars])#Must be a better way
	end
end


function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--dir", "-d"
            help     = "directory with nema simulations"
            arg_type = String
		"--pattern", "-p"
			help     = "Input file pattern to match"
			arg_type = String
			default  = "*.h5"
		"--odir", "-o"
            help     = "output directory"
            arg_type = String
		"--ofile", "-x"
            help     = "output file"
            arg_type = String
            default  = "evtdf.h5"
		"--filei", "-i"
	        help = "number of initial file in glob list"
	        default  = 1
			arg_type = Int
		"--filel", "-l"
		    help     = "number of last file in glob list"
		    default  = 1
			arg_type = Int
		"--loralgo", "-g"
			help     = "algorithm to use for LOR reconstruction "
			arg_type = String
			default  = "lor_kmeans"
		"--detconf", "-c"
			help     = "Detector configuration"
			arg_type = String
			default  = "default"
    end

    return parse_args(s)
end


if abspath(PROGRAM_FILE) == @__FILE__
	parsed_args = parse_commandline()
	println("Running makenema with arguments", parsed_args)
	makenema(parsed_args)
end
