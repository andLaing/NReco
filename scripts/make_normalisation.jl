if abspath(PROGRAM_FILE) == @__FILE__
	using Pkg
	Pkg.activate(normpath(joinpath(@__DIR__, "..")))
end

using ATools
using NReco

using ArgParse
using DataFrames
using Glob
using HDF5

function normalisation_histos(args::Dict{String, Any})
	indir     = args["dir"    ]
	f_pattern = args["pattern"]
	nbins     = args["nbins"  ]
	fov       = args["fov"    ]
	outfile   = args["out"    ]
	# Always assume FOV centred.
	bin_limits = [(-dim / 2, dim / 2) for dim in fov]
	println("Nbins = ", Int16.(nbins))
	h5open(outfile, "w") do h5out
		h5grp    = create_group(h5out, "fov_acceptance")
		gen_hist = zeros(Int, nbins)
		acc_hist = zeros(Int, nbins)

		for fn in glob(f_pattern, indir)
			primaries = readh5_todf(fn, "MC", "primaries")
			vertices  = readh5_todf(fn, "MC", "vertices" )

			origin_hist = ATools.hist3d(Matrix(primaries[!, [:x, :y, :z]]),
		 								nbins, bin_limits)
			gen_hist += origin_hist.weights

			# First version only rustpetalo 'first-vertex' equivalent.
			# Filter vertices for those with a track 1 and track 2 vertex in in LXe.
			function first_vertex(vrts::SubDataFrame)::Bool
				any(vrts.volume_id .== 0 .&& vrts.track_id .== 1) &&
				any(vrts.volume_id .== 0 .&& vrts.track_id .== 2)
			end
			vrt_evts  = getproperty.(keys(filter(first_vertex, groupby(vertices, :event_id))), :event_id)
			prim_filt = filter(row -> in(vrt_evts)(row.event_id), primaries)
			lxe_hist  = ATools.hist3d(Matrix(prim_filt[!, [:x, :y, :z]]),
									  nbins, bin_limits)
			acc_hist += lxe_hist.weights
		end
		h5grp["gen"] = gen_hist
		h5grp["acc"] = acc_hist
		println("Size check: ", size(acc_hist))
		# Test write an image file too.
		NReco.write_img(outfile[1:end-2] * "raw", nbins, fov,
						replace(Float32.(acc_hist ./ gen_hist), NaN => 0.0f0))
	end
	
end


function ArgParse.parse_item(::Type{Tuple{T, T, T}}, x::AbstractString) where T <: Real
	Tuple(parse.(T, split(x, ",")))
end


function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
		"--dir", "-d"
			help     = "directory with simulations."
			arg_type = String
			required = true
		"--pattern", "-p"
			help     = "input file pattern."
			arg_type = String
			default  = "MC*.h5"
		"--out", "-o"
			help     = "output file path."
			arg_type = String
			default  = "normalisation.h5"
		"--nbins", "-n"
			help     = "comma separated list of number of bins for fov dimensions."
			arg_type = Tuple{UInt16, UInt16, UInt16}
			default   = (UInt16(151), UInt16(151), UInt16(151))
		"--fov", "-m"
			help     = "comma separated list of FOV edge size."
			arg_type = Tuple{Float32, Float32, Float32}
			default  = (300.0f0, 300.0f0, 300.0f0)
	end
	parse_args(s)
end


if abspath(PROGRAM_FILE) == @__FILE__
	parsed_args = parse_commandline()
	println("Check1: ", parsed_args["nbins"], typeof(parsed_args["nbins"]))
	println("Check2: ", parsed_args["fov"], typeof(parsed_args["fov"]))
	normalisation_histos(parsed_args)
end