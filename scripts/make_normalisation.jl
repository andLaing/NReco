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


# Not sure how to write an N dimensional vector to file
# Brute force for now.
struct fov_hist_meta <: ATools.OutputDataset
	nbinx::Int16
	nbiny::Int16
	nbinz::Int16
	minx ::Float32
	maxx ::Float32
	miny ::Float32
	maxy ::Float32
	minz ::Float32
	maxz ::Float32
end


struct lor_hist_meta <: ATools.OutputDataset
	nbinr  ::Int16
	nbinphi::Int16
	nbinz  ::Int16
	nbinth ::Int16
	minr   ::Float32
	maxr   ::Float32
	minphi ::Float32
	maxphi ::Float32
	minz   ::Float32
	maxz   ::Float32
	minth  ::Float32
	maxth  ::Float32
end

function output_meta(data_type::Type{<: ATools.OutputDataset},
					 group    ::HDF5.Group                   ,
					 dset_name::String                       ,
					 values   ::ATools.OutputDataset         )
	out_dtype   = ATools.generate_hdf5_datatype(data_type)
	meta_dspace = dataspace([values])
	meta_dset   = create_dataset(group, dset_name, out_dtype, meta_dspace)
	write_dataset(meta_dset, out_dtype, [values]);
end


function first_vertex(vrts::SubDataFrame)::Bool
	any(vrts.volume_id .== 0 .&& vrts.track_id .== 1) &&
	any(vrts.volume_id .== 0 .&& vrts.track_id .== 2)
end


function filter_energy(gdf::SubDataFrame, minE::Float32)
	any(gdf.volume_id .== 0 .&& gdf.track_id .== 1 .&& gdf.pre_KE .> minE) &
	any(gdf.volume_id .== 0 .&& gdf.track_id .== 2 .&& gdf.pre_KE .> minE)
end


function normalisation_histos(args::Dict{String, Any})
	indir     = args["dir"    ]
	f_pattern = args["pattern"]
	nbins     = args["nbins"  ]
	fov       = args["fov"    ]
	outfile   = args["out"    ]
	lor_out   = args["lor"    ]
	# Always assume FOV centred.
	bin_limits = [(-dim / 2, dim / 2) for dim in fov]
	h5open(outfile, "w") do h5out
		h5grp    = create_group(h5out, "fov_acceptance")
		gen_hist = zeros(Int, nbins)
		acc_hist = zeros(Int, nbins)
		output_meta(fov_hist_meta, h5grp, "metadata",
					fov_hist_meta(nbins..., (bin_limits...)...))

		(lor_cols, lor_bins, lor_lim, lor_gen, lor_acc) = if lor_out
			lor_cols = [:r_lor, :phi_lor, :z_lor, :theta_lor]
			lor_lim  = [(0.0f0, sqrt(bin_limits[1][2]^2 + bin_limits[1][2]^2)),
						(-Float32(pi), Float32(pi)), bin_limits[3], (0.0f0, Float32(pi))]
			lor_bins = (Int(ceil(lor_lim[1][2] / 2)), 20, Int(nbins[3]), 10)
			lor_cols, lor_bins, lor_lim, zeros(Int, lor_bins), zeros(Int, lor_bins)
		else
			nothing, nothing, nothing, nothing, nothing
		end

		for fn in glob(f_pattern, indir)
			primaries = readh5_todf(fn, "MC", "primaries")
			vertices  = readh5_todf(fn, "MC", "vertices" )

			origin_hist = ATools.histNd(Matrix(primaries[!, [:x, :y, :z]]),
		 								nbins, bin_limits)
			gen_hist += origin_hist.weights
			if !isnothing(lor_cols)
				lor_space(x, y, z, vx, vy, vz) = NReco.lor_from_primary(vcat(x, y, z), vcat(vx, vy, vz))
				transform!(primaries, Not(:event_id) => ByRow(lor_space) => lor_cols)
				hist = ATools.histNd(Matrix(primaries[!, lor_cols]), lor_bins, lor_lim)
				lor_gen += hist.weights
			end


			# First version only rustpetalo 'first-vertex' equivalent.
			# Filter vertices for those with a track 1 and track 2 vertex in in LXe.
			filt_vrt  = filter(first_vertex, groupby(vertices, :event_id))
			filt_vrt  = filter(grp -> filter_energy(grp, args["eng"]), filt_vrt)
			vrt_evts  = getproperty.(keys(filt_vrt), :event_id)
			prim_filt = filter(row -> in(vrt_evts)(row.event_id), primaries)
			lxe_hist  = ATools.histNd(Matrix(prim_filt[!, [:x, :y, :z]]),
									  nbins, bin_limits)
			acc_hist += lxe_hist.weights
			if !isnothing(lor_cols)
				hist = ATools.histNd(Matrix(prim_filt[!, lor_cols]), lor_bins, lor_lim)
				lor_acc += hist.weights
			end
		end
		h5grp["gen"] = gen_hist
		h5grp["acc"] = acc_hist
		if !isnothing(lor_cols)
			lorgrp        = create_group(h5out, "lor_acceptance")
			output_meta(lor_hist_meta, lorgrp, "metadata",
						lor_hist_meta(lor_bins..., (lor_lim...)...))
			lorgrp["gen"] = lor_gen
			lorgrp["acc"] = lor_acc
		end

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
		"--eng", "-e"
			help     = "minimum true energy to consider."
			arg_type = Float32
			default  = 0.0f0
		"--lor", "-l"
			help     = "Option to save LOR space binning too."
			action   = :store_true
	end
	parse_args(s)
end


if abspath(PROGRAM_FILE) == @__FILE__
	parsed_args = parse_commandline()
	normalisation_histos(parsed_args)
end