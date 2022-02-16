<<<<<<< HEAD
using ATools

using Clustering
using DataFrames
using Distances
using HDF5
using LinearAlgebra
using Statistics
using StatsModels
=======
using DataFrames
using StatsModels
using Clustering
using Statistics
using LinearAlgebra
using HDF5
using ATools
>>>>>>> 4821f5f (First working  version of petalo output reader)

# Selection

"""
    primary_in_lxe(verticesdf::DataFrame)

Select primary photons in LXe.

Row │ event_id track_id parent_id x y z t moved pre_KE post_KE deposited process_id  volume_id
LXe           : volume_id = 0
Primary       : parent_id = 0
"""
function primary_in_lxe(verticesdf::DataFrame)
    vlxe      = select_by_column_value(verticesdf, "volume_id", 0)
    vlxephepr = select_by_column_value(vlxe, "parent_id", 0)
    return vlxephepr
end


"""
	xyz_dot(hitdf::DataFrame, simax::Hit)

Return the dot product between each SiPM in the event and the SiPM of max charge
"""
function xyz_dot(hitdf::DataFrame, simax::DataFrameRow)
	xyzmax = simax[[:x, :y, :z]]
	return dot.(eachrow(hitdf[!, [:x, :y, :z]]), Ref(xyzmax))
end


"""
	sipm_pos(sxyz::DataFrame, index::Integer)

Return the position of the SiPMs in the sipm_xyz database for a given index

"""
function sipm_pos(sxyz::DataFrame, index::Integer)
	return Array(select_by_index(sxyz, "sensor_id", index)[1,2:end])
end


"""
	find_xyz_sipm_qmax(hitdf::DataFrame)
Return the coordinates of the SiPM with maximum charge

"""
function find_xyz_sipm_qmax(hitdf::DataFrame)
	hitdf[argmax(hitdf.q), :]
end


"""
	sipmsel(hdf::DataFrame)

Return two data frames, separating the SiPMs in two groups depending on
the sign of the angle with the SiPM of max charge
"""
function sipmsel(hdf::DataFrame)
	simaxq = find_xyz_sipm_qmax(hdf)
	npr    = xyz_dot(hdf, simaxq)
	mask   = npr .> 0
	return hdf[mask, :], hdf[.!mask, :]
end


"""
	ksipmsel(hdf::DataFrame, ka::Vector{Int64})

Return two data frames, separating the SiPMs in two groups depending on
the value of vector ka (1 or 2)
"""
function ksipmsel(hdf::DataFrame, ka::Vector{Int64})
	return hdf[ka.==2, :], hdf[ka.==1, :]
end


# kmeans algorithm
"""
	get_hits_as_matrix(hitdf::DataFrame)

Given the dataframe hitdf (with fields x,y,z), return the
underlying matrix
"""
function get_hits_as_matrix(hitdf::DataFrame)
	f       = @formula(0 ~ x + y + z)
	f2      = apply_schema(f, schema(f, hitdf))
	_, pred = modelcols(f2, hitdf)
	return transpose(pred)
end

# lors

"""
	reassign_labels(b1::Hit, df1::DataFrame, b2::Hit, df2::DataFrame)
Check the phi of the two reconstruced positions and assign positive
phi to 1. If both are positve, no change is made.
"""
function reassign_labels(b1::Hit, df1::DataFrame, b2::Hit, df2::DataFrame)
	phi1 = atan(b1.y, b1.x)
	if phi1 >= 0.0
		return b1, df1, b2, df2
	end
<<<<<<< HEAD
	return b2, df2, b1, df1
=======
	return b2, df2, b1, df2
>>>>>>> 682fcf9 (Add reassignment of hemisphere labels)
end

"""
	lor_maxq(hitdf::DataFrame)

Compute lors using the SiPM of max charge (maxq algo) to divide the event
in two hemispheres and estimating the lor vertex from barycenter
"""
function lor_maxq(hitdf::DataFrame)
	hq2df, hq1df = sipmsel(hitdf)   # select using maxq
	b1 = baricenter(hq1df)          # baricenters
	b2 = baricenter(hq2df)
	return b1, b2, hq1df, hq2df
end


"""
	lor_kmeans(hitdf::DataFrame)

Compute lors using kmeans clustering.
Returns two estimations of vertices. One based in pure kmeans, the
other in barycenter.
"""
function lor_kmeans(hitdf::DataFrame)
	Mhits = transpose(Matrix(hitdf[!, [:x, :y, :z]]))
	kr    = kmeans(Mhits, 2, weights=hitdf.q)

	## Check that the cluster definitions are good.
	## Silhouette above 0.99 generally means k=2 a good model.
	## TODO Make minimum mean silhouette configuration variable?
	mean_sil = mean(silhouettes(kr, pairwise(SqEuclidean(), Mhits)))

	if mean_sil >= 0.99
		ka = assignments(kr)
		hq2df, hq1df = ksipmsel(hitdf, ka)   # select using kmeans list
		b1 = Hit(kr.centers[:, 1]..., sum(hitdf.q)) # baricenters
		b2 = Hit(kr.centers[:, 2]..., sum(hitdf.q))
		return b1, b2, hq1df, hq2df
	end
	return Hit(0, 0, 0, 0), Hit(0, 0, 0, 0), empty(hitdf), empty(hitdf)
end


"""
	baricenter(hdf::DataFrame)
	returns the barycenter of a cluster of hits
"""
function baricenter(hdf::DataFrame)
	xyzmeans = mean(Matrix(hdf[!, [:x, :y, :z]]), weights(hdf.q), dims=1)
	return Hit(xyzmeans[1], xyzmeans[2], xyzmeans[3], sum(hdf.q))
end


"""
	radial_correction(b::Hit, r::Float32, rsipm::Float32)

Take the estimated radius of interaction (r), and the radius of the sipm
and return the corrected positions
"""
function radial_correction(b::Hit, r::Float32)
	return ATools.calculate_xy(b.x, b.y, r)..., b.z
end


"""
	phistd(hitdf::DataFrame)

Compute the std deviation in phi weighted by charge, e.g:
Sqrt(1/Q Sum_i (phi_i - phi_mean) * qi )
"""
function phistd(hitdf::DataFrame)
	phi = fphi(hitdf)
	return std(phi, FrequencyWeights(hitdf.q), corrected=true)
end


"""
	transverse_angle(hitdf::DataFrame, adjust_cuadrants::Bool=true)
Calculate the transverse (phi) angle for each entry in the DataFrame
and adjust the negative values of the third cuadrant if requested.
"""
function transverse_angle(hitdf::DataFrame, adjust_cuadrants::Bool=true)
	phi = fphi(hitdf)
	if adjust_cuadrants
		quad3  = phi .< -pi / 2
		quad12 = any(phi .> 0)
		if quad12 && any(quad3)
			phi[quad3] .= 2 * pi .+ phi[quad3]
		end
	end
	return phi
end


"""
	xyzstd(hitdf::DataFrame)

Compute the std deviation in x weighted by charge, e.g:
Sqrt(1/Q Sum_i (phi_i - phi_mean) * qi )
"""
function xyzstd(hitdf::DataFrame, column::String="x")
	x = hitdf[!, column]
	return std(x, FrequencyWeights(hitdf.q), corrected=true)
end


"""
    filter_energies
"""
function filter_energies(df::DataFrame, qmin::Float32, qmax::Float32)
    interval = ATools.range_bound(qmin, qmax, ATools.OpenBound)
    filter(x -> interval.(x.q1) .& interval.(x.q2), df)
end

"""
Read petalo sim data and coerce into something resembling abracadabra
"""
function read_petalo(filename::String)
    (hits, parts, sns_xyz, sns_q, sns_tof) = h5open(filename) do h5in
        hits    = readh5_todf(h5in, "MC", "hits")
        parts   = readh5_todf(h5in, "MC", "particles")
        sns_pos = readh5_todf(h5in, "MC", "sns_positions")
        sns_q   = readh5_todf(h5in, "MC", "sns_response")
        sns_tof = readh5_todf(h5in, "MC", "tof_sns_response")
        return hits, parts, sns_pos, sns_q, sns_tof
    end

    ## start coercing
    processes = unique(parts.final_proc)
    volumes   = unique(hits.label)
    ## hits like vertices in abracadabra, need volume_id, parent_id, track_id, process_id
    hits = leftjoin(hits, parts[:, [:event_id, :particle_id, :final_proc]], on=[:event_id, :particle_id])
    transform!(hits, :label => ByRow(x -> findfirst(volumes .== x) - 1) => :volume_id,
        :particle_id => ByRow(x -> in((1, 2))(x) ? 0 : x - 1) => :parent_id,
        :final_proc  => ByRow(x -> findfirst(processes .== x)) => :process_id)
    rename!(hits, :time => :t, :particle_id => :track_id)
    transform!(hits, :event_id => ByRow(Int), renamecols=false)
    ## particles to act like primaries
    filter!(x -> x.primary == 1, parts)
    parts = parts[:, [:event_id, :initial_x, :initial_y, :initial_z,
        :initial_momentum_x, :initial_momentum_y, :initial_momentum_z]]
    rename!(parts, :initial_x => :x, :initial_y => :y, :initial_z => :z,
        :initial_momentum_x => :vx, :initial_momentum_y => :vy, :initial_momentum_z => :vz)
    transform!(parts, :event_id => ByRow(Int), renamecols=false)
    ## sns_q used as total_charge in abracadabra, sum the charge for each
    sns_q = combine(groupby(sns_q, [:event_id, :sensor_id]), :charge => sum => :charge)
    transform!(sns_q, :event_id => ByRow(Int), renamecols=false)
    ## tof_sns_response used as waveform. Need to convert time_bin to time (ns) and
    ## split multi-charge rows into charge rows. Also sensor_id negative for some reason
    transform!(sns_tof, :sensor_id => ByRow(x -> x < 0 ? -x : x), renamecols=false)
    transform!(sns_tof, :time_bin => (x -> x * 5e-3) => :time)
    ## Can't think of a better way to repeat rows according to charge:
    sns_tof = DataFrame(vcat([[x[[:event_id, :sensor_id, :time]] for i in 1:x.charge] for x in eachrow(sns_tof)]...))
    transform!(sns_tof, :event_id => ByRow(Int), renamecols=false)
    ## Get rid of the name in the sensor_xyz
    sns_xyz = sns_xyz[:, [:sensor_id, :x, :y, :z]]

    return ATools.PetaloDF(volumes, processes, sns_xyz, parts, hits, sns_q, sns_tof)
end




"""
    calibration_function
"""
function calibration_function(calibFunc::NReco.CalFunction, rmin::Real, rmax::Real)
    line_pars = h5open(calibFunc.cal_file) do h5cal
        bias   = read_attribute(h5cal[calibFunc.cal_grp], calibFunc.cal_std * "-bias" )
        lconst = read_attribute(h5cal[calibFunc.cal_grp], calibFunc.cal_std * "-const")
        llin   = read_attribute(h5cal[calibFunc.cal_grp], calibFunc.cal_std * "-lin"  )
        return bias, lconst, llin
    end
    cal_func = ATools.predict_interaction_radius(ATools.gpol1(collect(line_pars[2:end])),
        rmin, rmax, line_pars[1])
    return cal_func
end


"""
    calculate_interaction_radius!

"""
function calculate_interaction_radius!(df       ::DataFrame,
                                       predictor::Function,
                                       variable ::String,
                                       rmax     ::Union{Real, Nothing}=nothing)
    if variable == "cstd" && !(:cstd in propertynames(df))
        # We want to use the combined std which isn't necessarily saved in the H5
        comb_std(z, phi) = sqrt(z^2 + (rmax * phi)^2)
        transform!(df, [:zstd1, :phistd1] => ByRow(comb_std) => :cstd1,
            [:zstd2, :phistd2] => ByRow(comb_std) => :cstd2)
    end
    sym1 = Symbol(variable * "1")
    sym2 = Symbol(variable * "2")
    transform!(df, sym1 => predictor => :r1x, sym2 => predictor => :r2x)
    nothing
end
