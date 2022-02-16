using ATools

using DataFrames
using Distributions
using StatsBase


"""
	true_interaction(b1     ::Hit          , b2      ::Hit         ,
	                 volumes::Vector{Int64}, vertices::SubDataFrame)

Extract MC truth information about the LXe interactions.
Only considers vertices in volumes.
Returns the df row for the vertex closes to the barycentre in each case.
"""
function true_interaction(b1::Hit, b2::Hit,
	                      volumes::Vector{Int64},
						  vertices::SubDataFrame)::Union{Nothing, Tuple{DataFrameRow, DataFrameRow}}
	function dist_to_hit(h1::Hit)::Function
		hit_xyz = getfield.((h1,), [:x, :y, :z])
		function dist(x::Float32, y::Float32, z::Float32)::Float32
			dxyz(hit_xyz, vcat(x, y, z))
		end
		return dist
	end

	vrts = filter(df -> in(volumes).(df.volume_id), vertices)
	if isempty(vrts)
		@warn "Warning, no valid MC true vertex for current event."
		# TODO deal with these events correctly without getting weird rmin/rmax.
		return nothing
	end

	vrts = combine(df -> df[argmin(df.t), :], groupby(vrts, :track_id))
	
	transform!(vrts, [:x, :y, :z] => ByRow(dist_to_hit(b1)) => :dist1,
	                 [:x, :y, :z] => ByRow(dist_to_hit(b2)) => :dist2)

	cols = [:process_id, :x, :y, :z, :t, :pre_KE]
	## Allow them to be the same for now, should only happen
	## for bad reconstruction with no separation.
	row1 = argmin(vrts.dist1)
	row2 = argmin(vrts.dist2)

	return vrts[row1, cols], vrts[row2, cols]
end


"""
	check_valid_hits(hitdf::DataFrame, event::Int64)

Check the validity of the event hitmaps for continued processing
Returns true only of hitdf not empty and has at least 2 rows.
"""
function check_valid_hits(hitdf::Union{DataFrame, Nothing}, event::Int64)
	if hitdf === nothing
		@warn "Warning, hitdf evaluates to nothing for event = $event"
		return false
	end

	if nrow(hitdf) < 2
		@warn "Warning, hitdf is <2 for event = $event"
		return false
	end
	return true
end


"""
	check_valid_charge(h1::DataFrame, h2::DataFrame, dc::DetConf)

Check if the hemispheres have detected charge above threshold.
Requires that the sum of detected charge in each hemisphere
be in (dc.qmin, dc.qmax).
"""
function check_valid_charge(h1::DataFrame, h2::DataFrame, dc::DetConf)
	q1 = sum(h1.q)
	q2 = sum(h2.q)
	@info " total charge: q1 = $q1, q2 = $q2"
	if q1 < dc.qmin || q1 > dc.qmax
		@info "Warning, q1 is $q1"
		return false
	end

	if q2 < dc.qmin || q2 > dc.qmax
		@info "Warning, q2 is $q2"
		return false
	end
	return true
end


"""
	recovent(event     ::Integer      ,
	         dc        ::DetConf      ,
			 vertices  ::SubDataFrame ,
			 primaries ::SubDataFrame ,
			 sensor_xyz::DataFrame    ,
			 volumes   ::Vector{Int64},
			 waveform  ::SubDataFrame ,
			 lor_algo  ::Function     )

Return a dictionary with the variables characterising the event.
"""
function recovent(event     ::Integer,
		  dc        ::DetConf,
		  vertices  ::SubDataFrame,
		  primaries ::SubDataFrame,
		  sensor_xyz::DataFrame,
		  volumes   ::Vector{Int64},
		  waveform  ::SubDataFrame,
		  lor_algo  ::Function)
	q_bound     = ATools.range_bound(dc.qmin, dc.qmax, ATools.OpenBound)
	hemispheres = split_hemispheres(waveform, sensor_xyz,
	                                dc, q_bound, lor_algo)
	if isnothing(hemispheres)
		return nothing
	end
	b1, b2, ta1, ta2, hq1df, hq2df = hemispheres

	# find true position (and correlate with barycenter)
	mc_truth = true_interaction(b1, b2, volumes, vertices)
	if isnothing(mc_truth)
		return mc_truth
	end
	df1, df2 = mc_truth
	@info " True position in hemisphere 1" Array(df1[[:x, :y, :z]])
	@info " True position in hemisphere 2" Array(df2[[:x, :y, :z]])

	## Weights and phi positions
	int1_weights = FrequencyWeights(hq1df.q)
	phi1_values  = transverse_angle(hq1df)
	int2_weights = FrequencyWeights(hq2df.q)
	phi2_values  = transverse_angle(hq2df)

	true_rad1 = rxy(df1.x, df1.y)
	true_rad2 = rxy(df2.x, df2.y)
	# New (x,y) positions estimated from r1, r2
	xyz1 = radial_correction(b1, true_rad1)
	xyz2 = radial_correction(b2, true_rad2)

	ht1   = hq1df[argmin(hq1df.tmin), :]
	ht2   = hq2df[argmin(hq2df.tmin), :]
	tmin1 = ht1.tmin
	tmin2 = ht2.tmin

	@info " New (x,y,z) positions estimated from r1, r2 & r1q, r2q"
	@info " from r1:  x1 = $(xyz1[1]), y1=$(xyz1[2]), z1=$(xyz1[3])"
	@info " from r2:  x2 = $(xyz2[1]), y1=$(xyz2[2]), z1=$(xyz2[3])"

	@info " hit dataframe: size = $size(hitdf)"
	@debug first(hitdf, 5)

	reduced_event = ATools.EventParameters(
	event_id=event,
	phot1 = df1.process_id == 1,
	phot2 = df2.process_id == 1,
	# Primary particle origins
	xs = primaries.x[1],
	ys = primaries.y[1],
	zs = primaries.z[1],
	ux = primaries.vx[1],
	uy = primaries.vy[1],
	uz = primaries.vz[1],
	# barycenter on SiPM plane
	xr1    = b1.x,
	yr1    = b1.y,
	zr1    = b1.z,
	nsipm1 = nrow(hq1df),
	xr2    = b2.x,
	yr2    = b2.y,
	zr2    = b2.z,
	nsipm2 = nrow(hq1df),
	# total charge
	q1 = sum(hq1df.q),
	q2 = sum(hq2df.q),
	# Gamma energy at interaction
	E1 = df1.pre_KE,
	E2 = df2.pre_KE,
	# True interaction positions
	xt1 = df1.x,
	yt1 = df1.y,
	zt1 = df1.z,
	ti1 = df1.t,
	t1  = tmin1,
	tr1 = minimum(hq1df.trmin),
	xt2 = df2.x,
	yt2 = df2.y,
	zt2 = df2.z,
	ti2 = df2.t,
	t2  = tmin2,
	tr2 = minimum(hq2df.trmin),
	# find r1 and r2 (from True info)
	r1 = true_rad1,
	r2 = true_rad2,
	# Compute phistd and zstd of detected light
	phistd1   = std(phi1_values, int1_weights, corrected=true),
	zstd1     = std(hq1df.z    , int1_weights, corrected=true),
	widz1     = maximum(hq1df.z) - minimum(hq1df.z),
	widphi1   = maximum(phi1_values) - minimum(phi1_values),
	corrzphi1 = cor(hcat(hq1df.z, phi1_values), int1_weights)[1,2],
	phistd2   = std(phi2_values, int2_weights, corrected=true),
	zstd2     = std(hq2df.z    , int2_weights, corrected=true),
	widz2     = maximum(hq2df.z) - minimum(hq2df.z),
	widphi2   = maximum(phi2_values) - minimum(phi2_values),
	corrzphi2 = cor(hcat(hq2df.z, phi2_values), int2_weights)[1,2],
	# XYZ given true radius and barycentre
	x1 = xyz1[1],
	y1 = xyz1[2],
	z1 = xyz1[3],
	x2 = xyz2[1],
	y2 = xyz2[2],
	z2 = xyz2[3],
	# Average arrival time of first ntof SiPMs
	ta1 = ta1,
	ta2 = ta2,
	# Positions of SiPM with tmin
	xb1 = ht1.x[1],
	yb1 = ht1.y[1],
	zb1 = ht1.z[1],
	xb2 = ht2.x[1],
	yb2 = ht2.y[1],
	zb2 = ht2.z[1])

	return hq1df, hq2df, reduced_event
end


"""
	nema_analysis!(data_vec  ::Vector{ATools.EventParameters},
				   event     ::Integer,
				   dc        ::DetConf,
				   sensor_xyz::DataFrame,
				   volumes   ::Vector{Int64},
				   waveform  ::SubDataFrame,
				   primaries ::SubDataFrame,
				   vertices  ::SubDataFrame,
				   lor_algo  ::Function)

Fill the DataFrame for nema analysis
"""
function nema_dict!(data_vec      ::Vector{ATools.EventParameters}   ,
		    event     ::Integer     ,
		    dc        ::DetConf     ,
			sensor_xyz::DataFrame   ,
			volumes   ::Vector{Int64},
			waveform  ::SubDataFrame,
			primaries ::SubDataFrame,
			vertices  ::SubDataFrame,
			lor_algo  ::Function    )

	result = recovent(event, dc, vertices, primaries, sensor_xyz, volumes, waveform, lor_algo)

	if result !== nothing
		push!(data_vec, result[3])
	end
end


"""
	function nemareco(files    ::Vector{String},
					  dconf    ::DetConf,
		              file_i   ::Integer=1,
					  file_l   ::Integer=1,
					  phot     ::Bool=true
					  lor_algo ::Function=lor_maxq)

Return the evtdf DataFrame


   nsipm1  => Number of sipms in hemisphere 1 (nsipm2 for hemisphere 2)
   q1      => Total charge in hemisphere
   E1      => Energy of the gamma at interaction.
   r1      => True radius, gamma interaction point,
   phistd  =>  Phi standard deviations (weigthed by charge)
   zstd    => Z standard deviations (weigthed by charge)

   xs, ys, zs    => Position of the pair of gammas (in target),
   ux, uy, uz    => Direction vectors of gammas

   xt1, yt1, zt1 => True position of gammas,
   x1,   y1, z1  => Reco position of gammas (from barycenter and r1),
   xr1, yr1, zr  => Baricenter
   xb1, yb1, zb1 => Position of the sipm giving time stamp

   t1            => Time stamp of first photon,
   tr1           => Time stamp of first photon, after pdf and smearing
   ta1           => Average of time stamps, first five photons

"""
function nemareco(files    ::Vector{String},
				  dconf    ::DetConf,
	              file_i   ::Integer=1,
				  file_l   ::Integer=1,
				  lor_algo ::Function=lor_maxq)

	output_vector = ATools.EventParameters[]
	total_events  = 0

	for file in files[file_i:file_l]
		println("reading file = ", file)
		#pdf = read_abc(file)            # read file
		pdf = NReco.read_petalo(file)

		# count the total number of events generated.
		# This, like many other sections assumes primary generation of gammas.
		total_events += nrow(pdf.primaries)

		vertices  = groupby(filter(df -> df.parent_id .== 0, pdf.vertices), :event_id)
		primaries = groupby(pdf.primaries, :event_id)
		waveforms = groupby(pdf.waveform , :event_id)
		## For true vertex saving, can have detected charge from
		## gamma vertices in LXe and Steel_1 (naming stable?)
		volumes   = findall(pdf.volume_names .== "ACTIVE") .- 1

		for (event, wvf) in pairs(waveforms)
			nema_dict!(output_vector, event.event_id, dconf, pdf.sensor_xyz, volumes,
			wvf, primaries[values(event)], vertices[values(event)], lor_algo)
		end
	end
	return total_events, output_vector
end
