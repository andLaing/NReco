using DataFrames
using Distributions
using StatsBase
using ATools


"""
	true_xyz(b1::Hit, b2::Hit, df1::DataFrame, df2::DataFrame)

Return the sorted distance between baricenter and true
"""
function true_xyz(b1::Hit, df1::DataFrame, df2::DataFrame)
	# distance between baricenter and true
	d1 = dxyz([b1.x, b1.y, b1.z], [df1.x[1], df1.y[1], df1.z[1]])
	d2 = dxyz([b1.x, b1.y, b1.z], [df2.x[1], df2.y[1], df2.z[1]])

	if d2 < d1
		xt2 = [df1.x[1], df1.y[1], df1.z[1]]
		xt1 = [df2.x[1], df2.y[1], df2.z[1]]
	else
		xt1 = [df1.x[1], df1.y[1], df1.z[1]]
		xt2 = [df2.x[1], df2.y[1], df2.z[1]]
	end
	return xt1, xt2
end

"""
	recovent(event       ::Integer,
	dc          		 ::DetConf,
	df1         		 ::DataFrame,
	df2         		 ::DataFrame,
	primaries   		 ::SubDataFrame,
	sensor_xyz  		 ::DataFrame,
	waveform    		 ::SubDataFrame,
	lor_algo    		 ::Function)

Return a dictionary with the variables characterising the event.
"""
function recovent(event     ::Integer,
				  dc        ::DetConf,
				  df1       ::DataFrame,
				  df2       ::DataFrame,
				  primaries ::SubDataFrame,
				  sensor_xyz::DataFrame,
				  waveform  ::SubDataFrame,
				  lor_algo  ::Function)
	n3d = Dict()

	n3d[:phot1] =  df1.process_id[1] == 1
	n3d[:phot2] =  df2.process_id[1] == 1

	# Primary particle origins
	n3d[:xs] = primaries.x[1]
	n3d[:ys] = primaries.y[1]
	n3d[:zs] = primaries.z[1]
	n3d[:ux] = primaries.vx[1]
	n3d[:uy] = primaries.vy[1]
	n3d[:uz] = primaries.vz[1]

	# hit dataframe
	hitdf = recohits(event,sensor_xyz, waveform, dc.ecut, dc.pde, dc.sigma_tof)

	if hitdf === nothing
	@warn "Warning, hidtf evaluates to nothing for event = $event"
	return nothing
	end

	if nrow(hitdf) < 2
	@warn "Warning, hidtf is <2 for event = $event"
	return nothing
	end

	@info " hit dataframe: size = $size(hitdf)"
	@debug first(hitdf, 5)

	# reconstruct (x,y) : barycenter
	b1, b2, hq1df, hq2df = lor_algo(hitdf)

	n3d[:xr1]    = b1.x
	n3d[:yr1]    = b1.y
	n3d[:zr1]    = b1.z
	n3d[:nsipm1] = nrow(hq1df)

	n3d[:xr2]    = b2.x
	n3d[:yr2]    = b2.y
	n3d[:zr2]    = b2.z
	n3d[:nsipm2] = nrow(hq1df)

	@info " barycenters" b1 b2

	# total charge
	n3d[:q1] = sum(hq1df.q)
	n3d[:q2] = sum(hq2df.q)

	@info " total charge: q1 = $(n3d[:q1]), q2 = $(n3d[:q2])"
	if n3d[:q1] < dc.qmin || n3d[:q1] > dc.qmax
	@info "Warning, q1 is $(n3d[:q1]) for event $event"
	return nothing
	end
	if n3d[:q2] < dc.qmin || n3d[:q2] > dc.qmax
	@info "Warning, q2 is $(n3d[:q2]) for event $event"
	return nothing
	end

	# find true position (and correlate with barycenter)
	xt1, xt2 = true_xyz(b1, df1, df2)

	n3d[:xt1] = xt1[1]
	n3d[:yt1] = xt1[2]
	n3d[:zt1] = xt1[3]
	n3d[:t1]  = minimum(hq1df.tmin)
	n3d[:tr1] = minimum(hq1df.trmin)

	n3d[:xt2] = xt2[1]
	n3d[:yt2] = xt2[2]
	n3d[:zt2] = xt2[3]
	n3d[:t2]  = minimum(hq2df.tmin)
	n3d[:tr2] = minimum(hq2df.trmin)

	# find r1 and r2 (from True info)
	n3d[:r1] = rxy(xt1[1], xt1[2])
	n3d[:r2] = rxy(xt2[1], xt2[2])
	@info " True position in hemisphere 1" xt1
	@info " True position in hemisphere 1" xt2

	# Compute phistd and zstd1
	int1_weights    = FrequencyWeights(hq1df.q)
	phi_values      = fphi(hq1df)
	n3d[:phistd1]   = std(phi_values, int1_weights, corrected=true)
	n3d[:zstd1]     = std(hq1df.z   , int1_weights, corrected=true)
	n3d[:widz1]     = maximum(hq1df.z) - minimum(hq1df.z)
	n3d[:widphi1]   = maximum(phi_values) - minimum(phi_values)
	n3d[:corrzphi1] = cor(hcat(hq1df.z, phi_values), int1_weights)[1,2]
	int2_weights    = FrequencyWeights(hq2df.q)
	phi_values      = fphi(hq2df)
	n3d[:phistd2]   = std(phi_values, int2_weights, corrected=true)
	n3d[:zstd2]     = std(hq2df.z   , int2_weights, corrected=true)
	n3d[:widz2]     = maximum(hq2df.z) - minimum(hq2df.z)
	n3d[:widphi2]   = maximum(phi_values) - minimum(phi_values)
	n3d[:corrzphi2] = cor(hcat(hq2df.z, phi_values), int2_weights)[1,2]
	@info " phistd1 = $(n3d[:phistd1]), zstd1 = $(n3d[:zstd1])"
	@info " phistd2 = $(n3d[:phistd2]), zstd2 = $(n3d[:zstd2])"

	# New (x,y) positions estimated from r1, r2
	n3d[:x1], n3d[:y1], n3d[:z1] = radial_correction(b1, n3d[:r1])
	n3d[:x2], n3d[:y2], n3d[:z2] = radial_correction(b2, n3d[:r2])

	@info " New (x,y,z) positions estimated from r1, r2 & r1q, r2q"
	@info " from r1:  x1 = $(n3d[:x1]), y1=$(n3d[:y1]), z1=$(n3d[:z1])"
	@info " from r2:  x2 = $(n3d[:x2]), y1=$(n3d[:y2]), z1=$(n3d[:z2])"

	ntof1 = min(dc.ntof, nrow(hq1df))
	ntof2 = min(dc.ntof, nrow(hq2df))

	# sort reco times in ascending order
	t1s = sort(hq1df.trmin)
	t2s = sort(hq2df.trmin)
	# take average
	n3d[:ta1] = mean(t1s[1:ntof1])
	n3d[:ta2] = mean(t2s[1:ntof2])

	@info " time"
	@info " true:  t1 = $(n3d[:t1]), t2=$(n3d[:t2])"
	@info " smeared:  tr1 = $(n3d[:tr1]), tr2=$(n3d[:tr2])"
	@info " averaged:  ta1 = $(n3d[:ta1]), ta2=$(n3d[:ta2])"

	ht1  = select_by_column_value(hq1df, "tmin", n3d[:t1])
	ht2  = select_by_column_value(hq2df, "tmin", n3d[:t2])

	n3d[:xb1] = ht1.x[1]
	n3d[:yb1] = ht1.y[1]
	n3d[:zb1] = ht1.z[1]
	n3d[:xb2] = ht2.x[1]
	n3d[:yb2] = ht2.y[1]
	n3d[:zb2] = ht2.z[1]

	return hq1df, hq2df, n3d
end


"""
	nema_analysis!(     event       ::Integer,
						dc          ::DetConf,
						df1         ::DataFrame,
						df2         ::DataFrame,
						primaries   ::SubDataFrame,
						sensor_xyz  ::DataFrame,
						waveform    ::SubDataFrame,
						lor_algo    ::Function,
						n3d         ::Dict)

Fill the DataFrame for nema analysis
"""
function nema_dict!(n3df      ::DataFrame   ,
					event     ::Integer     ,
					dc        ::DetConf     ,
					df1       ::DataFrame   ,
					df2       ::DataFrame   ,
					primaries ::SubDataFrame,
					sensor_xyz::DataFrame   ,
					waveform  ::SubDataFrame,
					lor_algo  ::Function    )

	result = recovent(event, dc, df1, df2, primaries, sensor_xyz, waveform, lor_algo)

	if result !== nothing
		push!(n3df, result[3])
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

	# define data dictionary
	# TODO: It would be better to have a struct or similar with the keys and TYPES
	# This would protect the filling and maybe improve the push time.
	n3df = DataFrame(:phot1     => Bool[]   , :phot2     => Bool[]   ,
					 :nsipm1    => Int64[]  , :nsipm2    => Int64[]  ,
					 :q1        => Float32[], :q2        => Float32[],
					 :r1        => Float32[], :r2        => Float32[],
					 :phistd1   => Float32[], :phistd2   => Float32[],
					 :widphi1   => Float32[], :widphi2   => Float32[],
					 :zstd1     => Float32[], :zstd2     => Float32[],
					 :widz1     => Float32[], :widz2     => Float32[],
					 :corrzphi1 => Float32[], :corrzphi2 => Float32[],
					 :xt1       => Float32[], :xt2       => Float32[],
					 :yt1       => Float32[], :yt2       => Float32[],
					 :zt1       => Float32[], :zt2       => Float32[],
					 :t1        => Float32[], :t2        => Float32[],
					 :x1        => Float32[], :x2        => Float32[],
					 :y1        => Float32[], :y2        => Float32[],
					 :z1        => Float32[], :z2        => Float32[],
					 :xr1       => Float32[], :xr2       => Float32[],
					 :yr1       => Float32[], :yr2       => Float32[],
					 :zr1       => Float32[], :zr2       => Float32[],
					 :tr1       => Float32[], :tr2       => Float32[],
					 :xb1       => Float32[], :xb2       => Float32[],
					 :yb1       => Float32[], :yb2       => Float32[],
					 :zb1       => Float32[], :zb2       => Float32[],
					 :ta1       => Float32[], :ta2       => Float32[],
					 :xs        => Float32[], :ux        => Float32[],
					 :ys        => Float32[], :uy        => Float32[],
					 :zs        => Float32[], :uz        => Float32[])

	for file in files[file_i:file_l]               # loop on files
		println("reading file = ", file)
		pdf = read_abc(file)            # read file
		dfs = primary_in_lxe(pdf.vertices)       # primary photons in LXe

		## We are interested in events with two primary photons in LXe
		grp_dfs = filter(x -> any(x.track_id .== 1) && any(x.track_id .== 2),
							groupby(dfs, :event_id))
		primaries = groupby(pdf.primaries, :event_id)
		waveforms = groupby(pdf.waveform , :event_id)

		for (event, vdf) in pairs(grp_dfs)
			df1 = vdf[vdf.track_id .== 1, :]
            df2 = vdf[vdf.track_id .== 2, :]

			nema_dict!(n3df, event.event_id, dconf, df1, df2,
					   primaries[values(event)], pdf.sensor_xyz,
					   waveforms[values(event)], lor_algo)
    	end
	end
	return n3df
end
