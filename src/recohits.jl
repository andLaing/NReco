using DataFrames
using Distributions


"""
	recohits(event       ::Integer,
		 sensor_xyz  ::DataFrame,
		 waveform    ::DataFrame,
		 qcut        ::Float32,
		 pde         ::Float32,
		 sigma_tof   ::Float32)

For each event, take the sensor_xyz and waveform tables
and produce a hitdf (xyzqt) data frame, where each SiPM has information of position
(x,y,z), total charge in the SiPM and time (first photon).

"""
function recohits(event     ::Integer,
		  sensor_xyz::DataFrame,
		  waveform  ::SubDataFrame,
		  qcut      ::Float32,
		  pde       ::Float32,
		  sigma_tof ::Float32)

	wfm = waveform[rand(Float32, nrow(waveform)) .< pde, :]

	@debug "event =$event, waveform dataframe: size =$(size(wfm))"
	@debug first(wfm, 5)

	if nrow(wfm) == 0
		return nothing
	end

	# add a column of gaussian random numbers
	# representing the smearing of electronics and sensor
	tdist = Normal(0.0f0, sigma_tof)
	wfm[!, :dt] = rand(tdist, nrow(wfm))

	# add column of smeared times to true times
	transform!(wfm, [:time, :dt] => (+) => :mtime)

	@debug "waveform after prob cut and smeared time: size =$(size(wfm))"
	@debug first(wfm, 5)

	# For each sensor get tmin, tmin smeared and the number of photons.
	wfm_min = combine(groupby(wfm, :sensor_id), :time => minimum => :tmin,
		:mtime => minimum => :trmin, nrow => :q)

  	@debug "DF with time and q: size = $(size(wfm_min))"
	@debug first(wfm_min, 5)

	# cut on total charge (> qcut)
	filter!(x -> x.q .> qcut, wfm_min)

	# Get the positions of each of the selected sensors.
	xyzqt = leftjoin(wfm_min, sensor_xyz, on=:sensor_id)
	disallowmissing!(xyzqt, [:x, :y, :z])

	return xyzqt[!, [:x, :y, :z, :tmin, :trmin, :q]]
end


"""
	select_sensors(waveform::SubDataFrame, qcut     ::Float32,
	               pde     ::Float32     , sigma_tof::Float32)

Smear and select the SiPMs.
"""
function select_sensors(waveform::SubDataFrame, qcut     ::Float32,
	                    pde     ::Float32     , sigma_tof::Float32)
	# Filter according to probability and PDE (pass if prob < pde)
	wfm = waveform[rand(Float32, nrow(waveform)) .< pde, :]

	if nrow(wfm) == 0
		return nothing
	end

	# add a column of gaussian random numbers
	# representing the smearing of electronics and sensor
	tdist = Normal(0.0f0, sigma_tof)
	wfm[!, :dt] = rand(tdist, nrow(wfm))

	# add column of smeared times to true times
	transform!(wfm, [:time, :dt] => (+) => :mtime)
	
	 # cut on total charge (> qcut)
	return filter(grp -> nrow(grp) .> qcut, groupby(wfm, :sensor_id))
end


"""
	sensor_positions(waveform::GroupedDataFrame, sensor_xyz::DataFrame)

Get the minima for each sensor and the summed charge and combine with positions.
"""
function sensor_positions(waveform::GroupedDataFrame, sensor_xyz::DataFrame)
	wfm_min = combine(waveform, :time => minimum => :tmin,
	                  :mtime => minimum => :trmin, nrow => :q)
	# Get the positions of each of the selected sensors.
	xyzqt = leftjoin(wfm_min, sensor_xyz, on=:sensor_id)
	disallowmissing(xyzqt, [:x, :y, :z])
end


"""
	average_first_hits(waveform::GroupedDataFrame, sensors::Vector{Int64}, dc::DetConf)
calculate the average of the first ntof photons in a set of sensors.
"""
function average_first_hits(waveform::GroupedDataFrame, sensors::Vector{Int64}, dc::DetConf)
	nsens        = min(length(sensors), dc.ntof)
	h_sens       = in(sensors)
	sensor_times = combine(filter(grp -> h_sens(first(grp.sensor_id)), waveform), :mtime)
	sens_t       = sort(sensor_times.mtime)[1:nsens]
	mask         = sens_t .<= sens_t[1] + dc.sigma_tof * dc.nsigma
	mean(sens_t[mask])
end


"""
	split_hemispheres(waveform  ::SubDataFrame,
					  sensor_xyz::DataFrame   ,
					  dc        ::DetConf     ,
					  q_bound   ::Function    ,
					  lor_algo  ::Function    )

Reconstruct the SiPM hits (smearing etc) and split into hemispheres.
"""
function split_hemispheres(waveform  ::SubDataFrame,
	                       sensor_xyz::DataFrame   ,
						   dc        ::DetConf     ,
	                       q_bound   ::Function    ,
						   lor_algo  ::Function    )
	reco_wvf = select_sensors(waveform, dc.ecut, dc.pde, dc.sigma_tof)

	if isnothing(reco_wvf) || length(keys(reco_wvf)) < 2
		return nothing
	end

	## Get the minima with positions
	xyzqt = sensor_positions(reco_wvf, sensor_xyz)

	## Split into hemispheres
	b1, b2, hq1df, hq2df = lor_algo(xyzqt)

	if q_bound(b1.q) && q_bound(b2.q)
		## Assign numbers according to positive/negative phi.
		b1, hq1df, b2, hq2df = reassign_labels(b1, hq1df, b2, hq2df)
		
		## Get average start times for hemispheres
		ta1 = average_first_hits(reco_wvf, hq1df.sensor_id, dc)
		ta2 = average_first_hits(reco_wvf, hq2df.sensor_id, dc)
		return b1, b2, ta1, ta2, hq1df, hq2df
	end
	return nothing
end
