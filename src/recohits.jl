using DataFrames
using Distributions
using ATools


"""
	recohits(         event        ::Integer,
					  sensor_xyz  ::DataFrame,
					  waveform    ::DataFrame,
					  qcut        ::Float32,
					  pde         ::Float32,
					  sigma_tof   ::Float32)

For each event, take the sensor_xyz and waveform tables
and produce a hitdf (xyzqt) data frame, where each SiPM has information of position
(x,y,z), total charge in the SiPM and time (first photon).

"""
function recohits(event       ::Integer,
				  sensor_xyz  ::DataFrame,
				  waveform    ::SubDataFrame,
				  qcut        ::Float32,
				  pde         ::Float32,
				  sigma_tof   ::Float32)

  	# select the waveform of this event
	if nrow(waveform) == 0
		return nothing
	end

	## make a copy of the event.
	wfm = copy(waveform)

	# add a column with probability of surviving pdf cut (pass if prob < pde)
	wfm = wfm[rand(Float32, nrow(wfm)) .< pde, :]
	#wfm[!, "prob"] = rand(Float32, nrow(wfm))

	# add a column of charge
	# (each photon arriving to the SiPm in the waveform has charge of one)
	#wfm[!, "q"]    = Float32.(ones(nrow(wfm)))

	@debug "event =$event, waveform dataframe: size =$(size(wfm))"
	@debug first(wfm, 5)

	# SiPM pass PDE cut if prob < PDE
	#wfmc = select_by_column_value_lt(wfm, "prob", pde)
	if nrow(wfm) == 0
		return nothing
	end

	# add a column of gaussian random numbers
	# representing the smearing of electronics and sensor
	tdist = Normal(0.0, sigma_tof)
	wfm[!, :dt] = Float32.(rand(tdist, nrow(wfm)))

	# add column of smeared times to true times
	transform!(wfm, [:time, :dt] => (+) => :mtime)
	#wfmt = transform(wfmc, [:time, :dt] => (+) => :mtime)

	@debug "waveform after prob cut and smeared time: size =$(size(wfm))"
	@debug first(wfm, 5)

	# group by SiPMs and take minimum times
	# and count the entries to get total number of photons/pe
	wfm_min = combine(groupby(wfm, :sensor_id), :time => minimum => :tmin,
		:mtime => minimum => :trmin, nrow => :q)

	#@debug "waveform after grouping SiPMs and min time: size = $(size(wtmin))"
	#@debug first(wtmin, 5)

	# group by SiPMs and take minimum smeared time
	#wrtmin = combine(groupby(wfmt, :sensor_id), :mtime => minimum)
	#@debug "waveform after grouping SiPMs and min reco time: size = $(size(wrtmin))"
	#@debug first(wrtmin, 5)

	# group by SiPMs and compute the sum of the charge in the SiPMs
	#wtmq = combine(groupby(wfmt, :sensor_id), :q => sum)

	#@debug "waveform after grouping SiPMs and sum charge: size = $(size(wtmq))"
	#@debug first(wtmq, 5)

	#construct qt dataframe
	#wfmx  = DataFrame(sensor_id=wtmq.sensor_id,
	#				  tmin=wtmin.time_minimum,
	#				  trmin=wrtmin.mtime_minimum,
	#				  q=wtmq.q_sum)

  	@debug "DF with time and q: size = $(size(wfm_min))"
	@debug first(wfm_min, 5)

	# cut on total charge (> qcut)
	filter!(x -> x.q .> qcut, wfm_min)
	#qdft   = wfmx[wfmx.q .>qcut,:]

	xyzqt = leftjoin(wfm_min, sensor_xyz, on=:sensor_id)
	disallowmissing!(xyzqt, [:x, :y, :z])

	#select the id's of the sipms with charge above threshold
	#sids = qdft[!,:sensor_id]

	#compute positions of the SiPMs
	#pos   = sipm_pos.((sensor_xyz,),sids)
	#x = [p[1] for p in pos]
	#y = [p[2] for p in pos]
	#z = [p[3] for p in pos]

	# Construct hit data frame
	#xyzqt   = DataFrame(x=x,y=y,z=z,
	#				    tmin=qdft.tmin,
	#				    trmin=qdft.trmin,
	#				    q=qdft.q)

	return xyzqt[!, [:x, :y, :z, :tmin, :trmin, :q]]
end
