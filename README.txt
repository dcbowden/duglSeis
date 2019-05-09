
The matlab scripts here were designed for accessing, plotting and working with
seismic data. They are custom tailored for collaboration on a particular experiment,
from a seismic array in the Homestake Gold Mine, in South Dakota, US. 
(DUGL -> Deep Underground Seismic Lab). 

I do *NOT* claim these scripts will work out-of-the-box for any other data, but hope
others may find it useful to dissect or adapt for other purposes (for example, to
learn how to filter data, apply an instrument response, etc.). 

Many subroutines were taken or adapted from others, and I have tried to document 
this where needed. Notable pieces not included here are:
1) irisFetch.m    - Provided by IRIS to access their own webservers
		    (my own duglFetch.m calls were structured to mimic this)
2) TauP           - To estimate raypaths through the earth. My functions are just
		    a wrapper to the java functions provided by http://www.seis.sc.edu
		    
Other functions were adapted or inspired by SEIZMO by Garrett Euler,
and various scripts from Rob Porritt. 

All scripts assume a "sac-like" object structure, and are based on the readsac() and 
writesac() scripts, taken from Sandia Labs and with their own copyright rules.

I welcome any comments or bug reports, or tell me if I missed a proper citation! 
I do not plan to further develop or generalize these, however, and generally 
recommend migrating to ObsPy.

-----------------------------------------------------------------------------

All scripts work with data as "SAC" objects, meant to mimic the header style of
the program SAC. For example: 
s = duglFetch.Read('X6','1700','00','HHZ','2015 108 16:30:00','2015 109 06:29:00');
   s.DATA1     timeseries
   s.STLA      station lat
   s.STLO      station lon
   s.DELTA     timestep of samples
   ...etc...
All subsequent routines (to apply a filter, for example) are passed this sac object.

Last updated June 2018, Daniel Bowden, dbowden@caltech.edu, dcbowden@gmail.com

------------------------------------------------------------------------------

There are 4 main classes, each contain a number of functions:

duglFetch
	Init         (set up matlab path)
	Read        
	ReadEvent

duglPlot
	p1   	     (plot each on own axes)
	p2   	     (plot with shared axes)
	prs  	     (plot record section, sort by distance)
	power_spect 
	Map3D
	PlotTauP     (plot a circle Earth with raypaths)
	TimeTauP     (plot ray arrival times on timeseries)
	idTrace	     (identifies line under cursor)
	ForceSeconds (adjust x-axis)
	ResetTicks   (adjust x-axis)

duglProc
	bp   		(band-pass, wrapper for Filt)
	Filt 		(basic butterworth filter)
	applyResp 	(assumes pole-zero file present)
	getPZ		(sub-function to get poles and zeros)
	zpk2cmplx	(pole-zero to spectrum for convolution)
	TauP		(raypaths through the Earth)
	EstimateDelay   (estimate arrivals based on TauP)
	DistVicenty	(distance between lat/lon)
	power_spect	(get power spectrum. looks up instrument response unless provided)

duglSet
	No functions - hardcoded paths for data and metadata
	(this is how the scripts know station names, locations, etc.)

	!! Modify this if data is moved, or if you want to adapt these for
	   a different dataset.

There are also a couple other, more specific set of routines:

duglMod
	Used for propogating delays and attenuations around the 3D array
	Method is 'in development'

duglNoise
	Pure-matlab implementation of ambient noise cross correlations
	See "NoiseCC" folder for example calling scripts.

----------------------------------------------------
The duglFetch script was inspired by the irisFetch script, which connects
to IRIS servers. The example queries earthquake source information, but it 
can also fetch waveform data from any other station in the US (and much of the
world) with the irisFetch.Traces command. A script to convert the Traces
object to my "sac" object is included.

