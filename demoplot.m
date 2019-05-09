
%%%%%%
%%%%%%
%
%
%Three main classes, which each contain a number of functions:
%
%duglFetch
%	Init
%	Read
%	ReadEvent
%duglPlot
%	p1   	     (plot each on own axes)
%	p2   	     (plot with shared axes)
%	prs  	     (plot record section, sort by distance)
%	power_spect 
%	Map3D
%	PlotTauP     (plot a circle Earth with raypaths)
%	TimeTauP     (plot ray arrival times on timeseries)
%	idTrace	     (identifies line under cursor)
%	ForceSeconds (adjust x-axis)
%	ResetTicks   (adjust x-axis)
%duglProc
%	bp   		(band-pass, wrapper for Filt)
%	Filt 		(basic butterworth filter)
%	applyResp 	(assumes pole-zero file present)
%	getPZ		(sub-function to get poles and zeros)
%	zpk2cmplx	(pole-zero to spectrum for convolution)
%	TauP		(raypaths through the Earth)
%	EstimateDelay   (estimate arrivals based on TauP)
%	DistVicenty	(distance between lat/lon)
%
%Also, will need to set up local paths in "duglSet.m"
%


%%

%%
% IRISFETCH
% Provided by IRIS, allows for fetching of event data and waveforms from
% other seismic data.
%
% Add IRIS to path, or use duglFetch.Init
%
% Here, I'm querying an event I know roughly the time and magnitude:
%
%ev = irisFetch.Events('MinimumMagnitude',7.7,'startTime','2015-04-25 06:11:20','endTime','2015-04-25 06:13:26'); %nepal
ev = irisFetch.Events('MinimumMagnitude',6.1,'startTime','2015-04-24 13:50:16','endTime','2015-04-24 13:59:16') %canada



%%
% DUGLFETCH
% Inspiried by irisFetch format, will read specified station and component
% for a given timeperiod or event.


% If "ev" object is passed, that information will be added to the "s"
% objects.
% The "s" object holds header information meant to match SAC files
%


s = duglFetch.ReadEvent('X6','300|D4850','00','HHZ',ev,['+1h']);


%%
%
% DUGLPLOT 
% Various plotting routines. 


% Meant to mimic SAC style commands:
%  -- p1 to plot each trace on its own axes
%  -- p2 to plot traces on shared axes
%  
% Also:
%  -- power_spect
%  -- Map3D

figure(1)
duglPlot.p1(s);

figure(2)
duglPlot.Map3D;

%%
%
% DUGLPROC
% Tools for filtering, calculating spectrum, instrument response
%

figure(3)
scorr = duglProc.applyResp(s);
duglPlot.power_spect(scorr);
% 
% !
% Note that this spectrum contains an earthquake event, so looks different than
% background noise models
% 

figure(4)
% Read in all station's data
s = duglFetch.ReadEvent('X6','*','00','HHZ',ev,['+1h']);

%
% Bandpass from 10-20sec (0.05 to 0.1 Hz), 2 passes (acausal), 2 poles butterworth
%
s = duglProc.bp(s,.05,.1,2,2); 
duglPlot.p2(s);

%%
%
% TAUP
%
% The duglProc routine links to a TauP java library, which estimates
% arrival times and paths of seismic phases
%

figure(5)
taup = duglProc.TauP(s);
duglPlot.PlotTauP(taup);




% The Taup data (incident angle/dip and azimuth) to estimate delays on
% traces.
[lprop,tdif] = duglProc.EstimateDelay(s,taup);
figure(6)
hh = duglPlot.prs(lprop(:,1),taup(1),s,'scaling_factor',0.2);



figure(7)
hh = duglPlot.p2(s)
duglPlot.TimeTauP(taup)
% duglPlot.idTrace can be used after a line (data or arrival) has been clicked
% with Matlab's dataCursor


figure(8)
hh = duglPlot.p2(s,'timeshift',tdif,taup); 
xlim([4 5])
ylim([-5000 5000])


