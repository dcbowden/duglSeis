function s = iris2sac( trace )
%
% Personal wrapper to make readsac and irisfetch consistent.
%
% mytrace = irisFetch.Traces('CI','USC','*','BHZ','2012-09-07 07:03:09','2012-09-07 07:04:09','includePZ');
%
% s = readsac('usc.bhz.sac')
%


%
%                 network: 'CI'
%                 station: 'USC'
%                location: ''
%                 channel: 'BHZ'
%                 quality: 'M'
%                latitude: 34.0192
%               longitude: -118.2863
%               elevation: 58
%                   depth: 0
%                 azimuth: 0
%                     dip: -90
%             sensitivity: 627368000
%    sensitivityFrequency: 1
%              instrument: 'STS-2'
%        sensitivityUnits: 'M/S'
%                    data: [2400x1 double]
%             sampleCount: 2400
%              sampleRate: 40
%               startTime: 7.3512e+05
%                 endTime: 7.3512e+05
%                   sacpz: [1x1 struct]
%

s = []
for ii = 1:length(trace)
    s(ii).NPTS   = trace(ii).sampleCount;
    s(ii).B      = 0;
    s(ii).E      = (trace(ii).endTime - trace(ii).startTime)*24*3600; %seconds
    s(ii).IFTYPE = 'TIME SERIES FILE';
    s(ii).LEVEN  = 'TRUE';
    s(ii).DELTA  = 1/trace(ii).sampleRate;
    s(ii).IDEP   = trace(ii).sensitivityUnits;
    s(ii).KZDATE = datestr(trace(ii).startTime,'mmm dd, yyyy');
    s(ii).KZTIME = datestr(trace(ii).startTime,'HH:MM:SS.FFF');
    s(ii).KINST  = trace(ii).channel;
    s(ii).KSTNM  = trace(ii).station;
    s(ii).KNETWK = trace(ii).network;
    s(ii).STLA   = trace(ii).latitude;
    s(ii).STLO   = trace(ii).longitude;
    s(ii).STEL   = trace(ii).elevation;
    s(ii).DATA1  = trace(ii).data;
    s(ii).FILENAME = [datestr(trace(ii).startTime,'yyyymmddHHMMSSFFF'),...
                      '.',trace(ii).network,'.',trace(ii).station,'.',...
                      trace(ii).channel,'.sac'];
    %s(ii).       = trace(ii).
    % keep in this format for now
    s(ii).location = trace(ii).location;
    s(ii).quality = trace(ii).quality;
    s(ii).dip = trace(ii).dip;
    s(ii).instrument = trace(ii).instrument;
    s(ii).sacpz = trace(ii).sacpz;
    s(ii).sensitivityUnits = trace(ii).sensitivityUnits;
end
