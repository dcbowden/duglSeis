classdef duglFetch

% duglFetch
%
% s = duglFetch.Read('X6','1700','00','HHZ','2015 108 16:30:00','2015 109 06:29:00');
% s = duglFetch.Read('X6','1700','00','HHZ','2015-04-20 16:30:00','+1d +1h +4m +30s');
%
% 
% ev = irisFetch.Events('MinimumMagnitude',7.7,'startTime','2015-05-30 11:13:02','endTime','2015-05-30 11:33:02'); %japan
% s = duglFetch.ReadEvent('X6','300|D4850','00','HHZ',ev,'+1h');
%
% May need to run "duglFetch.Init" once to get java libraries set up on path



% Notes to self:
%   Will round msec in input to seconds
%   Uses matlab mapping toolbox

   properties (Constant = true)
	functions = {'Init','Read','ReadEvent',};
%        DATE_FORMATTER    = 'yyyy-mm-dd HH:MM:SS.FFF'; %default data format, in ms
   end %constant properties

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   methods(Static)
      %%%%%%%
      function Init
	javaaddpath([duglSet.IRIS_DIR,'/IRIS-WS-2.0.12.jar']);
      end %init

      %%%%%%%%
      function s = Read(network, stations_in, location, channel_in, startDate, endDate, varargin )
	% test = duglFetch.Read('X6','1700','00','HHZ','2015 108 16:30:00','2015 109 06:29:00');
	% test = duglFetch.Read('X6','1700','00','HHZ','2015-04-20 16:30:00','+1d +1h +4m +30s');

	getsacpz    = false;
	verbosity   = false;
	safeLocation = @(x) strrep(x,' ','-');
	% This is maybe more complete, but not sure if it's necessary. Slower
	%fixDelta = @(x) double(round(vpa(x)*1000)/1000);
	fixDelta = @(x) round(x*1000)/1000;

	args_in = [];
	for jj = 1:length(varargin)
		if(length(varargin{jj})>0)
			args_in{jj} = varargin{jj};
		end
	end
	flags = duglFetch.parse_args(args_in);
	
	if(length(strfind(endDate,'+'))>0)
		startInfo = duglFetch.makeDateInfo(startDate);
		endInfo   = duglFetch.plusDateInfo(startInfo,endDate);
	else
		startInfo    = duglFetch.makeDateInfo(startDate);
		endInfo      = duglFetch.makeDateInfo(endDate);
	end
	%startInfo
	%endInfo
	%location     = safeLocation(location);
	%if(startInfo.yyyy ~= endInfo.yyyy)
	%	disp('ERROR straddling year, not programmed for this')
	%	s = [];
	%	return
	%end
	[stations_pull,chan_pull] = duglFetch.parseStations(stations_in,channel_in);

	% HANDLE WILDCARDS IN STATION NAME?
	% HANDLE WILDCARDS IN CHANNEL?


	used_ista = 0;
	for ista = 1:length(stations_pull)
	    for ichan = 1:length(chan_pull)
		station = stations_pull{ista};
		chan = chan_pull{ichan};
		
		clear stmp;

		%% Read in Starting file
		startFileDir = [duglSet.DATA_DIR,'/',num2str(startInfo.yyyy),'/',num2str(startInfo.DOY,'%03u'),'/'];
		
		%% Old format, from antelope output
		%startFile = [startFileDir,'/',num2str(startInfo.yyyy),num2str(startInfo.DOY,'%03u'),'*',...
		%		'.',station,'.',chan];

		%% New naming convention:
		%% 2015/174/X6.A4100..HHZ.M.2015.174.000000.SAC
		%% 2015/174/X6.A4100..HHZ.Q.2015.174.000000000.SAC
		startFile = [startFileDir,'/X6.',station,'..',chan,'.Q.',...
		 		num2str(startInfo.yyyy),'.',num2str(startInfo.DOY,'%03u'),'*.SAC'];
		ssFiles = dir(startFile);
		if(length(ssFiles)>1)
			if(flags.verbose)
				disp(['Warning! More than one file matching: ',startFile]);
			end
		elseif(length(ssFiles)==0)
			continue;
		end

		for ii = 1:length(ssFiles)
			stmp(ii) = readsac([startFileDir,'/',ssFiles(ii).name]);
			stmp(ii).DELTA = fixDelta(stmp(ii).DELTA);
			stmp(ii).E = stmp(ii).B + stmp(ii).NPTS*stmp(ii).DELTA - stmp(ii).DELTA;;
		end

		%% Read in more files to span the time, if needed
		nstmp = length(stmp);
		if(endInfo.DOY ~= startInfo.DOY)
			%numExtraDays = endInfo.DOY-startInfo.DOY;
			numExtraDays = (endInfo.DOY+endInfo.yyyy*365) - (startInfo.DOY+startInfo.yyyy*365);
			for jj = 1:numExtraDays

				% special case: exactly 24 hours requested? Don't bother reading in next day
				if(jj==numExtraDays && endInfo.HH == 0 && endInfo.MM == 0 && endInfo.SS == 0)
					continue;
				end

				nextFileDir = [duglSet.DATA_DIR,'/',num2str(startInfo.yyyy),'/',num2str(startInfo.DOY+jj,'%03u'),'/'];
				%nextFile = [nextFileDir,'/',num2str(startInfo.yyyy),num2str(startInfo.DOY+jj,'%03u'),'*',...
				%	'.',station,'.',channel];
				nextFile = [nextFileDir,'X6.',station,'..',chan,'.Q.',...
						num2str(startInfo.yyyy),'.',num2str(startInfo.DOY+jj,'%03u'),'*.SAC'];


				% special case: roll into next year!
				if(startInfo.DOY+jj > 365)
					nextFileDir = [duglSet.DATA_DIR,'/',num2str(startInfo.yyyy+1),'/',num2str(startInfo.DOY+jj-365,'%03u'),'/'];
					nextFile = [nextFileDir,'X6.',station,'..',chan,'.Q.',...
						num2str(startInfo.yyyy+1),'.',num2str(startInfo.DOY+jj-365,'%03u'),'*.SAC']
				end

				ssFiles = dir(nextFile);
				if(length(ssFiles)>0)
					for ii = 1:length(ssFiles)
						stmp(nstmp+ii) = readsac([nextFileDir,'/',ssFiles(ii).name]);
						stmp(nstmp+ii).DELTA = fixDelta(stmp(nstmp+ii).DELTA);
						stmp(nstmp+ii).E = stmp(nstmp+ii).B + ...
						    stmp(nstmp+ii).NPTS*stmp(nstmp+ii).DELTA - stmp(nstmp+ii).DELTA;
						%stmp(jj+1) = readsac([nextFileDir,'/',sFile.name]);
					end
				end
				nstmp = length(stmp);
			end
%disp('PRE-COMBINE')
%disp(sprintf('%.8f %.8f',stmp.B))
%disp(sprintf('%.8f %.8f',stmp.E))
%disp(sprintf('%.8f %.8f',stmp.NPTS))
%disp(sprintf('%.12f %.12f',stmp.DELTA))
%disp(sprintf('%.8f %.8f',length(stmp(1).DATA1),length(stmp(1).DATA1)))
			return_separate_pieces = 0;
			if(return_separate_pieces)
				s=stmp;
				return
			end
			stmp = duglFetch.combinesac(stmp,flags);
		end
%disp('PRE-CUT')
%disp(sprintf('%.8f',stmp.B))
%disp(sprintf('%.8f',stmp.E))
%disp(sprintf('%.8f',stmp.NPTS))
%disp(sprintf('%.8f',length(stmp.DATA1)))


		if(length(stmp)>0)
			stmp = duglFetch.cutTime(stmp,startInfo,endInfo,flags);
%disp('POST-CUT')
%disp(sprintf('%.8f',stmp.B))
%disp(sprintf('%.8f',stmp.E))
%disp(sprintf('%.8f',stmp.NPTS))
%disp(sprintf('%.8f',length(stmp.DATA1)))
			if(length(stmp)>0)
				used_ista = used_ista+1;
				s(used_ista) = stmp;
			else
				%disp('Problem in cutTime! No data returned')
				%s(used_ista) = [];
			end
		else
			disp('No data found')
			s(used_ista) = [];
		end
	    end
	end
	if(exist('s')==0)
		disp('ERROR! No match found!')
		disp(['Requested: ',startInfo.datstr,' to ',endInfo.datstr])
		s = [];
	end
	%s = duglFetch.SetPositions(s);
      end % Fetch

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function s = ReadEvent(network, stations_in, location, channel, ev, endDate, varargin )
	% %Pull up more details about the event from IRIS server             
	% ev = irisFetch.Events('MinimumMagnitude',7.7,'startTime','2015-05-30 11:13:02','endTime','2015-05-30 11:33:02'); %japan
	% % Read in data from our array                                       
	% s = duglFetch.ReadEvent('X6','*','00','HHZ',ev,['+1h']);           

	args_in = [];
	for jj = 1:length(varargin)
		args_in{jj} = varargin{jj};
	end

	% Read in the data as normal
	s = duglFetch.Read(network,stations_in,location,channel,ev.PreferredTime,endDate,args_in);
	
	% Parse some additional SAC info:
	imagtyp = 57;
	switch ev.PreferredMagnitudeType
	  case 'MWW' %IMW (Moment Magnitude)
		imagtyp = 55;
	  otherwise 
	  	imagtyp = 57;
		display(['Not coded to interpret origin type:', ev.PreferredMagnitudeType]);
			%imb	52
			%ims	53
			%iml	54
			%imw	55
			%imd	56
			%imx	57
	end

	imagsrc = 71;
	switch ev.PreferredOrigin.Contributor
	  case 'NEIC COMCAT'
		imagsrc = 64; %iusgs
	  otherwise
		imagsrc = 71;
		display(['Not coded to interpret contributor:', ev.PreferredOrigin.Contributor]);
			%imb	52
			%ims	53
			%iml	54
			%imw	55
			%imd	56
			%imx	57
			%ineic	58
			%ipdeq	59
			%ipdew	60
			%ipde	61
			%iisc	62
			%ireb	63
			%iusgs	64
			%ibrk	65
			%icaltech	66
			%illnl	67
			%ievloc	68
			%ijsop	69
			%iuser	70
			%iunknown	71
	end

	ievtyp = 71;
	switch ev.Type
	  case 'earthquake' %iquake / ieq
		ievtyp = 40;
	  otherwise
		ievtyp = 71;
		display(['Not coded to interpret event type:', ev.Type]);
			%iunknown	71
			%inucl	37
			%ipren	38	nuclear pre-shot
			%ipostn	39	nuclear post-shot
			%iquake	40	earthquake
			%ipreq	41	foreshock
			%ipostq	42	aftershock
			%ichem	43	chemical explosion
			%iqb	72	quarry confirmed
			%iqb1	73	quarry designed shot
			%iqb2	74	quarry observed shot
			%iqbx	75	quarry single shot
			%iqmt	76	quarry-induced event
			%ieq	77	earthquake
			%ieq1	78	eq in swarm or aftershock seq
			%ieq2	79	felt earthquake
			%ime	80	marine explosion
			%iex	81	other explosion
			%inu	82	nuclear explosion
			%inc	83	nuclear cavity collapse
			%io_	84	other source
			%il	85	local unkown event
			%ir	86	regional unkown event
			%it	87	teleseism uknown event
			%iu	88	undetermined / conflicting
			%ieq3	89
			%ieq0	90
			%iex0	91
			%iqc	92
			%iqb0	93
			%igey	94
			%ilit	95
			%imet	96
			%iodor	97
	end
		

	if( isfield(ev,'Name') );
		kevnm = ev.Name;
	else
		kevnm = ev.PreferredOrigin.ContributorEventId;
	end
	%s = duglFetch.SetPositions(s); called already in "Read"

	% Add event info:
	for jj = 1:length(s)
		PubID = ev.PreferredOrigin.PublicId;
		idi = strfind(PubID,'originid=');
		if(length(idi)>0)
			id = PubID(idi+9:end);
		else
			id = '9999999';
		end
		s(jj).KEVNM   = id			              ;
		%s(jj).IEVREG  =                                          ;
		s(jj).EVLA    = ev.PreferredLatitude                      ;
		s(jj).EVLO    = ev.PreferredLongitude                     ;
		%s(jj).EVEL    = ev.                                      ;
		s(jj).EVDP    = ev.PreferredDepth                         ;
		s(jj).MAG     = ev.PreferredMagnitudeValue                ;
		s(jj).IMAGTYP = imagtyp                                   ;
		s(jj).IMAGSRC = imagsrc                                   ;
		s(jj).IEVTYP  = ievtyp                                    ;
		%s(jj).NEVID   = ev.PreferredOrigin.ContributorEventId     ;
		%s(jj).NORID   = ev.PreferredOrigin.ContributorOriginId    ;
		%s(jj).NWFID   =                                          ;
		%s(jj).KHOLE   =                                          ;

		% note! If the matlab mapping toolbox does not exist, can also
		%  use distance/azimuth from TauP (see duglProc for example)


		%% DCB 2/24/2017, fixing azimuth from station->source to source-> station
		%[arclen,az]  = distance(s(jj).STLA,s(jj).STLO,s(jj).EVLA,s(jj).EVLO);
		[arclen,az]  = distance(s(jj).EVLA,s(jj).EVLO,s(jj).STLA,s(jj).STLO);
		while (az<0)
			az = az+360;
		end
		baz = az+180;
		while(baz>360)
			baz = baz-360;
		end
		%s(jj).DIST    = deg2km(arclen)                               ;
		s(jj).DIST    = duglProc.DistVicenty(s(jj).STLA,s(jj).STLO,s(jj).EVLA,s(jj).EVLO)./1000;
		s(jj).AZ      = az                                         ;
		s(jj).BAZ     = baz                                         ;
		s(jj).GCARC   = arclen                                         ;
	end
	
      end % function ReadEvent

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%%%%%%
      function myDateInfo = makeDateInfo(dateInput)
	% Returns and object with fields for various date information
	% d1 = duglFetch.makeDateInfo('2015 001 00:00:00')
	% d1 = duglFetch.makeDateInfo('2015-01-01 00:00:00')

    	 myDateInfo.yyyy   = [];
    	 myDateInfo.DOY    = [];
    	 myDateInfo.mm     = [];
    	 myDateInfo.dd     = [];
    	 myDateInfo.HH     = [];
    	 myDateInfo.MM     = [];
    	 myDateInfo.SS     = [];
    	 myDateInfo.MS     = [];
    	 myDateInfo.datnum = [];
    	 myDateInfo.datstr = [];

	 if(ischar(dateInput)==0)
	 % Must be a datenum
		myDateInfo.datnum = dateInput;
		myDateInfo.datstr = datestr(dateInput,'yyyy-mm-dd HH:MM:SS.FFF');
		myVec = datevec(myDateInfo.datnum);
		myDateInfo.yyyy = myVec(1);	
		myDateInfo.mm = myVec(2);	
		myDateInfo.dd = myVec(3);	
		myDateInfo.HH = myVec(4);	
		myDateInfo.MM = myVec(5);	
		myDateInfo.SS = round(myVec(6));	
		%myDateInfo.SS = (myVec(6));	
		myDateInfo.MS = round(mod(myVec(6),1).*1000);
		% Calculate DOY
		myVec(:,4:end) = 0;
		datnum_date = datenum(myVec);
		myVec(:,2:end) = 0;
		datnum_year = datenum(myVec);
		doy = datnum_date - datnum_year;
		myDateInfo.DOY= doy;
	 elseif(dateInput(5) == ' ')
   	 %assume DOY input
		myDateInfo.yyyy = str2num(dateInput(1:4));
		myDateInfo.DOY = str2num(dateInput(6:8));
		myDateInfo.HH = str2num(dateInput(10:11));
		myDateInfo.MM = str2num(dateInput(13:14));
		myDateInfo.SS = str2num(dateInput(16:17));
		if(length(dateInput)>17)
			if(length(dateInput)==19)
				myDateInfo.MS = round(str2num(dateInput(19))).*100;
			elseif(length(dateInput)==20)
				myDateInfo.MS = round(str2num(dateInput(19:20))).*10;
			elseif(length(dateInput)==21)
				myDateInfo.MS = round(str2num(dateInput(19:21)));
			else
				disp(sprintf('Did not understand date: %s',dateInput))
				disp('yyyy-mm-dd HH:MM:SS.FFF or yyyy doy HH:MM:SS.FFF')
			end
		else
			myDateInfo.MS = 0;
		end
		%	myDateInfo.seconds = myDateInfo.HH*60*60 + ...
		%				myDateInfo.MM*60 + ...
		%				myDateInfo.SS + ...
		%				myDateInfo.MS./1000;
		% Formatting for nice string
		other = dateInput(9:end);
		d = horzcat(myDateInfo.yyyy,zeros(1,5)); 
		dateV = myDateInfo.DOY + datenum(d);
		beginning = datestr(datenum(dateV),'yyyy-mm-dd');
		myDateInfo.mm = str2num(beginning(6:7));
		myDateInfo.dd   = str2num(beginning(9:10));
		dateInput2 = [beginning,' ',other];
		myDateInfo.datnum = datenum(dateInput2);
		myDateInfo.datstr = datestr(myDateInfo.datnum,duglSet.DATE_FORMATTER);
	  else
		myDateInfo.datnum = datenum(dateInput);
		myDateInfo.datstr = datestr(myDateInfo.datnum,duglSet.DATE_FORMATTER);
		myVec = datevec(myDateInfo.datnum);
		myDateInfo.yyyy = myVec(1);	
		myDateInfo.mm = myVec(2);	
		myDateInfo.dd = myVec(3);	
		myDateInfo.HH = myVec(4);	
		myDateInfo.MM = myVec(5);	
		myDateInfo.SS = round(myVec(6));	
		%myDateInfo.SS = (myVec(6));	
		myDateInfo.MS = round(mod(myVec(6),1).*1000);
		% Calculate DOY
		myVec(:,4:end) = 0;
		datnum_date = datenum(myVec);
		myVec(:,2:end) = 0;
		datnum_year = datenum(myVec);
		doy = datnum_date - datnum_year;
		myDateInfo.DOY= doy;
	end
	month_names = [{'Jan'},{'Feb'},{'Mar'},{'Apr'},{'May'},{'Jun'},{'Jul'},{'Aug'},...
 	       	{'Sep'},{'Oct'},{'Nov'},{'Dec'}];
	myDateInfo.month = month_names{myDateInfo.mm};
      end
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%%%%%%
      function myDateInfo = plusDateInfo(dateInfo,modif)
	 % myDateInfo = duglFetch.plusDatinfo(dateInfo,'+1d +1h +4m +30s')
	 %
	 % Not programmed to handle milliseconds

    	 myDateInfo.yyyy   = [];
    	 myDateInfo.DOY    = [];
    	 myDateInfo.mm     = [];
    	 myDateInfo.dd     = [];
    	 myDateInfo.HH     = [];
    	 myDateInfo.MM     = [];
    	 myDateInfo.SS     = [];
    	 myDateInfo.MS     = [];
    	 myDateInfo.datnum = [];
    	 myDateInfo.datstr = [];

	 ds = '\+(?<day>\S+)d';
	 hs = '\+(?<hour>\S+)h';
	 ms = '\+(?<min>\S+)m';
	 ss = '\+(?<sec>\S+)s';
	 modifs = regexp(modif,[ds '|' hs '|' ms '|' ss],'names');

	 days  = horzcat(modifs.day); 
	 if(length(days)>0) 
		days=str2num(days); 
	 else 
	 	days=0;
	 end
	 hours  = horzcat(modifs.hour); 
	 if(length(hours)>0) 
		hours=str2num(hours); 
	 else 
	 	hours=0;
	 end
	 mins  = horzcat(modifs.min); 
	 if(length(mins)>0) 
		mins=str2num(mins); 
	 else 
	 	mins=0;
	 end
	 secs  = horzcat(modifs.sec); 
	 if(length(secs)>0) 
		secs=str2num(secs); 
	 else 
	 	secs=0;
	 end

	% Check for overflow
	while(secs >= 60)
		mins = mins+1;
		secs = secs-60;
	end
	while(mins >= 60)
		hours = hours+1;
		mins = mins-60;
	end
	while(hours >= 24)
		days = days+1;
		hours = hours-24;
	end


	 myDateInfo.yyyy = dateInfo.yyyy;
	 myDateInfo.DOY = dateInfo.DOY+days;
	 myDateInfo.HH  = dateInfo.HH+hours;
	 myDateInfo.MM  = dateInfo.MM+mins;
	 myDateInfo.SS  = dateInfo.SS+secs;
	 myDateInfo.MS = round(dateInfo.MS);

	if(mod(myDateInfo.SS,1)~=0)
		myDateInfo.MS = round(myDateInfo.MS + mod(myDateInfo.SS,1).*1000);
		myDateInfo.SS = round(myDateInfo.SS);
	end
	 while(myDateInfo.SS >= 60)
		myDateInfo.MM = myDateInfo.MM+1;
		myDateInfo.SS = myDateInfo.SS-60;
	end
	 while(myDateInfo.MM >= 60)
		myDateInfo.HH = myDateInfo.HH+1;
		myDateInfo.MM = myDateInfo.MM-60;
	end
	 while(myDateInfo.HH >= 24)
		myDateInfo.DOY = myDateInfo.DOY+1;
		myDateInfo.HH = myDateInfo.HH-24;
		% NOTE! Rolling over a full day leads to inconsistency
		% 1 day in datenum format is NOT equal to rolling up 24hr...
		% It seems Matlab accounts for leap seconds? 
	end
	 while(myDateInfo.DOY >= 366)
		myDateInfo.DOY = myDateInfo.DOY-365;
		myDateInfo.yyyy = myDateInfo.yyyy+1;
	 end

	 %myDateInfo.seconds = myDateInfo.HH*60*60 + ...
	 %       		myDateInfo.MM*60 + ...
	 %       		myDateInfo.SS;
	 temp_datnum = dateInfo.datnum+days+hours/24+mins/24/60+secs/24/60/60;
	 myVec = datevec(temp_datnum);
	 myDateInfo.mm = myVec(2);	
	 myDateInfo.dd = myVec(3);	
	 %myDateInfo.datnum = temp_datnum;
	% Recalculate b/c of rounding error / leap-second. (not sure which is the problem)
	% disp(sprintf('%4u-%02u-%02u %02u:%02u:%02u.%03u',...
	%			myDateInfo.yyyy,myDateInfo.mm,myDateInfo.dd,...
	%			myDateInfo.HH,myDateInfo.MM,myDateInfo.SS,myDateInfo.MS));

	% disp(sprintf('%4u-%02u-%02u %02u:%02u:%02u.%03u',...
	%			myDateInfo.yyyy,myDateInfo.mm,myDateInfo.dd,...
	%			myDateInfo.HH,myDateInfo.MM,myDateInfo.SS,myDateInfo.MS))
	 myDateInfo.datnum = datenum(sprintf('%4u-%02u-%02u %02u:%02u:%02u.%03u',...
				myDateInfo.yyyy,myDateInfo.mm,myDateInfo.dd,...
				myDateInfo.HH,myDateInfo.MM,myDateInfo.SS,myDateInfo.MS));

 	 myDateInfo.datstr = datestr(myDateInfo.datnum,duglSet.DATE_FORMATTER);

	 month_names = [{'Jan'},{'Feb'},{'Mar'},{'Apr'},{'May'},{'Jun'},{'Jul'},{'Aug'},...
 			{'Sep'},{'Oct'},{'Nov'},{'Dec'}];
	 myDateInfo.month = month_names{myDateInfo.mm};

      end % plus date info
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function sout = cutTime(s,startInfo,endInfo,varargin)

	if(length(varargin)>0)
		flags = varargin{1};
	else
		flags.lessOK = 0;
		flags.verbose= 1;
	end

	%%define nsec_begin
	nsec_begin =    (startInfo.DOY-s(1).NZJDAY)*24*60*60 + ...
			(startInfo.HH -s(1).NZHOUR)*60*60 + ...
			(startInfo.MM -s(1).NZMIN )*60 + ...
			(startInfo.SS -s(1).NZSEC ) - ...
			(s(1).NZMSEC / 1000);
	%%define nsec_end
	nsec_end = (endInfo.DOY-s(1).NZJDAY)*24*60*60 + ...
	           (endInfo.HH -s(1).NZHOUR)*60*60 + ...
	           (endInfo.MM -s(1).NZMIN )*60 + ...
	           (endInfo.SS -s(1).NZSEC ) - ...
		   (s(1).NZMSEC / 1000);

	%% Special case, straddling year
	if(startInfo.yyyy~=endInfo.yyyy)
		nsec_end = (endInfo.yyyy-s(1).NZYEAR)*24*60*60*365 + ...
		   (endInfo.DOY-s(1).NZJDAY)*24*60*60 + ...
	           (endInfo.HH -s(1).NZHOUR)*60*60 + ...
	           (endInfo.MM -s(1).NZMIN )*60 + ...
	           (endInfo.SS -s(1).NZSEC ) - ...
		   (s(1).NZMSEC / 1000);
	end
		
	%sout = s(1);
	%if(nsec_begin==0)
	%	nsec_begin = 1;
	%end


	% S contains some, probably longer, set of data: midnight to midnight
	% startInfo - endInfo requests some subset: 3am - 4am
	% nsec_begin - nsec_end is index of data to cut.

	%  |************************************|     Data
	%         |........|                          Request
	%   ------> nsec_begin
	%   ---------------> nsec_end


	
	if(nsec_end<0)
		% Requested time period ends before data available.
		if(flags.verbose)
			disp('ERROR cutTime! not in time range; request time period ends before data available')
		end
		sout = [];
		return
		%                     |*************|     Data
		%         |        |
		% nsec_be <------------
		%         nsec_end <---
	
	elseif(nsec_begin<0)
		% Default is to be strict, and not return waveform unless complete.
		% i.e.: if data requested (startInfo) is before data avaialble (s)

		%              |****************|     Data
		%         |     ...|
		% nsec_be <-----
		%              ----> nsec_end
		if(flags.lessOK==0)
			if(flags.verbose)
				disp('ERROR cutTime! not in time range; available data starts later than request')
				disp('Enable "lessOK" flag in request to allow for shorter traces to return')
				disp(s(1).FILENAME)
				disp(sprintf('Available:     %4u %3u %02u:%02u:%02u.%u + %fs',s(1).NZYEAR,s(1).NZJDAY,s(1).NZHOUR,s(1).NZMIN,round(s(1).NZSEC),s(1).NZMSEC,s(1).E))
				disp(sprintf('Request start: %4u %3u %02u:%02u:%02u.%u',startInfo.yyyy,startInfo.DOY,startInfo.HH,startInfo.MM,floor(startInfo.SS),round(mod(startInfo.SS,1)*100)))
				disp(sprintf('Request end:   %4u %3u %02u:%02u:%02u.%u',endInfo.yyyy,endInfo.DOY,endInfo.HH,endInfo.MM,floor(endInfo.SS),round(mod(endInfo.SS,1)*100)))
				disp(' ')
			end
			sout = [];
			return 
		else
			startInfo= duglFetch.plusDateInfo(startInfo, sprintf('+%fs',-nsec_begin));
			nsec_begin=0;
		     	if(flags.verbose)
				% update nsec_begin (basically ignore any front-end cutting)
				%startInfo
				%disp(sprintf('+%fs',-nsec_begin))
				disp('WARNING cutTime! Returning incomplete trace')
				disp(s(1).FILENAME)
				disp(' ')
		     	end
		end
	end

	% TEMP fix, 2 samples within error of sac output
	if(nsec_end/s(1).DELTA - length(s(1).DATA1) == 1)
		s(1).DATA1(end+1) = s(1).DATA1(end);
		if(flags.verbose)
			disp('Adding 1 sample!')
			disp(s(1).FILENAME)
		end
	end
	if(nsec_end/s(1).DELTA - length(s(1).DATA1) == 2)
		s(1).DATA1(end+2) = s(1).DATA1(end);
		if(flags.verbose)
			disp('Adding 2 samples!')
			disp(s(1).FILENAME)
		end
	end

	if(nsec_begin/s(1).DELTA > length(s(1).DATA1))
	      	if(flags.verbose)
			disp('ERROR cutTime! not in time range; request time period begins after data available ends')
			%nsec_begin/s(1).DELTA
			%nsec_end/s(1).DELTA
			%length(s(1).DATA1)
			%s.FILENAME

			%    |****************|                Data
			%                        |        |    Request
			
			disp(sprintf('length s(1).DATA1: %.8f', length(s(1).DATA1)))
			disp(sprintf('       NPTS      : %.8f', s(1).NPTS) )
			disp(sprintf('       s(1).E    : %.8f', s(1).E))
			disp(sprintf('length cutrequest: %.8f', round(nsec_end/s(1).DELTA) - round((nsec_begin)/s(1).DELTA)))
			disp(sprintf('nsec_end/delta   : %.8f', (nsec_end/s(1).DELTA)))
			disp(s(1).FILENAME)
			disp(' ')
		end
		sout = [];
		return 
	elseif(nsec_end/s(1).DELTA > length(s(1).DATA1))

		%    |****************|                Data
		%                |.....   |            Request

		if(flags.lessOK==0)
			if(flags.verbose)
				disp('ERROR cutTime! not in time range; data request goes beyond available data')
				disp('Enable "lessOK" flag in request to allow for shorter traces to return')
				disp(s(1).FILENAME)
				disp(sprintf('Available:     %4u %3u %02u:%02u:%02u.%u + %fs',s(1).NZYEAR,s(1).NZJDAY,s(1).NZHOUR,s(1).NZMIN,round(s(1).NZSEC),s(1).NZMSEC,s(1).E))
				disp(sprintf('Request start: %4u %3u %02u:%02u:%02u.%u',startInfo.yyyy,startInfo.DOY,startInfo.HH,startInfo.MM,floor(startInfo.SS),round(mod(startInfo.SS,1)*100)))
				disp(sprintf('Request end:   %4u %3u %02u:%02u:%02u.%u',endInfo.yyyy,endInfo.DOY,endInfo.HH,endInfo.MM,floor(endInfo.SS),round(mod(endInfo.SS,1)*100)))
				disp(' ')
		        end
			sout = [];
			return
		else
			nsec_end = floor(length(s(1).DATA1)*s(1).DELTA);
		      	if(flags.verbose)
				% update to go as far as we can
				disp('WARNING cutTime! Returning incomplete trace')
				disp(s(1).FILENAME)
				disp(' ')
		      	end
		end
	end

	%disp(sprintf('length s(1).DATA1: %.8f', length(s(1).DATA1)))
	%disp(sprintf('length cutrequest: %.8f', length( round((nsec_begin)/s(1).DELTA)+1 : round(nsec_end/s(1).DELTA))))
	%disp(sprintf('length cutrequest: %.8f', round(nsec_end/s(1).DELTA) - round((nsec_begin)/s(1).DELTA) ))
	%if( s(1).NPTS ~= length(round((nsec_begin)/s(1).DELTA)+1 : round(nsec_end/s(1).DELTA)))
	%round((nsec_begin)/s(1).DELTA)+1
	%round(nsec_end/s(1).DELTA)
	%length(s(1).DATA1)
	if( s(1).NPTS ~= round(nsec_end/s(1).DELTA) - round((nsec_begin)/s(1).DELTA) )
		s(1).DATA1 = s(1).DATA1( round((nsec_begin)/s(1).DELTA)+1 : round(nsec_end/s(1).DELTA) );
		%sout.DATA1 = s(1).DATA1( round((nsec_begin)/s(1).DELTA) : round(nsec_end/s(1).DELTA)-1 );
		%s(1).DATA1 = s(1).DATA1 - mean(s(1).DATA1);
	end
	sout = s(1);
	sout.NPTS = length(sout(1).DATA1);
	sout.E = nsec_end - nsec_begin; % number of seconds in the trace
	sout.E = sout.E-1*s(1).DELTA;   % because we start index at zero
	%sout.B = nsec_begin;
	%sout.E = sout.E - sout.B;
	sout.B = 0;			% because SAC file name, NZHOUR, etc. will contain start time

	% X6.D4850..HHE.Q.2015.001.134123680.SAC
	%    sta    cha   yyyy doy HHMMSSFFF
	sout.FILENAME = sprintf('X6.%s..%s.Q.%04u.%03u.%02u%02u%02u%03u.SAC',s(1).KSTNM,s(1).KCMPNM,...
			  startInfo.yyyy,startInfo.DOY,startInfo.HH,startInfo.MM,startInfo.SS,startInfo.MS);
	%	[num2str(startInfo.yyyy),num2str(startInfo.DOY,'%03u'),...
	%		num2str(startInfo.HH,'%02u'),num2str(startInfo.MM,'%02u'),...
	%		num2str(floor(startInfo.SS),'%02u'),'.',num2str(round(mod(startInfo.SS,1)*100),'%02u'),...
	%		'.',s(1).KSTNM,'.',s(1).KCMPNM];
	sout.NZJDAY = startInfo.DOY;
	sout.NZHOUR = startInfo.HH;
	sout.NZMIN = startInfo.MM;
	sout.NZSEC = floor(startInfo.SS);
	sout.NZMSEC = round(mod(startInfo.SS,1)*1000);
	
      end % cutTime
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
    end % methods static



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   methods(Static, Access=protected)      

      function flags = parse_args(args_in)
	nV = length(args_in);
	
	% List of flags
	flags.lessOK = 0;
	flags.verbose = 1;

	jj = 1;
	if( nV > 0)
	    while jj <= nV
		if(strcmp(args_in{jj},'lessOK'))
			flags.lessOK = 1;
		elseif(strcmp(args_in{jj},'quiet'))
			flags.verbose = 0;
		else
			disp(sprintf('Extra flag not understood: %s',args_in{jj}))
		end
		jj = jj+1;
	    end
	end


      end % 
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function [stations,components] = parseStations(stations_in,channel)
	ALL_stations = duglSet.Stations;
	ALL_components = duglSet.Components;
%	ALT_stations = duglSet.Alt_stations;
	stations = [];
	components = [];
	
	% Break down list of stations vs wildcard
	ii = 1;
	split = regexp(stations_in,'\|','split');
	for jj = 1:length(split);
		split{jj} = strrep(split{jj},'*','.');
		for kk = 1:length(ALL_stations)
			if( length(regexp(ALL_stations{kk},split{jj})) )
				stations{ii} = ALL_stations{kk};
				ii = ii+1;
			end
		end

%		for kk = 1:length(ALT_stations)
%			if( length(regexp(ALT_stations{kk},split{jj})) )
%				for kkk = 1:length(ALT_stations)
%					stations{ii} = stations
	end


	% Break down list of components vs wildcard
	ii = 1;
	split = regexp(channel,'\|','split');
	for jj = 1:length(split)
		split{jj} = strrep(split{jj},'*','.');
		for kk = 1:length(ALL_components)
			if( length(regexp(ALL_components{kk},split{jj})) )
				components{ii} = ALL_components{kk};
				ii = ii+1;
			end
		end
	end

	% Combine
	


      end %parseStations
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%%%%%%
      function sout = combinesac(s,varargin)

	if(length(varargin)>0)
		flags = varargin{1};
	else
		flags.lessOK = 0;
		flags.verbose= 1;
	end


	%s(1).E=round(length(s(1).DATA1)*s(1).DELTA)-1*s(1).DELTA; %recalc for precision
	%s(1).E=(length(s(1).DATA1)*s(1).DELTA)-1*s(1).DELTA; %recalc for precision
	
	% Combine multiple sac files into one continuous trace
	for jj = 2:length(s) 
		if(~strcmp(s(jj).KSTNM,s(1).KSTNM) || ~strcmp(s(jj).KCMPNM,s(1).KCMPNM))
			display(['Error! Multiple files/days read in belonging to',...
				 ' different stations/components']);
			return
		end
		s(jj).B = (s(jj).NZYEAR-s(1).NZYEAR)*24*60*60*365 +  (s(jj).NZJDAY-s(1).NZJDAY)*24*60*60 + ...
                                + (s(jj).NZHOUR-s(1).NZHOUR)*60*60 ...
				+ (s(jj).NZMIN-s(1).NZMIN)*60 + (s(jj).NZSEC-s(1).NZSEC) ...
				+ (s(jj).NZMSEC-s(1).NZMSEC)/1000;
		if(s(jj).B < s(jj-1).E)
			if(flags.verbose)
				disp(sprintf('%s',s(jj).FILENAME))
				disp(sprintf('%.8f',s(jj).B))
%figure
%duglPlot.p2(s,'tscale','sec','no_coloring')
				disp(sprintf('%s',s(jj-1).FILENAME))
				disp(sprintf('%.8f',s(jj-1).E))
				disp('WARNING! Sac file claiming it starts before neighbor ends');
			end
			%return
		end
		if( (s(jj).B - s(jj-1).E) > 10 )
			if(flags.verbose)
				disp('Warning! Significant gap in combining sac file');
				disp('   Returning only pieces before the gap');
			end
			sout = s(1);
			return
		end
		%s(jj).E=round(length(s(jj).DATA1)*s(jj).DELTA)-1*s(jj).DELTA; %recalc for precision
		%s(jj).E=(length(s(jj).DATA1)*s(jj).DELTA); %recalc for precision
		s(jj).E = s(jj).E + s(jj).B;


		%% +1 on sample because Matlab indexes from 1.
		%% (s.B is zero, but we can't call data(0))
		s(1).DATA1(round(s(jj).B/s(jj).DELTA)+1:round(s(jj).E/s(jj).DELTA)+1) = s(jj).DATA1;

		%s(1).E = s(jj).E-1/s(1).DELTA;
	end % for each separate file
	s(1).NPTS = length(s(1).DATA1);
	sout = s(1);

      end % combinesac
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function s=SetPositions(s);

        fid = fopen([duglSet.MAP_DIR,'/stations.txt']);
        stations_temp = textscan(fid,'%s %f %f %f %f %f %f');
	fclose(fid);
          % name lon lat depth(km) easting northing elvation(ft)
        stations = [stations_temp{2} stations_temp{3} stations_temp{4}];
        station_names = stations_temp{1};

        for kk = 1:length(s)
           for jj = 1:length(station_names)
                if(strmatch(s(kk).KSTNM,station_names(jj)))
			s(kk).STLA = stations(jj,2);
			s(kk).STLO = stations(jj,1);
			s(kk).STEL = stations(jj,3);
                end % if match
           end % foreach name
        end % foreach s
      end% function SetPositions

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end % methods static protected
end % classdef
