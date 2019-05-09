classdef duglNoise
   % for a given set of traces:
   % spec = preprocess(s);
   %disp('Fix this: spec.depth')
   %disp('Fix this: add xcorr and ifft steps')
   %
   %

   properties (Constant = true)
	functions = {'preprocess','do_smoothing','build_butterworth','run_PSD','rotate9','rotate4','pull9'};
   end %properties 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   methods(Static)
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function spec = preprocess(s,varargin)
	% s = duglFetch.Read('X6','800','00','HHZ','2015-07-28 01:00:00',['+1h']);
	% spec = duglNoise.preprocess(s);
	% spec = duglNoise.preprocess(s,'norm',0,'whit',0,'filt_order',0);
	% spec = duglNoise.preprocess(s,'norm',0,'whit',0,'filt_order',0,'amp_preserve',1);
	% 
	% amp_preserve can be 
	%    0 - nothing
	%    1 - within station
	%    2 - across all_


	%lowF = .1;
	%highF= 10; %python was 12

% Change 1/1/2018 for longer period, surface stations
lowF = 0.01;
highF = 8;

	proc_filt_order = 2;   % order, 0 for none
	proc_norm = 64;  % width, 0 for none
	proc_whit = 20;  % width, 0 for none
	amp_preserve = 0;

	iv = 1;
	while(iv <= length(varargin))
		switch(varargin{iv})
		  case 'lowF'
			lowF = varargin{iv+1};
			iv = iv+1;
		  case 'highF'
			highF = varargin{iv+1};
			iv = iv+1;
		  case 'filt_order'
			proc_filt_order = varargin{iv+1};
			iv = iv+1;
		  case 'norm'
			proc_norm = varargin{iv+1};
			iv = iv+1;
		  case 'whit'
			proc_whit = varargin{iv+1};
			iv = iv+1;
		  case 'amp_preserve'
			amp_preserve = varargin{iv+1};
			iv = iv+1;
		  otherwise
			disp(['Did not interpret: ',varargin{iv}]);
		  end
		iv = iv+1;
	end
	disp(sprintf('lowF: %f, highF: %f, order: %u, norm: %u, whit: %u, amp_pres: %u',lowF,highF,proc_filt_order,proc_norm,proc_whit,amp_preserve))
		
	if(amp_preserve == 1)
		% need to find the corresponding 3-components for each station
		for ii = 1:length(s)
			station_names{ii} = s(ii).KSTNM;
		end
		stas = unique(station_names);
		igood = 0;
		ibad = [];
		for ii = 1:length(stas)
			cmp3 = find(strcmp(stas{ii},station_names));
			if(length(cmp3)~=3)
				disp(sprintf('ERROR! Station: %s did not have 3 components present',stas{ii}))
				ibad = [ibad; ii];
				continue
			end
			igood = igood+1;
			stas_cmp3(igood,:) = cmp3;
			% stas_cmp3 is now a Nx3 matrix with each row associated Z,N,E components.
			%  refers to indices of s() objects
		end
		stas(ibad) = [];
		%stas
		%length(stas)
	end
		


	% Check length
	maxL = max([s.NPTS]);
	badL = find([s.NPTS]~=maxL);
	if(length(badL)>0)
		disp(['Dropping ',num2str(length(badL)),' incomplete traces'])
	end
	s(badL) = [];
	NFFT = 2^nextpow2(s(1).NPTS);

	% Check sample rate
	badSR = find([s.DELTA]~=s(1).DELTA);
	if(length(badSR)>0)
		disp('ERROR! Not yet programmed to handle different sample rates')
		return
	end

	
	for ii = 1:length(s)
		s(ii).DATA1 = s(ii).DATA1 - mean(s(ii).DATA1);
	end

	% taper
    	lwind=2000;
    	lwind=1000;
	taper = linspace(0,1,lwind)';
	for ii = 1:length(s)
		s(ii).DATA1(1:lwind) = s(ii).DATA1(1:lwind).*taper;
		s(ii).DATA1(end-lwind+1:end) = s(ii).DATA1(end-lwind+1:end).*flip(taper);
	end


	% filter
	if(proc_filt_order ~= 0);
		p = proc_filt_order;
		fnyq = (1/s(1).DELTA)/2;
		[B,A] = butter(p,[lowF highF]./fnyq);
		[H,W] = freqz(B,A,NFFT,'whole');
		for ii = 1:length(s)
			s(ii).DATA1 = filtfilt(B,A,s(ii).DATA1);
		end
	end

	% Time Domain Norm
	if(proc_norm ~= 0)
		%% setting width roughly equal to 1 sec
		%s_width = 2^(nextpow2(1/s(1).DELTA)-1);
		%s_width = 64; % to match python
		s_width = proc_norm; % from read-in arguments or default at top

		if(amp_preserve == 0)
			for ii = 1:length(s)
				env = abs(hilbert(s(ii).DATA1));
				env = duglNoise.do_smoothing(env,s_width);
				s(ii).DATA1 = s(ii).DATA1./env;
			end
		elseif(amp_preserve == 1)
			for ii = 1:length(stas)
				for jj = 1:3
					%disp(s(stas_cmp3(ii,jj)).FILENAME)
					env(jj,:) = abs(hilbert(s(stas_cmp3(ii,jj)).DATA1));
					env(jj,:) = duglNoise.do_smoothing(env(jj,:),s_width);
				end
				env3 = prctile(env,95);
				env3 = duglNoise.do_smoothing(env3,s_width)';
				for jj = 1:3
					s(stas_cmp3(ii,jj)).DATA1 = s(stas_cmp3(ii,jj)).DATA1./env3;
				end
			end
		elseif(amp_preserve == 2)
			for ii = 1:length(s)
				env(ii,:) = abs(hilbert(s(ii).DATA1));
				env(ii,:) = duglNoise.do_smoothing(env(ii,:),s_width);
			end
			env_all = prctile(env,95);
			env_all = duglNoise.do_smoothing(env_all,s_width)';
			for ii = 1:length(s)
				s(ii).DATA1 = s(ii).DATA1./env_all;
			end
		end
		
	end

	%%% FFT %%%
	for ii = 1:length(s)
		spec.data(ii,:) = fft(s(ii).DATA1,NFFT);
	end

	%%% Spectral Whiten %%%
	if(proc_whit ~= 0)
		%s_width = 40;
		%s_width = 20; % to match python
		s_width = proc_whit; % from read-in arguments or default
		spec_env = abs(spec.data);

		for ii = 1:size(spec,1)
		%for ii = 1:size(spec_env,1)
			spec_env(ii,:) = duglNoise.do_smoothing(spec_env(ii,:),s_width);
		end

		if(amp_preserve == 0)
			spec.data = spec.data./spec_env;
		elseif(amp_preserve == 1)
			for ii = 1:length(stas)
				spec_env3 = prctile(spec_env(stas_cmp3(ii,:),:),95);
				spec_env3 = duglNoise.do_smoothing(spec_env3,s_width);
				for jj = 1:3
					spec.data(stas_cmp3(ii,jj),:) = spec.data(stas_cmp3(ii,jj),:)./spec_env3;
				end
			end
		elseif(amp_preserve == 2)
			spec_env_all = prctile(spec_env,95);
			spec_env_all = duglNoise.do_smoothing(spec_env_all,s_width);
			for ii = 1:size(spec.data,1);
				spec.data(ii,:) = spec.data(ii,:)./spec_env_all;
			end
		end
	end
	
	%%% Re filter %%%
	if(proc_filt_order ~= 0);
		HH = abs(H').^2;
		spec.data = bsxfun(@times,spec.data,HH); % non-causal, zero lag, filtfilt
	end



	% Fill in other info
	for ii = 1:length(s)
		spec.station{ii} = s(ii).KSTNM;
		spec.lat(ii) = s(ii).STLA;
		spec.lon(ii) = s(ii).STLO;
		spec.depth(ii) = s(ii).STEL;
		spec.cmp{ii} = s(ii).KCMPNM;
	end
	spec.time = sprintf('%04i %03i %02i:%02i:%02i:%03i',...
			s(1).NZYEAR,s(1).NZJDAY,s(1).NZHOUR,s(1).NZMIN,s(1).NZSEC,s(1).NZMSEC);
	spec.datenum = datenum(s(1).NZYEAR,1,1,s(1).NZHOUR,s(1).NZMIN,s(1).NZSEC) + s(1).NZJDAY;
	spec.delta = s(1).DELTA;
	spec.NFFT =  NFFT;
	spec.npts = s(1).NPTS;
	spec.highF = highF;
	spec.lowF = lowF;
	%spec.proc = 'filt, norm, whit';
	spec.proc = sprintf('filt_order: %u, norm: %u, whit: %u, amp: %u',proc_filt_order,proc_norm,proc_whit,amp_preserve);

      end %preprocess
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function [this] = xcorr(spec,varargin);
	% will cross correlate every pair available in "spec"
	%  with no regard for ordering or missing stations.
	%  this is just the math.
	%
	% Not yet designed to calculate distances and record header info
	%  Must initStack separately
	%
	%   stack = duglNoise.initStack;
	%   spec = duglNoise.preprocess(s);
	%   temp = duglNoise.xcorr(spec);
	%   stack = duglNoise.addToStack(stack,temp);

	maxshift = 2048;
	maxshift = 1024;
	maxshift = 512;
	[M,N] = size(spec.data);
	NFFT = N;

	build_headers = 0;
        ii = 1;
        while ii < length(varargin)
           switch varargin{ii}
              case 'full_headers'
		% Would calculate station distance and azimuth each time
		% Not needed by default, since adding to a stack
		build_headers = 1;
	      otherwise 
		disp(['varargin not understood: ',varargin{ii}])
	   end
	   ii = ii+1;
	end

	x2_conj = conj(spec.data);
    
	last = 1;
	%for ii = 1:M-1 %excludes autocorrelation
	for ii = 1:M
	    	%jj = [ii+1:M]; %excludes autocorrelation
	    	jj = [ii:M];
 	    	%x1_conj = conj(spec.data(ii,:));
 	    	%temp = bsxfun(@times,spec.data(jj,:),x1_conj);
		%temp = bsxfun(@times,spec.data(ii,:),x2_conj);
            	%temp2 = ifft(temp,NFFT,2,'symmetric'); %too slow
		%temp2 = zeros(size(jj,2),size(temp,2));

	    	temp = bsxfun(@times,spec.data(ii,:),x2_conj(jj,:));
		temp2 = zeros(size(temp));
%save testing.mat temp temp2 spec x2_conj
%monkey
		for j = 1:length(jj)
		    temp2(j,:) = ifft(temp(j,:),NFFT,2,'symmetric');
		    %temp2(j,:) = ifft(temp(jj(j),:),NFFT,2,'symmetric');
		end

		kk = last:last-1+length(jj);
	    	ncf(kk,:) = [temp2(:,NFFT-maxshift+1:end) temp2(:,1:maxshift)];
		last = last+length(jj);
	end


	k = 1;
	for ii = 1:M
	    	jj = ii:M;
		for j = jj
			this.sta1{k} = [spec.station{ii},'.',spec.cmp{ii},'-',spec.station{j},'.',spec.cmp{j}];
			this.sta2{k} = [spec.station{j},'.',spec.cmp{j},'-',spec.station{ii},'.',spec.cmp{ii}];
			%disp([spec.station{ii},'.',spec.cmp{ii},'-',spec.station{j},'.',spec.cmp{j}])
			k = k+1;
		end
	end

	this.ncf = ncf;
	this.time = spec.time;
	this.datenum = spec.datenum;

	%this.delta  = spec.delta
	%this.NFFT   = spec.NFFT 
	%this.npts   = spec.npts 
	this.highF  = spec.highF;
	this.lowF   = spec.lowF ;
	this.proc   = spec.proc ;

	%disp('TODO! Build other headers (distance, locations, azimuths, etc)');

      end % xcorr
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function stack = initStack(varargin);
	% Initiates a stack object, with every possible station and component
	%
	%             ncf: [2628x1024 double]    noise correlation functions
	%         stacmp: {2628x1 cell}          name, i.e.: "B4850.HHN-D4850.HHN"
	%           loc1: [2628x3 double]        X,Y,Z of station1
	%           loc2: [2628x3 double]        X,Y,Z of station2 
	%           dist: [2628x1 double]        distance between the two
	%             az: [2628x1 double]        azimuth
	%            dip: [2628x1 double]        dip
	%           hits: [2628x1 double]        number of hits that went into this NCF
	%          times: {2628x1 cell}          list of every time that went into this NCF
	%       datenums: {2628x1 cell}          list of every time that went into this NCF
	%             ss: {24x1 cell}            each unique station
	%    ss_datenums: {24x1 cell}            each individual station's available times
	%      ss_lookup: [2628x2 double]        in the full ncf matrix, a lookup table of which stations
	%          notes: ''
	%          highF: []
	%           lowF: []
	%           proc: ''

	ALL_stations = duglSet.Stations;
	ALL_components = duglSet.Components;


	if(length(varargin)>0)
		ALL_stations = varargin{1};
	end
	if(length(varargin)>1)
		ALL_components = varargin{2};
	end
	ALL_stations

	sta_file = [duglSet.MAP_DIR,'/stations.txt'];
	%if(length(varargin)>0)
	%	sta_file = varargin{1};
	%end
	fid = fopen(sta_file);
	stations_temp = textscan(fid,'%s %f %f %f %f %f %f');
	fclose(fid);
	

	%M = length(ALL_stations)*length(ALL_components);
	k = 1;
	for ista = 1:length(ALL_stations)
	   for icmp = 1:length(ALL_components)
		ALL_stacmp{k} = [ALL_stations{ista},'.',ALL_components{icmp}];
		ALL_ss_index(k) = ista;

		% Get matching station location
		% Calculating distance will be redundant x3, but only calculated once...
		imatch = find(strcmp(stations_temp{1},ALL_stations{ista}));
		ALL_loc(k,:)  = [stations_temp{2}(imatch) stations_temp{3}(imatch) stations_temp{4}(imatch)];

		k = k+1;
	   end
	end
	%MM = M*(M-1)/2; % number of pairs % this one excludes autocorrelation
	M = length(ALL_stacmp);
	MM = M*(M-1)/2+M; % number of pairs

	maxshift = 2048;
	maxshift = 1024;
	maxshift = 512;
	stack.ncf = zeros(MM,maxshift*2);
	stack.stacmp = cell(MM,1);
	stack.loc1 = zeros(MM,3);
	stack.loc2 = zeros(MM,3);
	stack.dist = zeros(MM,1);
	stack.az   = zeros(MM,1);
	stack.dip  = zeros(MM,1);
	stack.hits = zeros(MM,1);
	stack.times = cell(MM,1);
	for ii = 1:MM
		stack.times{ii} = {};
	end
	stack.datenums = cell(MM,1);
	stack.ss = ALL_stations';
	stack.ss_datenums = cell(length(ALL_stations),1);
	stack.ss_lookup = zeros(MM,2);
	stack.notes = '';
	stack.highF = [];
	stack.lowF = [];
	stack.proc = '';

	k=1;
	for ii = 1:M
	    	for jj = ii:M
			% Name of pair
			stack.stacmp{k} = [ALL_stacmp{ii},'-',ALL_stacmp{jj}];
			stack.ss_lookup(k,:) = [ALL_ss_index(ii) ALL_ss_index(jj)];

			
			% Distance
			stack.loc1(k,:) = ALL_loc(ii,:);
			stack.loc2(k,:) = ALL_loc(jj,:);
			xdiff = duglProc.DistVicenty(ALL_loc(ii,2),ALL_loc(ii,1),...
						     ALL_loc(ii,2),ALL_loc(jj,1));
			ydiff = duglProc.DistVicenty(ALL_loc(ii,2),ALL_loc(ii,1),...
						     ALL_loc(jj,2),ALL_loc(ii,1));
			zdiff = ALL_loc(ii,3) - ALL_loc(jj,3);
			stack.dist(k) = sqrt(xdiff^2 + ydiff^2 + zdiff^2);


			% Azimuth and dip
			% Defined as direction from loc1 to loc2, ii to jj
			% Azimuth is degrees clockwise from north
			% Dip is degrees downward from horizontal
			%  (negative degrees for upward)
			az = azimuth(ALL_loc(ii,2),ALL_loc(ii,1),ALL_loc(jj,2),ALL_loc(jj,1));
			stack.az(k) = az;
			dip = atan2(zdiff,sqrt(xdiff^2 + ydiff^2))*180/pi;
			stack.dip(k) = dip;
			
			%disp([stack.stacmp{k},' ',num2str(stack.dist(k))])
			%disp([stack.stacmp{k}, ' ',num2str(az),' ',num2str(dip)])
			%plot3([ALL_loc(ii,1),ALL_loc(jj,1)],[ALL_loc(ii,2),ALL_loc(jj,2)],[ALL_loc(ii,3),ALL_loc(jj,3)],'k')


			k = k+1;
		end
	end
	stack.dist( isnan(stack.dist) ) = 0;
	%disp(num2str(k - 1));
	%disp(num2str(MM));
	
      end % initStack
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function filename = initHDF5(filename, stack);
	hdf5write(filename,'/ss',stack.ss)
	hdf5write(filename,'/stacmp',stack.stacmp,'writemode','append')
	hdf5write(filename,'/loc1',stack.loc1,'writemode','append')
	hdf5write(filename,'/loc2',stack.loc2,'writemode','append')
	hdf5write(filename,'/dist',stack.dist,'writemode','append')
	hdf5write(filename,'/az',stack.az,'writemode','append')
	hdf5write(filename,'/dip',stack.dip,'writemode','append')
	hdf5write(filename,'/ss_lookup',stack.ss_lookup,'writemode','append')
	hdf5write(filename,'/notes',stack.notes,'writemode','append')
	hdf5write(filename,'/highF',stack.highF,'writemode','append')
	hdf5write(filename,'/lowF',stack.lowF,'writemode','append')
	hdf5write(filename,'/proc',stack.proc,'writemode','append')
      end

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function filename = addHDF5(filename,stack,index)
	for ii = 1:length(stack.datenums)
		if(length(stack.datenums{ii})>0)
			this_datenum = stack.datenums{ii};
			this_time  = stack.times{ii};
			break
		end
	end
	if(length(this_time)==0)
		% Nothing was actually recorded for this file
		contintue
	end
	hdf5write(filename,['/',index,'/ncf'],stack.ncf,'writemode','append')
	hdf5write(filename,['/',index,'/hits'],stack.hits,'writemode','append')
	hdf5write(filename,['/',index,'/time'],this_time,'writemode','append')
	hdf5write(filename,['/',index,'/datenum'],this_datenum,'writemode','append')
	%hdf5write(filename,['/',index,'/times'],stack.times,'writemode','append')
	%hdf5write(filename,['/',index,'/datenums'],stack.datenums,'writemode','append')
	%hdf5write(filename,['/',index,'/ss_datenums'],stack.ss_datenums,'writemode','append')
      end

	
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function stack1 = mergeStack(stack1,stack2)
	%  stack1 = mergeStack(stack1,stack2) 
	%
	%   This will be faster than addStack, since it is assumed the 2
	%   stack objects are identical in length and order. 
	%   Output is a blind merge between the two.
	if(length(stack1.stacmp) == length(stack2.stacmp) && strcmp(stack1.stacmp{1},stack1.stacmp{1}))
	    % Can assume the two stacks are identical, quick merge
	    for ii = 1:length(stack1.stacmp)
		stack1.ncf(ii,:) = stack1.ncf(ii,:)+stack2.ncf(ii,:)./stack2.hits(ii);
		stack1.hits(ii) = stack1.hits(ii)+stack2.hits(ii);
	    end
% DCB 06/8/2017 Need to properly weight hits

	%DCB 04/10/2017, removing these trackers since stacking takes too long
	%    for ii = 1:length(stack1.stacmp)
	%	stack1.times{ii} = [stack1.times{ii} stack2.times{ii}];
	%	stack1.datenums{ii} = [stack1.datenums{ii} stack2.datenums{ii}];
	%    end
	%    for ii = 1:length(stack1.ss)
	%	stack1.ss_datenums{ii} = [stack1.ss_datenums{ii} stack2.ss_datenums{ii}];
	%    end

	else
		disp('ERROR! stacks not compatible. Not coded to handle this case yet')
	end
	if(length(stack1.proc)==0)
		stack1.highF = stack2.highF;  % assumes no user-error of mixing different proc or freq 
		stack1.lowF  = stack2.lowF;
		stack1.proc  = stack2.proc;
	elseif(strcmp(stack1.proc,stack2.proc)==0)
		disp('ERROR! stacks contain different processing flags!')
	end
      end % mergeStack
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function stack = addStack(stack,this)
	% stack = addStack(stack,temp)
	%
	% temp in this case is a quickly created stack object from duglNoise.xcorr
	% and does not necessarily contain all stations, header info, or use correct ordering.
	%
	num_stacked = 0;
	temp_ss_hit = zeros(size(stack.ss));
	for ista = 1:length(this.sta1)
	   imatch1 = find(strcmp(this.sta1{ista},stack.stacmp));
	   imatch2 = find(strcmp(this.sta2{ista},stack.stacmp));
	   if(length(imatch1)==1)
		stack.ncf(imatch1,:)    = stack.ncf(imatch1,:) + this.ncf(ista,:);
		stack.hits(imatch1)     = stack.hits(imatch1)+1;
		stack.times{imatch1}    = [stack.times{imatch1},{this.time}];
		stack.datenums{imatch1} = [stack.datenums{imatch1},this.datenum];
		temp_ss_hit( stack.ss_lookup(imatch1,1) ) = 1;
		temp_ss_hit( stack.ss_lookup(imatch1,2) ) = 1;
		num_stacked = num_stacked+1;

	   elseif(length(imatch2)==1)
		stack.ncf(imatch2,:)    = stack.ncf(imatch2,:) + flip(this.ncf(ista,:));
		stack.hits(imatch2)     = stack.hits(imatch2)+1;
		stack.times{imatch2}    = [stack.times{imatch2},{this.time}];
		stack.datenums{imatch2} = [stack.datenums{imatch2},this.datenum];
		temp_ss_hit( stack.ss_lookup(imatch2,1) ) = 1;
		temp_ss_hit( stack.ss_lookup(imatch2,2) ) = 1;
		num_stacked = num_stacked+1;
	   else
		disp(['No Match! ',this.sta{ista}])
	   end
	end
	ss_hit = find(temp_ss_hit)';
	for ihit = ss_hit
		stack.ss_datenums{ ihit } = [stack.ss_datenums{ihit},this.datenum];
	end
	disp([num2str(num_stacked),' added to global stack'])

	stack.highF = this.highF;  % assumes no user-error of mixing different proc or freq 
	stack.lowF  = this.lowF;
	stack.proc  = this.proc;

      end % addStack
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %function stack = addDay(stack,this)
      %  
      %end % addDay
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function s = stack2sac(stack,str)
	use = duglNoise.findStationSubset(stack,str);
	disp(['Found ',num2str(length(use)),' matches'])
	%stack
	%use
	s = sacstruct(length(use));
	delta = .04;
	maxshift = 2048;
	maxshift = 1024;
	maxshift = 512;
	for i = 1:length(use);
		ii = use(i);
		spl = regexp(stack.stacmp{ii},'-','split');
		spl1 = regexp(spl{1},'\.','split');
		spl2 = regexp(spl{2},'\.','split');

		s(i).FILENAME = stack.stacmp{ii};
		s(i).KCMPNM   = [spl1{2}(3) spl2{2}(3)];
		s(i).DATA1    = stack.ncf(ii,:)';
		s(i).DELTA    = delta;
		s(i).O        = 0;
		s(i).NPTS     = length(s(i).DATA1);
		s(i).B	      = -delta*maxshift;
		s(i).E        = delta*maxshift;
		s(i).USER0    = stack.hits(ii);
     		s(i).NVHDR    = 6;
     		s(i).IFTYPE   = 'ITIME';
     		s(i).KSTNM    = spl1{1};
     		s(i).EVLA     = stack.loc1(ii,2);
     		s(i).EVLO     = stack.loc1(ii,1);
		s(i).EVEL     = stack.loc1(ii,3)./1000;
		%s(i).EVDP     = -stack.loc1(ii,3)./1000;
     		s(i).STLA     = stack.loc2(ii,2);
     		s(i).STLO     = stack.loc2(ii,1);
		%s(i).STDP     = -stack.loc2(ii,3)./1000;
		s(i).STEL     = stack.loc2(ii,3)./1000;
     		s(i).USER6    =-(stack.dist(ii))./1000;
     		s(i).DIST     =(stack.dist(ii))./1000;
     		s(i).LCALDA   =0;
     		s(i).LEVEN    =1;
		s(i).AZ       = stack.az(ii);
		s(i).USER1    = stack.az(ii);
		s(i).USER2    = stack.dip(ii);
				
		
		% LATERAL Distance
		ldist = duglProc.DistVicenty( s(i).EVLA, s(i).EVLO, s(i).STLA, s(i).STLO );
		ldist( isnan(ldist) ) = 0;
		s(i).USER3  = ldist./1000;

		% VERTICAL Distance
		vd  = abs(s(i).EVEL - s(i).STEL);
		s(i).USER4 = vd./1000;


     		%writesac(s);  
	end
	disp('Trimming traces')
	cut = [-10 10];
	for i = 1:length(use)
		s(i).B = cut(1);
		s(i).E = cut(2);
		npts = s(i).NPTS;
		b = npts/2 + cut(1)/s(1).DELTA;
		e = npts/2 + cut(2)/s(1).DELTA;
		s(i).DATA1 = s(i).DATA1(b:e);
		s(i).NPTS = length(s(i).DATA1);
	end
      end %stack2sac
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function addPaths(stack,varargin)
	use_only_hit = 0; % uses all in stack by default
	use_colors = 0;
	jj = 1;
	while(jj <= length(varargin))
		switch varargin{jj}
		   case 'hit'
			use_only_hit = 1;
			use = find(stack.hits);
		   case 'subset'
			jj = jj+1;
			str = varargin{jj};
			use = duglNoise.findStationSubset(stack,str);
			disp(['Found ',num2str(length(use)),' matches'])
			use_only_hit = 1;
		   case 'vel'
			use_colors = 1;
			jj = jj+1;
			cc = varargin{jj};
			if(length(cc) ~= length(stack.stacmp))
				disp('Error! Vector of colors/velocities not matching length of raypath object')
				return 
			end
	           otherwise 
		     	disp(['varargin not understood: ',varargin{ii}])
		end
		jj = jj + 1;
	end
	

	disp('Adding paths to current figure')
	if(use_only_hit)
		disp(['Using ',num2str(length(use)),' paths'])
		for iuse = 1:length(use)
			ii = use(iuse);
			%disp([stack.stacmp{ii},' ',num2str(stack.loc1(ii,1)),' ',num2str(stack.loc2(ii,1))])
			if(use_colors==0)
				plot3([stack.loc1(ii,1),stack.loc2(ii,1)],[stack.loc1(ii,2),stack.loc2(ii,2)],[stack.loc1(ii,3),stack.loc2(ii,3)],'k');
			else
				plot3([stack.loc1(ii,1),stack.loc2(ii,1)],[stack.loc1(ii,2),stack.loc2(ii,2)],[stack.loc1(ii,3),stack.loc2(ii,3)],'k');
			end
		end
	else
		% Just use all paths listed in stack
		disp(['Using ',num2str(length(stack.stacmp)),' paths'])
		for ii = 1:length(stack.stacmp)
			if(use_colors==0)
				plot3([stack.loc1(ii,1),stack.loc2(ii,1)],[stack.loc1(ii,2),stack.loc2(ii,2)],[stack.loc1(ii,3),stack.loc2(ii,3)],'k');
			else
				plot3([stack.loc1(ii,1),stack.loc2(ii,1)],[stack.loc1(ii,2),stack.loc2(ii,2)],[stack.loc1(ii,3),stack.loc2(ii,3)],'color',[cc(ii,:)]);
			end
		end
	end
      end	
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function smoothed = do_smoothing(d,s_width)
	d = abs(d);
	win = ones(s_width*2,1)./(s_width*2);
	smoothed = conv(d,win,'same');
      end %do_smoothing
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function H = build_butterworth(s,lowF,highF,n,p)
	% H = build_butterworth(s,lowF,highF,n,p)
	%
	% To apply:
	%  
	% sp = fft(s(1).DATA1,NFFT);
	% sp2 = (sp.*hb);
	% d2 = real(ifft(sp2,NFFT));
	% d2(nsamp+1:end) = [];
	%
	% Results identical to "filter(B,A,data)"
	%
	% To get a-causal, zero-phase "filtfilt"
	% sp2 = time-reverse[ sp.*H ] .*H
	%     = conj(sp).*conj(H).*H
	%     = sp.*H.*conj(H)
	%     = sp.*(abs(H).^2)
	if(nargin < 3)
		disp('ERROR! Usage:')
		disp('s = bp(s,low,high,n,p)')
		return
	elseif (nargin < 4)
		n=2;
		p=2;
	end
	for ii = 1:length(s)
		if(s(ii).NPTS~=s(1).NPTS)
			disp('ERROR! must all be same length')
			return
		end
	end
	nsamp = s(1).NPTS;
	NFFT = 2^nextpow2(nsamp);
	fnyq = (1./s(1).DELTA)./2;
	%[Z,P,K] = butter(p,[lowF highF]./fnyq);
	[B,A] = butter(p,[lowF highF]./fnyq);
	[H,W] = freqz(B,A,NFFT,'whole');
      end %build_butterworth
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function [pwdB2,FF2] = run_PSD(s)
	for jj = 1:length(s)
		NFFT = 2^nextpow2(s(jj).NPTS);
		data = s(jj).DATA1;
		data = data./1e9;
		[pw,FF] = pwelch(data,[],[],NFFT,1/s(jj).DELTA,'onesided');
		pwdB = 10*log10(pw);


		% Down sample...
		pwdB(1) = [];
		FF(1) = [];
		if(length(FF)~=524288)
			disp('Size not right for hardcoded downsample')
			continue
		end
		dsamp = 2^10;
		pwdB_resh = reshape(pwdB,dsamp,512);
		pwdB2(jj,:) = mean(pwdB_resh);
		FF_resh = reshape(FF,dsamp,512);
		FF2(jj,:) = mean(FF_resh);
	end
      end %colect_PSD
% DCB todo: connect this to other PSD scripts?
% finish plot_PSD?
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function [sout] = rotate9(s);
	if(length(s)~=9)
		disp('ERROR! Need to give exactly 9 components')
		return
	end
	az = s(1).USER1;  % should all have the same 
	inc = 90 - s(1).USER2; % user2 is "geology dip"
	disp(sprintf('az = %f, inc = %f ',az,inc))
	spl = regexp(s(1).FILENAME,'-','split');
	spl1 = regexp(spl{1},'\.','split');
	spl2 = regexp(spl{2},'\.','split');
	if(any(90-[s.USER2]~=inc))
		disp('ERROR! All must be of same pair, different angles given')
		return
	end
		
	for ii = 1:length(s)
		switch s(ii).KCMPNM
		  case 'ZZ'; ZZ=s(ii);
		  case 'ZE'; ZE=s(ii);
		  case 'ZN'; ZN=s(ii);
		  case 'EZ'; EZ=s(ii);
		  case 'EE'; EE=s(ii);
		  case 'EN'; EN=s(ii);
		  case 'NZ'; NZ=s(ii);
		  case 'NE'; NE=s(ii);
		  case 'NN'; NN=s(ii);
		end
	end
	
	% rotation matrix defined from: https://service.iris.edu/irisws/rotation/docs/1/help/
	% where: 
	%   - [L Q T] = M * [Z N E]
	%   - inc is incidence angle away from vertical
	%   - azimuth is clockwise from north (like a compass)

	M = [cosd(inc)  -sind(inc)*sind(az)   -sind(inc)*cosd(az);
		sind(inc)   cosd(inc)*sind(az)   cosd(inc)*cosd(az);
		0         -cosd(az)          sind(az)];
	out = kron(M,M)*[ZZ.DATA1 ZE.DATA1 ZN.DATA1 EZ.DATA1 EE.DATA1 EN.DATA1 NZ.DATA1 NE.DATA1 NN.DATA1]';

	% Rotated into LQT, where L is P-wave, Q is SV, T is SH

	% copy over other misc header info
	LL=ZZ; LQ=ZZ; LT=ZZ; QL=ZZ; QQ=ZZ; QT=ZZ; TL=ZZ; TQ=ZZ; TT=ZZ; 

	LL.DATA1 = out(1,:); LL.FILENAME = [spl1{1},'.HHL-',spl2{1},'.HHL']; LL.KCMPNM = 'LL'; 
	LQ.DATA1 = out(2,:); LQ.FILENAME = [spl1{1},'.HHL-',spl2{1},'.HHQ']; LQ.KCMPNM = 'LQ';
	LT.DATA1 = out(3,:); LT.FILENAME = [spl1{1},'.HHL-',spl2{1},'.HHT']; LT.KCMPNM = 'LT'; 
	QL.DATA1 = out(4,:); QL.FILENAME = [spl1{1},'.HHQ-',spl2{1},'.HHL']; QL.KCMPNM = 'QL'; 
	QQ.DATA1 = out(5,:); QQ.FILENAME = [spl1{1},'.HHQ-',spl2{1},'.HHQ']; QQ.KCMPNM = 'QQ'; 
	QT.DATA1 = out(6,:); QT.FILENAME = [spl1{1},'.HHQ-',spl2{1},'.HHT']; QT.KCMPNM = 'QT'; 
	TL.DATA1 = out(7,:); TL.FILENAME = [spl1{1},'.HHT-',spl2{1},'.HHL']; TL.KCMPNM = 'TL'; 
	TQ.DATA1 = out(8,:); TQ.FILENAME = [spl1{1},'.HHT-',spl2{1},'.HHQ']; TQ.KCMPNM = 'TQ'; 
	TT.DATA1 = out(9,:); TT.FILENAME = [spl1{1},'.HHT-',spl2{1},'.HHT']; TT.KCMPNM = 'TT'; 

	sout = [LL LQ LT QL QQ QT TL TQ TT];
      end
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function [stack2] = rotateStack(stack);
	
	station_lookups_to_rotate = unique(stack.ss_lookup,'rows');
	disp(sprintf('Rotating %u station pairs',size(station_lookups_to_rotate,1)))
	
	% Copy over generic information to new stack array
	stack2.ss = stack.ss;
	stack2.ss_datenums = stack.ss_datenums;
	stack2.notes = stack.notes;
	stack2.highF = stack.highF;
	stack2.lowF  = stack.lowF ;
	stack2.proc  = stack.proc ;
	ist = 0;

	
	for ii = 1:size(station_lookups_to_rotate,1)
		if(station_lookups_to_rotate(ii,1)==station_lookups_to_rotate(ii,2))
			continue
		end
		icmp = find(sum(bsxfun(@eq, stack.ss_lookup, station_lookups_to_rotate(ii,:)),2)==2);
		%disp(sprintf('%u',length(icmp)))
		if(length(icmp)<9)
			disp(sprintf('ERROR! %s does not have all 9 components',stack.stacmp(icmp(1))))
		end
		if(any(stack.az(icmp)~=stack.az(icmp(1))))
			disp(sprintf('ERROR! %s az not consistent with other components',stack.stacmp(icmp(1))))
		end

		% Now have the 9 components needed for rotation
		IN = zeros(length(stack.ncf(icmp(1),:)),9);
		az = stack.az(icmp(1));
		inc = 90 - stack.dip(icmp(1));

		% sort into a properly ordered array
		for ii = 1:9
			spl = regexp(stack.stacmp{icmp(ii)},'-','split');
			spl1 = regexp(spl{1},'\.','split');
			spl2 = regexp(spl{2},'\.','split');
			this_cmp = [spl1{2}(3) spl2{2}(3)];
			switch this_cmp
			  case 'ZZ'; IN(:,1) = stack.ncf(icmp(ii),:); %disp(['ZZ = ',stack.stacmp{icmp(ii)}])
			  case 'ZE'; IN(:,2) = stack.ncf(icmp(ii),:); %disp(['ZE = ',stack.stacmp{icmp(ii)}])
			  case 'ZN'; IN(:,3) = stack.ncf(icmp(ii),:); %disp(['ZN = ',stack.stacmp{icmp(ii)}])
			  case 'EZ'; IN(:,4) = stack.ncf(icmp(ii),:); %disp(['EZ = ',stack.stacmp{icmp(ii)}])
			  case 'EE'; IN(:,5) = stack.ncf(icmp(ii),:); %disp(['EE = ',stack.stacmp{icmp(ii)}])
			  case 'EN'; IN(:,6) = stack.ncf(icmp(ii),:); %disp(['EN = ',stack.stacmp{icmp(ii)}])
			  case 'NZ'; IN(:,7) = stack.ncf(icmp(ii),:); %disp(['NZ = ',stack.stacmp{icmp(ii)}])
			  case 'NE'; IN(:,8) = stack.ncf(icmp(ii),:); %disp(['NE = ',stack.stacmp{icmp(ii)}])
			  case 'NN'; IN(:,9) = stack.ncf(icmp(ii),:); %disp(['NN = ',stack.stacmp{icmp(ii)}])
			end
		end
		M = [cosd(inc)  -sind(inc)*sind(az)   -sind(inc)*cosd(az);
			sind(inc)   cosd(inc)*sind(az)   cosd(inc)*cosd(az);
			0         -cosd(az)          sind(az)];
		out = kron(M,M)*IN';

		% Now re-adjust names and fill output stack
		for ii = 1:9
			stack2.ncf(ist+ii,:) = out(ii,:);
		end
		stack2.stacmp{ist+1} = [spl1{1},'.HHL-',spl2{1},'.HHL'];
		stack2.stacmp{ist+2} = [spl1{1},'.HHL-',spl2{1},'.HHQ'];
		stack2.stacmp{ist+3} = [spl1{1},'.HHL-',spl2{1},'.HHT'];
		stack2.stacmp{ist+4} = [spl1{1},'.HHQ-',spl2{1},'.HHL'];
		stack2.stacmp{ist+5} = [spl1{1},'.HHQ-',spl2{1},'.HHQ'];
		stack2.stacmp{ist+6} = [spl1{1},'.HHQ-',spl2{1},'.HHT'];
		stack2.stacmp{ist+7} = [spl1{1},'.HHT-',spl2{1},'.HHL'];
		stack2.stacmp{ist+8} = [spl1{1},'.HHT-',spl2{1},'.HHQ'];
		stack2.stacmp{ist+9} = [spl1{1},'.HHT-',spl2{1},'.HHT'];
	
		stack2.loc1(ist+1:ist+9,:) = bsxfun(@times,ones(9,3),stack.loc1(icmp(1),:));
		stack2.loc2(ist+1:ist+9,:) = bsxfun(@times,ones(9,3),stack.loc2(icmp(1),:));
		stack2.dist(ist+1:ist+9) = ones(9,1).*stack.dist(icmp(1));
		stack2.dip(ist+1:ist+9) = ones(9,1).*stack.dip(icmp(1));
		stack2.az(ist+1:ist+9) = ones(9,1).*stack.az(icmp(1));
		% times / datenums / hits?
		stack2.hits(ist+1:ist+9) = ones(9,1);
		
		ist = ist+9;
	end
		


      end %rotateStack
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%









      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function [RR,RT,TR,TT] = rotate4(EE,EN,NE,NN,baz);
      %  i = taup(1).incidentangle;
      %  baz = taup(1).azimuth+180;
%	M = [cosd(i)  -sind(i)*sind(baz)   -sind(i)*cosd(baz);
%		sind(i)   cosd(i)*sind(baz)   cosd(i)*cosd(baz);
%		0         -cosd(baz)          sind(baz)];
	M = [ cosd(baz) sind(baz);
		-sind(baz) cosd(baz)];
disp('Fix this: baz or az?')
disp('Fix this: finish kronicker command')
	
      end
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function [EE,EN,NE,NN] = pull4(str)
      % [EE,EN,NE,NN] = pull4(str)
      % str can be any regexp. for 2 stations ie: str = './final/*800*300*';
	spl = regexp(str,'\*','split')
      %   spl =  './final/'    '800'    'YATES'    ''

	disp([spl{1:end-2},'*HHE*',spl{end-1},'*HHE*',spl{end}])
	dir([spl{1:end-2},'*HHE*',spl{end-1},'*HHE*',spl{end}])
	disp([spl{1:end-2},'*HHE*',spl{end-1},'*HHN*',spl{end}])
	dir([spl{1:end-2},'*HHE*',spl{end-1},'*HHN*',spl{end}])
	disp([spl{1:end-2},'*HHN*',spl{end-1},'*HHE*',spl{end}])
	dir([spl{1:end-2},'*HHN*',spl{end-1},'*HHE*',spl{end}])
	disp([spl{1:end-2},'*HHN*',spl{end-1},'*HHN*',spl{end}])
	dir([spl{1:end-2},'*HHN*',spl{end-1},'*HHN*',spl{end}])
	
	EE = readsac([spl{1:end-2},'*HHE*',spl{end-1},'*HHE*',spl{end}]);
	EN = readsac([spl{1:end-2},'*HHE*',spl{end-1},'*HHN*',spl{end}]);
	NE = readsac([spl{1:end-2},'*HHN*',spl{end-1},'*HHE*',spl{end}]);
	NN = readsac([spl{1:end-2},'*HHN*',spl{end-1},'*HHN*',spl{end}]);


	%s = readsac(str);
	%if(length(s)<4)
	%	disp('ERROR!, Not at least 4 components pulled')
	%	disp('Only found: ')
	%	for ii = 1:length(s)
	%		disp(s(ii).FILENAME)
	%	end
	%	return
	%end
	%for ii = 1:length(s)
	%   switch s(ii).KCMPNM
	%       case 'EE'
	%	EE = s(ii);
	%       case 'EN'
	%	EN = s(ii);
	%       case 'NE'
	%	NE = s(ii);
	%       case 'NN'
	%	NN = s(ii);
	%       otherwise
	%	disp(['No match! ',s(ii).KSTNM,' ',s(ii).FILENAME])
        %    end
	%end
      end % pull4
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function [ZZ,ZE,ZN,EZ,EE,EN,NZ,NE,NN] = pull9(str)
      % [ZZ,ZE,ZN,EZ,EE,EN,NZ,NE,NN] = duglNoise.pull9(str)
      % str can be any regexp. ie: str = './final/*800*300*';
	s = readsac(str);
	if(length(s)~=9)
		disp('ERROR!, Not exactly 9 components pulled')
		disp('Only found: ')
		for ii = 1:length(s)
			disp(s(ii).FILENAME)
		end
		return
	end
	for ii = 1:9
	   switch s(ii).KCMPNM
	       case 'ZZ'
		ZZ = s(ii);
	       case 'ZE'
		ZE = s(ii);
	       case 'ZN'
		ZN = s(ii);
	       case 'EZ'
		EN = s(ii);
	       case 'EZ'
		EZ = s(ii);
	       case 'EE'
		EZ = s(ii);
	       case 'EN'
		EN = s(ii);
	       case 'NZ'
		NZ = s(ii);
	       case 'NE'
		NE = s(ii);
	       case 'NN'
		NN = s(ii);
	       otherwise
		disp(['No match! ',s(ii).KSTNM,' ',s(ii).FILENAME])
            end
	end
      end % pull9
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function stack = pullSubsetStack(stack,str);
	use = duglNoise.findStationSubset(stack,str);
	disp(['Found ',num2str(length(use)),' matches'])

	if(isfield(stack,'ncf'))
		stack.ncf       = stack.ncf(use,:);
	end
	if(isfield(stack,'loc1'))
		stack.loc1      = stack.loc1(use,:);
	end
	if(isfield(stack,'loc2'))
		stack.loc2      = stack.loc2(use,:);
	end
	if(isfield(stack,'dist'))
		stack.dist      = stack.dist(use);
	end
	if(isfield(stack,'az'))
		stack.az        = stack.az(use);
	end
	if(isfield(stack,'dip'))
		stack.dip       = stack.dip(use);
	end
	if(isfield(stack,'hits'))
		stack.hits      = stack.hits(use);
	end
	if(isfield(stack,'ss_lookup'))
		stack.ss_lookup = stack.ss_lookup(use,:);
	end
	if(isfield(stack,'stacmp'))
		stack.stacmp    = stack.stacmp(use);
	end
	if(isfield(stack,'times'))
		stack.times     = stack.times(use);
	end
	if(isfield(stack,'datenums'))
		stack.datenums  = stack.datenums(use);
	end
	if(isfield(stack,'phvel'))
		stack.phvel  = stack.phvel(use);
	end
	if(isfield(stack,'grvel'))
		stack.grvel  = stack.grvel(use);
	end

      end % pullSubsetStack
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function stack = stackBP(stack,lowF,highF,n,p);
	if(nargin < 3)
		disp('ERROR! Usage:')
		disp('stack = duglNoise.stackBP(stack,low,high,n,p)')
		return
	elseif (nargin < 4)
		n=2;
		p=2;
	end

	
	nsamp = size(stack.ncf,2);
	taperl = round(nsamp*.10);
	taper = linspace(0,1,taperl)';
	stack.ncf(:,1:taperl) = bsxfun(@times,stack.ncf(:,1:taperl),taper');
	stack.ncf(:,nsamp-taperl+1:end) = bsxfun(@times,stack.ncf(:,nsamp-taperl+1:end),flip(taper)');
	
	sps = 1/.04;
	fnyq = sps/2;
	[B,A] = butter(p,[lowF highF]./fnyq);
	if n==1; stack.ncf=filter(B,A,stack.ncf')';
	elseif n==2; stack.ncf=filtfilt(B,A,stack.ncf')';
	else disp('ERROR! Number of passes can be only 1 or 2');
	end

	stack.notes = [stack.notes,' BP:',num2str(lowF),'-',num2str(highF),' ;'];
      end % stackBP
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function nprs(stack,str,varargin)
	% duglNoise.prs(stack,'D4850&HHZ.*HHZ')
	%  search string can be multiple, subsequent searches
	%  separated by &
	args_in = [];
	flags.clear_before_plotting = 1;
	flags.scaling_factor = 12;
	maxshift = 2048;
	maxshift = 1024;
	maxshift = 512;
	delta = .04;

	jj = 1;
	while(jj <= length(varargin))
		switch varargin{jj}
		   case 'noclear'
			flags.clear_before_plotting = 0;
		   case 'scaling_factor'
			jj = jj+1;
			flags.scaling_factor = varargin{jj};
		end
		jj = jj+1;
	end
			
		
	if(flags.clear_before_plotting)
		clf
	end

	% Pull some subset out of the total stack
	use = duglNoise.findStationSubset(stack,str);
	disp(['Found ',num2str(length(use)),' matches'])


	%Need to scale the seismogram's amplitude
	%  First, loop through to find max y-axis (distance)
	maxd = max(stack.dist(use));

	%Determin scaling
	scaling=maxd/flags.scaling_factor;
	
	
	%Loop through and plot data
	for  ii = use
		tt = [-maxshift:maxshift-1].*delta;
		norm = max(stack.ncf(ii,:))/scaling;
		hh = plot(tt,stack.ncf(ii,:)./norm-stack.dist(ii),'k');
		set(hh,'DisplayName',stack.stacmp{ii})
		set(hh,'LineWidth',.1)
		ax = axis;
		hh = text(ax(2)*1.05,-stack.dist(ii),stack.stacmp{ii});
		hold on
	end
	xlabel('Time (s)')
	ylabel('Distance (km)')
	siz = get(gca,'Position')
	siz(3) = siz(3)*.8;
	set(gca,'Position',[siz])

      end % nprs
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function plotSpec(spec)
        % spec = duglNoise.preprocess(s);
	% plotSpec(spec)
	for ii = 1:size(spec.data,1)
		ff = 1/spec.delta/2 * linspace(0,1,spec.NFFT/2+1);
		NFFT = spec.NFFT;
		spect = spec.data(ii,:);
		size(spect)
		size(spect(1:NFFT/2+1))
		size(ff)
		max(abs(spect))
		loglog(ff, abs(spect(1:NFFT/2+1)) ./ max(abs(spect)))
		hold on
	end
      end % plotSpec
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function use = findStationSubset(stack,str)
	use_final = 1:length(stack.stacmp);
	spl = regexp(str,'&','split');
	for ispl = 1:length(spl)
		str2 = spl{ispl};
		%disp(['Pulling out: ',str2])
		%%str2(str2=='*')='.';

		use = use_final;
		use_final = [];

		if(regexp(str2,'dip'))
			value = str2num(str2(5:end));
			operator = str2(4);
			switch operator
			   case '<'
				good = find(abs([stack.dip])<value);
			   case '>'
				good = find(abs([stack.dip])>value);
			   otherwise
				good = [];
				disp(sprintf('Not understood: %s',str2))
			   end
			   disp(sprintf('Found %u with dip%s%u',length(good),operator,value))
			use_final = [use_final,good(:)'];
		else % Must be a station name?
			for ii = use
				if(regexp(stack.stacmp{ii},str2))
					%disp(stack.stacmp{ii});
					use_final = [use_final,ii];
				end
			end
		end
	end
	use = use_final;
      end % findStationSubset
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function stack = load_Yao_picks(filenamepath)
	% stack = duglNoise.load_Yao_picks(filenamepath);
	%
	% To load in saved dispersion picks from Yao's EGF processing.
	% Could be modified to some other code output...
	% What's important here is re-reading and re-absorbing to my "stack" object conventions

	%filenamepath = '/net/gauss/data1/dbowden/dugl_analysis/NoiseCC/matl_10hz_amp1/OUTPUT_SAC_covarianceweighted_GOOD/EGF_Analysis/EGFs/';

	sta_name = duglSet.Stations;
	sta_loc = duglSet.Locations;

	stack = [];
	icount = 0;
	Disp_files = dir([filenamepath,'/*Disp*']);
	for ii = 1:length(Disp_files)
			

		fid = fopen([filenamepath,Disp_files(ii).name]);
		tline1 = fgetl(fid);
		tline2 = fgetl(fid);
		C = textscan(fid,'%f %f %f %u','HeaderLines',0); % headerlines=0 because fgetl already read out twice
		fclose(fid);
	
		spl = regexp(Disp_files(ii).name,'\.|-|_','split');
		%CDisp.T.stackcn.1700-C2000_0.05-5s_TT.dat
		%GDisp.EGFcn.LHS-D4100_0.05-5s_ZZ.dat
		if(strcmp(spl{1},'CDisp'))
			sta1 = spl{4};
			sta2 = spl{5};
			cmp = spl{9};
		elseif(strcmp(spl{1},'GDisp'))
			sta1 = spl{3};
			sta2 = spl{4};
			cmp = spl{8};
		end
		stacmp = [sta1,'.HH',cmp(1),'-',sta2,'.HH',cmp(2)];

		if(ii>1)
			icount = length(stack.stacmp);
			ifile = icount+1; % default, add to end of stack
			for jj = 1:icount
				if(strcmp(stack.stacmp(jj),stacmp))
					ifile = jj;
				end
			end
		else
			ifile=1;
		end

	
		[loc1] = sscanf(tline1,'%f %f');
		stack.loc1(ifile,1) = loc1(1)-360;
		stack.loc1(ifile,2) = loc1(2);
		stack.loc1(ifile,3) = 0;
		[loc2] = sscanf(tline2,'%f %f');
		stack.loc2(ifile,1) = loc2(1)-360;
		stack.loc2(ifile,2) = loc2(2);
		stack.loc2(ifile,3) = 0;
		stack.stacmp{ifile,1} = stacmp;

		%disp(sprintf('%s    %s',stack.stacmp{ifile},stacmp))

		% Lookup and re-write 
		i1 = find(strcmp(sta_name,sta1));
		loc1 = sta_loc(i1,:);
		stack.loc1(ifile,:) = loc1;

		i2 = find(strcmp(sta_name,sta2));
		loc2 = sta_loc(i2,:);
		stack.loc2(ifile,:) = loc2;

		% Distance
		xdiff = duglProc.DistVicenty(loc1(2),loc1(1),...
					     loc1(2),loc2(1));
		ydiff = duglProc.DistVicenty(loc1(2),loc1(1),...
					     loc2(2),loc1(1));
		zdiff = loc1(3) - loc2(3);
		stack.dist(ifile) = sqrt(xdiff^2 + ydiff^2 + zdiff^2);
		% Azimuth and dip
		% Defined as direction from loc1 to loc2, ii to jj
		% Azimuth is degrees clockwise from north
		% Dip is degrees downward from horizontal
		%  (negative degrees for upward)
		az = azimuth(loc1(2),loc1(1),loc2(2),loc2(1));
		stack.az(ifile) = az;
		dip = atan2(zdiff,sqrt(xdiff^2 + ydiff^2))*180/pi;
		stack.dip(ifile) = dip;

		% "d1" is how Yao's code calculated distance
		%  I prefer only using lateral distance, and using Vicenty's for more accuracy
		d1 = deg2km(distance(loc1(2),loc1(1),loc2(2),loc2(1))).*1000;
		d1 = sqrt( zdiff^2 + d1^2 );

		dist_lateral = sqrt(xdiff^2 + ydiff^2);
		dist_ratio  = dist_lateral / d1;
		C{2} = C{2}.*dist_ratio;  % correct the velocities accordingly
		
		
		vec = [C{1} C{2}];
		vec(vec(:,2)==0,:) = [];
		if(strcmp(spl{1},'CDisp'))
			stack.phvel{ifile,1} = vec;
		elseif(strcmp(spl{1},'GDisp'))
			stack.grvel{ifile,1} = vec;
		end
	end
	sl = length(stack.stacmp);
	if(length(stack.phvel) < sl)
		stack.phvel{sl} = [];
	end
	sl = length(stack.stacmp);
	if(length(stack.grvel) < sl)
		stack.grvel{sl} = [];
	end


      end % load_Yao_picks
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   end

end
