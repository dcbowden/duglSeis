classdef duglProc

% Note! Station coords in text file different than on SAC files

   properties (Constant = true)
	options = {'applyResp','getPZ','zpk2cmplx','TauP','EstimateDelay'};
   end %constant properties

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   methods(Static)

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function s = bp(s,lowF,highF,n,p);
	% sout = bp(s,lowF,highF,n,p);
	%
	% Wrapper for duglProc.Filt
	% Syntax meant to mimic SAC command for filtering
	%
	% Includes a linear taper, which "Filt" does not
	% Includes labeling, which "Filt" does not
	if(nargin < 3)
		disp('ERROR! Usage:')
		disp('s = bp(s,low,high,n,p)')
		return
	elseif (nargin < 4)
		n=2;
		p=2;
	end
	
	num = max(size(s));

	for i = 1:num(1)
		nsamp = length(s(i).DATA1);
	        taperl = round(nsamp*.10);
	        taper = linspace(0,1,taperl)';
        	s(i).DATA1(1:taperl) = s(i).DATA1(1:taperl).*taper;
        	s(i).DATA1(nsamp-taperl+1:end) = s(i).DATA1(nsamp-taperl+1:end) .* flip(taper);
	
		sps = 1/s(i).DELTA;
		data = duglProc.Filt(s(i).DATA1,lowF, highF,n,p,sps);
		s(i).DATA1 = data;
		s(i).misc.bpL = lowF;
		s(i).misc.bpH = highF;
	end
      end % function bp
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function output = Filt(data,lowF,highF,n,p,sps)
	fnyq = sps/2;
	[B,A] = butter(p,[lowF highF]./fnyq);
	if n==1; output=filter(B,A,data);
	elseif n==2; output=filtfilt(B,A,data);
	else disp('ERROR! Number of passes can be only 1 or 2');
	end
      end % function Filt
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function sout = downSamp(s,ds_factor)
	%if(mod(length(s(1).DATA1),ds_factor) ~= 0)
	%	disp('ERROR, better to use integer divisor as factor!')
	%	return
	%end
	sout = s;
	for ii = 1:length(s)
		sout(ii).DELTA = sout(ii).DELTA.*ds_factor;

		% Matlab's downsample just cherry-picks every n'th sample
		%sout(ii).DATA1 = downsample(s(ii).DATA1,ds_factor);

		% Matlab's decimate does a proper job of filtering
		sout(ii).DATA1 = decimate(s(ii).DATA1,ds_factor);
		%  Note: python script used hamming window of order 31, matlab's default is 
		%   Chebyshev order 8
		
		sout(ii).NPTS = length(sout(ii).DATA1);
		sout(ii).E=round(length(sout(ii).DATA1)*sout(ii).DELTA)-1*sout(ii).DELTA; %recalc for precision
	end
      end % function downSamp
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %function sout = attachResp(s,varargin)
      %  % Instead of applying instrument response, just attach to the s object
      %end % function attachResp
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function sout = applyResp(s,varargin)
	return_type = 'VEL';
	for ii = 1:length(varargin)
		if(strcmp(varargin{ii},'ACC'))
			return_type = 'ACC';
		elseif(strcmp(varargin{ii},'DIS'))
			return_type = 'DIS';
		elseif(strcmp(varargin{ii},'VEL'))
			return_type = 'VEL';
		else
			disp(['Not understood: ',varargin{ii}])
		end
		disp(return_type)
	end
	
        sout = s; %Copy the object structure.
        for jj = 1:length(s)
		if(s(jj).RESP0 == 1)
			disp(['Already applied response to s(',num2str(jj),')'])
			disp('Skipping...')
			return
		end
			
		sout(jj).RESP0 = 1;

		% Remove trend
		s(jj).DATA1 = detrend(s(jj).DATA1);

		%% Remove mean (not needed if remove trend)
		%s(jj).DATA1 = s(jj).DATA1 - mean(s(jj).DATA1);
		
		% Taper
		nsamp = length(s(jj).DATA1);
		taperl = nsamp.*0.05;
		taper = hann(taperl*2);
		taper = taper(1:taperl);
        	s(jj).DATA1(1:taperl) = s(jj).DATA1(1:taperl).*taper;
        	s(jj).DATA1(nsamp-taperl+1:end) = s(jj).DATA1(nsamp-taperl+1:end) .* flip(taper);


		NFFT = 2^nextpow2(s(jj).NPTS);
		ff = 1/s(jj).DELTA/2 * linspace(0,1,NFFT/2+1);

		% Calculate transfer function "h"
		[zz,pp,const,pz_valid] = duglProc.getPZ(s(jj));
		switch(return_type)
			case 'ACC'
			   H_resp = duglProc.zpk2cmplx(ff,zz,pp,const,0).';
			case 'VEL'
			   H_resp = duglProc.zpk2cmplx(ff,zz,pp,const,1).';
			case 'DIS'
			   H_resp = duglProc.zpk2cmplx(ff,zz,pp,const,2).';
		end
		%h = h.*1e9;

		% Set up a filter
		p = 2;
		lowF = 7e-4;
		highF = 30;
		fnyq = 1/s(jj).DELTA/2;
		[B,A] = butter(p,[lowF highF]./fnyq);
		[H_filter,W] = freqz(B,A,NFFT,'whole');
		% Butter and Freqz make the same thing, just pole-zero vs freq-domain
		
	
		% FFT the data
		spect = fft(s(jj).DATA1,NFFT);


		%% Plot (for testing only)
		%figure
		%semilogx(ff,abs(spect(1:NFFT/2+1)) ./ max(abs(spect)))
		%hold on
		%semilogx(ff,abs(H_filter(1:NFFT/2+1)))
		%semilogx(ff,abs(H_resp)./mean(abs(H_resp)));
		
		% Apply filter
		spect = spect.*(abs(H_filter).^2);
			% spect = spect.*H_filter;
			% is a causal filter
			% this should give same results as "filter(B,A,data)"
			% 
			% To get anti-causal, zero phase, like "filtfilt":
			% sp = time-reverse[ sp.*H ] .*H
			%    = conj(sp).*conj(H).*H
			%    = sp.*H.*conj(H)
			%    = sp.*(abs(H).^2)

		spect(1:NFFT/2+1) = spect(1:NFFT/2+1).*H_resp(1:NFFT/2+1);
			% Note, the H_resp is defined for 1/2 spectra
			% while H_filter defined for whole... fix this for consistency?

		% FFT back, after re-adjusting for only 1/2
		spect(NFFT/2+2:end) = conj(spect(NFFT/2:-1:2)); 
		spect(1) = 0;
		spect(NFFT/2+1)=abs(spect(NFFT/2+1));
		d2 = real(ifft(spect,NFFT,'symmetric'));
		d2(s(jj).NPTS+1:end) = [];
			% (is this similar to fftshift?)
			% (the 'real' and 'imag' parts should be same as Hilbert)

		sout(jj).DATA1 = d2;
		%% To match SAC files, put in nm/s
		%% sout(jj).DATA1 = d2.*1e9;
		sout(jj).IDEP = 07;
		sout(jj).RESP0 = 1;
		sout(jj).FILENAME = [sout(jj).FILENAME,'.R'];
	end
      end %function applyResp
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function [zz,pp,const,pz_valid,varargout] = getPZ(s)
      % [zz,pp,const,pz_valid] = duglProc.getPZ(s)
      % [zz,pp,const,pz_valid,azimuth] = duglProc.getPZ(s)
	if(length(s) > 1)
		disp('Error! getPZ only set up to take 1 sac file at a time')
		return
	end
	sac_start = datenum(s.NZYEAR,0,0,s.NZHOUR,s.NZMIN,s.NZSEC)+s.NZJDAY;
	sac_end   = sac_start+round(s.NPTS*s.DELTA)/60/60/24;

	%disp(datestr(sac_start))
	%disp(datestr(sac_end))
	% List all possible PZ files for this station/component
	pz_files = dir([duglSet.RESP_DIR,'/*',s.KSTNM,'_',s.KCMPNM,'*']);
	for kk = 1:length(pz_files) % search for the correct time range
		pzfi = regexp(pz_files(kk).name,'__(?<y1>\d\d\d\d)\.(?<doy1>\d\d\d)\.(?<h1>\d\d)\.(?<m1>\d\d).+_(?<y2>\d\d\d\d)\.(?<doy2>\d\d\d)\.(?<h2>\d\d)\.(?<m2>\d\d)','names');
		pzfi.start = datenum(str2num(pzfi.y1),0,0,str2num(pzfi.h1),str2num(pzfi.m1),0) + str2num(pzfi.doy1);
		pzfi.end   = datenum(str2num(pzfi.y2),0,0,str2num(pzfi.h2),str2num(pzfi.m2),0) + str2num(pzfi.doy2);
		%disp(datestr(pzfi.start))
		%disp(datestr(pzfi.end))

		if(sac_start > pzfi.start && sac_end < pzfi.end)
			pz_valid = pzfi;
			this_pz_file = [duglSet.RESP_DIR,'/',pz_files(kk).name];
%disp(this_pz_file)
			pz_fid = fopen(this_pz_file);
    			%initialize some flags and various variables
    			buff = 0;
    			nzeros = 0;
    			npoles = 0;
    			zero_scan_flag = -1;
    			pole_scan_flag = -1;
    			pole_count_flag = 1;
    			zero_count_flag = 1;
    			zz=0;
    			pp=0;
			azimuth = 0;
    			% loop over the entire file
    			while (buff ~= -1)
    			    % read the next line in the file
    			    buff = fgets(pz_fid);
    			    % check to make sure it isn't the end of the file
    			    if (buff ~= -1)
    			        tmp = sscanf(buff, '%s" "%d');
    			        if length(strfind(buff,'ZEROS'))>0
    			            zero_scan_flag = 1;
    			            pole_scan_flag = 0;
    			            tmp = sscanf(buff, '%s %d');
    			            nzeros = tmp(6);
    			    	%disp(['Num Zeros: ',num2str(nzeros)]);
    			        elseif length(strfind(buff,'POLES'))>0
    			            pole_scan_flag = 1;
    			            zero_scan_flag = 0;
    			            tmp = sscanf(buff,'%s %d');
    			            npoles = tmp(6);
    			    	%disp(['Num Poles: ',num2str(npoles)]);
    			        elseif length(strfind(buff,'CONSTANT'))>0
    			            pole_scan_flag = 0;
    			            tmp = sscanf(buff, '%s %e');
    			            const = tmp(9);
    			    	%disp(['Constant: ',num2str(const)]);
				%* AZIMUTH     (deg) : 90.0 
				elseif length(strfind(buff,'AZIMUTH'))>0
				    spl = regexp(buff,' : ','split');
				    azimuth = str2num(spl{2});
    			        elseif (zero_scan_flag == 1)
    			            if (zero_count_flag < nzeros+1)
    			                tmp = sscanf(buff,'%f %f');
    			                zero_count_flag = zero_count_flag + 1;
    			                if (zero_count_flag > 1)
    			                    j=zero_count_flag-1;
    			                    zz(j) = tmp(1) + tmp(2)*i;
    			                end
    			            end
    			        elseif (pole_scan_flag == 1) 
    			            if (pole_count_flag < npoles+1)
    			                tmp = sscanf(buff, '%f %f');
    			                pole_count_flag = pole_count_flag + 1;
    			                if (pole_count_flag > 1) 
    			                    j=pole_count_flag-1;
    			                    pp(j) = tmp(1) + tmp(2)*i;
    			                end
    			            end
    			        end
    			    end
    			end
    			
    			% fill in missing poles and zeros
    			if (length(zz) < nzeros)
    			    for j=length(zz)+1:nzeros
    			        zz(j) = 0+0i;
    			    end
    			end
    			if (length(pp) < npoles)
    			    for j=length(pp)+1:npoles
    			        pp(j) = 0+0i;
    			    end
    			end
    			
			varargout{1} = azimuth;
			%disp(num2str(azimuth))
    			
    			fclose(pz_fid);
			return
		end %if pz in time-range
	end % foreach possible pz_file

	%disp('ERROR! No matching pz file found!')
	disp(sprintf('!No PZ file %s %s %s',datestr(sac_start),datestr(sac_end),s.FILENAME))
	zz = 0;
	pp = 0;
	const = 0;
	pz_valid = 0;
	varargout{1} = 90; % azimuth default
	return
      end %function getPZ
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function h = zpk2cmplx(f,z,p,const,wpow);
	% Inspired by Seizmo scripts by Garrett Euler
	% Default wpow = 0;  displacement
	% Change from Seizmo: for wpow = 1 == velocity:
	%	The constant extracted from SACPZ is A0*SENS
	%	(A0 is scaling for pole-zero transform)
	%	(SENS is volts <--> m/s)
	%	To get velocity transfer function (not displacement), 
	%	we divide h by i*2*pi*f, and same for A0, but NOT SENS.
	%
	%	Maybe there is a smarter way?
	z =cplxpair(z(:));
	p =cplxpair(p(:));

	nf = numel(f);
	nz = numel(z);
	np=numel(p)  ;
	% Calculate A0 ourselves, based on pole-zeros
	%   A0 scales the transfer function such that:
	%   h( f==1Hz ) = 1;
	f0 = 1; % Assuming for all our seismometers this is true
	w0=complex(0,2*pi*f0.');
	hn = abs(prod(w0(ones(nz,1),:) - z(:,ones(1,1)),1)...
		./prod(w0(ones(np,1),:) - p(:,ones(1,1)),1));
	A0 = 2*pi*f0/hn;
	SENS = const / A0;
	
	% for every frequncy, n:
	% w[n] = 2*pi*i*f  (omega)
	% h[n] = A0 * ((w[n]-z0) * (w[n]-z1) * ...) /
	%              ((w[n]-p0) * (w[n]-p1) * ...))
	w = complex(0,2*pi*f(:).');
	%h = prod(w(ones(nz,1),:) - z(:,ones(nf,1)),1)...
	h = A0.*prod(w(ones(nz,1),:) - z(:,ones(nf,1)),1)...
		./prod(w(ones(np,1),:) - p(:,ones(nf,1)),1);
	
	% This was coded for an example where output of P-Z file is units vel.
	% since #p = #z+1, prod(z)/prod(p) leaves one extra omega on top.
	% wpow=1 removes the 1 omega, returning response h to default units
	% wpow=2 divides by an extra omega --> integrate once to displacement
	% wpow=0 leaves *w --> derivative once to acc
	h = h./w.^wpow;
	if(wpow>0) 
		h(f==0)=0; 
	end
	%h = h./SENS.*1e9;
	h = h./SENS;
      end %function zpk2cmplx
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function [pwdB_sm,T_sm] = power_spect(s,varargin);
       % [pwdB_sm,T_sm] = duglProc.power_spect(s);
       % [pwdB_sm,T_sm] = duglProc.power_spect(s,H_resp);
       %
       % Option to use pre-calculated H_resp, acc normally
       % H_resp must be of the same NFFT as used here...

	% Assume the only varargin in pre-calculated H_resp
	if(length(varargin)>1)
		H_resp = varargin{1};
	end

	data = s(1).DATA1;
	if(s(1).RESP0==1)
		disp('This already has instrument response removed!')
		disp('We want to apply after, skipping.')
		return
	end
	data = detrend(data);

	nsecs_per_win = 800;
	npts_per_win = nsecs_per_win./s.DELTA;

	%window = hann(npts_per_win);
	window = tukeywin(npts_per_win,.2);
	[pw,FF] = pwelch(data,window,numel(window)*.75,[],1/s.DELTA,'onesided');


	% Apply instrument response
	if(exist('H_resp')==0)
		[zz,pp,const,pz_valid] = duglProc.getPZ(s);
		if(length(zz)==1 && const==0)
			disp('No instrument response, no power spectrum...')
			pwdB_sm = 0;
			T_sm = 0;
			return
		end
		H_resp = duglProc.zpk2cmplx(FF,zz,pp,const,0).';
	end
	%pw = pw.*h./1e9;
	pw = pw.*abs(H_resp.*conj(H_resp));

	% Convert to decibel
	pwdB = 10.*log10(pw);

	%%% Smooth more like McNamara
	%% By Period
	TT = 1./FF;
	T_oct_left(1) = 1/50;  % 50 Hz highest kept
	ii = 1;
	% nsecs_per_win/10 longest period kept

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% full octave width at 1/8 octave intervals 
	%  each step 2^(1/8)
	%while T_oct_left(ii)*2^(1/8)*2^(8/8) < nsecs_per_win/5
	%	ii = ii+1;
	%	T_oct_left(ii) = T_oct_left(ii-1)*2^(1/8);
	%end
	%T_oct_right = 2^(8/8).*T_oct_left;
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% half octave width at 1/8 octave intervals 
	while T_oct_left(ii)*2^(1/8)*2^(4/8) < nsecs_per_win/5
		ii = ii+1;
		T_oct_left(ii) = T_oct_left(ii-1)*2^(1/8);
	end
	T_oct_right = 2^(4/8).*T_oct_left;
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% half octave width at 1/16 octave intervals 
	%while T_oct_left(ii)*2^(1/16)*2^(4/8) < nsecs_per_win/5
	%	ii = ii+1;
	%	T_oct_left(ii) = T_oct_left(ii-1)*2^(1/16);
	%end
	%T_oct_right = 2^(4/8).*T_oct_left;
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	T_oct_center = sqrt(T_oct_left.*T_oct_right);
	pwdB_sm = zeros(length(T_oct_center),1);
	for ii = 1:length(T_oct_center)
		inRange = find(sum([TT>T_oct_left(ii)   TT<T_oct_right(ii)],2)==2);
		pwdB_sm(ii) = mean(pwdB(inRange));
	end
	T_sm = T_oct_center;


	%	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%	% Smooth once over
	%	clear T_oct_left T_oct_right
	%	T_oct_left(1) = 1/50;  % 50 Hz highest kept
	%	ii = 1;
	%	% full octave width at 1/8 octave intervals 
	%	%  each step 2^(1/8)
	%	while T_oct_left(ii)*2^(1/8)*2^(8/8) < nsecs_per_win/5
	%		ii = ii+1;
	%		T_oct_left(ii) = T_oct_left(ii-1)*2^(1/8);
	%	end
	%	T_oct_right = 2^(8/8).*T_oct_left;
	%	T_oct_center = sqrt(T_oct_left.*T_oct_right);
	%	pwdB_sm2 = zeros(length(T_oct_center),1);
	%	for ii = 1:length(T_oct_center)
	%		inRange = find(sum([T_sm'>T_oct_left(ii)   T_sm'<T_oct_right(ii)],2)==2);
	%		pwdB_sm2(ii) = mean(pwdB_sm(inRange));
	%	end
	%	T_sm2 = T_oct_center;





	%	if(length(varargout)==2)
	%		%% By Frequency
	%		F_oct_left(1) = 1/(nsecs_per_win/5);
	%		ii = 1;
	%		while F_oct_left(ii)*2^(1/8)*2 < 50  % 50Hz highest possible
	%			ii = ii+1;
	%			F_oct_left(ii) = F_oct_left(ii-1)*2^(1/8);
	%		end
	%		F_oct_right = 2.*F_oct_left;
	%		F_oct_center = sqrt(F_oct_left.*F_oct_right);
	%		pwdB_smF = zeros(length(F_oct_center),1);
	%		for ii = 1:length(F_oct_center)
	%			inRange = find(sum([FF>F_oct_left(ii)   FF<F_oct_right(ii)],2)==2);
	%			pwdB_smF(ii) = mean(pwdB(inRange));
	%		end
	%		F_sm = F_oct_center;

	%	end


      end %function 
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function taup = TauP(varargin);
	%     - TauP toolkit developed by:
	%        H. Philip Crotwell, Thomas J. Owens, Jeroen Ritsema
	%        Department of Geological Sciences
	%        University of South Carolina
	%        http://www.seis.sc.edu
	%        crotwell@seis.sc.edu
	%        citation: Crotwell, H. P., T. J. Owens, and J. Ritsema (1999). The TauP Toolkit: Flexible seismic travel-time and ray-path utilities, Seismological Research Letters 70, 154–160.
	% 
	%
	%	- Original matlab interface by Qin Li in 2002
	%	- Snippets from SEIZMO package by Garrett Euler and from Rob Porritt
	%
	%	taup = TauP(s);
	%		(where s is SAC object)
	%
	%	taup = TauP(source);
	%		(where source is object with:
	%			source.depth;
	%		        source.srclat;
	%			source.srclon; )
	%
	%	taup = duglProc.TauP(s,'P');
	%		(for a specific phase)
	%
	%	taup will be object with variables like:
	%	    taup(ii).modelname
	%	    taup(ii).event    (location)
	%	    taup(ii).station  (location)
	%	    taup(ii).depth    (event depth)
	%	    taup(ii).distance
	%	    taup(ii).mindistance
	%	    taup(ii).phase    (name)
	%	    taup(ii).puristphase
	%	    taup(ii).time
	%	    taup(ii).rayparameter
	%	    taup(ii).path     (contains many points of the ray's path)

	model = 'iasp91'; % 1D earth model, one example
	%center_of_array = [44.3488 -103.7570 ]; %switched lat and lon from duglSet?
	center_of_array = duglSet.center_of_array;
	stlat = center_of_array(2);
	stlon = center_of_array(1);

	passed_struct = varargin{1};
	if(length(varargin)>1)
		phases = varargin{2};
	else
		phases = 'ttbasic';
	end

	if( isfield(passed_struct,'EVLA') ); % seismogram passed
		s(1) = passed_struct(1);
		depth =  s(1).EVDP;
		srclat = s(1).EVLA;
		srclon = s(1).EVLO;
	elseif( isfield(passed_struct,'PreferredLatitude') ); % event passed
		ev = passed_struct;
         	srclat = ev.PreferredLatitude;
         	srclon = ev.PreferredLongitude;
         	depth  = ev.PreferredDepth;
	elseif( isfield(passed_struct,'srclat'));
		depth = passed_struct.depth;
		srclat = passed_struct.srclat;
		srclon = passed_struct.srclon;
	else
		display('ERROR! Could not understand TauP input')
		return 
	end
	% else: load other structure? FRAME?

	% try adding *.jar if no TauP class exists
	if(~exist('edu.sc.seis.TauP.MatTauP_Path','class'))
	    fs=filesep;
	    jars = dir([duglSet.TAUP_DIR,'/lib/*.jar']);
	    
	    for i=1:numel(jars)
	        if(~ismember([duglSet.TAUP_DIR,'/lib/',jars(i).name],javaclasspath))
	            javaaddpath([duglSet.TAUP_DIR,'/lib/',jars(i).name]);
	        end
	    end
	end
	sc=javaObject('edu.sc.seis.TauP.SphericalCoords');
        %deg=sc.distance(ev(1),ev(2),st(1),st(2));
	deg = sc.distance(srclat,srclon,stlat,stlon);
	az  = sc.azimuth(srclat,srclon,stlat,stlon);
	
	% create path object for velocity model
	tpobj=javaObject('edu.sc.seis.TauP.MatTauP_Path',model);
	% calculate paths for the specified depth, phase & distance
	tpobj.setSourceDepth(depth);
	tpobj.setPhaseNames(phases);
	tpobj.calculate(deg);
	tpobj.setEv(srclat,srclon);
	tpobj.setAz(az);
	tpobj.pathInterpolate;
	narr=tpobj.getNumArrivals;
	R2D=180/pi;
	
	% struct output
	taup(1:narr,1)=struct('modelname',[],'event',[],'station',[],'depth',[],...
	    'distance',[],'mindistance',[],'phase',[],'puristphase',[],...
	    'time',[],'rayparameter',[],'path',[],'azimuth',[],'incidentangle',[]);
	for ii=1:narr
	    % get phase info
	    arr=tpobj.getMatArrival(ii-1);
	    taup(ii).modelname=model;
	    taup(ii).event=[srclat srclon];
	    taup(ii).station=[stlat stlon];
	    taup(ii).depth=arr.getSourceDepth;
	    taup(ii).distance=arr.getDistDeg;
	    taup(ii).mindistance=arr.getModuloDistDeg;
	    taup(ii).phase=char(arr.getName);
	    taup(ii).puristphase=char(arr.getPuristName);
	    taup(ii).time=arr.getTime;
	    taup(ii).rayparameter=arr.getRayParam/R2D;
	
	    % extract path info
	    pts=arr.getMatPath;
	    taup(ii).path.depth=pts.getDepths;
	    taup(ii).path.distance=pts.getDistances;
	    taup(ii).path.time=pts.getTimes;
	    taup(ii).path.latitude=pts.getLats;
	    taup(ii).path.longitude=pts.getLons;

	    taup(ii).azimuth = az;

	    % incident angle.. could get from ray parameter also P = Rsin(theta)/V
	    inc_dep = taup(ii).path.depth(end)-taup(ii).path.depth(end-5);
	    R = 6371; 
	    inc_dis = (taup(ii).path.distance(end)-taup(ii).path.distance(end-5))*R*sind(1);
	    inc = -atan2(inc_dep,inc_dis)*180/pi;
	    taup(ii).incidentangle = inc;
	end
      end %function TauP
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function [lprop,tdif] = EstimateDelay(varargin)
	% [lprop,tdif] = EstimateDelay(az,inc)
	% [lprop,tdif] = EstimateDelay(az,inc,s)
	% [lprop,tdif] = EstimateDelay(az,inc,s,taup)
	% [lprop,tdif] = EstimateDelay(s,taup)

	default_v = 5;

	have_tau = 0;
	have_s = 0;
	nums = [];
	for ii = 1:length(varargin)
		if(isfield(varargin{ii},'path'))
			taup = varargin{ii};
			have_tau = 1;
		elseif(isfield(varargin{ii},'EVLA'))
			s = varargin{ii};
			have_s = 1;
		else
			nums = [nums,varargin{ii}];
		end
	end

	% Check that we have an az and inc angle to use
	if(length(nums) == 2 || length(nums)==3)
		az = nums(1);
		inc = nums(2);
	elseif(length(nums)>3)
		disp('ERROR! Can not understand input numbers')
		nums
		return 
	elseif(~have_tau && have_s) %User not defined, and no taup given
		%Run taup ourselves, assuming we want first P
		taup = duglProc.TauP(s,'P'); 
		[n,m] = min([taup.time]);
		inc = taup(m).incidentangle;
		az = s(1).AZ;
	elseif(have_tau)
		inc = [taup.incidentangle];
		az = [taup.azimuth];
	else
		disp('ERROR! Not enough input information given')
		return
	end

	if(length(nums) == 1)
		v = nums(1);
	elseif(length(nums)==3)
		v = nums(3);
	else
		v = default_v;
	end
	disp(['Using v = ',num2str(v)])

	% Decide which stations to calculate delays for
	%	fid = fopen([duglSet.MAP_DIR,'/stations.txt']);
	%	stations_temp = textscan(fid,'%s %f %f %f %f %f %f');
	%	fclose(fid);
	%	  % name lon lat depth(km) easting northing elvation(ft)
	%	stations = [stations_temp{2} stations_temp{3} stations_temp{4}];
	%	station_names = stations_temp{1};
	station_names = duglSet.Stations;
	stations = duglSet.Locations;

	% if s given, only calculate those delays 
	% (order of lprop should correspond to s)
	if(have_s)
		st2 = []; % temporary stations file
		stnm2 = [];
		for kk = 1:length(s)
		   for jj = 1:length(station_names)
			if(strmatch(s(kk).KSTNM,station_names(jj)))
				st2 = [st2; stations(jj,:)];
				stnm2 = [stnm2; {station_names{jj}}];
			end
		   end
		end
		stations = st2;
		station_names = stnm2;
	end

	% Calculate delay times relative to an arbitrary origin (center of array).
	%center_of_array = [-103.757047 44.348816 877.761077 ];
	center_of_array = duglSet.center_of_array;

	% Wavefronts traveling with some theta (x-y, 360deg from E), 
	%  and some phi up from x-y axis (90 - incident angle)

	theta = 90-az; %geologist defines azimuth as degrees clockwise from north
		       %math defines theta as degrees counter-clockwise from east
	phi = 90-inc; %seismologist defines incoming angle as dip down from horizontal
		      %math defines angle down from vertical
	%
	% We will define all station's angles alpha and beta similarly
	for ii = 1:length(station_names)
		% center_of_array( x , y )
		% stations(index, [ x,y,z ])

		%% Could also use raw feet, but would need to convert surface stations,
		%% Testing shows they agree to within ~0.2m
		%%
		%% DuglProc.DistVicenty( y1 , x1, y2, x2 )
		xdif = duglProc.DistVicenty(center_of_array(2),center_of_array(1),...
					    center_of_array(2),stations(ii,1));
		ydif = duglProc.DistVicenty(center_of_array(2),center_of_array(1),...
						stations(ii,2),center_of_array(1));
		%% correct to be +/-
		if(center_of_array(1)>stations(ii,1))
			xdif = -xdif;
		end
		if(center_of_array(2)>stations(ii,2))
			ydif = -ydif;
		end
		zdif = stations(ii,3) - center_of_array(3);
		pos(ii,:) = [xdif,ydif,zdif];
		% Could also use raw feet to calculate cartesian distances, 
		%  testing shows vicenty's algorithm agrees to within ~0.2m
		L = sqrt(xdif^2 + ydif^2 + zdif^2);
		alpha = 180/pi*atan2(ydif,xdif);
		beta = 180/pi*acos(zdif/L);
		%disp(sprintf('%s %f %f %f %f %f',station_names{ii},xdif,ydif,zdif,alpha,beta))
		lprop(ii,:) = L*cosd(alpha)*sind(beta).*cosd(theta).*sind(phi)...
				+ L*sind(alpha)*sind(beta).*sind(theta).*sind(phi)...
				+ L*cosd(beta).*cosd(phi);
		tdif(ii,:) = lprop(ii,:)./1000 ./ v;
	end
	%disp(sprintf('%f %f',theta,phi))
	%for ii = 1:length(station_names)
	%	disp(sprintf('%s %f',station_names{ii},lprop(ii)))
	%end
      end %EstimateDelay
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function [lprop,tdif] = distFromStation(s,name)
	%[lprop,tdif] = duglProc.distFromStation(s,name)
	ista = 0;
	for ii = 1:length(s)
		if(strcmp(s(ii).KSTNM,name))
			ista = ii;
		end
	end
	if(ista==0)
		disp('No matching station in "s"')
		return
	end

	%	% Decide which stations to calculate delays for
	%	fid = fopen([duglSet.MAP_DIR,'/stations.txt']);
	%	stations_temp = textscan(fid,'%s %f %f %f %f %f %f');
	%	fclose(fid);
	%	  % name lon lat depth(km) easting northing elvation(ft)
	%	stations = [stations_temp{2} stations_temp{3} stations_temp{4}];
	%	station_names = stations_temp{1};
	station_names = duglSet.Stations;
	stations = duglSet.Locations;

	% (order of lprop should correspond to s)
	st2 = []; % temporary stations file
	stnm2 = [];
	for kk = 1:length(s)
	   for jj = 1:length(station_names)
		if(strmatch(s(kk).KSTNM,station_names(jj)))
			st2 = [st2; stations(jj,:)];
			stnm2 = [stnm2; {station_names{jj}}];
		end
		if(strmatch(s(ista).KSTNM,station_names(jj)))
			st1 = stations(jj,:);
			stnm1 = station_names{jj};
		end
	   end
	end

	for ii = 1:length(st2)
		xdif = duglProc.DistVicenty(st1(2),st1(1),...
					    st1(2),st2(ii,1));
		ydif = duglProc.DistVicenty(st1(2),st1(1),...
					    st2(ii,2),st1(1));
		zdif = st2(ii,3) - st1(3);
		lprop(ii) = sqrt(xdif^2 + ydif^2 + zdif^2);
	end
	lprop( isnan(lprop) ) = 0;
	lprop = lprop';
	v = 5;
	tdif = lprop./1000 ./ v;
		
	


	
      end %distFromStation
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function dis = DistVicenty(lat1,lon1,lat2,lon2)
	%  dis = DistVicenty(lat1,lon1,lat2,lon2)
	%
	%  Original algorithm source:
	%  T. Vincenty, "Direct and Inverse Solutions of Geodesics on the Ellipsoid
	%  with Application of Nested Equations", Survey Review, vol. 23, no. 176,
	%  April 1975, pp 88-93
	%  Adapted by Michael Kleder, 2004.
	%
	if abs(lat1)>90 | abs(lat2)>90
	    error('Input latitudes must be between -90 and 90 degrees, inclusive.')
	end
	% Supply WGS84 earth ellipsoid axis lengths in meters:
	a = 6378137; % definitionally
	b = 6356752.31424518; % computed from WGS84 earth flattening coefficient definition
	% convert inputs in degrees to radians:
	lat1 = lat1 * 0.0174532925199433;
	lon1 = lon1 * 0.0174532925199433;
	lat2 = lat2 * 0.0174532925199433;
	lon2 = lon2 * 0.0174532925199433;
	% correct for errors at exact poles by adjusting 0.6 millimeters:
	if abs(pi/2-abs(lat1)) < 1e-10;
	    lat1 = sign(lat1)*(pi/2-(1e-10));
	end
	if abs(pi/2-abs(lat2)) < 1e-10;
	    lat2 = sign(lat2)*(pi/2-(1e-10));
	end
	f = (a-b)/a;
	U1 = atan((1-f)*tan(lat1));
	U2 = atan((1-f)*tan(lat2));
	lon1 = mod(lon1,2*pi);
	lon2 = mod(lon2,2*pi);
	L = abs(lon2-lon1);
	if L > pi
	    L = 2*pi - L;
	end
	lambda = L;
	lambdaold = 0;
	itercount = 0;
	while ~itercount | abs(lambda-lambdaold) > 1e-12  % force at least one execution
	    itercount = itercount+1;
	    if itercount > 50
	        warning('Points are essentially antipodal. Precision may be reduced slightly.');
	        lambda = pi;
	        break
	    end
	    lambdaold = lambda;
	    sinsigma = sqrt((cos(U2)*sin(lambda))^2+(cos(U1)*...
	        sin(U2)-sin(U1)*cos(U2)*cos(lambda))^2);
	    cossigma = sin(U1)*sin(U2)+cos(U1)*cos(U2)*cos(lambda);
	    sigma = atan2(sinsigma,cossigma);
	    alpha = asin(cos(U1)*cos(U2)*sin(lambda)/sin(sigma));
	    cos2sigmam = cos(sigma)-2*sin(U1)*sin(U2)/cos(alpha)^2;
	    C = f/16*cos(alpha)^2*(4+f*(4-3*cos(alpha)^2));
	    lambda = L+(1-C)*f*sin(alpha)*(sigma+C*sin(sigma)*...
	        (cos2sigmam+C*cos(sigma)*(-1+2*cos2sigmam^2)));
	    % correct for convergence failure in the case of essentially antipodal points
	    if lambda > pi
	        warning('Points are essentially antipodal. Precision may be reduced slightly.');
	        lambda = pi;
	        break
	    end
	end
	u2 = cos(alpha)^2*(a^2-b^2)/b^2;
	A = 1+u2/16384*(4096+u2*(-768+u2*(320-175*u2)));
	B = u2/1024*(256+u2*(-128+u2*(74-47*u2)));
	deltasigma = B*sin(sigma)*(cos2sigmam+B/4*(cos(sigma)*(-1+2*cos2sigmam^2)...
	    -B/6*cos2sigmam*(-3+4*sin(sigma)^2)*(-3+4*cos2sigmam^2)));
	dis = b*A*(sigma-deltasigma);
      end % DistanceVicenty
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function [L,Q,T] = rotateLQT(Z,E,N,i,baz);
      %  i = taup(1).incidentangle;
      %  baz = taup(1).azimuth+180;
      %  [L,Q,T] = duglProc.rotateLQT(Z,E,N,i,baz);
	M = [cosd(i)  -sind(i)*sind(baz)   -sind(i)*cosd(baz);
		sind(i)   cosd(i)*sind(baz)   cosd(i)*cosd(baz);
		0         -cosd(baz)          sind(baz)];
	L = Z; Q = E; T = N; % just to copy other header info
	out = M*[Z.DATA1,E.DATA1,N.DATA1]';
	L.DATA1 = out(1,:);
	Q.DATA1 = out(2,:);
	T.DATA1 = out(3,:);

	L.KCMPNM = 'HHL';
	rep = strfind(L.FILENAME,'HHZ');
	if(length(rep)>0)
		L.FILENAME(rep:rep+2) = 'HHL';
	end

	Q.KCMPNM = 'HHQ';
	rep = strfind(Q.FILENAME,'HHE');
	if(length(rep)>0)
		Q.FILENAME(rep:rep+2) = 'HHQ';
	end

	T.KCMPNM = 'HHT';
	rep = strfind(T.FILENAME,'HHN');
	if(length(rep)>0)
		T.FILENAME(rep:rep+2) = 'HHT';
	end
	
      end %rotate LQT
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function [T,R] = rotateTR(N,E,varargin)

	if(length(varargin)>0)
		baz = varargin{1}
	else
		baz = N.BAZ;
	end

	if(length(varargin)>1)
		disp('Too many arguments! Not understood')
		return
	end

	T=N;R=E; % just to copy over other header info

	M = [ -cosd(baz) +sind(baz);
		-sind(baz) -cosd(baz)];

	out = M*[E.DATA1,N.DATA1]';
	T.DATA1 = out(1,:)';
	R.DATA1 = out(2,:)';
	

	channel_first2 = N.KCMPNM(1:2);

	T.KCMPNM = [channel_first2,'T'];
	rep = strfind(T.FILENAME,[channel_first2,'N']);
	if(length(rep)>0)
		T.FILENAME(rep:rep+2) = [channel_first2,'T'];
	end

	R.KCMPNM = [channel_first2,'R'];
	rep = strfind(R.FILENAME,[channel_first2,'E']);
	if(length(rep)>0)
		R.FILENAME(rep:rep+2) = [channel_first2,'R'];
	end
      end %rotate TR
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function [N,E] = rotateNE(s1,s2,varargin)
	%[N,E] = duglFetch.rotateNE(s1,s2, azimuth)
	% Note that to pull BH1,BH2 back to BHN,BHE
	%  we need to rotate *backwards* by -azimuth
	%
	% Full example:
	% tr(1) = irisFetch.Traces('IU','RSSD','00','BH1',st,en);
	% tr(2) = irisFetch.Traces('IU','RSSD','00','BH2',st,en);
	% tr(3) = irisFetch.Traces('IU','RSSD','00','BHZ',st,en);
	% s = iris2sac(tr);
	% ch = irisFetch.Channels('channel','IU','RSSD','00','BH1');
	% pull_st = datenum(st);
	% pull_en = datenum(en);
	% for ii = 1:length(ch)
	% 	this_st = datenum(ch(ii).StartDate);
	% 	this_en = datenum(ch(ii).EndDate);
	% 	if( (pull_st > this_st) && (pull_en < this_en) )
	% 		use_epoch = ii;
	% 	end
	% end
	% az = ch(use_epoch).Azimuth;
	% Z2 = s(3);
	% [N2,E2]  = duglProc.rotateNE(s(1),s(2),-az);
	%
	% Check that you're doing it right by comparing to:
	% https://service.iris.edu/irisws/rotation/docs/1/builder/ 

	if(length(varargin)>0)
		az = varargin{1};
	else
		az = s1.AZ;
	end

	if(length(varargin)>1)
		disp('Too many arguments! Not understood')
		return
	end

	N=s1;E=s2; % just to copy over other header info

	M = [ cosd(az) sind(az);
		-sind(az) cosd(az)];

	out = M*[s1.DATA1,s2.DATA1]';
	N.DATA1 = out(1,:)';
	E.DATA1 = out(2,:)';
	

	channel_first2 = s1.KCMPNM(1:2);

	N.KCMPNM = [channel_first2,'N'];
	N.KINST = [channel_first2,'N'];
	rep = strfind(N.FILENAME,[channel_first2,'N']);
	if(length(rep)>0)
		N.FILENAME(rep:rep+2) = [channel_first2,'N'];
	end

	E.KCMPNM = [channel_first2,'E'];
	E.KINST = [channel_first2,'E'];
	rep = strfind(E.FILENAME,[channel_first2,'E']);
	if(length(rep)>0)
		E.FILENAME(rep:rep+2) = [channel_first2,'E'];
	end
      end %rotate NE
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    end % methods static
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %methods(Static, Access=protected)      
   %end % methods static protected
end % classdef
