classdef duglPlot

   properties (Constant = true)
	options = {'p1,p2,power_spect'};
	default_depth_coloring        = 1;
	default_auto_scale_traces     = 0;
	default_clear_before_plotting = 1;
	default_scaling_factor        = 10;
	default_points                = 0;
	default_tscale		      = 0;  % 0 is 'auto'
	default_xlab		      = 'Time'; % auto
	default_spect_use_freq	      = 0;
	default_dist_def              = 'abs';
	default_cut1		      = 0;
	default_cut2                  = 0;
   end %constant properties


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   methods(Static)

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function hh = p1(varargin)
	%s_base = assignin('base','name','s')
	
	args_in = [];
	for jj = 1:length(varargin)
		args_in{jj} = varargin{jj};
	end
	[s,flags] = duglPlot.parse_args(args_in);
	[i_longest,xlab,tscale]  = duglPlot.find_longest(s);
	if(flags.tscale == 0)
		% Do nothing, let "find_longest" determine best
		%  (ie: seconds vs minutes vs hours)
	else
		% overwrite tscale with whatever was passed in flags
		tscale = flags.tscale;
		xlab = flags.xlab;
	end
	if(flags.clear_before_plotting)
		clf
	end
	if(flags.depth_coloring)
		cs = duglPlot.getDepthColoring(s,'plot');
	else
		cs = lines;
	end
	n = length(s);
	for jj = 1:length(s)
		hs(jj) = subplot(n,1,jj);
		tt = linspace(s(jj).B,s(jj).E,length(s(jj).DATA1)).*tscale;
		hh = plot(tt,s(jj).DATA1,'color',cs(jj,:));
		hold on
		%if(flags.depth_coloring)
		%	hh = plot(tt,s(jj).DATA1,'color',cs(jj,:));
		%	%hold on  %held in getDepthColoring
		%else
		%	hh = plot(tt,s(jj).DATA1,'k');
		%	hold on
		%end
		set(hh,'DisplayName',s(jj).FILENAME)
		set(hh,'LineWidth',.1)
		mytitle = s(jj).FILENAME;
		mytitle(mytitle=='_') = ' ';
		title(mytitle)
		hold on
	end
	xlabel(xlab)
	linkaxes(hs,'x');
	%if(flags.auto_scale_traces==0)
	%	linkaxes(hs,'xy');
	%end
      end %function p1
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function h = p2(varargin)
	% figure_handle = duglPlot.p2(s)
	% figure_handle = duglPlot.p2(s,'noclear','scale_trace','no_coloring')
	% figure_handle = duglPlot.p2(s,'tscale','min');
	args_in = [];
	for jj = 1:length(varargin)
		args_in{jj} = varargin{jj};
	end
	[s,flags] = duglPlot.parse_args(args_in);
	[i_longest,xlab,tscale]  = duglPlot.find_longest(s);
	[t0,tend,t1] = duglPlot.find_tindex(s);
	if(flags.tscale == 0)
		% Do nothing, let "find_longest" determine best
		%  (ie: seconds vs minutes vs hours)
	else
		% overwrite tscale with whatever was passed in flags
		tscale = flags.tscale;
		xlab = flags.xlab;
	end

	if(flags.clear_before_plotting)
		disp('Clearing')
		clf
	end
	n = length(s);
	if(length(s)>25)
		disp('Warning! Too many plots overlain')
	end
	if(flags.depth_coloring)
		cs = duglPlot.getDepthColoring(s,'plot');
	else
		cs = lines;
	end
	if(isfield(flags,'tdif'))
		tshift = flags.tdif(:,1)*tscale;
	else
		tshift = zeros(length(s),1);
	end
	legs = [];
	

	for jj = 1:length(s)
%disp(sprintf('t1      = %.8f',t1(jj)))
		tt = linspace(s(jj).B+t1(jj),s(jj).E+t1(jj),length(s(jj).DATA1)).*tscale;
		%tt = linspace(s(jj).B+t1(jj),s(jj).E+t1(jj),length(s(jj).DATA1));
%disp(sprintf('tt(1)   = %.8f',tt(1)))
%disp(sprintf('s(jj).B = %.8f',s(jj).B))
%disp(sprintf('tt(end) = %.8f',tt(end)))
%disp(sprintf('s(jj).E = %.8f',s(jj).E))
%disp(sprintf('tshift  = %.8f',tshift(jj)))
%save tt.mat tt
		if(flags.auto_scale_traces)
			data = s(jj).DATA1;
			data = data - mean(data);
			thisMax = max(abs(data));
			data = data./thisMax;
		else
			data = s(jj).DATA1;
		end
		if(flags.depth_coloring)
			hh = plot(tt+tshift(jj),data,'color',cs(jj,:));
			%hold on  %held in getDepthColoring
		else
			%hh = plot(tt+tshift(jj),data);
			hh = plot(tt+tshift(jj),data,'color',cs(jj,:));
			hold on
			mytitle = s(jj).FILENAME;
			mytitle(mytitle=='_') = ' ';
			legs = [legs,{mytitle}];
		end
		set(hh,'DisplayName',s(jj).FILENAME)
		set(hh,'LineWidth',.1)
		xlabel(xlab)
	end
	%xlim([t0-t0 tend-t0].*tscale)
	if(~flags.depth_coloring)
		legend(legs)
	end
	il = i_longest;
	mytitle = [num2str(s(il).NZYEAR),'-',num2str(s(il).NZJDAY),' ',...
			num2str(s(il).NZHOUR,'%02.f'),':',num2str(s(il).NZMIN,'%02.f'),':',...
			num2str(s(il).NZSEC,'%02.f'),':',num2str(s(il).NZMSEC,'%02.f')];
	if(isfield(s(1),'misc'))
		if(isfield(s(1).misc,'bpL'))
			mytitle = [mytitle,'  BP: ',sprintf('%0.3f %0.3f',s(1).misc.bpL,s(1).misc.bpH)];
		end
	end
	if(isfield(flags,'tdif'))
		mytitle = [mytitle,'   Aligned by phase: ',flags.taup(1).phase];
	end
	title(mytitle);
	h = gcf;
      end %function
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function h = prs(lprop,taup,varargin)
	% h = prs(lprop,taup,s, ...)
	%
	% "plot record section"
	% will plot each trace, sorted by distance from event along wavefront
	args_in = [];
	for jj = 1:length(varargin)
		args_in{jj} = varargin{jj};
	end
	[s,flags] = duglPlot.parse_args(args_in);
	[i_longest,xlab,tscale]  = duglPlot.find_longest(s);
	if(flags.clear_before_plotting)
		clf
	end
	if(flags.depth_coloring)
		cs = duglPlot.getDepthColoring(s,'plot');
	end
	if(size(lprop,2)~=1 || length(taup)>1)
		disp('WARNING! size(lprop)>1 || length(taup)>1!')
		disp('Using first taup and lprop only')
		lprop = lprop(:,1);
		taup = taup(1);
	end
	n = length(s);
	maxd = max((lprop))-min((lprop));
	scaling = maxd/flags.scaling_factor;
	for ii = 1:n
		tt = linspace(s(ii).B,s(ii).E,length(s(ii).DATA1)).*tscale;
		dis = lprop(ii); %meters;
		norm = max(s(ii).DATA1)/scaling;
		if(flags.depth_coloring)
			hh = plot(tt,s(ii).DATA1/norm+dis,'color',cs(ii,:));
			hold on  %held in getDepthColoring
		else
			hh = plot(tt,s(ii).DATA1/norm+dis);
			hold on
		end
		set(hh,'DisplayName',s(ii).FILENAME);
		set(hh,'LineWidth',.1);
	end
	h = gcf;
	xlabel(xlab);
	ylabel('Distance [m] along wavefront')
	il = i_longest;
	mytitle = [num2str(s(il).NZYEAR),'-',num2str(s(il).NZJDAY),' ',...
			num2str(s(il).NZHOUR),':',num2str(s(il).NZMIN),':',...
			num2str(s(il).NZSEC),':',num2str(s(il).NZMSEC),'  '...
			'Distance by phase: ',taup.phase];
	if(isfield(s(1),'misc'))
		if(isfield(s(1).misc,'bpL'))
			mytitle = [mytitle,'  BP: ',sprintf('%0.3f %0.3f',s(1).misc.bpL,s(1).misc.bpH)];
		end
	end
	title(mytitle)
	grid 
	grid minor
      end % prs
		
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function nprs(varargin)
	% duglPlot.nprs(s)
	% duglPlot.nprs(s,'scaling_factor',30,'lateral')
	% duglPlot.nprs(s,'scaling_factor',30,'lateral','cut',-5,5)


	args_in = [];
	for jj = 1:length(varargin)
		args_in{jj} = varargin{jj};
	end
	[s,flags] = duglPlot.parse_args(args_in);




	maxshift = 2048;
	maxshift = 1024;
	maxshift = 512;
	delta = .02;

	%args_in = [];
	%flags.clear_before_plotting = 1;
	%flags.scaling_factor = 12;
	%jj = 1;
	%while(jj <= length(varargin))
	%	switch varargin{jj}
	%	   case 'noclear'
	%		flags.clear_before_plotting = 0;
	%	   case 'scaling_factor'
	%		jj = jj+1;
	%		flags.scaling_factor = varargin{jj};
	%	end
	%	jj = jj+1;
	%end
	if(flags.clear_before_plotting)
		clf
		cc = [0 0 0];
	else
		cc = [lines(1)];
	end

	if(strcmp(flags.dist_def,'abs'))
		dists = [s.DIST];
	elseif(strcmp(flags.dist_def,'lateral'))
		dists = [s.USER3];
	elseif(strcmp(flags.dist_def,'vertical'))
		dists = [s.USER4];
	end


	% Plot the traces
	maxd = max(abs(dists));
	scaling=maxd/flags.scaling_factor;
	for ii = 1:length(s)
		tt = linspace(s(ii).B,s(ii).E,length(s(ii).DATA1));
		norm = max(s(ii).DATA1)/scaling;
		hh = plot(tt,s(ii).DATA1/norm-dists(ii),'color',[cc(1,:)]);
		set(hh,'DisplayName',s(ii).FILENAME)
		set(hh,'LineWidth',.1)
		ax = axis;
		hold on
	end


	% Plot location of envelope max or abs max
	%	for ii = 1:length(s)
	%		tt = linspace(s(ii).B,s(ii).E,length(s(ii).DATA1));
	%		pos_half = [round(length(s(ii).DATA1)/2)+2:length(s(ii).DATA1)];
	%		%env = abs(hilbert(s(ii).DATA1(pos_half)));
	%		env = (s(ii).DATA1(pos_half));
	%		[n,m] = max(env);
	%		if(m==1)
	%			continue
	%		end
	%		norm = max(s(ii).DATA1)/scaling;
	%		plot(tt(pos_half(m)),s(ii).DATA1(pos_half(m))/norm-dists(ii),'r*');
	%		hold on
	%		hh = plot(tt,s(ii).DATA1/norm-dists(ii),'color',[lines(1)]);
	%	end
		

	% Set x-limits 
	xlabel('Time (s)')
	if(flags.cut1~=0 || flags.cut2~=0)
		x1 = flags.cut1;	
		x2 = flags.cut2;
		xlim([x1 x2])
	else
		x1 = ax(1);
		x2 = ax(2);
	end

	% Add labels
	for ii = 1:length(s)
		hh = text(x2*1.05,-dists(ii),s(ii).FILENAME);
	end

	% Y-axis
	if(strcmp(flags.dist_def,'abs'))
		ylabel('Distance (km)')
	elseif(strcmp(flags.dist_def,'lateral'))
		ylabel('Lateral Distance (km)')
	elseif(strcmp(flags.dist_def,'vertical'))
		ylabel('Vertical Distance (km)')
	end

	siz = get(gca,'Position');
	siz(3) = siz(3)*.8;
	set(gca,'Position',[siz])
      end %nprs
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function plot_PSD(pwdB2,FF2,varargin)
% DCB copy in
      end % plot_PSD
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function obj_h = idTrace;
	obj_h = gco;
	type_h = get(obj_h,'type');
	ax = axis;
	if type_h ~= 'line' % Wrong selection
	    disp('The selected object is not a line. Please retry.')
	    x_out=[];y_out=[];z_out=[];
	    return
	end
	xx = get(obj_h,'XData');
	yy = get(obj_h,'YData');
	nn = get(obj_h,'DisplayName');
	display(nn)
      end %function

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function fig_handle = power_spect(varargin)
	args_in = [];
	for jj = 1:length(varargin)
		args_in{jj} = varargin{jj};
	end
	[s,flags] = duglPlot.parse_args(args_in);
	if(flags.clear_before_plotting)
		disp('Clearing')
		clf
	end

	if(flags.depth_coloring)
		cs = duglPlot.getDepthColoring(s,'semilogx');
	end

	for jj = 1:length(s)

		NFFT = 2^nextpow2(s(jj).NPTS);
		%% McNamara:
		%% truncate time series down to nearest 2^N
		%% NFFT = s(jj).NPTS/2+1;

		if(s(jj).RESP0~=1)
			disp('Applying instrument response')
			s(jj) = duglProc.applyResp(s(jj),'ACC');
			%s(jj).DATA1 = s(jj).DATA1.*1.589459;
		end
		disp('Calculating spectra')
		data = s(jj).DATA1;
		%data = data./1e9;
		[pw,FF] = pwelch(data,[],[],NFFT,1/s(jj).DELTA,'onesided');
		%% 'power'?
		

%		% Should be possible to apply instrument response on-the-fly
%		% directly to output of pwelch. But pwelch is only real, and 
%		% instrument response correction doesn't come out right...
%		% Ignoring this for now
%		% Update: See duglProc.power_spect for this
%		if(s(jj).RESP0~=1)
%			% Calculate and apply our own instrument response
%			[zz,pp,const,pz_valid] = duglProc.getPZ(s(jj));
%			h = duglProc.zpk2cmplx(FF,zz,pp,const,1).';
%			%pw = pw.*h./1e9;
%			pw = pw.*h;
%		end
		pwdB = 10.*log10(pw);
		%semilogx(FF,pwdB,'color',[set3(ii,:)])
		if(flags.spect_use_freq == 1)
			if(flags.depth_coloring)
				semilogx(FF,pwdB,'color',cs(jj,:));
			else
				semilogx(FF,pwdB)
			end
		else
			if(flags.depth_coloring)
				semilogx(1./FF,pwdB,'color',cs(jj,:));
			else
				semilogx(1./FF,pwdB)
			end
		end
			
		hold on
		leg_list{jj} = s(jj).FILENAME;
	end

	max(FF)
	duglPlot.updateLegend(leg_list,flags);
	load noisemodel.mat
	% nlnm( period, model )   low noise 
	% nhnm( period, model )   high noise
	if(flags.spect_use_freq == 1)
		hh1 = semilogx(1./nlnm(:,1),nlnm(:,2),'c');
		hh2 = semilogx(1./nhnm(:,1),nhnm(:,2),'c');
		xlabel('Frequency')
	else
		hh1 = semilogx(nlnm(:,1),nlnm(:,2),'c');
		hh2 = semilogx(nhnm(:,1),nhnm(:,2),'c');
		xlabel('Period')
	end
	set(hh1,'DisplayName','NoiseModel');
	set(hh2,'DisplayName','NoiseModel');
	title('Power Spectra')
	ylabel('Power [10*log10( units^2 / Hz](dB)')
	% "units" could be [m], [m/s], or [m/s^2]
      end %function

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function [pw,FF] = power_spect2(varargin)
	% get PSD given input SAC file
	args_in = [];
	for jj = 1:length(varargin)
		args_in{jj} = varargin{jj};
	end
	[s,flags] = duglPlot.parse_args(args_in);
	if(length(s)>1)
		disp('ERROR in power_spect2: Only 1 trace allowed per plot')
		return
	end
	data = s(1).DATA1;
	if(s(1).RESP0==1)
		disp('This already has instrument response removed!')
		disp('We want to apply after, skipping.')
		return
	end

	nsecs_per_win = 800;
	npts_per_win = nsecs_per_win./s.DELTA;

	%window = hann(npts_per_win);
	window = tukeywin(npts_per_win,.2);
	[pw,FF] = pwelch(data,window,numel(window)*.75,[],1/s.DELTA,'onesided');

	%	% Do it the McNamara way
	%	n_psd = floor( length(s.DATA1)/npts_per_win*4 )
	%	for ii = 1:n_psd
	%		i_start = 1 + (ii-1)*npts_per_win*1/4;
	%		i_end = ii*length(window)*1/4;
	%		data2 = data(i_start:i_end);
	%		data2 = detrend(data2);
	%		taperl = .1*length(data2);
	%		taper = flip((cos([1:taperl]./taperl*pi)+1)./2)';
	%		
	%		data2(1:taperl) = data2(1:taperl).*taper;
	%		data2(length(data2)-taperl+1:end) = data2(length(data2)-taperl+1:end).*flip(taper);
	%		NFFT = length(data2)/2+1;
	%		sp = fft(data2,NFFT)
	%		% Not finished.... 3/2/2017. One sided FFT?
	%	end
		
	

	% Apply instrument response
	[zz,pp,const,pz_valid] = duglProc.getPZ(s);
	H_resp = duglProc.zpk2cmplx(FF,zz,pp,const,0).';
	%pw = pw.*h./1e9;
	pw = pw.*abs(H_resp.*conj(H_resp));

	% Convert to decibel
	pwdB = 10.*log10(pw);


	semilogx(1./FF,pwdB)
	hold on

	% Smooth more like McNamara
	TT = 1./FF;
	T_oct_left(1) = 1/50;  % 40 Hz highest kept
	ii = 1;
	% nsecs_per_win/10 longest period kept
	while T_oct_left(ii)*2^(1/8)*2 < nsecs_per_win/5
		ii = ii+1;
		T_oct_left(ii) = T_oct_left(ii-1)*2^(1/8);
	end
	T_oct_right = 2.*T_oct_left;
	T_oct_center = sqrt(T_oct_left.*T_oct_right);
	pwdB_sm = zeros(length(T_oct_center),1);
	for ii = 1:length(T_oct_center)
		inRange = find(sum([TT>T_oct_left(ii)   TT<T_oct_right(ii)],2)==2);
		pwdB_sm(ii) = mean(pwdB(inRange));
	end
	T_sm = T_oct_center;

	semilogx(T_sm,pwdB_sm)
	hold on
	%semilogx(1./T_sm,pwdB_sm)
	
	load noisemodel.mat
	% nlnm( period, model )   low noise 
	% nhnm( period, model )   high noise
	hh1 = semilogx(nlnm(:,1),nlnm(:,2),'c');
	hh2 = semilogx(nhnm(:,1),nhnm(:,2),'c');
	xlabel('Period')
	xlim([1e-2 1e2])

      end
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function fig_handle = power_spect_color(varargin)
	args_in = [];
	for jj = 1:length(varargin)
		args_in{jj} = varargin{jj};
	end
	[s,flags] = duglPlot.parse_args(args_in);
	if(length(s)>1)
		disp('ERROR in power_spect_color: Only 1 trace allowed per plot')
		return
	end
	if(flags.clear_before_plotting)
		disp('Clearing')
		clf
	end

	if(s(1).RESP0~=1)
		disp('Applying instrument response')
		s(1) = duglProc.applyResp(s(1));
		%s(1).DATA1 = s(1).DATA1.*1.589459;
	end
	data = s(1).DATA1;
	%data = data./1e9;

	nsecs_per_win = 800;
	npts_per_win = nsecs_per_win./s.DELTA;
	freqRange = [1/nsecs_per_win/10  20];

	window = hann(npts_per_win);
	%n_psd = ceil( length(s.DATA1)/npts_per_win );
	n_psd = floor( length(s.DATA1)/npts_per_win*4); % for 75% overlap, 4 windows per complete
	%spect = zeros( npts_per_win/2+1 , n_psd);

	[pw,FF] = pwelch(data,window,numel(window)*.75,[],1/s(jj).DELTA,'onesided');
	


	%fff = zeros( npts_per_win/2+1 , n_psd);
	for ii = 1:n_psd
		i_start = 1 + (ii-1)*npts_per_win*1/4;
		%i_end = ii*length(window)*1/4;   % DCB Modifying 7.11.2018, window only 1/4 long, not 1/4 overlap!
		i_end = i_start + npts_per_win; 
		%disp(sprintf('%u %u',i_start,i_end))

		[pw,FF] = pwelch(data(i_start:i_end),window,[],numel(window),1/s(jj).DELTA,'onesided');
		%[pw,FF] = pwelch(data,window,0,numel(window),1/s(jj).DELTA);
		% Convert to m from m/s.
		%p(2:end) = p(2:end) ./ ((2*pi*f(2:end)) .^ 2);
		%p = sqrt(p);

		pwdB = 10.*log10(pw);
	
		spect(:,ii) = pwdB;
		%fff(:,ii) = FF;
	end


	% Apply instrument response

	% Smooth more like McNamara
	TT = 1./FF;
	perRange = 1./freqRange;
	T_oct_left(1) = 1/40;
	ii = 1;
	while T_oct_left(ii)*2^(1/8)*2 < nsecs_per_win/10
		ii = ii+1;
		T_oct_left(ii) = T_oct_left(ii-1)*2^(1/8);
	end
	T_oct_right = 2.*T_oct_left;
	T_oct_center = sqrt(T_oct_left.*T_oct_right);
	for ii = 1:length(T_oct_center)
		inRange = find(sum([TT>T_oct_left(ii)   TT<T_oct_right(ii)],2)==2);
		spect_sm(ii) = mean(spect(inRange));
	end

	% Resample
%		ds = 80;
%		tmp2 = zeros( length(FF) , 1);
%		spect3 = zeros( floor(length(FF)/ds), n_psd);
%		F_log2 = logspace(log10(FF(2)),log10(FF(end)),length(FF));
%		F_log3 = logspace(log10(FF(2)),log10(FF(end)),floor(length(FF)./ds));
%		for ii = 1:n_psd
%			tmp2 = smooth(interp1(FF,spect(:,ii),F_log2),ds);
%			spect3(:,ii) = interp1(F_log2,tmp2,F_log3);
%		end


	%p_bins = logspace(log10min(spect(:))
	%p_bins = linspace(-225,-50,200);
    	%h = hist(spect3',p_bins);

	


    	% Set up grid.
    	[x,y] = meshgrid(F_log3,p_bins);

	% Normalize histogram in each frequency bin and get percentiles.
	ll = []; ul = []; ml = [];
	for mm=1:size(h,2)
	  h(:,mm) = h(:,mm)/sum(h(:,mm));
	  ll(mm) = prctile(spect3(mm,:),5);
	  ul(mm) = prctile(spect3(mm,:),95);
	  ml(mm) = mean(spect3(mm,:));
	end
	
	semilogx(F_log3,ll,'k-')
	hold on
	p1 = pcolor(x,y,h);
	shading interp
	hold on
	semilogx(F_log3,ll,'k-','linewidth',1.5)
	semilogx(F_log3,ul,'k-','linewidth',1.5)

	load noisemodel.mat
	% nlnm( period, model )   low noise 
	% nhnm( period, model )   high noise
	%if(flags.spect_use_freq == 1)
	hh1 = semilogx(1./nlnm(:,1),nlnm(:,2),'LineStyle','--','Color',[0.4 0.4 0.4],'Linewidth',2);
	hh2 = semilogx(1./nhnm(:,1),nhnm(:,2),'LineStyle','--','Color',[0.4 0.4 0.4],'Linewidth',2);
	xlabel('Frequency')
	xlim([10^-2 25])
	ylim([-225 -50])
    	set(gcf,'renderer','painters');
    	cax = caxis;
    	caxis([0 cax(2)./4])
    	%cmap = parula(1000);
    	%cmap(1,:) = [1 1 1];
    	%colormap(cmap);
	ylabel('Power [10*log_1_0(m^2/s^4/Hz](dB)')
      end  % power_spect_color


      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function spectro_db(sta,cha,varargin)
	%duglPlot.spectro_db('ROSS','HHZ');
	%duglPlot.spectro_db('ROSS','HHZ','2015 190 00:00:00','+10d');
	
	%% Default t_start and t_end

	% This is one year period where most stations were complete
	%	t_start = duglFetch.makeDateInfo('2015 182 00:00:00');
	%	t_end   = duglFetch.makeDateInfo('2016 182 00:00:00');

	% This is fullest extent of 24-station coverage
	% (currently pre-processed in database)
	t_start = duglFetch.makeDateInfo('2015/05/20 00:00:00');
	t_end   = duglFetch.makeDateInfo('2016/09/24 00:00:00');

	if(length(varargin)==2)
		t_start = duglFetch.makeDateInfo(varargin{1});
		if(length(strfind(varargin{2},'+'))>0)
			t_end   = duglFetch.plusDateInfo(t_start,varargin{2});
		else
			t_end = duglFetch.makeDateInfo(varargin{2});
		end
	elseif(length(varargin)==0)
		disp('Using default, full time range')
	else
		disp('Problem with input, not plotting!')
		return
	end
	disp(sprintf('%s-%s  from %s to %s',sta,cha,t_start.datstr,t_end.datstr));

	load([duglSet.SPECT_DIR,'/check_times.mat'])
	%load([duglSet.SPECT_DIR,'/',sta,'_',cha,'_pwsp_filled.mat'])
	%load([duglSet.SPECT_DIR,'/',sta,'_',cha,'_pwsp.mat'])
	load([duglSet.SPECT_DIR,'/',sta,'_',cha,'_P.mat'])

	testRange = [P.times>t_start.datnum P.times<t_end.datnum];
	inRange = find(sum(testRange,2)==2);


	clf
	hold on

	T = P.T;
	times = P.times(inRange);
	pwdB = P.pwdB(inRange,:);

	pcolor(times,T,pwdB')
	shading flat
	hold on
	set(gca,'YScale','log')
	datetick	
	ylim([10^-1 10^2])
	xlim([times(1) times(end)])
	caxis([-160 -110])
	title([P.sta,' ',P.cha])
	ylabel('Period')
	ax = axis;
	chsize(14)
	h = colorbar;
	ylabel(h,'Power [10*log_1_0(m^2/s^4/Hz](dB)');
	
	% compare times in database to check_times
	for ii = 1:length(check_times)
	        if(check_times(ii)>t_start.datnum && check_times(ii)<t_end.datnum)
			in_catalog = find( abs(times - check_times(ii))<0.000001 );
			if(length(in_catalog) == 0) % Then not in catalog
				%disp(datestr(check_times(ii),duglSet.DATE_FORMATTER))
				plot([check_times(ii) check_times(ii)],[ax(3:4)],'w')
			end
		end
	end
	
	% Pull set times from *.mat file with saved spectral info
	%
	%for ii = 1:length(times)
	%  if(pwdB(ii,end)>-120)
	%	plot([times(ii) times(ii)],[ax(3:4)],'w')
	%  end
	%end
	
      end  % spectro_db
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function [varargout] =  psd_db(sta,cha,varargin)
	%duglPlot.psd_db('ROSS','HHZ');
	%duglPlot.psd_db('ROSS','HHZ','2015 190 00:00:00','+10d');
	% [ff, ml, ll, ul, times] = duglPlot.psd_db('ROSS','HHZ','clean_glitch','freq')
	
	%% Default t_start and t_end

	% This is one year period where most stations were complete
	%	t_start = duglFetch.makeDateInfo('2015 182 00:00:00');
	%	t_end   = duglFetch.makeDateInfo('2016 182 00:00:00');

	% This is fullest extent of 24-station coverage
	% (currently pre-processed in database)
	t_start = duglFetch.makeDateInfo('2015/05/20 00:00:00');
	t_end   = duglFetch.makeDateInfo('2016/09/24 00:00:00');

	% Other defaults
	clean_glitch = 0;
	plot_period = 0;


	nV = length(varargin);
	jj = 1;
	while(jj<=nV)
		if(strcmp(varargin{jj}(1:3),'201'))
			t_start = duglFetch.makeDateInfo(varargin{jj});
			jj = jj + 1;
			if(length(strfind(varargin{jj},'+'))>0)
				t_end   = duglFetch.plusDateInfo(t_start,varargin{jj});
			else
				t_end   = duglFetch.makeDateInfo(varargin{jj});
			end
		elseif(strcmp(varargin{jj},'freq'))
			plot_period = 0;
		elseif(strcmp(varargin{jj},'per'))
			plot_period = 1;
		elseif(strcmp(varargin{jj},'clean_glitch'))
			clean_glitch = 1;
		else
			disp('Problem with input, not plotting!')
			return
		end
		jj = jj + 1;
	end

	disp(sprintf('%s-%s  from %s to %s',sta,cha,t_start.datstr,t_end.datstr));

	load([duglSet.SPECT_DIR,'/check_times.mat'])
	load([duglSet.SPECT_DIR,'/',sta,'_',cha,'_P.mat'])

	testRange = [P.times>t_start.datnum P.times<t_end.datnum];
	inRange = find(sum(testRange,2)==2);

	T = P.T;
	pwdB = P.pwdB(inRange,:);
	times = P.times(inRange);

	% Clean for periods with some kind of recenter or glitch
	% (shows up as long period anomaly)
	bad = [];
	for jj = 1:size(pwdB,1)
	  if(pwdB(jj,end)>-120)
		bad = [bad; jj];
	  end
	end
	disp(sprintf('Found %.2f%% with station problems / glitch',length(bad)/size(P.pwdB,1)*100))
	if(length(bad)>0 && clean_glitch)
		disp(sprintf('Removing %.2f%% because of glitch',length(bad)/size(P.pwdB,1)*100))
		pwdB(bad,:) = [];
		times(bad) = [];
	else
		disp('Not removing anything');
	end

	p_bins = -200:1:-50;

	% Plot 1 thing in semilog to initialize axis
	semilogx(1,-150,'k-')
	hold on

	if( plot_period )
		xout = T;

		h = hist(pwdB,p_bins);

		[x,y] = meshgrid(T,p_bins);
		% Normalize histogram in each frequency bin and get percentiles.
		ll = []; ul = []; ml = [];
		for mm=1:size(h,2)
		  h(:,mm) = h(:,mm)/sum(h(:,mm));
		  ll(mm) = prctile(P.pwdB(:,mm),5);
		  ul(mm) = prctile(P.pwdB(:,mm),95);
		  ml(mm) = mean(P.pwdB(:,mm));
		end
		x_halfstep = mean(diff(log10(P.T)))*0.5;
		xx = 10.^(log10(x)-x_halfstep);
		yy = y - mean(diff(p_bins))*0.5;
		p1 = pcolor(xx,yy,h);

		%semilogx(P.T,ll,'k-','linewidth',1.5)
		%semilogx(P.T,ul,'k-','linewidth',1.5)

		load noisemodel.mat
		hh1 = semilogx(nlnm(:,1),nlnm(:,2),'c');
		hh2 = semilogx(nhnm(:,1),nhnm(:,2),'c');
		xlabel('Period')
		xlim([5e-2 1e2])

		for lp = [.1 .2 .3 .4 .5 .6 .7 .8 .9 1 2 3 4 5 6 7 8 9 10 20 30 40 50 60 70 80 90 100]
			plot([lp lp],[-200 -50],'k-','linewidth',.1)
		end
		for lp = [ -200: 10 : -80]
			plot([5e-2 10^2],[lp lp],'k-','linewidth',.1)
		end
	else

		ff = flip(1./T);
		xout = ff;
		pwdB_ff = fliplr(pwdB);

		h = hist(pwdB_ff,p_bins);
		
		[x,y] = meshgrid(ff,p_bins);
		% Normalize histogram in each frequency bin and get percentiles.
		ll = []; ul = []; ml = [];
		for mm=1:size(h,2)
		  h(:,mm) = h(:,mm)/sum(h(:,mm));
		  this_ff = pwdB_ff(:,mm);
		  this_ff( isinf(this_ff) ) = [];
		  ll(mm) = prctile(this_ff,5);
		  ul(mm) = prctile(this_ff,95);
		  ml(mm) =    mean(this_ff);
		end
		x_halfstep = mean(diff(log10(ff)))*0.5;
		xx = 10.^(log10(x)-x_halfstep);
		yy = y - mean(diff(p_bins))*0.5;
		p1 = pcolor(xx,yy,h);

		%semilogx(ff,ll,'k-','linewidth',1.5)
		%semilogx(ff,ul,'k-','linewidth',1.5)
		semilogx(ff,ml,'k-','linewidth',1.5)

		load noisemodel.mat
		hh1 = semilogx(1./nlnm(:,1),nlnm(:,2),'c');
		hh2 = semilogx(1./nhnm(:,1),nhnm(:,2),'c');
		xlabel('Freq')
		%xlim([5e-2 1e2])
		xlim([1e-2 20])
		

		for lf = [.01 .02 .03 .04 .05 .06 .07 .08 .09 .1 .2 .3 .4 .5 .6 .7 .8 .9 1 2 3 4 5 6 7 8 9 10 20]
			plot([lf lf],[-200 -50],'k-','linewidth',.1)
		end
		for lp = [ -200: 10 : -80]
			plot([1e-2 20],[lp lp],'k-','linewidth',.1)
		end
	end
	shading flat
	ylim([-200 -80])
	%colorbar
	title([P.sta,' ',P.cha])
	ylabel('Power [10*log_1_0(m^2/s^4/Hz](dB)')
	caxis([0 .15])
	chsize(14)

	varargout{1} = xout;
	varargout{2} = ml;
	varargout{3} = ll;
	varargout{4} = ul;
	varargout{5} = times;


      end  % psd_db
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function fig_handle = single_spectro(varargin)
	% duglPlot.single_spectro(s)
	% 
	% hardcoded to divide data into 500 windows
	% Thanks to Tanner for help with spectrogram func.

	args_in = [];
	for jj = 1:length(varargin)
		args_in{jj} = varargin{jj};
	end
	[s,flags] = duglPlot.parse_args(args_in);
	if(flags.clear_before_plotting)
		disp('Clearing')
		clf
	end

	if(flags.depth_coloring)
		cs = duglPlot.getDepthColoring(s,'semilogx');
	end
	if(length(s) > 1)
		disp('ERROR, can only do one spectrogram at a time')
		disp('ex:  >> duglPlot.single_spectro(s(1));')
		return
	end

	NFFT = 2^nextpow2(s.NPTS);
	if(s.RESP0~=1)
		disp('Applying instrument response')
		s = duglProc.applyResp(s);
	end
	data = s.DATA1./1e9;

	nwins = 500;
	npts_per_win = s.NPTS/nwins;
	NFFT = min(2^14, 2^nextpow2(npts_per_win));

	% Set up hann window with specified duration.
	window = hann(npts_per_win);

	% Get psd spectrogram, freq array, t array.
	% We use 50% overlap between time segments.
	[sp,f,t,p] = spectrogram(data,window,numel(window)/2,NFFT,1./s.DELTA);
	% sp --> spectro
	% f --> freq
	% t --> times sampled
	% p --> spectro but power units

	% remove zero frequency
	f(1) = []; sp(1,:) = []; p(1,:) = [];
	pwdB = 10.*log10(p);
	
	
	% Plotting
	ax(1) = subplot(4,1,1:3);  % 75% for spectogram
	hold on
	if(flags.spect_use_freq == 1)
		pcolor(t,f,pwdB);
		ylim([1./s.DELTA/NFFT 1./s.DELTA./2])
	else
		pcolor(t,1./f,pwdB);
		ylim([.1 100])
	end
	shading flat
	colorbar
	%caxis([-255 -120])
	xlim([0 s.NPTS*s.DELTA])
	ylabel('Freq (Hz)')
	
	ax(2) = subplot(4,1,4);
	tt = linspace( 0, s.NPTS*s.DELTA, s.NPTS);
	plot(tt,s.DATA1);
	xlabel('Seconds')
	xlim([0 s.NPTS*s.DELTA])
	
	% set subplot positions to line up
	% make a colorbar, lock in the position, then delete
	colorbar
	drawnow
	pos = get(gca,'Position');
	colorbar('delete');
	drawnow
	set(gca,'Position',pos);
      end %single_spectro
	
	
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function Map3D(varargin)
	% duglPlot.Map3D
	% duglPlot.Map3D('noclear', 'nocoloring', 'points','full')
	% 
	% 'points' adds rasterized points to map instead of just colors
	% 	   easier to see topography, but much larger file

	args_in = [];
	for jj = 1:length(varargin)
		args_in{jj} = varargin{jj};
	end
	[s,flags] = duglPlot.parse_args(args_in);
	flags.points		    = duglPlot.default_points;
	flags.fullmap		    = 0;  % Will use a zoomed-in range by default
	flags.surf 		    = 1; % include surface mesh
	for ii = 1:length(flags.other)
	    switch flags.other{ii}
	      case 'points'
		flags.points = 1;
	      case 'full'
		flags.fullmap = 1;
	      case 'nosurf'
		flags.surf = 0;
	      otherwise 
		disp(['Still not recognized: ',flags.other{ii}])
	    end
	end
	
	if(flags.clear_before_plotting)
		disp('Clearing')
		clf
	end

	% Precalculated topography
	%    M2 zoomed in a little, M3 zoomed in alot
	if(flags.fullmap==0)  % default
		load([duglSet.MAP_DIR,'/M3.mat']);
		M = M3;
	else
		load([duglSet.MAP_DIR,'/M2.mat']);
		M = M2;
		%load([duglSet.MAP_DIR,'/M3.mat']);
		%M = M3;
	end


	% Load stations file
	fid = fopen([duglSet.MAP_DIR,'/stations.txt']);
	stations_temp = textscan(fid,'%s %f %f %f %f %f %f');
	fclose(fid);
	% name lon lat depth(km) easting northing elvation(ft)
	stations = [stations_temp{2} stations_temp{3} stations_temp{4}];
	station_names = stations_temp{1};

	% Always remove RSSD
	stations(end,:) = []
	station_names(end) = [];
	

	if(flags.fullmap==0) % default
		% For now, removing "Far-away" sites.
		skip = [];
		for ii = 1:length(station_names)
			if(strcmp(station_names(ii),'TPK') ||...
			   strcmp(station_names(ii),'SHL') ||...
			   strcmp(station_names(ii),'DEAD') ||...
			   strcmp(station_names(ii),'RSSD'))
				skip = [skip,ii];
			end
		end
		stations(skip,:) = [];
		station_names(skip) = [];
	end
	
	
	fid = fopen([duglSet.MAP_DIR,'/OtherTools/shafts.txt']);
	shafts_temp = textscan(fid,'%s %f %f %f %f %f %f');
	fclose(fid);
	shafts = [shafts_temp{2} shafts_temp{3} shafts_temp{4}];
	shaft_names = shafts_temp{1};
	Yates_shaft_line.X = [shafts(1,1) shafts(1,1)];
	Yates_shaft_line.Y = [shafts(1,2) shafts(1,2)];
	Yates_shaft_line.Z = [shafts(1,3) 0];
	Ross_shaft_line.X = [shafts(2,1) shafts(2,1)];
	Ross_shaft_line.Y = [shafts(2,2) shafts(2,2)];
	Ross_shaft_line.Z = [shafts(2,3) 0];

	if(flags.depth_coloring)
		for ii = 1:length(station_names)
			stmp(ii).KSTNM = station_names{ii};
		end
		cs = duglPlot.getDepthColoring(stmp,'no_legend');
	end


	%Mf = flipud(M);
	[XXX,YYY] = meshgrid(XX,YY);
	if(flags.surf)
		surf(XXX,YYY,M)
		hold on
	end

	if(flags.points)
		nn = size(M,1)*size(M,2);
		XXXX = reshape(XXX,[1 nn]);
		YYYY = reshape(YYY,[1 nn]);
		MMMM = reshape(M,[1 nn]);
		scatter3(XXXX(1:4:nn),YYYY(1:4:nn),MMMM(1:4:nn),4,'k.')
		hold on
	end
	hold on
	colormap(demcmap(M));
	shading flat
	view(6.5,20)

	size(stations);
	size(cs);
	if(flags.depth_coloring)
		scatter3(stations(:,1),stations(:,2),stations(:,3),18,cs,'filled')
		%scatter3(stations(:,1),stations(:,2),stations(:,3),18,'k','markertype','o','filled')
	else
		scatter3(stations(:,1),stations(:,2),stations(:,3),18,'ro','filled')
	end
	text(stations(:,1),stations(:,2),stations(:,3),station_names)
	%plot3(Ross_shaft_line.X,  Ross_shaft_line.Y, Ross_shaft_line.Z,'r')
	%plot3(Yates_shaft_line.X,  Yates_shaft_line.Y, Yates_shaft_line.Z,'r')

	daspect([1 1 100000])


	% Compass:
	compass = [-103.74 44.35 1000];
	plot3([compass(1) compass(1)+.002],[compass(2) compass(2)],[compass(3) compass(3)],'k','linewidth',2)
	text(compass(1)+.002,compass(2),compass(3),'E')
	plot3([compass(1) compass(1)],[compass(2) compass(2)+.002],[compass(3) compass(3)],'k','linewidth',2)
	text(compass(1),compass(2)+.002,compass(3),'N')
	plot3([compass(1) compass(1)],[compass(2) compass(2)],[compass(3) compass(3)+200],'k','linewidth',2)
	text(compass(1),compass(2),compass(3)+200,'Z')
	
      end %function
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function ForceSec
	ha = findall(gcf,'type', 'axes');
	%if(length(ha) > 1)
	%	tgca = get(ha(1));
	%else
	%	tgca = gca;
	%end
	mul = 1;
	if(strcmp(ha(1).XLabel.String,'Time [min]'))
		mul = 60;
		set(ha(1),'Tag','min')
	elseif(strcmp(ha(1).XLabel.String,'Time [hour]'))
		mul = 60*60;
		set(ha(1),'Tag','hour')
	elseif(strcmp(ha(1).XLabel.String,'Time [day]'))
		mul = 60*60*24;
		set(ha(1),'Tag','day')
	elseif(strcmp(ha(1).XLabel.String,'Time [sec]'))
		set(ha(1),'Tag','sec')
		display('WARNING! Already in Seconds')
	else
		display('Could not determine current x-axis')
	end
	xticks = ha(1).XTickLabel;
	for jj = 1:length(xticks)
		ha(1).XTickLabel{jj} = num2str( str2num(xticks{jj})*mul );
	end
	ha(1).XLabel.String = 'Time [sec]';
	if(length(ha) > 1)
		for ii = 2:length(ha)
			for jj = 1:length(xticks)
				ha(ii).XTickLabel{jj} = num2str( str2num(xticks{jj})*mul );
				str2num(xticks{jj})*mul
			end
		end
	end
			
			
      end % function ForceSec

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function ResetTicks
	switch get(gca,'Tag')
	   case 'min'
		xlabel('Time [min]')
	   case 'hour'
		xlabel('Time [hour]')
	   case 'day'
		xlabel('Time [day]')
	   case 'sec'
		xlabel('Time [sec]')
	end
	tgca = gca;
	tgca.XTickLabelMode = 'auto';
      end % function ResetTicks

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function PlotTauP(taup);
	[s,flags] = duglPlot.parse_args([]);
	if(flags.clear_before_plotting)
		disp('Clearing')
		clf
	else
		figure
	end
	ax = gca;


	% plot grid
	[cx,cy]=duglPlot.circle(6871,360);
	plot(ax,cx,cy,'color',[0.2 0.2 0.2],'linewidth',1,...
	    'tag','plot_tauppath_grid');
	
	set(ax,'position',[0.025 0.05 0.95 0.9]); 
	 hold on
	
	[cx1,cy1]=duglPlot.circle(6871,180);
	[cx2,cy2]=duglPlot.circle(6771,180);
	plot(ax,[cx1; cx2],[cy1; cy2],'color',[0.2 0.2 0.2],'linewidth',1,...
	    'tag','plot_tauppath_grid');
	[cx1,cy1]=duglPlot.circle(6871,36);
	[cx2,cy2]=duglPlot.circle(6671,36);
	plot(ax,[cx1; cx2],[cy1; cy2],'color',[0.2 0.2 0.2],'linewidth',1,...
	    'tag','plot_tauppath_grid');
	[cx,cy]=duglPlot.circle(6371,4);
	plot(ax,cx(1:2:3),cy(1:2:3),'color',[0.2 0.2 0.2],'linewidth',1,...
	    'tag','plot_tauppath_grid');
	plot(ax,cx(2:2:4),cy(2:2:4),'color',[0.2 0.2 0.2],'linewidth',1,...
	    'tag','plot_tauppath_grid');
	
	% plot major discontinuities
	for i=[6371 3480 1220]
	    [cx,cy]=duglPlot.circle(i,360);
	    plot(ax,cx,cy,'color','k','linewidth',2,...
	        'tag','plot_tauppath_major_discon');
	end
	
	% plot minor discontinuities
	for i=[5961 5711 3780]
	    [cx,cy]=duglPlot.circle(i,360);
	    plot(ax,cx,cy,'color',[0.5 0.5 0.5],'linewidth',1,...
	        'tag','plot_tauppath_minor_discon');
	end
	axis(ax,'off','equal');
		
	cc = parula(length(taup));
	% loop over phases
	for i=1:length(taup)
	    evangle = -taup(1).distance;
	    % plot phase path
	    cx=(6371-taup(i).path.depth)...
	        .*sin((evangle+taup(i).path.distance)/180*pi);
	    cy=(6371-taup(i).path.depth)...
	        .*cos((evangle+taup(i).path.distance)/180*pi);
	    %plot(ax,cx,cy,...
	    %    'displayname',taup(i).phase,'tag','tauppath_raypath');
	    if(i==1)
	      plot(ax,cx,cy,...
	        'displayname',taup(i).phase,'tag','tauppath_raypath',...
		'color',[cc(i,:)],'linewidth',2);
	    else
	      plot(ax,cx,cy,...
	        'displayname',taup(i).phase,'tag','tauppath_raypath',...
		'color',[cc(i,:)]);
	    end
	end
	    lh=legend(ax,[ flipud(findobj(ax,'tag','tauppath_raypath'))'],...
	        get([ flipud(findobj(ax,'tag','tauppath_raypath'))'],...
	        'displayname'),'location','northwestoutside');
	    set(lh,'color','none','fontsize',6,'interpreter','none');
	
      end %function PlotTauP
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function TimeTauP(taup)
	% Adds indicators on time-series (assumes current figure)
	
	tgca = gca;
	% taup.time is in seconds
	if(strcmp(tgca.XLabel.String,'Time [min]'))
		tscale = 60;
	elseif(strcmp(tgca.XLabel.String,'Time [hour]'))
		tscale = 60*60;
	elseif(strcmp(tgca.XLabel.String,'Time [day]'))
		tscale = 60*60*24;
	elseif(strcmp(tgca.XLabel.String,'Time [sec]'))
		tscale = 1;
	else
		disp('Warning! Unrecognized time interval on plot')
		disp('  Assuming x-axis is in seconds')
		tscale = 1;
	end

	ax = axis;
	for kk = 1:length(taup)
		arrival_time = taup(kk).time/tscale;
		hh = plot([arrival_time arrival_time],[ax(3:4)],'k--');
		set(hh,'DisplayName',taup(kk).phase);
		set(hh,'LineWidth',.1);
		ht = text(arrival_time,ax(4).*.9,taup(kk).phase),
	end
		
	
	
      end %function TimeTauP
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function [s,flags] = parse_args(args_in)
	nV = length(args_in);
	
	% List of flags
	flags.depth_coloring        = duglPlot.default_depth_coloring;
	flags.auto_scale_traces     = duglPlot.default_auto_scale_traces;
	flags.clear_before_plotting = duglPlot.default_clear_before_plotting;
	flags.scaling_factor        = duglPlot.default_scaling_factor;
	flags.tscale 		    = duglPlot.default_tscale;
	flags.spect_use_freq	    = duglPlot.default_spect_use_freq;
	flags.dist_def              = duglPlot.default_dist_def;
	flags.cut1		    = duglPlot.default_cut1;
	flags.cut2		    = duglPlot.default_cut2;
	flags.other		    = [];

	jj = 1;
	% First arg, if anything, must be SAC data
	if( nV>0 && isfield(args_in{1},'KSTNM'))
	     s = args_in{1};
	     jj = 2;
	else
	     if(evalin('base',['exist(''s'')']))
   	     	s = evalin('base','s');
	     	display('Using workspace variable "s" for data')
	     else
		display('No sac data found')
		s = [];
	     end
	end

	% More?
	if( nV>0 )
	     while jj<=nV
		% Coloring?
   		if(strcmp(args_in{jj},'depth_coloring'))
			jj = jj+1;
			flags.depth_coloring = args_in{jj}
		elseif(strcmp(args_in{jj},'no_coloring'))
			flags.depth_coloring = 0;
		elseif(strcmp(args_in{jj},'yes_coloring'))
			flags.depth_coloring = 1;

		% amp Scaling?
		elseif(strcmp(args_in{jj},'scale_trace'))
			flags.auto_scale_traces = 1;
		% time scaling
		elseif(strcmp(args_in{jj},'tscale'))
			jj = jj+1;
			switch args_in{jj}
			  case 'sec'
				flags.tscale = 1;
				flags.xlab = 'Time [sec]';
			  case 'min'
				flags.tscale = 1/60;
				flags.xlab = 'Time [min]';
			  case 'hour'
				flags.tscale = 1/60/60;
				flags.xlab = 'Time [hour]';
			end
		% Clear?
		elseif(strcmp(args_in{jj},'clear_before_plotting'))
			jj = jj+1;
			flags.clear_before_plotting = args_in{jj};
		elseif(strcmp(args_in{jj},'noclear') || strcmp(args_in{jj},'noclf'))
			flags.clear_before_plotting = 0;

		%scaling for record section?
		elseif(strcmp(args_in{jj},'scaling_factor'))
			jj = jj+1;
			flags.scaling_factor = args_in{jj};

		%distance definition for record section?
		elseif(strcmp(args_in{jj},'lateral'))
			flags.dist_def = 'lateral';
		elseif(strcmp(args_in{jj},'vertical'))
			flags.dist_def = 'vertical';
		elseif(strcmp(args_in{jj},'absolute'))
			flags.dist_def = 'abs';

		%preset x-axis cut-time?
		elseif(strcmp(args_in{jj},'cut'))
			jj = jj+1;
			flags.cut1 = args_in{jj};
			jj = jj+1;
			flags.cut2 = args_in{jj};

		% shift calculated?
		elseif(strcmp(args_in{jj},'timeshift'))
			jj = jj+1;
			flags.tdif = args_in{jj};
			jj = jj+1;
			flags.taup = args_in{jj};
		elseif(strcmp(args_in{jj},'period'))
			flags.spect_use_freq = 0;
		elseif(strcmp(args_in{jj},'freq') || strcmp(args_in{jj},'frequency'))
			flags.spect_use_freq = 1;
		else
			disp(sprintf('Not understood flag: %s',args_in{jj}))
			flags.other = [flags.other,{args_in{jj}}];
		end
		jj = jj+1;
	    end
	end
      end %function parse_args
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function [i_longest,xlab,tscale] = find_longest(s)
      % finds the longest of traces in object "s", and returns
      % so we know whether appropriate to plot units of minute/second/hour
	[n i_longest] = max([s.NPTS]);
	tt = linspace(s(i_longest).B,s(i_longest).E,length(s(i_longest).DATA1));
	tscale = 1;
	xlab = 'Time [sec]';
	if(tt(end)>1e3)
		tscale = 1/60;
		xlab = 'Time [min]';
		tt = tt./60;
		if(tt(end)>400)
			xlab = 'Time [hours]';
			tscale = 1/60/60;
		end
	end
      end %function
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function [t0,tend,t1] = find_tindex(s)
		if(isnan(s(1).NZJDAY))
			t1 = [s.B];
			t2 = [s.E];
			t0 = min(t1);
			tend = max(t2);
			%t1 = t1-t0;
			t1 = zeros(size(t2));
		else
			t1 =    [s.NZJDAY].*24*60*60 + ...
				[s.NZHOUR].*60*60 + ...
				[s.NZMIN] .*60 + ...
				[s.NZSEC] + ...
				[s.NZMSEC] ./ 1000;
			t2 = t1+[s.E];

			t0 = min(t1);
			tend = max(t2);
			t1 = t1-t0;
		end
      end
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function cs = getDepthColoring(s,plottype);
	 c0 = [0 0 0];
	 %c300 = [203    24    29]./256; % red
	 %c800 = [251   106    74]./256;
	 %c1700 = [106    81   163 ]./256; %purple
	 %c2000 = [158   154   200 ]./256;
	 %c4100 = [33   113   181 ]./256; %blue
	 %c4850 = [ 35   139    69]./256; %green


	%% Previous colors
	% c300 = [ 212 0 0   ] ./256; %red 
	% c800 = [ 255 85 0  ] ./256; %orange
	% c1700 =[ 255 0 255 ] ./256; %pink
	% c2000 =[ 85 85 255 ] ./256; %blue
	% c4100 =[ 85 170 0  ] ./256; %dark green
	% c4850 =[ 0 255 128 ] ./256; %light green

	%% From colorblind colorpallete
    	 c300 = [1.0000         0         0];
    	 c800 = [1.0000    0.6445         0];
    	 c1700= [0.7977    0.1953    0.7969];
    	 c2000= [1.0000         0    1.0000];
    	 c4100= [     0    0.5000    1.0000];
    	 c4850= [     0    0.8000    1.0000];





	 cs = zeros(length(s),3);
	 for ii = 1:length(s)
	 	if(length(regexp(s(ii).KSTNM,'300')))
	 		cs(ii,:) = c300;
	 	elseif(length(regexp(s(ii).KSTNM,'800')))
	 		cs(ii,:) = c800;
	 	elseif(length(regexp(s(ii).KSTNM,'800')))
	 		cs(ii,:) = c800;
	 	elseif(length(regexp(s(ii).KSTNM,'1700')))
	 		cs(ii,:) = c1700;
	 	elseif(length(regexp(s(ii).KSTNM,'2000')))
	 		cs(ii,:) = c2000;
	 	elseif(length(regexp(s(ii).KSTNM,'4100')))
	 		cs(ii,:) = c4100;
	 	elseif(length(regexp(s(ii).KSTNM,'4850')))
	 		cs(ii,:) = c4850;
	 	else
	 		cs(ii,:) = c0;
	 	end
	 end
	 if(strcmp(plottype,'no_legend'))
		%nothing
	 else
	     if(strcmp(plottype,'plot'))
	     	plot(0,0,'color',c0,'linewidth',2)
	     	hold on
	     	plot(0,0,'color',c300,'linewidth',2)
	     	plot(0,0,'color',c800,'linewidth',2)
	     	plot(0,0,'color',c1700,'linewidth',2)
	     	plot(0,0,'color',c2000,'linewidth',2)
	     	plot(0,0,'color',c4100,'linewidth',2)
	     	plot(0,0,'color',c4850,'linewidth',2)
	     elseif(strcmp(plottype,'semilogx'))
	     	semilogx(0,0,'color',c0,'linewidth',2)
	     	hold on
	     	semilogx(0,0,'color',c300,'linewidth',2)
	     	semilogx(0,0,'color',c800,'linewidth',2)
	     	semilogx(0,0,'color',c1700,'linewidth',2)
	     	semilogx(0,0,'color',c2000,'linewidth',2)
	     	semilogx(0,0,'color',c4100,'linewidth',2)
	     	semilogx(0,0,'color',c4850,'linewidth',2)
	     end
	     legend('0','300','800','1700','2000','4100','4850')
	  end

      end %function
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function updateLegend(leg_list,flags)
	if(~flags.depth_coloring)
	   if(~flags.clear_before_plotting)
		prevleg = get(legend);
		leg_list = {prevleg.String{:},leg_list{:}};
		hh = findobj(gcf,'Type','Line');
		for kk = 1:length(hh)
			if(strcmp(hh(kk).DisplayName,'NoiseModel'))
				delete(hh(kk));
			end
		end
	   end
	end
	legend(leg_list)
      end % function updadeLegend
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function [cx,cy] = circle(R,nsteps);
	ang=0:360/nsteps:360;
	cx=R*sind(ang);
	cy=R*cosd(ang);
      end


   end % methods static
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   methods(Static, Access=protected)      
	
   end % methods static protected
end % classdef
