classdef duglSet

   properties (Constant = true)
	functions = {'none'};
        DATE_FORMATTER    = 'yyyy-mm-dd HH:MM:SS.FFF'; %default data format, in ms
	%DATA_DIR = '/net/gauss/data2/dbowden/Homestake_SAC/';
	DATA_DIR = '/Volumes/DBEXT3openB/Homestake_SAC';
	RESP_DIR = '/Volumes/DBEXT3openB/Homestake_SAC/Response';
        %SPECT_DIR    = '/net/gauss/data1/dbowden/dugl_analysis/Power_spect_tukey_OLD/'; 
        %SPECT_DIR    = '/net/gauss/data1/dbowden/dugl_analysis/Power_spect_tukey_halfoctave/'; 
	SPECT_DIR = '/Volumes/DBEXT3openB/Homestake_analysis/Power_spect_tukey_halfoctave';
        %MAP_DIR    = '/net/gauss/data1/dbowden/ligo/Homestake_sacmat/Map/';
	%MAP_DIR = '/Volumes/DBEXT3openA/dugl_analysis/Map';
	MAP_DIR = '/Users/danielbowden/Documents/MATLAB/Homestake_MAP/Map/';
        TAUP_DIR    = '/Users/danielbowden/Documents/MATLAB/matTaup/';
        IRIS_DIR    = '/Users/danielbowden/Documents/MATLAB';
	Components = {'HHZ','HHE','HHN'};
	center_of_array = [-103.757047 44.348816 877.761077 ];  % arbitrary, but needs to be consistent for some scripts
	surface_height = 1604;




	Stations = {'ROSS','YATES','ORO','WTP','SHL','DEAD','LHS','RRDG','TPK','300','800','1700','A2000','B2000','C2000','D2000','E2000','A4100','C4100','D4100','A4850','B4850','C4850','D4850','RSSD'};
	Locations = [	-103.7576 	44.3450 	1637.6  	%ROSS	
			-103.7515 	44.3522 	1634.3  	%YATES 	
			-103.7521 	44.3435 	1544.1  	%ORO 	
			-103.7425 	44.3537 	1561.9  	%WTP 	
			-103.7095 	44.3169 	1768.3  	%SHL     
			-103.7531 	44.3826 	1518.5  	%DEAD    
			-103.7748 	44.3473 	1658.0  	%LHS     
			-103.7654 	44.3595 	1677.0  	%RRDG    
			-103.7978 	44.3411 	1698.0  	%TPK     
			-103.756874	44.346445	1505.1024	%300	
			-103.758864	44.347394	1349.9592 	%800	
			-103.758421	44.351723	1073.84088	%1700	
			-103.762288	44.351133	983.25432 	%A2000	
			-103.760497	44.349168	983.339664	%B2000	
			-103.763709	44.351615	983.111064	%C2000	
			-103.7689	44.353313	983.028768	%D2000	
			-103.771554	44.356405	983.144592	%E2000	
			-103.753484	44.344886	342.494616	%A4100	
			-103.75068	44.35120 	342.326976	%C4100	
			-103.751083	44.34358 	342.41232 	%D4100	
			-103.762542	44.340474	115.165632	%A4850	
			-103.75813	44.346273	114.87912 	%B4850	
			-103.752818	44.346315	114.616992	%C4850	
			-103.750509	44.353012	115.18392 	%D4850	
			-104.035900	44.121200	2090 	];	%RSSD	








   end %constant properties
end %classdef

        %DATA_DIR    = '/net/gauss/data1/dbowden/dugl/sacdb/sac/';                         % realtime, but incomplete
        %DATA_DIR    = '/net/gauss/data1/dbowden/dugl_complete/HomestakeData/sacdb/sac';    % includes bailers, not real time
	%DATA_DIR = '/net/gauss/data1/dbowden/dugl_complete/Homestake_Fetch';
        %RESP_DIR    = '/net/gauss/data1/dbowden/dugl/Response_files/'; 
        %RESP_DIR    = '/net/gauss/data1/dbowden/dugl_complete/Homestake_Fetch/Response'; 
