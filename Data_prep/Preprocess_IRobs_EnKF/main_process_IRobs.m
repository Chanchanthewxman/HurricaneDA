% Main algorithm producing GOES observations in EnKF-readable format & other features
% Author: Zhu (Judy) Yao. 2022.

% -------------------------------------------------------------------------------------
% Set up control parameters
% -------------------------------------------------------------------------------------
control = struct;
% --- domain type
control.fix_domain = true; % if domain if fixed or movable
% Assume: fixed domain --> best-track data is not needed; Movable domain: best-track data is needed.
if ~control.fix_domain
	control.bestrack_dir = '../../Preprocess_Obs/raw_Obs/Bestrack/'; % directory where best-track files are
end  
% --- Path
control.obs_dir = '/expanse/lustre/projects/pen116/cpruett/Tools_WRF_EnKF/BERYL/Data/Preprocess_Obs/raw_Obs/GOES_IR/'; % directory where raw GOES files are
control.obs_collect_dir = '/expanse/lustre/projects/pen116/cpruett/Tools_WRF_EnKF/BERYL/Data/Preprocess_Obs/raw_Obs/Collected_IR/'; %
control.geogrid_dir = '/expanse/lustre/projects/pen116/cpruett/Tools_WRF_EnKF/BERYL/Data/Preprocess_Domain/';

control.output_dir = '../../Tools_WRF_EnKF/BERYL/Data/Preprocess_Obs/toEnKFobs/GOES_IR/'; % directory where this algorithm outputs
% --- Storm information
control.storm_phase = {'BERYL',};  
control.period = {{'202406260000','202407091200'},};
% --- Satellite information
control.favCH = [8,];
control.facWL = {'6.2um', };
% --- WRF setup
control.domain = 'd01'; % WRF domain number
control.dx = 6; % WRF resolution: # km
% --- Other
control.filter_reso = [24]; %[36;24]; % data thinning distance in KM; SCL (successive covariance localization) is used.
control.roi_oh = {[100,100]}; %{[200,0]; [30,30]}; % ROI in KM [other variables, hydrometeors]
% the purpose of above two lines: 1) thin the obs in a 18-by-18 km box, with a 200-km radius of inﬂuence (ROI) for non-hydro variables and 0-km ROI for hydrometeors, and 2) thin the obs in a 12-by-12 km box with a 30-km ROI for all variables.

control.obsError = 3; % observation error
control.Sat_alt = 35000; % satellite altitude in KM
% -------------------------------------------------------------------------------------

dnow = datetime(now, 'ConvertFrom', 'datenum');
disp(['........Running the program at ', char(dnow), ' .............']); %datetime(now,'InputFormat','dd-MM-yyyy HH:mm:SS')]);

tStart_all = tic; % Start the time watch for the whole program
% -------------------------------------------------------------------------------------
% Loop through each storm object 
% -------------------------------------------------------------------------------------
for istorm = 1:length(control.storm_phase)

    % -------------------------------------------------------------------------------------
    % --- Make subdirectory for output
    % -------------------------------------------------------------------------------------
    if ~exist([control.output_dir,],'dir')
        [~, msg, ~] = mkdir(control.output_dir);
        if isempty(msg)
            disp(['Successfully created a subdirectory in ',control.output_dir]);
        else
            error('Error: ',msg);
        end
    end
    
    % ---------------------------------------------------------------------------------------------------------------
    % --- Generate hours of interest 
    % --- if movable domain: additionally generate best-track locations between synoptic times using linear interpolation
    % ---------------------------------------------------------------------------------------------------------------

	if control.fix_domain

		start_str = control.period{istorm}{1};
		end_str   = control.period{istorm}{2};
		% Convert to datetime format
		start_time = datetime(start_str, 'InputFormat', 'yyyyMMddHHmm');
		end_time   = datetime(end_str,   'InputFormat', 'yyyyMMddHHmm');
		% Generate time vector with hourly intervals
		time_vec = start_time:hours(1):end_time;
		% Convert datetime array back to string in 'yyyymmddHHMM' format
		% Generate DAtimes
		time_strs = datestr(time_vec, 'yyyymmddHHMM');

        % Collect DAtimes
        DAtimes_str = strings(size(time_strs,1),1);
        for ibr = 1:size(time_strs,1)
            DAtimes_str(ibr) = convertCharsToStrings(time_strs(ibr,:));
        end

	else 
	    [bestrack_str,start_time,end_time] = Hourly_Bestrack(istorm, control); % (cell)

		% Write hourly best-track data to a file
		filename = strcat(control.output_dir,'/bestrack_perHour_',start_time,'_',end_time);
	    disp(['Output the hourly best-track location and time: ',filename]);
		formatSpec = '%12s%12.3f%12.3f\n';
	    fileID = fopen(filename,'w');
		for itime = 1:size(bestrack_str,1)
	        fprintf(fileID, formatSpec, ...
		        bestrack_str{itime,1}, bestrack_str{itime,3}(1), bestrack_str{itime,3}(2));
	    end
		fclose(fileID);

		% Collect DAtimes
	    DAtimes_str = strings(size(bestrack_str,1),1);
		for ibr = 1:size(bestrack_str,1)
			DAtimes_str(ibr) = convertCharsToStrings(bestrack_str{ibr,1});
	    end

	end

    % -------------------------------------------------------------------------------------
    % --- Collect GOES Obs files within the period of interest into a directory
    % -------------------------------------------------------------------------------------
    disp('Collecting useful GOES obs files for this study......');
    Collect_GOESR(istorm, DAtimes_str, control);

    % -------------------------------------------------------------------------------------
    % --- Loop through each DAtime and output obs
    % -------------------------------------------------------------------------------------
    disp('Processing to EnKF readable obs......');
    Tb_dir = [control.obs_collect_dir,'*'];
    Tb_files = strsplit(ls(Tb_dir));
    Tb_files = Tb_files(~cellfun('isempty',Tb_files)); % get rid of annoying empty cell    

    for it = 1:size(DAtimes_str)
        files = [];
        for iTb = 1:length(Tb_files)
            Tb_file = Tb_files{iTb};
            if contains(Tb_file, DAtimes_str{it,1})
                files = [files, string(Tb_file)];
            else
                continue;
            end
        end
        
        DAtime = DAtimes_str(it);
        disp(['Dealing with the file at DAtime: ',DAtime]);
        if length(files) == 1
            SingleCH_write(istorm, DAtime, files, control);                
        %else
			% Muitiple-channel preprocessing is not developed. 
            %MultiCh_write(istorm, DAtime, files, control);
        end
   end
end

% -------------------------------------------------------------------------------------
% Diagnose elapsed time
% -------------------------------------------------------------------------------------
tEnd_all = toc(tStart_all);
disp(['Took ', num2str(tEnd_all), ' seconds to finish the whole workflow!']);
































