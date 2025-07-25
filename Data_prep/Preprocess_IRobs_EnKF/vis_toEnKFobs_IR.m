% Validate the GOES-R obs that is readable by EnKF by plotting

% -------------------------------------------------------------------------------------
% Set up control parameters
% -------------------------------------------------------------------------------------
control = struct;
% ----Path
control.obs_dir = '../../Preprocess_Obs/toEnKFobs/GOESR_IR/';
control.output_dir = '../../Preprocess_Obs/Visual/toEnKFobs/GOESR_IR/';
control.geogrid_dir = '../../Preprocess_Domain/';
% ---Storm information
control.storm_phase = 'IRMA';  
% --- WRF setup
control.domain = 'd03';
control.dx = 3; % WRF resolution: 3 km
% ---Satellite information
control.satName = "GOES16";
control.favCH = "8";
control.facWL = {'6.2um', };
% ---Other
control.filter_reso = [18;12];
% 18-by-18 km box with a 200-km radius of inﬂuence for non-hydro variables; 12-by-12 km box with a 30-km radius of inﬂuence for all variables
control.roi_oh = {[200,0], [30,30]}; % roi [other variables, hydrometeors]


dnow = datetime(now, 'ConvertFrom', 'datenum');
disp(['........Plotting at ', char(dnow), ' .............']);

% -------------- More set up --------------
tStart_all = tic; % Start the time watch for the whole program
obs_dir = [control.obs_dir,control.storm_phase,'/*_so'];
bestrack_dir = [control.obs_dir,control.storm_phase,'/bestrack_perHour'];
obs_files = strsplit(ls(obs_dir));
obs_files = obs_files(~cellfun('isempty',obs_files));

if ~exist([control.output_dir,control.storm_phase],'dir')
    [~, msg, ~] = mkdir(control.output_dir,control.storm_phase);
     if isempty(msg)
        disp(['Successfully created a subdirectory in ',control.output_dir,' for ',control.storm_phase]);
     else
        error('Error: ',msg);
     end
end

% -------------- Read raw best-track file ------------
fid = fopen(bestrack_dir);
bt_record = textscan(fid,'%s','delimiter','');
fclose(fid);
bt_str_all = string(bt_record{1}(:));

loc_storm = strings(length(bt_str_all), 3);
for ir = 1:length(bt_str_all)
    bt_str_per = strsplit(bt_str_all(ir));
    loc_storm(ir,:) = bt_str_per; % time,lat,lon
end

% -------------- Iterate through so files ------------
for iso = 1:length(obs_files)
    DA_time = loc_storm(iso,1);
    disp(['DA time is ', char(DA_time), ' ......']);
    so_file = obs_files{iso};
    
    disp(['Reading obs file: ',so_file,' ......']);
    % - Read the so file
    fid = fopen(so_file);
    obs_record = textscan(fid,'%s','delimiter','');
    fclose(fid);
    obs_str_all = string(obs_record{1}(:));
    len_record = length(obs_record{1}(:)); % how many records there are per so file
    obs_4 = zeros(len_record,4);
    for ir = 1:len_record
        obs_str_per = strsplit(obs_str_all(ir));
        obs_4(ir,:) =  str2double(obs_str_per(1,4:7)); % latitude, longitude, Tb value, ROI for hydro
    end

    % - Separate obs with large ROI from small ROI
    idx_largeROI =  obs_4(:,4) == 0;
    idx_smallROI = obs_4(:,4) == 30;
    obs = cell(length(control.roi_oh),1);
    obs{1,1} = obs_4(idx_largeROI,1:3);
    obs{2,1} = obs_4(idx_smallROI,1:3);

    % - Define map boundaries using the geo_em.d03.nc generated by WPS geogrid.exe
    geo_file =  [control.geogrid_dir,control.storm_phase,'/',char(DA_time),'/geo_em.',control.domain,'.nc'];
    disp(['Reading ', geo_file, '......']);
    xlat_m = ncread(geo_file,'XLAT_M');
    xlon_m = ncread(geo_file,'XLONG_M');
    min_xlat = double(min(xlat_m,[],'all')-0.5);
    max_xlat = double(max(xlat_m,[],'all')+0.5);
    min_xlon = double(min(xlon_m,[],'all')-0.5);
    max_xlon = double(max(xlon_m,[],'all')+0.5);
   
    % ---------- Plot the figure -----------
    figure; hFig=gcf; set(hFig, 'Position', [0 0 800 800]);
    % scatter Tbs on a projected map
    m_proj('mercator','lon',[min_xlon max_xlon],'lat',[min_xlat max_xlat]);
    for iroi = 1:length(control.roi_oh)
        if iroi == 1
            H = m_scatter(obs{iroi,1}(:,2), obs{iroi,1}(:,1),50,obs{iroi,1}(:,3),'^','filled');
            hold on;
        else
            H = m_scatter(obs{iroi,1}(:,2), obs{iroi,1}(:,1),50,obs{iroi,1}(:,3),'v','filled');
            hold on;
        end
    end 
    lat_bt = str2double(loc_storm(iso,2)); lon_bt = str2double(loc_storm(iso,3));
    m_scatter(lon_bt,lat_bt,50,0, '*'); % 0 is represented by black color in this colormap
    hold on;
    % use the customized colormap
    colormap(IR_colormap(0.5)); caxis([185 325]); 
    cb = colorbar;
    set(cb,'Fontsize', 23);
    cb.Label.String = 'Brightness Temperature (K)';
    % add coastline 
    m_coast('color','k');
    % grid lines 
    lon_range = round(min_xlon:2:max_xlon);
    lat_range = round(min_xlat:2:max_xlat);
    m_grid('xtick',lon_range,'ytick',lat_range,'tickdir','out','fontsize',22);
    xlh = xlabel(['DA time: ', DA_time], 'Fontsize',22,'fontweight','bold');
    xlh.Position(2) = xlh.Position(2) - 0.01;  % move the label 0.02 data-units further down 
    % title
    title_char1 = control.storm_phase + ": " + control.satName + " CH"+ control.favCH;
    %title_char2 = "Filtered Resolution: " + num2str(control.filter_reso(iroi))  + " KM";
    title_char =  title_char1; %[title_char1, title_char2];
    title(title_char,'Fontsize',20)

    %save_dir = [control.output_dir+string(control.storm_phase)+'/'+string(DA_time)+'_'+num2str(control.filter_reso(iroi))+'km'+'.png'];
    save_dir = [control.output_dir+string(control.storm_phase)+'/'+control.storm_phase+'_'+string(DA_time)+'.png'];
    saveas(gcf,save_dir);       
    disp(['Saving ', char(save_dir), ' ......']); 

end





