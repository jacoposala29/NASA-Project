clear all
close all

set(0, 'DefaultFigureVisible', 'on')

%%%%%%% group of cases to compare
dir_in = {...
%     '/Users/jacoposala/Downloads/MHWs_data_ouput_Aug16/ECCO_daily_avg_box_NEP_1992_2017/' ...
%     '/Users/jacoposala/Downloads/MHWs_data_ouput_Aug16/OISST_daily_avg_box_NEP_1992_2017/' ...
%     '/Users/jacoposala/Downloads/MHWs_data_ouput_Aug16/ECCO_daily_avg_box_SWP_1992_2017/' ...
    '/Users/jacoposala/Downloads/MHWs_data_ouput_Aug16/ECCO_daily_avg_box_TASMAN_1992_2017/' ...
%     '/Users/jacoposala/Downloads/MHWs_data_ouput_Aug16/OISST_daily_avg_box_SWP_1992_2017/' ...
    '/Users/jacoposala/Downloads/MHWs_data_ouput_Aug16/OISST_daily_avg_box_TASMAN_1992_2017/' ...
    };
file_in = {...
    'ECCOv4r4_heat_daily_box_1992_2017_prcnt90_smooth_noTrend_minLen_5tsteps_maxGap_2tsteps_withAVE.mat' ...
%     'ECCOv4r4_heat_daily_box_1992_2017_prcnt90_smooth_noTrend_minLen_5tsteps_maxGap_2tsteps_withAVE.mat' ...
%     'ECCOv4r4_heat_daily_box_1992_2017_prcnt90_smooth_noTrend_minLen_5tsteps_maxGap_2tsteps_withAVE.mat' ...
    'oisst_v2_1992_2017_prcnt90_smooth_noTrend_minLen_5tsteps_maxGap_2tsteps_withAVE.mat'...
%     'oisst_v2_1992_2017_prcnt90_smooth_noTrend_minLen_5tsteps_maxGap_2tsteps_withAVE.mat' ...
%     'oisst_v2_1992_2017_prcnt90_smooth_noTrend_minLen_5tsteps_maxGap_2tsteps_withAVE.mat' ...
    };
ix = 1;
iy = 1;

% date_start = datenum(2013,11,1); % NEP
date_start = datenum(2012,3,1); % TASMAN
date_start = datenum(2012,3,1); % TASMAN

% date_end   = datenum(2015,6,1); % NEP
date_end   = datenum(2016,7,1); % TASMAN

%%%%%%% another group of cases to compare
% dir_in = {'/Users/jacoposala/Downloads/MHWs_data_ouput_Aug16/ECCO_monthly_heat_zlev01_1993_2016/' ...
%     '/Users/jacoposala/Downloads/MHWs_data_ouput_Aug16/ECCO_monthly_heat_zlev01_1993_2016/' ...
%     '/Users/jacoposala/Downloads/MHWs_data_ouput_Aug16/ECCO_monthly_heat_zlev01_1992_2017/' ...
%     '/Users/jacoposala/Downloads/MHWs_data_ouput_Aug16/ECCO_monthly_heat_zlev01_1992_2017/' ...
%     };
% file_in = {...
%     'ECCOv4r4_heat_zlev01_1993_2016_prcnt90_noTrend_minLen_1tsteps_maxGap_0tsteps_withAVE.mat' ...
%     'ECCOv4r4_heat_zlev01_1993_2016_prcnt90_noTrend_minLen_1tsteps_maxGap_2tsteps_withAVE.mat' ...
%     'ECCOv4r4_heat_zlev01_1992_2017_prcnt90_noTrend_minLen_1tsteps_maxGap_0tsteps_withAVE.mat' ...
%     'ECCOv4r4_heat_zlev01_1992_2017_prcnt90_noTrend_minLen_1tsteps_maxGap_2tsteps_withAVE.mat' ...
%     };
% 
% ix = 5548; %1;
% iy = 9; %1;
% 

% for OISST 
% ix=211;
% iy=90

% date_start = datenum(2000,1,1);
% date_end   = datenum(2016,12,31);

% %%%%%%% another group of cases to compare
% dir_in = {'/Users/jacoposala/Downloads/MHWs_data_ouput_Aug16/ECCO_daily_avg_box_NEP_OHC_1992_2017/' ...
%     '/Users/jacoposala/Downloads/MHWs_data_ouput_Aug16/ECCO_daily_avg_box_NEP_OHC_1992_2017/' ...
%     '/Users/jacoposala/Downloads/MHWs_data_ouput_Aug16/ECCO_daily_avg_box_NEP_OHC_1992_2017/' ...
%     };
% file_in = {...
%     'ECCOv4r4_heat_daily_box_1992_2017_prcnt90_smooth_noTrend_minLen_5tsteps_maxGap_0tsteps_withAVE.mat' ...
%     'ECCOv4r4_heat_daily_box_1992_2017_prcnt90_smooth_noTrend_minLen_5tsteps_maxGap_2tsteps_withAVE.mat' ...
%     'ECCOv4r4_heat_daily_box_1992_2017_prcnt90_smooth_noTrend_minLen_5tsteps_maxGap_3tsteps_withAVE.mat' ...
%     };
% 
% date_start = datenum(2013,11,1);
% date_end   = datenum(2016,1,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% starting here: %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fnt_sz = 18;
markers = {'*' 'o','s','d','<','>'};
cols    = {'k' 'r' 'b' 'm' 'c' 'y'};
styles  = {'-' '--' '-' '--' '-' '--'};


data_in = {};

for i=1:length(file_in)
    data_in{i} = load([dir_in{i} file_in{i}]);
end


% percentile smooth (if used)
bfr_title = 'Percentile used (may be smooth, if smoothing was applied)';
figure('color','w','Units', 'inches','position',[0 0 21 8])
for i=1:length(file_in)
    plot(...
        squeeze(data_in{i}.find_MHWs_info.data_percentile3d(ix,iy,:)),...
        "Color",cols{i},"LineStyle",styles{i},'LineWidth',length(file_in)+5-i)
    hold on
    
end
title(bfr_title)
set(gca,'FontSize',fnt_sz,'xlim',[1 sum(year(data_in{i}.find_MHWs_info.data_datenum)==year(data_in{i}.find_MHWs_info.data_datenum(1)))])

% percentile
if isfield(data_in{1}.find_MHWs_info,'data_percentile3d_not_smooth')
    bfr_title = 'Percentile (before smoothing, if smoothing was applied)';
    figure('color','w','Units', 'inches','position',[0 0 21 8])
    for i=1:length(file_in)
        plot(...
            squeeze(data_in{i}.find_MHWs_info.data_percentile3d_not_smooth(ix,iy,:)),...
            "Color",cols{i},"LineStyle",styles{i},'LineWidth',length(file_in)+5-i)
        hold on

    end
    title(bfr_title)
    set(gca,'FontSize',fnt_sz,'xlim',[1 sum(year(data_in{i}.find_MHWs_info.data_datenum)==year(data_in{i}.find_MHWs_info.data_datenum(1)))])
end

% close all
% events tstep_part_of_event_msk
figure('color','w','Units', 'inches','position',[0 0 21 16])
tiledlayout(length(file_in),1)
for i=1:length(file_in)

    nexttile()
    bfr_title = file_in{i};

    clear bfr_msk
    
    bar(data_in{i}.find_MHWs_info.data_datenum,...
        squeeze(data_in{i}.find_MHWs_info.data_used4MHWs(ix,iy,:)),...
        "faceColor",[.85 .85 .85],"LineStyle",'-')
    hold on

    bfr_msk = squeeze(data_in{i}.find_MHWs_info.tstep_part_of_event_msk(ix,iy,:));
    bar(data_in{i}.find_MHWs_info.data_datenum(bfr_msk),...
        squeeze(data_in{i}.find_MHWs_info.data_used4MHWs(ix,iy,bfr_msk)),...
        "faceColor",'y')
    
    bfr_msk = squeeze(data_in{i}.find_MHWs_info.data_used4MHWs(ix,iy,:))>...
        squeeze(data_in{i}.find_MHWs_info.data_percentile3d(ix,iy,:));
    bar(data_in{i}.find_MHWs_info.data_datenum(bfr_msk),...
        squeeze(data_in{i}.find_MHWs_info.data_used4MHWs(ix,iy,bfr_msk)),...
        "faceColor",[1 .8 0])

    bfr_msk = squeeze(data_in{i}.find_MHWs_info.start_tstep_msk(ix,iy,:));
    xline(data_in{i}.find_MHWs_info.data_datenum(bfr_msk),'r')

    bfr_msk = squeeze(data_in{i}.find_MHWs_info.end_tstep(ix,iy, ...
        ~isnan(squeeze(data_in{i}.find_MHWs_info.end_tstep(ix,iy,:)))));
    xline(data_in{i}.find_MHWs_info.data_datenum(bfr_msk),'k',"LineStyle",'-')

    bfr_msk = squeeze(data_in{i}.find_MHWs_info.peak_tstep_msk(ix,iy,:));
    % xline(data_in{i}.find_MHWs_info.data_datenum(bfr_msk),'r',...
    %     "LineStyle",'--','LineWidth',2)
    plot(data_in{i}.find_MHWs_info.data_datenum(bfr_msk),...
        squeeze(data_in{i}.find_MHWs_info.data_used4MHWs(ix,iy,bfr_msk)),'*r')

    plot(data_in{i}.find_MHWs_info.data_datenum(:),...
        squeeze(data_in{i}.find_MHWs_info.data_percentile3d(ix,iy,:)),...
        "Color",'r',"LineStyle",'-','LineWidth',1)
    
    title(bfr_title,'interpreter','none')
    set(gca,'FontSize',fnt_sz,'xlim',[date_start date_end])

    ylabel({'Data to' 'define MHWs'})
    datetick('x','mm/dd/yy','keeplimits')

end


