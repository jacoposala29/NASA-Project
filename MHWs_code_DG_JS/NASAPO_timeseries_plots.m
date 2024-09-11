clear all
close all 

set(0, 'DefaultFigureVisible', 'on')

dir_path = '/Users/dogi7244/Work/DATA/Sala/MHWs_chapter2/OISST_box_test_June19/output/figures/';
fname    = 'oisst_v2_1993_2016_prcnt90_noTrend_minLen_1tsteps_withAVE';

data_in = load([dir_path fname]);

% let's select indices of the point of interest
ind_dim1 = 1;
ind_dim2 = 1;

data_in_1d = {};

fields_in = fieldnames(data_in.find_MHWs_info);

data_in_1d.time_dnum = [datenum(data_in.find_MHWs_info.years(1),1,1): ...
    datenum(data_in.find_MHWs_info.years(end),12,31)];

for i=1:length(fields_in)
    clear bfr
    bfr = eval(['data_in.find_MHWs_info.' fields_in{i}]);
    if length(size(bfr))==3
        eval(['data_in_1d.' fields_in{i} ' = squeeze(bfr(ind_dim1,ind_dim2,:));'])
    else
        eval(['data_in_1d.' fields_in{i} ' = bfr;'])
    end
end

data_in_1d.data_noTrend = data_in_1d.data_seasonal + data_in_1d.data_used4MHWs;
data_in_1d.data_seasonal_anom = data_in_1d.data_seasonal - mean(data_in_1d.data_seasonal);
data_in_1d.data_percentile3d_plus_seas = data_in_1d.data_seasonal + data_in_1d.data_percentile3d;
data_in_1d.data_percentile3d_not_smooth_plus_seas = data_in_1d.data_seasonal + data_in_1d.data_percentile3d_not_smooth;

fields2plot = {{'data_seasonal' 'data_percentile3d_plus_seas' ...
    'data_percentile3d_not_smooth_plus_seas' 'data_noTrend'} ...
    {'data_seasonal' 'data_percentile3d_plus_seas' 'data_percentile3d_not_smooth_plus_seas'} ...
    {'data_used4MHWs' 'data_percentile3d' 'data_percentile3d_not_smooth'}};

fields2plot_cols = {{'b' 'r' [.5 .5 .5] 'k'} {'b' 'r' [.5 .5 .5]} {'k' 'r' [.5 .5 .5]}};

x_lim_select_ALL = { ...
    [data_in_1d.time_dnum(1) data_in_1d.time_dnum(end)+1] ...
    [datenum(2012,3,1) datenum(2016,8,31)]...
    };
% [datenum(2015,3,1) datenum(2016,6,30)]...
% [datenum(2008,3,1) datenum(2011,3,31)]...
    
for ilim=1:length(x_lim_select_ALL)
    figure('color','w','Units', 'inches','position',[0 0 1200 800])
    tiledlayout(length(fields2plot),1)
    for i=1:length(fields2plot)
        nexttile()
        for j=1:length(fields2plot{i})
            clear bfr
            bfr = eval(['data_in_1d.' fields2plot{i}{j}]);
            
            plot(data_in_1d.time_dnum,bfr,'LineWidth',2,'color',fields2plot_cols{i}{j})

            grid on
            set(gca,'fontsize',26,'xlim',x_lim_select_ALL{ilim})
            datetick('x','mm/dd/yy','keeplimits','keepticks') 
            hold on
        end
    end
end

