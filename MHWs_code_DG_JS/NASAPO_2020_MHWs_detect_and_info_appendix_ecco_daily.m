function data_out = NASAPO_2020_MHWs_detect_and_info_appendix_ecco_daily(input_params)%years,zlev,...
%     tag_field,tag_file)

% please note that tag_file in the daily files is also used to read the
% variable in the file

% input_vars.field     = input_params.field;%tag_field;
% user_path_ecco = ['/Users/dgiglio/Work/DATA/ECCOv4r4/nc_files_JS/'  ...
%     input_vars.field '_global_tiles/daily/'];
% user_path_ecco = '/scratch/alpine/jasa1084/inputs_MHW_project_blanca/ECCO_daily/k0/1992_2018/';

% input_vars.year_start = 2004;
input_vars.year_start_fname = input_params.year_start_fname;%1992;

% input_vars.year_end   = 2018;
input_vars.year_end_fname = input_params.year_end_fname;%2018;

% input_vars.num_days   = 5113;
input_vars.num_days   = input_params.num_days;%9495;

% tag_years_metadata = '1993_2017';
% tag_years_metadata = '1993_2017';

% input_vars.ecco_name= 'ECCOv4r4';
input_vars.ecco_name= input_params.ecco_name;%'ECCOv4r4';

% input_vars.num_of_tiles = 13;
input_vars.num_of_tiles = input_params.num_of_tiles;% 13;

% delta days before/after to calculate percentile - should be the same for
% OISST daily
data_out.delta_days_4_percentile = input_params.delta_days_4_percentile;%5;
data_out.num_gap_tsteps_in_event = input_params.num_gap_tsteps_in_event;

% expected size of the variables (for some reason, G_forcing has different
% dimensions from other variables)
size1 = [90    90    input_params.num_of_tiles input_params.num_days];

size2 = [];

%%% adjust for OHC
if ~isempty(input_params.zlev)
    data_out.case_tag = [input_vars.ecco_name '_' input_params.tag_field '_zlev' num2str(input_params.zlev,'%02d') '_daily'];
else
    data_out.case_tag = [input_params.ecco_name '_' input_params.tag_field '_daily'];
end
    
input_vars.fpath     = input_params.user_path_ecco;

switch input_params.tag_field
    case 'heat'
        %         input_vars.var_names = {'DATA' ...
        %             'G_total' 'G_forcing' ...
        %             'G_advection_conv' 'adv_vConv' ...
        %             'G_advection_conv_mer' 'G_advection_conv_zon' ...
        %             'G_diffusion_conv' 'dif_vConv' ...
        %              }; % 'adv_meridH'
        input_vars.var_names = {'DATA' ...
            'G_total' 'G_forcing' ...
            'G_advection' 'adv_vConv' ...
            'G_diffusion' 'dif_vConv' ...
             }; % 'adv_meridH'
        data_out.cax_percentile = input_params.cax_percentile;%[0 6];

    otherwise
        to complete
        input_vars.var_names = {'DATA' ...
            'G_total' 'G_forcing' ...
            'G_advection_conv' 'adv_vConv' ...
            'G_diffusion_conv' 'dif_vConv' ...
             }; % 'adv_meridH'
        data_out.cax_percentile = input_params.cax_percentile;%[0 1];
end
%input_vars.ecco_other_vars = {'area' 'XC_lon' 'YC_lat' 'Z_depth' 'vol'};
input_vars.ecco_other_vars = {'area' 'XC_lon' 'YC_lat' 'Z_depth'};

data_out.var_names = input_vars.var_names;
%
data_out.flag_dim_lon_lat_time = input_params.flag_dim_lon_lat_time;% false;
data_out.flag_monthly_timeseries = input_params.flag_monthly_timeseries;%false;
data_out.flag_daily_timeseries   = input_params.flag_daily_timeseries;%true;

%

% load other vars
% vars_all = input_vars.ecco_other_vars;
% vars_all{end+1} = input_vars.var_names{1};

vars_all = horzcat(input_vars.ecco_other_vars,input_vars.var_names);

time_flag = true;
for i=1:length(vars_all)
    clear bfr
    %     switch tag_field
    %         case 'salt'
    %             bfr = read_vars_in_ncfile([input_vars.fpath input_vars.ecco_name ...
    %                 '_' tag_field '_' vars_all{i} '_' ...
    %                 num2str(input_vars.year_start) '_' num2str(input_vars.year_end+1) '.nc']);
    %         case 'heat'
    %             bfr = read_vars_in_ncfile([input_vars.fpath input_vars.ecco_name '_' vars_all{i} '_' ...
    %                 num2str(input_vars.year_start) '_' num2str(input_vars.year_end+1) '.nc']);
    %     end
    
    switch vars_all{i}
        case {'XC_lon', 'YC_lat', 'Z_depth','vol'}
            
            bfr = read_vars_in_ncfile([input_vars.fpath input_params.ecco_name ...
                '_' vars_all{i} '_' ...
                input_params.tag_years_metadata '.nc']);
            
            eval([vars_all{i} ' = bfr.' vars_all{i} ';'])
        case {'area'}
            bfr = read_vars_in_ncfile([input_vars.fpath input_params.ecco_name ...
                '_' vars_all{i} '_' ...
                input_params.tag_years_metadata '.nc']);
            
            eval([vars_all{i} '_bfr = bfr.' vars_all{i} ';'])
        otherwise
            if ~isempty(input_params.zlev)
                bfr = read_vars_in_ncfile([input_vars.fpath input_params.ecco_name ...
                    '_' vars_all{i} input_params.tag_file '_zlev' num2str(input_params.zlev) ...
                    '_' num2str(input_params.year_start_fname) '_' num2str(input_params.year_end_fname) '.nc']);
            else
                bfr = read_vars_in_ncfile([input_vars.fpath input_params.ecco_name ...
                    '_' vars_all{i} input_params.tag_file ...
                    '_' num2str(input_params.year_start_fname) '_' num2str(input_params.year_end_fname) '.nc']);
            end
            
            eval(['data = bfr.' vars_all{i} input_params.tag_file ';'])
            
            %             % convert to kelvin
            %             if contains(tag_file,'zint') && strcmp(vars_all{i},'DATA')
            %                 data = data + 273.15;
            %             end
            
            if sum(size(data)-size1~=0)==0
                disp([vars_all{i} ': correct size'])
            elseif sum(size(data)-size2~=0)==0
                code here if to permute
                %data = permute(data,[1 2 4 3 5]);
            else
                check_what_wrong_with_input_dimensions
            end
            data_out.size_source = size(data);

            data = reshape(data,[size(data,1)*size(data,2) size(data,3) size(data,4)]);
            
            % save time variable
            if time_flag
                % time
                switch bfr.time_units
                    case 'hours since 2004-01-01 12:00:00'
                        time_datenum = datenum(2004,1,1,12,0,0) + bfr.time./24;
%                     case 'hours since 1992-01-02 12:00:00'
                    case 'days since 1992-01-02 12:00:00'
                        time_datenum = datenum(1992,1,2,12,0,0) + bfr.time./24;
                    case 'hours since 1992-01-02 12:00:00'
                        time_datenum = datenum(1992,1,2,12,0,0) + bfr.time./24;
                    otherwise
                        case2code
                end
                msk_years    = year(time_datenum)>=min(input_params.years) & ...
                    year(time_datenum)<=max(input_params.years);

                data_out.data_datenum = time_datenum(msk_years);

                time_flag = false;
            end
            
            % select data in the years of interest
            data                    = data(:,:,msk_years);


            % except for DATA, remove seasonal cycle and trend
            if ~strcmp(vars_all{i},'DATA')
                % remove seasonal cycle for a daily timeseries
                data = NASAPO_2020_MHWs_detect_and_info_appendix_rm_seas_daily(data, data_out.data_datenum, ...
                    input_params.delta_days_4_percentile);
                
                % detrend
                data = reshape(detrend(reshape(data,size(data,1)*size(data,2),[])')',...
                    size(data,1),size(data,2),[]);
            end
            
            %
            eval([vars_all{i} ' = data;'])
    end
end

data_out.size_source(4) = sum(msk_years);


data_out.xc          = repmat(XC_lon,[1 1 1 data_out.size_source(4)]);
% convert longitudes to be from 20 to less than 380
data_out.xc(data_out.xc<20) = data_out.xc(data_out.xc<20) + 360;

data_out.yc          = repmat(YC_lat,[1 1 1 data_out.size_source(4)]);
data_out.area4d      = repmat(area_bfr,[1 1 1 data_out.size_source(4)]);
data_out.z_source    = Z_depth;

% if ~contains(tag_file,'zint') %isempty(tag_file)
%     data_out.Z_depth     = Z_depth(zlev);
%     data_out.zlev        = zlev;
% end

% data_out.Z_depth     = Z_depth(zlev);
% data_out.zlev        = zlev;


data_out.data = DATA;

n = 0;
for i=1:length(input_vars.var_names)
    if ~strcmp(input_vars.var_names{i},'DATA')
%         eval(['data_out.' input_vars.var_names{i} ' = ' ...
%             input_vars.var_names{i}  '(:,:,msk_years);'])
        eval(['data_out.' input_vars.var_names{i} ' = ' ...
            input_vars.var_names{i}  ';'])
        n = n + 1;
        data_out.data2ave_ALL_from_data_in{n} = input_vars.var_names{i};
    end
end



% data_out.x = squeeze(data_in.lon);
% data_out.y = squeeze(data_in.lat);
data_out.maplonlimit = [min(data_out.xc(:)) max(data_out.xc(:))];

% [data_out.X,data_out.Y] = ndgrid(data_out.x,data_out.y);

% define indices to look at in the timeseries plot at the end
data_out.i_all = input_params.i_all;%[6586 5548]; % 156.5W (around 204 east), 35.1757S
data_out.j_all = input_params.j_all;%[9 9];

% caxis for percentile plot


% variables to average over MHW event-relevant timesteps (e.g. over all the
% timesteps for an event, over the final ones only...) 
% NOTE: see also data_out.data2ave_ALL_from_data_in above
data_out.data2ave_ALL_from_find_MHWs_info = {'data_used4MHWs'};

return
