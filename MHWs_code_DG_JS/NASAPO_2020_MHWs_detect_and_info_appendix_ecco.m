function data_out = NASAPO_2020_MHWs_detect_and_info_appendix_ecco(input_params)%years,zlev,...
    %tag_field,tag_file)


input_vars.year_start_filename = input_params.year_start_filename;%1993;
input_vars.year_end_filename   = input_params.year_end_filename;%2016;
input_vars.ecco_name= input_params.ecco_name;%'ECCOv4r4';
input_vars.num_of_tiles = input_params.num_of_tiles;%13;
input_vars.ecco_monthly_num_of_levs = input_params.ecco_monthly_num_of_levs;%50;
input_vars.fpath     = input_params.user_path_ecco;
input_vars.ecco_other_vars = {'area' 'XC_lon' 'YC_lat' 'Z_depth' 'vol'};


% expected size of the variables (for some reason, G_forcing has different
% dimensions from other variables)
size1 = [input_params.ecco_num_of_points_in_tiles      input_params.ecco_num_of_points_in_tiles      ...
    input_vars.num_of_tiles   input_vars.ecco_monthly_num_of_levs    ...
    12*(input_params.year_end_in_file-input_params.year_start_in_file+1)];

size2 = [input_params.ecco_num_of_points_in_tiles      input_params.ecco_num_of_points_in_tiles      ...
    input_vars.ecco_monthly_num_of_levs    input_vars.num_of_tiles   ...
    12*(input_params.year_end_in_file-input_params.year_start_in_file+1)];


if isempty(input_params.tag_file)
    data_out.case_tag = [input_vars.ecco_name '_' input_params.tag_field '_zlev' ...
        num2str(input_params.zlev,'%02d')];
    
elseif (isempty(input_params.zlev) && contains(input_params.tag_file,'zint')) || ...
        contains(input_params.tag_file,'zint_k')
    data_out.case_tag = [input_vars.ecco_name '_' input_params.tag_field input_params.tag_file];
else
    checkwhatwrong
end
    

switch input_params.tag_field
    case {'heat','ohc_to50m'}
        if isfield(input_params, 'flag_use_only_DATA') && ...
                input_params.flag_use_only_DATA
            input_vars.var_names = {'DATA'};
        else
            input_vars.var_names = {'DATA' ...
                'G_total' 'G_forcing' ...
                'G_advection_conv' 'adv_vConv' ...
                'G_diffusion_conv' 'dif_vConv' ...
                }; % 'adv_meridH' 'G_advection_conv_mer' 'G_advection_conv_zon' ...
                
        end
        data_out.cax_percentile = input_params.cax_percentile;%[0 6];

    otherwise
        input_vars.var_names = {'DATA' ...
            'G_total' 'G_forcing' ...
            'G_advection_conv' 'adv_vConv' ...
            'G_diffusion_conv' 'dif_vConv' ...
            }; % 'adv_meridH'
        data_out.cax_percentile = input_params.cax_percentile;%[0 1];
end

data_out.var_names = input_vars.var_names;
data_out.flag_dim_lon_lat_time = input_params.flag_dim_lon_lat_time;%false;
data_out.flag_monthly_timeseries = input_params.flag_monthly_timeseries;%true;
data_out.num_gap_tsteps_in_event = input_params.num_gap_tsteps_in_event;

vars_all = horzcat(input_vars.ecco_other_vars,input_vars.var_names);

time_flag = true;

for i=1:length(vars_all)
    clear bfr bfr_tag_years
    switch vars_all{i}
        case {'XC_lon', 'YC_lat', 'Z_depth','vol', 'area'}
            bfr_tag_years = input_params.tag_years_metadata;
        otherwise
            bfr_tag_years = [num2str(input_vars.year_start_filename) '_' num2str(input_vars.year_end_filename)];
    end
    switch input_params.tag_field
        case {'salt','ohc_to50m'}
            bfr = read_vars_in_ncfile([input_vars.fpath input_vars.ecco_name ...
                '_' input_params.tag_field '_' vars_all{i} '_' ...
                bfr_tag_years '.nc']);
        case 'heat'
            bfr = read_vars_in_ncfile([input_vars.fpath input_vars.ecco_name '_' vars_all{i} '_' ...
                bfr_tag_years '.nc']);
    end
    
    switch vars_all{i}
        case {'XC_lon', 'YC_lat', 'Z_depth','vol'}
            eval([vars_all{i} ' = bfr.' vars_all{i} ';'])
        case {'area'}
            eval([vars_all{i} '_bfr = bfr.' vars_all{i} ';'])
        otherwise
            eval(['data = bfr.' vars_all{i} ';'])
            % convert to kelvin
            if contains(input_params.tag_file,'zint') && strcmp(vars_all{i},'DATA') && ...
                    ~contains(input_params.tag_file,'zint_k')
                data = data + 273.15;
            end
%             disp(size(data))
%             disp('-----')
%             disp(size(size1))
            if sum(size(data)-size1~=0)==0
                % disp('sum(size(data)-size1~=0)==0')
                
                if isempty(input_params.tag_file) || contains(input_params.tag_file,'zint_k')
                    data = squeeze(data(:,:,:,input_params.zlev,:));
                    
                else
                    %%%%%%%%%% for xint
                    if ~exist('DZ_depth','var')
                        if sum(size(vol)==[input_params.ecco_num_of_points_in_tiles      input_params.ecco_num_of_points_in_tiles      input_vars.ecco_monthly_num_of_levs    13])==4
                            
                        elseif sum(size(vol)==[input_vars.ecco_monthly_num_of_levs input_params.ecco_num_of_points_in_tiles      input_params.ecco_num_of_points_in_tiles      13])==4
                            vol = permute(vol,[2 3 1 4]);
                        else
                            checkwhatwrongwithsize
                        end
                        
                        DZ_depth = repmat(...
                            permute(vol,[1 2 4 3]) ./ ...
                            repmat(area_bfr,[1 1 1 size(vol,3)]),...
                            [1 1 1 1 size(data,5)]);
                        
                        bfrk = strsplit(strrep(strrep(input_params.tag_file,'_zint_k',''),'k',''),'_');
                        k1   = str2num(bfrk{1});
                        k2   = str2num(bfrk{2});
                        
                        data_out.layer_thickness = ...
                            sum(DZ_depth(80,80,1,k1:k2,1),4);
                        
                        data_out.depth_at_k1 = Z_depth(k1);
                        data_out.depth_at_k2 = Z_depth(k2);
                        
                        disp(['>>>>>>>> Layer thickness:' ...
                            num2str(data_out.layer_thickness) ' m'])
                        
                        
                    end
                    
                    %%%%%%%%%%
                    data = squeeze(nansum(data(:,:,:,k1:k2,:).*...
                        DZ_depth(:,:,:,k1:k2,:),4));
                    %%%%%%%%%%
                end
            elseif sum(size(data)-size2~=0)==0
                % disp('elseif')
                
                if isempty(input_params.tag_file) || contains(input_params.tag_file,'zint_k')
                    data = squeeze(data(:,:,input_params.zlev,:,:));
                else
                    %%%%%%%%%%
                    data = permute(data,[1 2 4 3 5]);
                    data = squeeze(nansum(data(:,:,:,k1:k2,:).*...
                        DZ_depth(:,:,:,k1:k2,:),4));
                    %%%%%%%%%%
                end
            else
                % check_what_wrong_with_input_dimensions
                disp('Input dimensions do not match the expected size1 or size2.');
                disp(['Size of data: ', mat2str(size(data))]);
                disp(['Expected size1: ', mat2str(size1)]);
                disp(['Expected size2: ', mat2str(size2)]);
            end
            data_out.size_source = size(data);

            data = reshape(data,[size(data,1)*size(data,2) size(data,3) size(data,4)]);
            
            % save time variable
            if time_flag
                % time
                bfr_years = input_params.year_start_in_file:input_params.year_end_in_file;
                time_datenum = nan(length(bfr_years)*12,1);
                n = 0;
                for iyr=1:length(bfr_years)
                    for mm=1:12
                        n = n+1;
                        time_datenum(n) = datenum(bfr_years(iyr),mm,15);
                    end
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
                data = NASAPO_2020_MHWs_detect_and_info_appendix_rm_seas_monthly(data);
                
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

if isempty(input_params.tag_file)
    data_out.Z_depth     = Z_depth(input_params.zlev);
    data_out.zlev        = input_params.zlev;
end

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
data_out.i_all = input_params.i_all;%[5548]; % 156.5W (around 204 east), 35.1757S
data_out.j_all = input_params.j_all;%[9];

% caxis for percentile plot

% variables to average over MHW event-relevant timesteps (e.g. over all the
% timesteps for an event, over the final ones only...) 
% NOTE: see also data_out.data2ave_ALL_from_data_in above
data_out.data2ave_ALL_from_find_MHWs_info = {'data_used4MHWs'};

return
