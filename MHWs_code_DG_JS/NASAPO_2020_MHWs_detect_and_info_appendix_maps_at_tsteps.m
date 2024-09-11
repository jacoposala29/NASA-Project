function NASAPO_2020_MHWs_detect_and_info_appendix_maps_at_tsteps(data_in,find_MHWs_info,vars2plot,...
    vars2plot_cax,vars2plot_title,tag_fig_out,fig_path)

scatter_size = 6;
    
for i=1:length(data_in.data_datenum)
    close all
    figure;set(gcf,'color','w','position',[0       0      1800        1200]);
    bfr_fpath = [];
    
    for ivar=1:length(vars2plot)
        nexttile()
        
        bfr = eval(vars2plot{ivar});
        bfr_add = '';
        
        if ivar==1
            bfr_add = [', ' datestr(data_in.data_datenum(i),'mmm-yyyy')];
        elseif ivar==length(vars2plot)
            bfr_fpath = [fig_path '/monthly_maps/'];
            bfr_add = '';
        end
        
        if ~data_in.flag_dim_lon_lat_time
            % ecco_xc,ecco_yc,ecco_data should all have the same
            bfr_ecco_data   = reshape(bfr,data_in.size_source);
            bfr_ecco_xc     = data_in.xc(:,:,:,1);
            bfr_ecco_yc     = data_in.yc(:,:,:,1);
            
            bfr_ecco_data   = bfr_ecco_data(:,:,:,i);
            
            % put on a regular grid
            data_on_grid = NASAPO_2020_MHWs_detect_and_info_appendix_ecco_to_regular_grid(...
                bfr_ecco_xc,bfr_ecco_yc,bfr_ecco_data);
            
            bfr2pl   = data_on_grid.data;
            bfr2pl_X = data_on_grid.X;
            bfr2pl_Y = data_on_grid.Y;
        else
            bfr2pl   = bfr(:,:,i);
            bfr2pl_X = data_in.X;
            bfr2pl_Y = data_in.Y;
        end
        
        
        switch vars2plot{ivar}
            case {'find_MHWs_info.data_used4MHWs','data_in.G_total',...
                    'data_in.G_forcing','data_in.G_diffusion_conv','data_in.G_advection_conv'}
                NASAPO_2020_MHWs_detect_and_info_appendix_plot_map(...
                    bfr2pl_X,bfr2pl_Y,bfr2pl,...
                    data_in.maplonlimit,vars2plot_cax{ivar},...
                    [vars2plot_title{ivar} bfr_add],...
                    [data_in.case_tag '_minLen' num2str(find_MHWs_info.delta_tstep) 'tsteps_' ...
                    datestr(data_in.data_datenum(i),'yyyymm') tag_fig_out],...
                    bfr_fpath);
            otherwise
                NASAPO_2020_MHWs_detect_and_info_appendix_plot_map_scatter(...
                    bfr2pl_X,bfr2pl_Y,bfr2pl,...
                    data_in.maplonlimit,vars2plot_cax{ivar},...
                    [vars2plot_title{ivar} bfr_add],scatter_size,...
                    [data_in.case_tag '_minLen' num2str(find_MHWs_info.delta_tstep) 'tsteps_' ...
                    datestr(data_in.data_datenum(i),'yyyymm') tag_fig_out],...
                    bfr_fpath);
        end
        
    end
    
end
return