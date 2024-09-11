function NASAPO_2020_MHWs_detect_and_info_appendix_plot_map(X,Y,d,...
    maplonlimit,cax,fig_title,fig_tag,fig_path)

axesm('robinson', 'maplonlimit', maplonlimit, 'frame', 'on');

pcolorm(Y',X',d');shading flat;colorbar;
geoshow('landareas.shp','FaceColor',[.85 .85 .85])

set(gca,'fontsize',28)
if ~isempty(fig_title)
    if ~contains(fig_title,'\')
        title(fig_title,'interpreter','none')
    else
        title(fig_title)
    end
end

colormap('Parula')
if ~isempty(cax)
    caxis(cax)
    
    if abs(min(cax))==max(cax)
        colormap([winter(16);[.9 .9 .9];flipud(autumn(16))])
    end
    
    if min(cax)==0 && max(cax)==ceil(max(cax))
        colormap(jet(max(cax)*4+1))
    end
end

axis off

if ~isempty(fig_path)
    print('-dpng',[fig_path fig_tag '.png'])
end

return