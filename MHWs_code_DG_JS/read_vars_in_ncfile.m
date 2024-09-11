function data = read_vars_in_ncfile(filename)
% read all the variables in the nc file and store them in the structure
% data: offsets and scalefactors need to be applied afterwards by the user

f = filename;

ncid = netcdf.open(f,'NC_NOWRITE');

% get info about the file (i.e. var names)
finfo = ncinfo(f);

% select attributes of interest
att_names = {'scale_factor' 'scale_factor' 'add_offset' '_FillValue' ...
    'missing_value' 'units' 'long_name' 'standard_name' 'calendar'};

for j=1:length(finfo(1).Variables(:))
    
    % check if variable name looks OK and make changes if needed
    var_name = getfield(finfo(1).Variables(j),'Name');
    % replace plus sign
    var_name2save = strrep(var_name,'+','_plus_');
    % replace -
    var_name2save = strrep(var_name2save,'-','_');
    % add a char at the start if the name starts with a number
    if ~isnan(str2double(var_name2save(1))) && ~strcmp(var_name2save(1),'i')
        % disp(['>> ' var_name])
        var_name2save = ['bfr_' var_name2save];
    end
    
    % if there are "__" at the beginning and end, we will remove them
    if length(var_name2save)>1 && strcmp(var_name2save(1:2),'__')
        var_name2save = var_name2save(3:end);
    end
    
    if length(var_name2save)>1 && strcmp(var_name2save(end-1:end),'__')
        var_name2save = var_name2save(1:end-2);
    end
        
    % var id
    varid = netcdf.inqVarID(ncid,var_name);
    
    % get var
    try 
        
%         eval(['data.' var_name2save ...
%             ' = squeeze(netcdf.getVar(ncid,varid,''double''));']) %
%             changed on Nov29 to make it work with ECCO monthly ohc
        eval(['data.' var_name2save ...
            ' = (netcdf.getVar(ncid,varid,''double''));'])
    catch
       try
            eval(['data.' var_name2save ...
                ' = squeeze(netcdf.getVar(ncid,varid,''char''));'])
        catch
            display(['Could not read variable ' var_name])
        end
    end
    
    % try and see if there are attributes you may need
    for ii=1:length(att_names)
        
        try
            eval(['data.' var_name2save ...
                '_' att_names{ii} ' = squeeze(netcdf.getAtt(ncid,varid,''' att_names{ii} '''));'])
        catch
            
        end
    end
end

netcdf.close(ncid)

return