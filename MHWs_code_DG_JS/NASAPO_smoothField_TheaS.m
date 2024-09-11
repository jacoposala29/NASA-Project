function smoothedField = NASAPO_smoothField_TheaS(field, smoothType, smoothWindow)

switch smoothType
    case 'loess'
        smoothedField = smoothdata(field,'loess',smoothWindow);
    case 'movmean'
        if mod(smoothWindow, 2) > 0 % if window size is odd
            smoothedField = movmean(field,smoothWindow,'Endpoints','discard');
        else
            smoothWindow = smoothWindow/2;
            iMonthPlot = 1;
            for iMonth = 1+smoothWindow:length(field)-smoothWindow
                fieldWindow = field(iMonth-smoothWindow:iMonth+smoothWindow);
                fieldWindow(1) = fieldWindow(1)/2;
                fieldWindow(end) = fieldWindow(end)/2;
                smoothedField(iMonthPlot) = mean(fieldWindow);
                iMonthPlot = iMonthPlot + 1;
            end
        end
end

end