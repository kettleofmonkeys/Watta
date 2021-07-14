function [found, iceLayer] = post_find_ice_surf_newest( startPond, endPond)
        global g_dif2Top g_dif2Bottom g_dif2Last g_DepthVal g_heights ...
        g_top g_bottom g_next g_heights g_Flat_Change ...
        g_TopSmooth_bin g_Flat_Change_bin g_BottomSmooth_bin_fix g_x_dist g_BottomSmooth_bin ...
        g_iceLayer g_dist g_heights_orig

    found = false;
%     fRegion = g_Flat_Change_bin(startPond:endPond);
    fRegion = g_Flat_Change(startPond:endPond);

    fRegion(1:startPond-1) = NaN;
    iceLayer(1:startPond-1) = NaN;

    if( endPond < size(fRegion,1)+1)
        fRegion(endPond+1:end) = NaN;
        iceLayer(endPond+1:end) = NaN;
    end    
    %this is used to pull the point where the 2nd derivative can 
    %distinguish ice from lake cover. Typically therere are two peaks
%    iceLayer = zeros(size(g_Flat_Change_bin)); %zeros(size(g_heights_orig));
%    iceLayer(1:startPond-1) = NaN;
%    iceLayer(endPond+1:end) = NaN;
    
    diff2_is = fRegion;
    sCheck = nansum(~isnan(diff2_is));
    if(sCheck < 10)
        iceLayer = NaN;
        found = false;
        return;
    end
    
%    diff2(depthEst > 1.5) = NaN;    %only the shallow estimates for ice
    pdSix = fitdist(diff2_is,'Kernel');
    minVal = min(diff2_is);
    maxVal = max(diff2_is);
    if( (maxVal - minVal) < 0.2)
        found = false;
        return;        
    end
    %xDft = minVal-1:0.05:maxVal+1;
    xDft = minVal-1:0.02:maxVal+1;

    ySix = pdf(pdSix,xDft);
    diff1 = gradient(ySix) ./ gradient(xDft);
    diff2 = gradient(diff1) ./ gradient(xDft);
    [peakDiffs, diff2_line, prob2,  iceL, foundN] = find_peaks_of_curve_simple_Tr(diff1, diff2, xDft, ySix );
    
    % if ic
    % breakPoint = (peakDiffs(1) - peakDiffs(2)) / 2;
%    peakDiffs(diff2_line > 0) = [];
%    if( size(peakDiffs) 
%    diff2_line(diff2_line > 0) = [];
    %{
    aDf = ySix;
    xTest = xDft;
    aDf(xTest > peakDiffs(1)) = NaN;
    aDf(xTest < peakDiffs(2)) = NaN;
    xDifOfMin_ind = find(aDf == min(aDf));
    breakPoint = xTest(xDifOfMin_ind);
    %}
    peak2ndD_lake = peakDiffs(1);
    peak2ndD_ice = peakDiffs(2);
%    peak2ndD_ice = max(peakDiffs(1), peakDiffs(2));
%    peak2ndD_lake = min(peakDiffs(1), peakDiffs(2));
    if( foundN && (peak2ndD_ice - peak2ndD_lake) > 0.1 )
%        breakPt = 
        breakPt = iceL.height;%(peak2ndD_lake + peak2ndD_ice) ./ 2;  
        iceLayer = fRegion >  breakPt;
%        iceLayer = rIceLayer;
        found = true;
    else
        iceLayer = NaN;
        found = false;
    end
    t = 2;

end
