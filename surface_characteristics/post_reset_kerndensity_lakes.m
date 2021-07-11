function op = post_reset_kerndensity_lakes(subheights, top_filt, val_St, photon_count, sizeN)
    %%% detect bottom in given segment defined by val_St (start photon) and
    %%% sizN (number of photons in segment. 
    
    %%% create arrays for recalculation
    re_notEnoughInformation = zeros(sizeN);
    re_bottom_orig = nan(sizeN);
    re_next_orig = nan(sizeN);
    re_bottom = nan(sizeN);
    re_next = nan(sizeN);
    re_isNarrowLayer = nan(sizeN);
    re_isBottomResolutionDependent = nan(sizeN);
    re_isNextBottom = nan(sizeN);
    re_interpolated = nan(sizeN);
    re_bottom_prob = nan(sizeN);
    re_next_prob = nan(sizeN);
    re_dif2Bottom = nan(sizeN);
    re_dif2Last = nan(sizeN);
    
    re_bottom_lg = nan(sizeN);
    re_next_lg = nan(sizeN);
    re_dif2Bottom_lg = nan(sizeN);
    re_dif2Last_lg = nan(sizeN);
    szPht = sizeN(1); 

    %%% loop through all photons in segment to detect bottom
    for iM = 1: szPht       %%% if using parellel computing, change 'for' to 'parfor'
        midPoint_index = iM;
        midPoint_index = midPoint_index+val_St-1;
        %choose a larger window (where there will inherently be NaN values
        wideBerth = 3000;           
        maxStart = max(midPoint_index - wideBerth,1);
        maxEnd = min(midPoint_index + wideBerth, size(subheights,1));
        s_heights = subheights(maxStart:maxEnd); %create a subarray 
        topValue = top_filt(iM+val_St-1);
        inVl = midPoint_index - maxStart + 1;       %this is relative to the beginning of the subsection
        fVal = set_single_depth_value_under_lyr( inVl, photon_count, topValue, s_heights );
        
        %%
        re_notEnoughInformation(iM) = fVal.notEnoughInformation;
        if( re_notEnoughInformation(iM) == 0)
            re_isNarrowLayer(iM) = fVal.isNarrowLayer;    
            re_isBottomResolutionDependent(iM) = fVal.isBottomResolutionDependent;
            re_isNextBottom(iM) = fVal.isNextBottom;
            re_interpolated(iM) = fVal.interpolated;
            re_bottom_lg(iM) = fVal.bottom_lg;
            re_next_lg(iM) = fVal.next_lg;
            re_dif2Bottom_lg(iM) = fVal.dif2Bottom_lg;
            re_dif2Last_lg(iM) = fVal.dif2Last_lg;
            re_bottom(iM) = fVal.bottom;
            re_next(iM) = fVal.next;
            re_bottom_orig(iM) = fVal.bottom_orig;
            re_next_orig(iM) = fVal.next_orig;
            re_bottom_prob(iM) = fVal.bottom_prob;
            re_next_prob(iM) = fVal.next_prob;
            re_dif2Bottom(iM) = fVal.dif2Bottom;
            re_dif2Last(iM) = fVal.dif2Last; 
        end
    end
    op.re_notEnoughInformation = re_notEnoughInformation;
    op.re_bottom_orig = re_bottom_orig; 
    op.re_next_orig = re_next_orig; 
    op.re_bottom = re_bottom; 
    op.re_next = re_next; 
    op.re_isNarrowLayer = re_isNarrowLayer; 
    op.re_isBottomResolutionDependent = re_isBottomResolutionDependent;
    op.re_isNextBottom = re_isNextBottom; 
    op.re_interpolated = re_interpolated;
    op.re_bottom_prob = re_bottom_prob;
    op.re_next_prob = re_next_prob;
    op.re_dif2Bottom = re_dif2Bottom;
    op.re_dif2Last = re_dif2Last;
    op.re_bottom_lg = re_bottom_lg;
    op.re_next_lg = re_next_lg;
    op.re_dif2Bottom_lg = re_dif2Bottom_lg;
    op.re_dif2Last_lg = re_dif2Last_lg;

end
