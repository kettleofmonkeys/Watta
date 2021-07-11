function op = set_single_depth_value( iM, photon_count  )
    global g_isNarrowLayer g_isChangeSlope ...
        g_isBottomResolutionDependent g_isNextBottom g_interpolated ...
        g_top g_bottom g_next g_subsurf...
        g_top_orig g_bottom_orig g_next_orig ...
        g_top_prob g_bottom_prob g_next_prob ...
        g_dif2Top g_dif2Bottom g_dif2Last ...
        g_top_lg g_bottom_lg g_next_lg g_dif2Top_lg g_dif2Bottom_lg g_dif2Last_lg ...
        g_notEnoughInformation
    % This function simply calls set_singe_depth_value, the primary function where we determine the elevation of the water/ice surface and bottom, and stores its output in the global variables 
    fVal = set_single_depth_value_sub( iM, photon_count  );
    g_notEnoughInformation(iM) = fVal.notEnoughInformation;
    g_isNarrowLayer(iM) = fVal.isNarrowLayer;    
    g_isBottomResolutionDependent(iM) = fVal.isBottomResolutionDependent;
    g_isNextBottom(iM) = fVal.isNextBottom;
    g_interpolated(iM) = fVal.interpolated;
    g_top_lg(iM) = fVal.top_lg;
    g_bottom_lg(iM) = fVal.bottom_lg;
    g_next_lg(iM) = fVal.next_lg;
    g_dif2Top_lg(iM) = fVal.dif2Top_lg;
    g_dif2Bottom_lg(iM) = fVal.dif2Bottom_lg;
    g_dif2Last_lg(iM) = fVal.dif2Last_lg;
    g_subsurf(iM) = fVal.subsurf;
    
    g_top(iM) = fVal.top;
    g_bottom(iM) = fVal.bottom;
    g_next(iM) = fVal.next;
    g_top_orig(iM) = fVal.top_orig;
    g_bottom_orig(iM) = fVal.bottom_orig;
    g_next_orig(iM) = fVal.next_orig;
    g_top_prob(iM) = fVal.top_prob;
    g_bottom_prob(iM) = fVal.bottom_prob;
    g_next_prob(iM) = fVal.next_prob;
    g_dif2Top(iM) = fVal.dif2Top;
    g_dif2Bottom(iM) = fVal.dif2Bottom;
    g_dif2Last(iM) = fVal.dif2Last; 

end