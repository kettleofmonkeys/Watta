function op = set_single_depth_value_re( iM, photon_count  )
    global g_re_isNarrowLayer ...
        g_re_isBottomResolutionDependent g_re_isNextBottom g_re_interpolated ...
        g_re_top g_re_bottom g_re_next g_re_subsurf...
        g_re_top_orig g_re_bottom_orig g_re_next_orig ...
        g_re_top_prob g_re_bottom_prob g_re_next_prob ...
        g_re_dif2Top g_re_dif2Bottom g_re_dif2Last ...
        g_re_top_lg g_re_bottom_lg g_re_next_lg g_re_dif2Top_lg g_re_dif2Bottom_lg g_re_dif2Last_lg ...
        g_re_notEnoughInformation

    % primary function where we determine a depth
    fVal = set_single_depth_value_sub( iM, photon_count  );
    g_re_notEnoughInformation(iM) = fVal.notEnoughInformation;
    g_re_isNarrowLayer(iM) = fVal.isNarrowLayer;    
    g_re_isBottomResolutionDependent(iM) = fVal.isBottomResolutionDependent;
    g_re_isNextBottom(iM) = fVal.isNextBottom;
    g_re_interpolated(iM) = fVal.interpolated;
    g_re_top_lg(iM) = fVal.top_lg;
    g_re_bottom_lg(iM) = fVal.bottom_lg;
    g_re_next_lg(iM) = fVal.next_lg;
    g_re_dif2Top_lg(iM) = fVal.dif2Top_lg;
    g_re_dif2Bottom_lg(iM) = fVal.dif2Bottom_lg;
    g_re_dif2Last_lg(iM) = fVal.dif2Last_lg;
    g_re_subsurf(iM) = fVal.subsurf;
    
    g_re_top(iM) = fVal.top;
    g_re_bottom(iM) = fVal.bottom;
    g_re_next(iM) = fVal.next;
    g_re_top_orig(iM) = fVal.top_orig;
    g_re_bottom_orig(iM) = fVal.bottom_orig;
    g_re_next_orig(iM) = fVal.next_orig;
    g_re_top_prob(iM) = fVal.top_prob;
    g_re_bottom_prob(iM) = fVal.bottom_prob;
    g_re_next_prob(iM) = fVal.next_prob;
    g_re_dif2Top(iM) = fVal.dif2Top;
    g_re_dif2Bottom(iM) = fVal.dif2Bottom;
    g_re_dif2Last(iM) = fVal.dif2Last; 

end