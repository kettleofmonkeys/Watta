function op = sub_recalculation( startPoint, endPoint )
    % recalculate the bottom within each presumed lake with a smaller
    % window of 30 photons (or as defined in g_photon_count_recalc) 
    global g_recalc g_St g_Ed g_re_top_new_interp g_top_new_interp ...
        g_re_subsurf_new g_subsurf_new

    g_recalc = true;
    g_St = startPoint;
    g_Ed = endPoint;
    de_calc_depth_sub();    
    %note that no outlier removal is happening here
    create_smoothed_pond_mn_re();   %this might really be redundant
    
    g_re_top_new_interp = g_top_new_interp(g_St:g_Ed);
    g_re_subsurf_new = g_subsurf_new(g_St:g_Ed);
    


end

 function op = de_calc_depth_sub()

    global g_photon_count g_re_lats g_re_lons g_heights g_subsetOfHeights ...
        g_re_isNarrowLayer  g_lats g_lons g_dist ...
        g_re_isBottomResolutionDependent g_re_isNextBottom g_re_interpolated ...
        g_re_top g_re_bottom g_re_next g_subsurf_new...
        g_re_top_orig g_re_bottom_orig g_re_next_orig ...
        g_re_top_prob g_re_bottom_prob g_re_next_prob ...
        g_re_dif2Top g_re_dif2Bottom g_re_dif2Last ...
        g_re_top_lg g_re_bottom_lg g_re_next_lg g_re_dif2Top_lg g_re_dif2Bottom_lg g_re_dif2Last_lg ...
        g_re_notEnoughInformation  ...        
        g_St g_Ed g_re_dist g_top_new_interp g_top_filt
    size(g_heights)
	g_St
	g_Ed
	
    photon_set = g_heights(g_St:g_Ed);
    g_re_lats = g_lats(g_St:g_Ed);
    g_re_lons = g_lons(g_St:g_Ed);
    g_re_dist = g_dist(g_St:g_Ed);

    g_subsetOfHeights = photon_set;
    
    photon_count = g_photon_count;
    top_filt = g_top_filt;
    val_St = g_St;
    subheights = g_heights;
    subheights(subheights >= g_subsurf_new) = NaN;  %%% remove top photons, so we can focus on detecting the bottom    
    sizeN = size(photon_set);
    
    top_filt = top_filt;
    val_St = val_St;
    photon_count = photon_count;
    sizeN = sizeN;
    resRel = post_reset_kerndensity_lakes(subheights, top_filt, val_St, photon_count, sizeN);  %%% this performs the kernel density estimate again and estimates the bototm
    
    %%% output
    g_re_notEnoughInformation =    resRel.re_notEnoughInformation;
    g_re_bottom_orig = resRel.re_bottom_orig;
    g_re_next_orig = resRel.re_next_orig;
    g_re_bottom = resRel.re_bottom;
    g_re_next = resRel.re_next;
    g_re_isNarrowLayer = resRel.re_isNarrowLayer;
    g_re_isBottomResolutionDependent = resRel.re_isBottomResolutionDependent;
    g_re_isNextBottom = resRel.re_isNextBottom;
    g_re_interpolated = resRel.re_interpolated;
    g_re_bottom_prob = resRel.re_bottom_prob;
    g_re_next_prob = resRel.re_next_prob;
    g_re_dif2Bottom = resRel.re_dif2Bottom;
    g_re_dif2Last = resRel.re_dif2Last;
    g_re_bottom_lg = resRel.re_bottom_lg;
    g_re_next_lg = resRel.re_next_lg;
    g_re_dif2Bottom_lg = resRel.re_dif2Bottom_lg;
    g_re_dif2Last_lg = resRel.re_dif2Last_lg;

end
