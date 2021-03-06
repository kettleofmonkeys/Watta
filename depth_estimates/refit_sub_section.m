
function op = refit_sub_section(st,ed)

    global g_heights g_St g_Ed ... 
        g_isNarrowLayer g_re_isNarrowLayer ...
        g_isBottomResolutionDependent g_re_isBottomResolutionDependent ...
        g_isNextBottom g_re_isNextBottom ...
        g_interpolated g_re_interpolated ...
        g_top g_re_top ...
        g_bottom g_re_bottom ...
        g_next g_re_next ...
        g_subsurf g_re_subsurf...
        g_top_orig g_re_top_orig ...
        g_bottom_orig g_re_bottom_orig ...
        g_next_orig g_re_next_orig ...
        g_top_prob g_re_top_prob ...
        g_bottom_prob g_re_bottom_prob ...
        g_next_prob g_re_next_prob ...
        g_dif2Top g_re_dif2Top ...
        g_dif2Bottom g_re_dif2Bottom ...
        g_dif2Last g_re_dif2Last ...
        g_top_lg g_re_top_lg ...
        g_bottom_lg g_re_bottom_lg ... 
        g_next_lg g_re_next_lg  ....
        g_dif2Top_lg g_re_dif2Top_lg ...
        g_dif2Bottom_lg g_re_dif2Bottom_lg ...
        g_dif2Last_lg g_re_dif2Last_lg ...
        g_notEnoughInformation g_re_notEnoughInformation ...
        g_surfaceStrength g_re_surfaceStrength ...        
        g_TopSmooth g_re_TopSmooth ...
        g_BottomSmooth g_re_BottomSmooth ...
        g_top_density_ph g_re_top_density_ph...
        g_unreliablePhDens_surf g_re_unreliablePhDens_surf ...
        g_impossibleHeight g_re_impossibleHeight ...
        g_weight_dist_top g_re_weight_dist_top ...
        g_threshold_spreadDens g_re_threshold_spreadDens ...
        g_top_new  g_re_top_new ...
        g_top_new_interp g_re_top_new_interp  ...
        g_bottom_new g_re_bottom_new ...
        g_next_new g_re_next_new ...
        g_altered_top g_re_altered_top ...
        g_subsurf_new g_re_subsurf_new ...
        g_acceptable_heights g_re_acceptable_heights

        g_isNarrowLayer(st:ed) =  g_re_isNarrowLayer;
        g_isBottomResolutionDependent(st:ed) =   g_re_isBottomResolutionDependent;
        g_isNextBottom(st:ed) =   g_re_isNextBottom;
        g_interpolated(st:ed) =   g_re_interpolated;
        g_bottom(st:ed) =   g_re_bottom;
        g_next(st:ed) =   g_re_next;
        g_bottom_orig(st:ed) =   g_re_bottom_orig;
        g_next_orig(st:ed) =   g_re_next_orig;
        g_bottom_prob(st:ed) =   g_re_bottom_prob;
        g_next_prob(st:ed) =   g_re_next_prob;
        g_dif2Bottom(st:ed) =   g_re_dif2Bottom;
        g_dif2Last(st:ed) =   g_re_dif2Last;
        g_bottom_lg(st:ed) =   g_re_bottom_lg; 
        g_next_lg(st:ed) =   g_re_next_lg;
        g_dif2Bottom_lg(st:ed) =   g_re_dif2Bottom_lg;
        g_dif2Last_lg(st:ed) =   g_re_dif2Last_lg;
        g_notEnoughInformation(st:ed) =   g_re_notEnoughInformation;

end

