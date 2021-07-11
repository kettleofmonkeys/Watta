function op = post_produce_final_estimates_bob()
    %interpret lake depth estimates from lake surface and bottom, correcting for refraction 
    global g_lakesurface_final g_subsurf_ice g_surface_ice g_bottom_final ...
        g_depth_all g_depth_real g_depth_corr_all g_depth_corr_real ...
        g_dist g_x_dist g_BottomSmooth g_BottomSmooth_bin_fix g_Flat_Change ...
        g_bottom_final g_bottom_final_real  g_bottom_refit_edges_bin ...
        g_iceLayer g_interpolated g_next_new ...
        g_heights_orig g_next g_bottom g_surface_ice_bin g_iceLayer_bin g_Flat_Change_bin ...
        g_subsurf_ice_bin g_Next_bin g_surfaceType g_lakebottom_final_splined g_top_new_interp g_lakes
    
    g_surfaceType = g_interpolated == 1;          % this will be altered. Intermediate values
    g_lakesurface_final = g_top_new_interp;
    g_lakesurface_final(g_lakes==0)=nan;
    
    g_bottom_final = g_BottomSmooth;
    g_bottom_final_real = g_bottom_final;
    g_bottom_final_real(g_surfaceType == 1) = NaN; %limits to where bottoms were actually found without interpolation
    
    g_iceLayer(isnan(g_iceLayer)) = 0;
    
    %%% compute lake depth, apply correction for refraction (g_depth_corr_all)
    g_depth_all = g_lakesurface_final - g_bottom_final;
    g_depth_corr_all = g_bottom_final + 0.25416*g_depth_all;
    
    %%% compute lake depth, apply correction for refraction, without
    %%% interpolated parts
    g_depth_real = g_lakesurface_final - g_bottom_final_real;
    g_depth_corr_real = g_bottom_final_real + 0.25416*g_depth_real;

    %%% surface ice (only where ice layer was detected) 
    g_surface_ice = g_lakesurface_final;
    g_surface_ice(~g_iceLayer) = NaN;

    %%% subsurface ice (> bottom)
    g_subsurf_ice = g_next_new;
    g_subsurf_ice(g_next <= g_bottom_final) = NaN; 
    g_subsurf_ice(isnan(g_lakesurface_final)) = NaN;
end

