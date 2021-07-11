function op = post_detect_ice_layer_bob(lake_edges )
    %find subsurface ice layers in each lake    
    global  g_iceLayer g_Flat_Change
    
    g_iceLayer = nan(size(g_Flat_Change));
    for i = 1:size(lake_edges.idx_breaks_left)
        leftT = lake_edges.idx_breaks_left(i);
        rightT = lake_edges.idx_breaks_right(i);
        [found, iceLayer] = post_find_ice_surf_newest( leftT, rightT);
        if( found )
            g_iceLayer(leftT:rightT) = iceLayer; 
        end
    end
    
    g_iceLayer(isnan(g_iceLayer)) = 0;
    
end
