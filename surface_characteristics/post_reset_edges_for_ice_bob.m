function op = post_reset_edges_for_ice_bob(lake_edges )
    % Where subsurface ice layers are near the edges of the lake, 
    % we assume  these are linked to the lake itself and create a 
    % smoothed transition from the bottom to the surface   
    global g_dist g_subsurf_ice g_bottom_new g_Flat_Change g_next_new ...
        g_interpolated g_bottom_refit_edges g_heights_orig g_subsurf_new  

    g_subsurf_ice = g_next_new;
    g_subsurf_ice(g_next_new <= g_bottom_new) = NaN;
    g_subsurf_ice(g_Flat_Change==0) = NaN;

    idx_breaks_left = lake_edges.idx_breaks_left;
    idx_breaks_right = lake_edges.idx_breaks_right;

    tAdjustEdges = g_bottom_new;
    
    for i = 1:size(idx_breaks_left)
        length_lake=abs(g_dist(idx_breaks_left(i))-g_dist(idx_breaks_right(i)));
        buffer_bottom_max = min(110,length_lake*0.25);      %%% only consider the edges of the lake, either 110 m or 25% of the lake width
        buffer_bottom_min = -1*buffer_bottom_max;           %%% same, for right side
        buffer_bottom_sections = 10;                        %%% look for subsurface ice in steps of 10 m
    
        foundEdge = false;
        mNmDpth = g_subsurf_new(idx_breaks_left(i));         %start with a "depth" at the top
        buffer_bottom=0; 
        while buffer_bottom<buffer_bottom_max && foundEdge == false   %%% search until end of buffer is reached, or an edge is found
            idx_nan_left=find(g_dist>=  g_dist(idx_breaks_left(i))+buffer_bottom...
                & g_dist<=g_dist(idx_breaks_left(i))+buffer_bottom+buffer_bottom_sections );% indices of photons within the segment of width buffer_bottom_sections            
            if(~isempty(idx_nan_left))                
                mNmDpth_n = nanmean(g_subsurf_ice(idx_nan_left));   %%% if photons found in the segment, compute mean elevation of subsurf_ice photons
                if(~isnan(mNmDpth_n))
                    if( mNmDpth_n >= mNmDpth)   %If this average is higher than the last one -> edge found
                        foundEdge = true;    
                    else                
                        mNmDpth = mNmDpth_n;      %if it's lower than the previous iteration, reset and keep trucking
                        buffer_bottom = buffer_bottom + buffer_bottom_sections;
                    end
                else   %no subsurf ice found. we could either stop (if the sections in buffer_bottom_sections are wide enough
                       %or we coudld keep searching. I'm going to stop.
                    foundEdge = true; % actually no subsurface ice was found..   
                end                
            else   
                buffer_bottom = buffer_bottom+buffer_bottom_sections;   %%% no photons in segment, go to next one
            end            
        end     
        idx_nan_left=find(g_dist>=  g_dist(idx_breaks_left(i))...
                & g_dist<=g_dist(idx_breaks_left(i))+buffer_bottom ); %%% indices of photons belonging to left edge
        
        %%% same for right edge    
        foundEdge = false;
        mNmDpth = g_subsurf_new(idx_breaks_right(i));         
        buffer_bottom=0; 
        while buffer_bottom>buffer_bottom_min && foundEdge == false
            idx_nan_right=find(g_dist<=  g_dist(idx_breaks_right(i))- buffer_bottom ...
                & g_dist>=g_dist(idx_breaks_right(i))-buffer_bottom - buffer_bottom_sections );% ...
            if(~isempty(idx_nan_right))                
                mNmDpth_n = nanmean(g_subsurf_ice(idx_nan_right));
                if(~isnan(mNmDpth_n))
                    if( mNmDpth_n >= mNmDpth)
                        foundEdge = true;   
                    else                
                        mNmDpth = mNmDpth_n;      
                        buffer_bottom = buffer_bottom + buffer_bottom_sections;
                    end
                else  
                    foundEdge = true;    
                end                
            else   
                buffer_bottom = buffer_bottom+buffer_bottom_sections;
            end            
        end
        idx_nan_right=find(g_dist<=  g_dist(idx_breaks_right(i))...
                & g_dist>=g_dist(idx_breaks_right(i))-buffer_bottom ); %%% indices of photons belonging to right edge
        

     
        tLeft = max(g_subsurf_ice(idx_nan_left), tAdjustEdges(idx_nan_left));
        tRight = max(g_subsurf_ice(idx_nan_right), tAdjustEdges(idx_nan_right));

        tAdjustEdges(idx_nan_left) = tLeft;%g_subsurf_ice_bin(idx_nan_left);
        tAdjustEdges(idx_nan_right) = tRight;%g_subsurf_ice_bin(idx_nan_right);
 
    end

    g_bottom_refit_edges = tAdjustEdges;

    
end