    function op = post_find_edges_of_lake_newest_bob()
    % defines individual lakes (connected surfaces) from flat lakes
    global g_heights g_dist g_index_change_slope g_lake_score g_classes g_Flat_Change g_top_new_interp g_top_filt
    %%% refine the breakpoints of the lakes considered (defined by
    %%% g_classes) so that they better match the edges of these lakes 

    fin_edges = nan(size(g_heights));
    
    %%% find segments of lake surface, assume that segments separated by less
    %%% than X meter are still connected
    [idx_seg,~]=find(g_classes==g_lake_score);  %%% only use the lakes with the classes defined in g_classes; returns the index of the photons which belong to these lakes    
    if (length(idx_seg)>0)
    %%% ascending/descending distance (g_dist)
    if (g_dist(1)>g_dist(end))   % descending
     [idx_breaks_right,~,idx_change_class]=intersect(idx_seg,g_index_change_slope); % common points in idx_seg and g_index_change_slope, which contains the locations of the breakpoints (based on slope now), including first and last photon
     if (idx_change_class(end)==length(g_index_change_slope))
      idx_breaks_left=g_index_change_slope(idx_change_class(1:end-1)+1);
      idx_breaks_right=idx_breaks_left(1:end-1);
     else
      idx_breaks_left=g_index_change_slope(idx_change_class(1:end)+1);
     end
    else
     [idx_breaks_left,~,idx_change_class]=intersect(idx_seg,g_index_change_slope);
     if (idx_change_class(end)==length(g_index_change_slope))
      idx_breaks_right=g_index_change_slope(idx_change_class(1:end-1)+1);
      idx_breaks_left=idx_breaks_left(1:end-1);
     else
       idx_breaks_right=g_index_change_slope(idx_change_class+1);
     end
    end
    
    
max_sep=17;     %%% merge lakes separated by less than 17 m
sep_surf=g_dist(idx_breaks_left(2:(end)))-g_dist(idx_breaks_right(1:(end-1))); %%% distance betweenlakes
idx_connect=find(sep_surf<max_sep);
idx_breaks_right(idx_connect)=[];   %%% merge
idx_breaks_left(idx_connect+1)=[];  %%% merge




%%% refine edges (iteratively, if n_iter > 1)
g_Flat_Change_post=g_Flat_Change;
idx_breaks_left_new=idx_breaks_left;
idx_breaks_right_new=idx_breaks_right;
perc_lim_top=90; %%% used to define upper part of lake surface elevation
perc_lim_top_out=98;   %%% futher outlier removal for lake surface

n_iter=1
for j = 1:n_iter
 if (j>1)
     idx_breaks_left=idx_breaks_left_new;
     idx_breaks_right=idx_breaks_right_new
   for i=1:length(idx_breaks_left) 
        %%% left
        [i , length(idx_breaks_left)]
        idx_seg=idx_breaks_left(i):idx_breaks_right(i);
        ph_seg=g_top_filt(idx_breaks_left(i):idx_breaks_right(i));  %%% top photons of lake
        n_ph_seg=length(idx_seg);
        max_gap=min(100,round(n_ph_seg*0.25)); %% max number of photons for gaps to be connected

        %%% find edges which are actually not part of lake
        lim_top_out=prctile(ph_seg,perc_lim_top_out);               %%% find rough estimate of elevation representative for top (excluding potential outliers
        lim_top=prctile(ph_seg(ph_seg<=lim_top_out),perc_lim_top);  %%% define threshold for edge refinement 

        idx_edges=find(g_top_filt(idx_breaks_left(i):idx_breaks_right(i))>lim_top);  %%% identify photons above threshold
        idx_gap=find(diff(idx_edges)>max_gap);                                       %%% gaps between photons above threshold
        %%% if one of first 10 photons is above limit, this is probably not
        %%% lake surface and should be removed together with neigbouring
        %%% photons if gap doesn't exceed the limit
        if (sum(idx_edges==1:10)>0)   
            if (~isempty(idx_gap))
             idx_breaks_left_new(i)=idx_seg(idx_edges(idx_gap(1)));     %%% stop removing at photon where too large a gap occurs
            else
                if (idx_edges(end)+1<length(idx_seg))                   %%% makes sure lake isn't too small
                    idx_breaks_left_new(i)=idx_seg(idx_edges(end)+1);
                end
            end
            g_Flat_Change_post(idx_breaks_left(i):idx_breaks_left_new(i)-1)=NaN;
        end
        %%% right side, if one of last 10 photons is above limit, remove 
        if (sum(idx_edges==(n_ph_seg-9):n_ph_seg)>0)  
            if (~isempty(idx_gap))
                idx_breaks_right_new(i)=idx_seg(idx_edges(idx_gap(end)+1));
            else
                if (idx_edges(1)-1>1)
                    idx_breaks_right_new(i)=idx_seg(idx_edges(1)-1);
                end
            end
            g_Flat_Change_post(idx_breaks_right_new(i)+1:idx_breaks_right(i))=NaN;
        end
        idx_seg=idx_breaks_left_new(i):idx_breaks_right_new(i);
        %%% remove remaining extreme positive outliers in g_top_new_interp and g_top_filt. Replace with median
        %%% top value in neighbourhood
        idx_out=find(g_top_new_interp(idx_seg)>lim_top_out);
        if (length(idx_out)<0.80*length(idx_seg))           
         g_top_new_interp(idx_seg(idx_out))=nan;
          for (k = 1:length(idx_out))
              l=1;
              while isnan(g_top_new_interp(idx_seg(idx_out(k))))
               g_top_new_interp(idx_seg(idx_out(k)))=nanmedian(g_top_new_interp(idx_seg(max(idx_out(k)-200*l,1)):idx_seg(min(idx_out(k)+200*l,length(idx_seg)))));
                l=l+1;
              end
          end
        end
        idx_out=find(g_top_filt(idx_seg)>lim_top_out);
         if (length(idx_out)<0.80*length(idx_seg)) 
             g_top_filt(idx_seg(idx_out))=nan;
                for (k = 1:length(idx_out))
                  l=1;
                   while isnan(g_top_filt(idx_seg(idx_out(k))))
                       g_top_filt(idx_seg(idx_out(k)))=nanmedian(g_top_filt(idx_seg(max(idx_out(k)-200*l,1)):idx_seg(min(idx_out(k)+200*l,length(idx_seg)))));
                         l=l+1;
                   end
                end
          end
end
 end  
    %%% update variables
    g_Flat_Change=g_Flat_Change_post;
    idx_breaks_left = idx_breaks_left_new;
    idx_breaks_right = idx_breaks_right_new;
    op.idx_breaks_left = idx_breaks_left_new;
    op.idx_breaks_right = idx_breaks_right_new;
    %%% set bottom values not in lake to NaN
    

end
else
    op.idx_breaks_left=[];  
    op.idx_breaks_right=[]; 
end
end
