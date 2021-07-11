
function op = post_final_smoothing(idx_left,idx_right)
    global g_re_top g_re_bottom  g_subsetOfHeights g_photon_count g_re_TopSmooth g_re_BottomSmooth ...
    g_re_dist g_heights g_Flat_Change g_St g_Ed g_top_new_interp g_dist g_re_top g_re_bottom lake_ends g_subsetOfHeights g_photon_count g_re_TopSmooth g_re_BottomSmooth g_top_new_interp ...
    g_BottomSmooth g_BottomSmooth_low lake_edges g_bottom_refit_edges g_re_dist g_heights g_Flat_Change g_St g_Ed g_top_new_interp g_top_filt g_dist g_re_top g_re_bottom lake_ends g_subsetOfHeights g_photon_count g_re_TopSmooth g_subsurf_ice...
    
    %%% index of segment 
    idx_seg=(idx_left:idx_right);
    if (~isempty(idx_seg))
        l_seg=abs(range(g_dist(idx_seg))); %%% length of segment
        n_seg=length(idx_seg);             %%% number of photons within segment 
        %%% add max. 2000 surface photons left/right of lake edges for smoothing
        idx_surf_l=max(1,idx_left-2000):max(1,idx_left-1);
        idx_surf_r=min(length(g_top_new_interp),idx_right+1):min(length(g_top_new_interp),idx_right+2000);
        surf_l=g_top_new_interp(idx_surf_l);  %%% surface photons left of lake
        surf_r=g_top_new_interp(idx_surf_r);  %%% surface photons right of lake
        dist_l=g_dist(idx_surf_l);      %%% coordinates of left surfae photons
        dist_r=g_dist(idx_surf_r);      %%% coordinates of right surface photons
    
        %%% start iterative outlier removal
        buffer_bottom_max = min(120,l_seg*0.25);   %%% accept larger outlier in this part near the edges
        copy_bottom_refit_edges=g_bottom_refit_edges(idx_left:idx_right); %%% create copy
	n_iter=5;        
        for j=1:n_iter    
            n_buff_edge=min(round(n_seg*0.10),100);                       %%% use g_top_new_interp for smooting in buffer of 100 photons, or 10% of lake photons
            %%% combine surface photons left of lake, left buffer surface photons,
            %%% bottom photons, right buffer surface photons and surface
            %%% photons right of lake
            subsetOfHeights_bot=[surf_l; g_top_new_interp(idx_seg(1:n_buff_edge));copy_bottom_refit_edges((n_buff_edge+1):(end-n_buff_edge));g_top_new_interp(idx_seg((end-(n_buff_edge-1)):end)); surf_r];
            %%% first create very smooth bottom                
            sZ_bot = size(subsetOfHeights_bot,1);    %%% number of photons to be smoothed
            xV_bot = 1:sZ_bot;                       %%% index of photons to be smoothed
            %%% define smoothing span
            smoothingValue_bot_smooth = (75*10)/length(idx_seg);
            %%% apply smoothing
            subsetOfHeights_bot_smooth = smooth(xV_bot,subsetOfHeights_bot,smoothingValue_bot_smooth,'rloess');
            subsetOfHeights_bot_smooth=subsetOfHeights_bot_smooth((length(idx_surf_l)+1):end-length(idx_surf_r));
            %%% remove bottom points > g_top
            idx_over=find(subsetOfHeights_bot_smooth>g_top_new_interp(idx_seg));
             if (~isempty(idx_over))
                subsetOfHeights_bot_smooth(idx_over)=g_top_new_interp(idx_seg(idx_over));
            end
            %%% set first and and point to g_top_new_interp
            subsetOfHeights_bot_smooth(1)=g_top_new_interp(idx_seg(1));
            subsetOfHeights_bot_smooth(end)=g_top_new_interp(idx_seg(end)); 
            %%% find outliers, based on ratio of difference between original bottom
            %%% and smoothed version, and local depth. Threshold for
            %%% rejection depends on # iteration
            outl=(g_bottom_refit_edges(idx_seg)-subsetOfHeights_bot_smooth)./(g_top_new_interp(idx_seg)-subsetOfHeights_bot_smooth);
            idx_out=find(abs(outl)>max(0.5-(j*0.05),0.25));
            if (~isempty(idx_out))
              idx_out_seg=idx_seg(idx_out);  %%% index of outliers in global array
              %%% near edges, accept outliers within 0.5-1 m (value depends
              %%% on # iteration)
              keep_left=find(g_dist(idx_out_seg)<min(g_dist(idx_seg(1))+buffer_bottom_max,g_dist(end))&abs(g_bottom_refit_edges(idx_out_seg)-subsetOfHeights_bot_smooth(idx_out))<min(0.5+(n_iter-j)*.25,1));   %%% if near edges, don't reject these outliers
              keep_right=find(max(g_dist(idx_seg(end))-buffer_bottom_max,g_dist(1))<g_dist(idx_out_seg)&abs(g_bottom_refit_edges(idx_out_seg)-subsetOfHeights_bot_smooth(idx_out))<min(0.5+(n_iter-j)*.25,1));
              if (~isempty(keep_right))     
                [~,pos_keep_right]=intersect(idx_out,keep_right);   %%% index pf photons to keep within idx_out
                idx_out(pos_keep_right)=[];                           %%% remove from idx_out
              end
              if (~isempty(keep_left))
                [~,pos_keep_left]=intersect(idx_out,keep_left);   %%% index pf photons to keep within idx_out
                idx_out(pos_keep_left)=[];                          %%% remove from idx_out
              end
 %           idx_out_seg=idx_seg(idx_out);                          %%% remove outliers
             copy_bottom_refit_edges(idx_out)=nan;
           end
            %%% in copy_bottom_refit_edges,replace photons if subsurf_ice is
            %%% closer to smoothed estimate
            %%% in first and last 10%, use subsurf_ice if higher than copy_bottom_refit_edges
            diff_sub_bot=g_subsurf_ice(idx_seg)-subsetOfHeights_bot_smooth;
            idx_repl=find(abs(diff_sub_bot)<1&abs(diff_sub_bot)<abs(copy_bottom_refit_edges-subsetOfHeights_bot_smooth));
            idx_repl_left=find(g_dist(idx_seg)<(g_dist(idx_seg(1))+0.10*l_seg)&0<diff_sub_bot&diff_sub_bot<1);
            idx_repl_right=find(g_dist(idx_seg)>g_dist(idx_seg(end))-0.10*l_seg&0<diff_sub_bot&diff_sub_bot<1);
            copy_bottom_refit_edges(idx_repl)=g_subsurf_ice(idx_seg(idx_repl));
            copy_bottom_refit_edges(idx_repl_left)=g_subsurf_ice(idx_seg(idx_repl_left));
            copy_bottom_refit_edges(idx_repl_right)=g_subsurf_ice(idx_seg(idx_repl_right));
        end        
      %%% final smoothing after outlier removal
      subsetOfHeights_bot=[surf_l; g_top_new_interp(idx_seg(1:n_buff_edge));copy_bottom_refit_edges((n_buff_edge+1):(end-n_buff_edge));g_top_new_interp(idx_seg((end-(n_buff_edge-1)):end)); surf_r];
      %%% start smoothing
      sZ_bot = size(subsetOfHeights_bot,1); %%% number of photons to be smoothed
      xV_bot = 1:sZ_bot;                      %%% index of photons to be smoothed
      %%%use 2 smoothing spans    
      smoothingValue_bot3 = (g_photon_count*3)/length(idx_seg);   
      smoothingValue_bot5 = (g_photon_count*5)/length(idx_seg);
      %%% apply smoothing
      g_re_BottomSmooth3 = smooth(xV_bot,subsetOfHeights_bot,smoothingValue_bot3,'rloess');
      g_re_BottomSmooth5 = smooth(xV_bot,subsetOfHeights_bot,smoothingValue_bot5,'rloess');
      %%% remove surface photons left and right of lake
      g_re_BottomSmooth3=g_re_BottomSmooth3((length(idx_surf_l)+1):end-length(idx_surf_r));
      g_re_BottomSmooth5=g_re_BottomSmooth5((length(idx_surf_l)+1):end-length(idx_surf_r));
      %%% remove bottom points > g_top
      idx_over=find(g_re_BottomSmooth3>g_top_new_interp(idx_seg));
      if (~isempty(idx_over))
            g_re_BottomSmooth3(idx_over)=g_top_new_interp(idx_seg(idx_over));
      end
      idx_over=find(g_re_BottomSmooth5>g_top_new_interp(idx_seg));
      if (~isempty(idx_over))
          g_re_BottomSmooth5(idx_over)=g_top_new_interp(idx_seg(idx_over));
      end 
      %%% set first and and point to g_top_new_interp
      g_re_BottomSmooth3(1)=g_top_new_interp(idx_seg(1));
      g_re_BottomSmooth3(end)=g_top_new_interp(idx_seg(end));
      g_re_BottomSmooth5(1)=g_top_new_interp(idx_seg(1));
      g_re_BottomSmooth5(end)=g_top_new_interp(idx_seg(end)); 

      %%% export smoothed bottom
      g_BottomSmooth(idx_left:idx_right) =  g_re_BottomSmooth3;
      if (isempty(isfinite(g_BottomSmooth_low)))        %%% create g_BottomSmooth_low if it doesn't exist yet
                g_BottomSmooth_low=g_BottomSmooth;
      end
        g_BottomSmooth_low(idx_left:idx_right) = g_re_BottomSmooth5;
    else
        g_BottomSmooth_low = g_BottomSmooth; 
    end
end
