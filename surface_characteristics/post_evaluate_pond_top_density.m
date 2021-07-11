function op = post_evaluate_pond_top_density()
    % Compare the difference between the local photon density at any given pond top estimate 
    % and that of the background rate/density around the mean surface, as well as  the distance 
    % from the photon height to the local mean height of surrounding pho tons, to arrive at a weighted 
    % likelihood that this is a lake surface vs. background noise
    % Assign g_surfaceType to 1 (surface alone) 2(lake) 3(ice layer)
    global   g_dist g_top g_heights_orig g_top_density_ph...
       g_surf_dens  g_background_rate g_binSurf g_unreliablePhDens_surf ...
       g_impossibleHeight g_weight_dist_top g_weight_dens_top g_threshold_spreadDens
   
%%% clean any remaining extreme outliers in g_heigths
  perc_lim_low=1;  %%% used to define upper part of lake surface elevation
  perc_lim_hi=99;  %%% futher outlier removal for lake surface
  length_track=range(g_dist);
  length_seg=2e3;               %%% look for outliers in 2 km segments
  n_seg_track=ceil(length_track/length_seg);
  for i = 1:n_seg_track
    idx_ph_seg=find(g_dist>=(i-1)*length_seg&g_dist<i*length_seg);
    %%% find percentiles in segment
    lim_hi=prctile(g_heights_orig(idx_ph_seg),perc_lim_hi);
    lim_low=prctile(g_heights_orig(idx_ph_seg),perc_lim_low);
    idx_hi=find(g_heights_orig(idx_ph_seg)>lim_hi);
    idx_lo=find(g_heights_orig(idx_ph_seg)<lim_low);
    g_heights_orig(idx_ph_seg(idx_hi))=nan;
    g_heights_orig(idx_ph_seg(idx_lo))=nan;
  end
  
  
  %%% define output arrays
  lHeights = size(g_heights_orig,1);
  g_unreliablePhDens_surf = zeros(size(g_heights_orig));
  g_impossibleHeight = zeros(size(g_heights_orig));
  g_weight_dist_top = nan(size(g_heights_orig));
  g_weight_dens_top = nan(size(g_heights_orig));
  g_top_density_ph = nan(size(g_heights_orig));
    
  %%% loop through all photons 
    wideBerth = 75*2;     %%% compute photon density in bins of 75*2 photons around center photon      
    range_top = 2;        %%% compute photon density in 2 m vertical range around g_top and compare this to background photon density
    g_threshold_spreadDens = 0.075;   %%% minimum required difference between surface photon density and background density
    for iT = 1:size(g_top,1)  
        curVal_top = g_top(iT);
        startOrig = max(iT - wideBerth,1);
        endOrig = min(iT + wideBerth, lHeights);
        curVal_td = post_calculate_phot_density(curVal_top-range_top, curVal_top+range_top, startOrig, endOrig);   %%% computes photon density around top
        g_top_density_ph(iT) = curVal_td;
        spreadDens = g_surf_dens(iT) - g_background_rate(iT);       %%% difference between surface photon density and background photon density
        g_unreliablePhDens_surf(iT) = spreadDens< g_threshold_spreadDens; %if the difference is essentially negligible
        rank_ph_dens = max(0, ((curVal_td - g_background_rate(iT))/spreadDens)) * 4; %the higher the rate, the higher the likelihood its real
        rank_ph_dens = min(5,rank_ph_dens);
        
        if( curVal_top > (g_binSurf(iT)- 10) && curVal_top < (g_binSurf(iT)+ 10))
        %if the current top value is within a 10 m rangeof the peak in the histogram, then it's given high rank
            rank_distance = 5;
        else 
        % else it gets a zero rank
            if( curVal_top > g_binSurf(iT)+10)
                maxV = g_binSurf(iT)+20;
                distV = (maxV - curVal_top) / 10;
                if( distV < 1)
                    g_impossibleHeight(iT) = 1;
                    rank_distance = 0;
                else
                    rank_distance = distV * 4;   
                end
            else
                minV = g_binSurf(iT)-10;
                distV = (curVal_top - minV) / 10;
                if( distV < 1)
                    g_impossibleHeight(iT) = 1;
                    rank_distance = 0;
                else
                    rank_distance = distV * 4;   
                end

            end
            
        end 
        g_weight_dist_top(iT) = rank_distance;
        g_weight_dens_top(iT) = rank_ph_dens;        
    end
end