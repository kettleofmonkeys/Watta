function op = post_evaluate_pond_bottom_density_bob()
    evaluate_bottom_pond_density(); %%% Compare the difference between the local photon density at any given pond lake bottom estimate and that of the background rate/density around the mean surface and produce a weighted likelihood that  this is a lake bottom 
    calculate_new_bottom();  %%% Interpolate a new lake bottom when the rating at the photon location is too  unreliable.
end

function op = calculate_new_bottom()
    % given weights as determined in evaluate_bottom_pond_density, interpolate 
    % a new lake bottom when the rating at the photon location is too unreliable
    global g_weight_dens_bottom g_impossibleDepth g_bottom_revision_required ...
    g_bottom_new g_heights_orig g_BottomSmooth g_photon_count
    
    %%% remove unreliable bottom estimate
    g_bottom_revision_required = g_weight_dens_bottom < 0.3;   %kill these, but bottoms can be really light on density
    g_bottom_revision_required(g_impossibleDepth == 1) = 1;
    g_bottom_new(g_bottom_revision_required == 1) = NaN;
    
    
    allCount = size(g_bottom_new,1);
    if sum(~isnan(g_bottom_new)) < 2
        return;         %%% not enough photons remaining for interpolation
    end
    
    aX_all = 1:allCount;
    aX_or = aX_all;

    %%% remove nans
    aX_or(isnan(g_bottom_new)) = [];
    or_Bot = g_bottom_new;
    or_Bot(isnan(or_Bot)) = [];    
    
    allCount = size(or_Bot,1);

    if sum(~isnan(or_Bot)) < 2
        return;
    end
    g_bottom_new= interp1(aX_or, or_Bot, aX_all);
    g_bottom_new = g_bottom_new';

end

function op = evaluate_bottom_pond_density()
    % Compare the difference between the local photon density at any 
    % given pond lake bottom estimate and that of the background rate/density
    % around the  mean surface, as well as the distance from the photon height 
    % to the local mean height of surrounding photons, to arrive at a weighted
    % likelihood that this is a lake bottom
    
    global  g_bottom_new g_heights_orig g_bottom_density_ph...
       g_surf_dens  g_background_rate g_binSurf g_unreliablePhDens_bottom ...
       g_weight_dens_bottom g_impossibleDepth g_top_new_interp g_subsurf_new ...
       g_threshold_spreadDens g_bottom
   
    lHeights = size(g_heights_orig,1);         %%% total number of photons
    %%% define arrays
    g_unreliablePhDens_bottom = zeros(size(g_heights_orig));
    g_impossibleDepth = zeros(size(g_heights_orig));
    g_weight_dens_bottom = nan(size(g_heights_orig));
    g_bottom_density_ph = nan(size(g_heights_orig));
    lSz = size(g_bottom_new,1);
    
    wideBerth = 75*2;                         % use 75 photons left/right of centre photon          
    for iT = 1:lSz
        curVal_bot = g_bottom_new(iT);
        aTop_Val = g_top_new_interp(iT);
        %%% compare bottom and top value
        if( curVal_bot >= aTop_Val + 1)       %this means it was miscalculated horribly
            g_impossibleDepth(iT) = 1;
            continue;
        end

        %%% calculate photon density in 2m - band around estimated bottom 
        startOrig = max(iT - wideBerth,1);
        endOrig = min(iT + wideBerth, lHeights);
        abs_top = repmat(curVal_bot+2, 1, endOrig - startOrig+1);
        topOfDensity_est = min(abs_top, g_subsurf_new(startOrig:endOrig)'); %%% only use photons below the lake surface
        curVal_td = post_calculate_phot_density(curVal_bot-2, topOfDensity_est', startOrig, endOrig);  %%% calculate photon density around bottom photon
      
        g_bottom_density_ph(iT) = curVal_td;
        spreadDens = g_surf_dens(iT) - g_background_rate(iT);   %%% difference between photon density around top and background photon density
        
        %%% compare to surface and background density in wide neighbourhood
        minT = max(1, iT-3000); maxT = min(size(g_surf_dens,1),iT+3000);
        threshV = g_surf_dens(minT:maxT) - g_background_rate(minT:maxT);
        threshold_spreadDens = max(0, nanmean(threshV)-2*nanstd(threshV));  %%% define threshold as mean - 2 standard deviations of difference in surface and background density in wide neighbourhood
        g_unreliablePhDens_bottom(iT) = spreadDens < threshold_spreadDens || threshold_spreadDens <0.01 ; %if the difference in density is small  -> unreliable bottom estimate

        %%% assign a rank, stored in g_weight_dens_bottom
        rank_ph_dens = max(0, ((curVal_td - g_background_rate(iT))/spreadDens)) * 4; %the higher the rate, the higher the likelihood its real
        rank_ph_dens = min(5,rank_ph_dens);
        g_weight_dens_bottom(iT) = rank_ph_dens;
            
    end
    
    
    
    
    
end