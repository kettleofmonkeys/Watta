function op = post_calculate_lake_prob()
% falgs lakes where the probability is too low (for example when the photon density is insucient or the bottom is too shallow
    global g_weight_dens_bottom g_impossibleDepth ...
        g_impossibleHeight g_weight_dist_top g_weight_dens_top ...
        g_altered_top g_high_Slope_region g_stdDevBot ...
        g_Determine_Lake g_bottom_new g_top_new_interp g_heights_orig ...
        g_topTest1 g_topTest2 ...
        g_botTest1 g_botTest2 g_botTest3 g_botTest4
    
    calculate_bot_stdDev();
    %indicators from the establishment of the top
    %false top, only surface actually detected. very strong indicator
    topOnlySurfDetected = g_altered_top == 1 | g_altered_top == 2;
    

    %this means that it was an outlier, but the next bottom was not
    %relevant (slightly less of an indicator of a NOT lake)
    topBadlyDetected = g_altered_top == 3 | g_altered_top == 4;
    
    
    %now the bottom
    bot_PhDensityLow = g_weight_dens_bottom < 3;    % low photon density bottom
    bot_higherThanTop = g_impossibleDepth == 1;     % means bottom calc overshot the top
    bot_tooShallow = (g_top_new_interp - g_bottom_new) < 1.5; % shallow lake

    stdDevTooHigh = g_stdDevBot > 5;    %bottom is meaningless.
    
    %each of these should indicate NOT LAKE
    topTest1 = g_top_new_interp;
    topTest1(~topOnlySurfDetected) = NaN;
    topTest2 = g_top_new_interp;
    topTest2(~topBadlyDetected) = NaN;
    botTest1 = g_bottom_new;
    botTest1(~bot_PhDensityLow)  = NaN;
    botTest2 = g_bottom_new;
    botTest2(~bot_higherThanTop)  = NaN;
    botTest3 = g_bottom_new;
    botTest3(~bot_tooShallow)  = NaN;
    botTest4 = g_bottom_new;
    botTest4(~stdDevTooHigh)  = NaN;
    
    g_topTest1 = topTest1;
    g_topTest2 = topTest2;
    g_botTest1 = botTest1;
    g_botTest2 = botTest2;
    g_botTest3 = botTest3;
    g_botTest4 = botTest4;
   
end

    function op = calculate_bot_stdDev()
    global g_photon_count g_stdDevBot g_BottomSmooth
    
    lHeights = size(g_BottomSmooth,1);          %%% size of array g_BottomSmooth
    g_stdDevBot = nan(size(g_BottomSmooth));    %%% output
    for iT = 1:lHeights                         %%% loop through photons and calculate standard deviation in window 3*g_photon_count left and right  
        wideBerth = g_photon_count*3;           
        startOrig = max(iT - wideBerth,1);
        endOrig = min(iT + wideBerth, lHeights);
        g_stdDevBot(iT) = nanstd(g_BottomSmooth(startOrig:endOrig));
    end
    
end
