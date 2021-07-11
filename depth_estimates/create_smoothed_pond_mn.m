    function op = create_smoothed_pond_mn()
    global g_topOutliers g_bottomOutliers
    
    %This series of functions corrects for any overshots or 
    %undershots (where the surface/bottom are unnecessarily low/high) 
        
    % We first detect the outliers  g_topOutliers 
    % for the top and g_bottomOutliers for the bottom
    % and then correct the ones we find
    % When we address the bottom anomalies, we make sure not to reintroduce
    % a case where the top is anomalously high
    % 
    % The input for each of these is 
    %(1) a Threshold value: the number of
    % Median Absolute Deviations (MAD) within a window that an outlier must exceed
    % (2) a WindowScale - which multiplies the number of photons on each
    % side that should be considered when determining whether it's an
    % anomaly
    % Note that there's no particular reason either of these values have to
    % be whole numbers.
    % This is done iteratively, widening the window on a second step, to catch cases where
    % anomalies are in a run, such that a very wide window is required to
    % wash out the anomaly
    
    create_dist();  %so we can calculate distances
    
    disp('correct overshot top');
    g_topOutliers = calculateOutliers_movmed( 1, 4, 5 );   %%% look for outliers in the top layer (g_top)
    correct_overshots_on_top( 4, 5);                       %%% correct these outliers

    disp('correct undershots bottom');
    g_bottomOutliers = calculateOutliers_movmed( 0, 2, 2 );%%% look for outliers in the bottom layer (g_bottom)
    correct_undershots_on_bottom(2,2);                     %%% correct these outliers
    
    disp('correct BIG overshots top');                     %%% look for any large outliers that remain in the top layer (g_top)
    g_topOutliers = calculateOutliers_movmed( 1, 5, 6 );   %%% correct these large outliers
    correct_overshots_on_top( 5, 6);

    disp('correct BIG undershots bottom');                 %%% look for any large outliers that remain in the bottom layer (g_bottom)
    g_bottomOutliers = calculateOutliers_movmed( 0, 8, 2 );%%% correct these large outliers
    correct_undershots_on_bottom(8,2);
    
    disp('Interpolate missing values');
    interpolate_missing_values();                          % interpolate missing values in bottom layer (g_bottom)
    disp('Smooth top and bottom layers');
    final_smoothing();                                     % apply a smoothing to the top (g_top) and bottom layer (g_bottom). Output is saved in g_TopSmooth and g_BottomSmooth, respectively

end

function outlierList = calculateOutliers_movmed( topOrBottom, ThresholdFactor, WindowScale )
    global g_top g_bottom g_photon_count 
    %%% identify outliers that are more than (ThresholdFactor x median
    %%% absolute deviations [MAD]) away from the median of the window of
    %%% size g_photon_count*WindowScale*2
    if(topOrBottom == 1)                                   %%% top layer (g_top) 
        outlierList = isoutlier(g_top,'movmedian',g_photon_count*WindowScale*2,'ThresholdFactor',ThresholdFactor);
    else
        outlierList = isoutlier(g_bottom,'movmedian',g_photon_count*WindowScale*2,'ThresholdFactor',ThresholdFactor);
    end
    
end

function isOutlier = calculateOutliers_local(topOrBottom, ThresholdFactor, WindowScale, curValX, curValY)
    global g_top g_bottom  g_outlierTop g_outlierBottom g_photon_count ...
          g_heights
    %%%   Similar to calculateOutliers_movmed, but for one specific point in g_top or g_bottom
    %%%   Output is a logical value (true/false)
    szLst = size(g_heights,1);                              %%% number of total photons
    sStart = max(curValX - g_photon_count*WindowScale,1);   %%% first photon in range (g_photon_count*WindowScale) around center photon
    sEnd = min(curValX + g_photon_count*WindowScale, szLst);%%% last photon in range (g_photon_count*WindowScale) around center photon

    if(topOrBottom == 1)                                    %%% looking at top or bottom?                 
        valuesTest = g_top;
        valuesTest(g_outlierTop == 1) = NaN;                %%% don't use outliers
    else
        valuesTest = g_bottom;
        valuesTest(g_outlierBottom == 1) = NaN;             %%% don't use outliers
    end
    valuesTest(curValX) = curValY;                          %%% add back the value of the outlier considered                     
    TF = isoutlier(valuesTest(sStart:sEnd),'median','ThresholdFactor', ThresholdFactor);    %%% returns TRUE or FALSE for all points in array valuesTest

    testLoc = curValX - sStart+1;                           %%% index in valuesTest of the point considered 
    isOutlier = TF(testLoc);                                %%% is the point considered identified as an outlier (T/F)
    
end

function op = correct_overshots_on_top( ThresholdFactor, WindowScale )
    global g_top g_bottom g_next g_heights g_dif2Top g_dif2Bottom g_dif2Last...
           g_top_prob g_bottom_prob g_next_prob ...
        g_photon_count g_topOutliers 
    %this function corrects the outliers in the top that were identified in
    %calculateOutliers_movmed by re-estimating the top at that location
    %with a larger number of photons
   
    szLst = size(g_heights,1);                                  %%% number of photons
    for indexT = 2:szLst-1                                      %%% loop through all photons except first and last
        overShotTop = g_topOutliers(indexT) > 0;                %%% check if top photon is an outlier identified in calculateOutliers_movmed, if not -> continue
        if( ~overShotTop)
            continue;
        end
        
        correctedVal = false;                                   %%% if top photon is outlier -> correct by widening the range of photons (g_photon_count) used in the estimation of the surface/bottom 
        scaleFactor = 1;
        while (correctedVal == false && scaleFactor < 5)        %%% if no better solution found after 4 iterations,exit
            curVal_Top = g_top(indexT);
            overShotTop = calculateOutliers_local(1, ThresholdFactor, WindowScale, indexT, curVal_Top);     %%% is top photon still an outlier?
            if( overShotTop == 1 )                                                                          %%% if so, re-estimate top at this location with a larger number of photons
                scaleFactor = scaleFactor + 1;                                                              
                set_single_depth_value( indexT, g_photon_count*scaleFactor  );
            else                                                                                            %%% if not, stop loop
                correctedVal = true;
            end
        end
    end
end

function op = correct_undershots_on_bottom( ThresholdFactor, WindowScale )
    global g_top g_bottom g_next g_heights g_dif2Top g_dif2Bottom g_dif2Last...
        g_photon_count g_bottomOutliers g_top_prob g_bottom_prob g_next_prob ...
        g_interpolated
    %This function corrects the outliers in the bottom that were identified in
    %calculateOutliers_movmed by re-estimating the bottom at that location
    %with a larger number of photons
   
    szLst = size(g_heights,1);                                  %%% number of photons
    for indexT = 2:szLst-1                                      %%% loop through all photons except first and last
        underShotBottom = g_bottomOutliers(indexT) > 0;         %%% check if bottom photon is an outlier identified in calculateOutliers_movmed, if not -> continue
        if( ~underShotBottom )
            continue;
        end
        if( g_interpolated(indexT))                             %%% check if photon 
            continue;
        end

        correctedVal = false;                                    %%% if bottom photon is outlier -> correct by widening the range of photons (g_photon_count) used in the estimation of the surface/bottom 
        scaleFactor = 1;
        while (correctedVal == false && scaleFactor < 5)         %%% if no better solution found after 4 iterations,exit
            curVal_Bottom = g_bottom(indexT);
            underShotBottom = calculateOutliers_local(2, ThresholdFactor, WindowScale, indexT, curVal_Bottom);  %%% is bottom photon still an outlier?

            if( underShotBottom == 1 )                           %%% if so, re-estimate bottom at this location with a larger number of photons
                scaleFactor = scaleFactor + 1;
                bandwidth = 0.1;                                 %%% first check if widening the photon search window does not result in an outlier for the top
                [heights, diff2_line, probabilities, isLake] = subfn_pull_subset_for_estimate_il(indexT, g_photon_count*scaleFactor, bandwidth);
                overShot_Top = calculateOutliers_local( 1, ThresholdFactor, WindowScale, indexT, heights(1));
                if(~overShot_Top)                                %%% if new top estimate is not an outlier, re-estimate
                    set_single_depth_value( indexT, g_photon_count*scaleFactor);
                end
            else
                correctedVal = true;
            end
        end
    end
end

function op = interpolate_missing_values()
    global g_bottom g_heights
    %%% this function interpolates missing values in g_bottom
    allCount = size(g_heights,1);

    if allCount < 2                     %%% less than two points available, interpolation not possible
        return;
    end
    
    aX_all = 1:allCount;                %%% index of all photons (incl. nan)
    aX_or = aX_all;                     %%% create copy of aX_all
    aX_or(isnan(g_bottom)) = [];        %%% remove nan values in copy of aX_all
    or_Bottom = g_bottom;               %%% create copy of bottom estimates
    or_Bottom(isnan(or_Bottom)) = [];   %%% remove nan values in bottom estimates

    allCount = size(or_Bottom,1);

    if allCount < 2                     %%% less than two valid points available remaining, interpolation not possible
        return;
    end

    g_bottom = interp1(aX_or, or_Bottom, aX_all);   %%% do the interpolation
    g_bottom = g_bottom';

end

function op = final_smoothing()
    global g_top g_bottom  g_heights g_photon_count g_TopSmooth g_BottomSmooth ...
    % This function smooths the bottom and top using the rloess method with
    % a span of g_photon_count*5/(total number of photons)
    
    sZ = size(g_heights,1);
    xV = 1:sZ;

    smoothingValue = (g_photon_count*5)/sZ;

    g_TopSmooth = smooth(xV,g_top,smoothingValue,'rloess');
    g_BottomSmooth = smooth(xV,g_bottom,smoothingValue,'rloess');
end
