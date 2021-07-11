function op = create_smoothed_pond_mn_re()
    global g_re_topOutliers g_re_bottomOutliers g_re_top g_top_new_interp g_St g_Ed
    
    %%% this is basically a copy of create_smoothed_pond_mn, but redefined
    %%% for the bottom recalculation
    
    g_re_top = g_top_new_interp(g_St:g_Ed);
    
    %This series of functions corrects for any overshots or 
    %undershots (where the surface/bottom are unnecessarily low/high) in
    %the selection of photons in g_bottom between index g_St and g_Ed
        
    % We first detect the outliers g_bottomOutliers for the bottom
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

    disp('correct undershots bottom');
    g_re_bottomOutliers = calculateOutliers_movmed( 0, 2, 2 );
    correct_undershots_on_bottom(2,2);
    
    disp('correct sBIG undershots bottom');
    g_re_bottomOutliers = calculateOutliers_movmed( 0, 8, 2 );
    correct_undershots_on_bottom(8,2);
    disp('Interpolate missing values');
    interpolate_missing_values();
end




function outlierList = calculateOutliers_movmed( topOrBottom, ThresholdFactor, WindowScale )
    global g_re_top g_re_bottom g_photon_count 

    if(topOrBottom == 1)
        outlierList = isoutlier(g_re_top,'movmedian',g_photon_count*WindowScale*2,'ThresholdFactor',ThresholdFactor);
        %for a more rigorous (read: annoyingly dense) set of outliers, we
        %can artificially create a window and then use hte "gesd" function
        %instead of movmedian. This is sensitive to outliers overlapping
    else
        outlierList = isoutlier(g_re_bottom,'movmedian',g_photon_count*WindowScale*2,'ThresholdFactor',ThresholdFactor);
    end
    
end

function isOutlier = calculateOutliers_local(topOrBottom, ThresholdFactor, WindowScale, curValX, curValY)
    global g_re_top g_re_bottom  g_re_outlierTop g_re_outlierBottom g_photon_count ...
          g_subsetOfHeights

    szLst = size(g_subsetOfHeights,1);
    sStart = max(curValX - g_photon_count*WindowScale,1);
    sEnd = min(curValX + g_photon_count*WindowScale, szLst);

    if(topOrBottom == 1)
        valuesTest = g_re_top;
    else
        valuesTest = g_re_bottom;
        valuesTest(g_re_outlierBottom == 1) = NaN;
    end
    valuesTest(curValX) = curValY;
    TF = isoutlier(valuesTest(sStart:sEnd),'median','ThresholdFactor', ThresholdFactor);

    testLoc = curValX - sStart+1;
    
    isOutlier = TF(testLoc);
    
end

function op = correct_undershots_on_bottom( ThresholdFactor, WindowScale )
    global  g_re_bottom  g_subsetOfHeights ...
        g_photon_count g_re_bottomOutliers ...
        g_re_interpolated g_re_top ...
        g_re_isNarrowLayer g_re_isBottomResolutionDependent ...
        g_re_isNextBottom g_re_bottom_lg ...
        g_re_next_lg g_re_dif2Bottom_lg ...
        g_re_dif2Last_lg g_re_next ...
        g_re_bottom_orig g_re_next_orig ...
        g_re_bottom_prob g_re_next_prob ...
        g_re_dif2Bottom g_re_dif2Last g_re_notEnoughInformation ...
        g_St g_heights g_top_filt
    
    %this function finds any instance where the pond top exceeds the wider mean
    %with outliers removed
    top_filt = g_top_filt;

    subheights = g_heights;
    szLst = size(g_subsetOfHeights,1);
    for indexT = 2:szLst-1 %g_photon_count*2+1:szLst-g_photon_count*2
        underShotBottom = g_re_bottomOutliers(indexT) > 0;
        if( ~underShotBottom )
            continue;
        end
        if( g_re_interpolated(indexT))
            continue;
        end

        correctedVal = false;
        scaleFactor = 1;
        while (correctedVal == false && scaleFactor < 5)
            curVal_Bottom = g_re_bottom(indexT);
            underShotBottom = calculateOutliers_local(2, ThresholdFactor, WindowScale, indexT, curVal_Bottom);

            if( underShotBottom == 1 ) 
                
                scaleFactor = scaleFactor + 1;
                bandwidth = 0.1;
                
                midPoint_index = indexT;
                midPoint_index = midPoint_index+g_St-1;
                %choose a larger window (where there will inherently be NaN values
                wideBerth = 3000;           
                maxStart = max(midPoint_index - wideBerth,1);
                maxEnd = min(midPoint_index + wideBerth, size(subheights,1));
                s_heights = subheights(maxStart:maxEnd); %create a subarray 
                inVl = midPoint_index - maxStart + 1;       %this is relative to the beginning of the subsection

                topValue = g_top_filt(indexT+g_St-1);
                [heights, diff2_line, probabilities, isLake] = subfn_pull_subset_for_estimate_il_sub(inVl, g_photon_count*scaleFactor, bandwidth, s_heights);

                overShot_Top = calculateOutliers_local( 1, ThresholdFactor, WindowScale, indexT, heights(1));

                if(~overShot_Top)
                    fVal = set_single_depth_value_under_lyr( inVl, g_photon_count*scaleFactor, topValue, s_heights);
                    g_re_notEnoughInformation(indexT) = fVal.notEnoughInformation;
                    if( g_re_notEnoughInformation(indexT) == 0)
                        g_re_isNarrowLayer(indexT) = fVal.isNarrowLayer;    
                        g_re_isBottomResolutionDependent(indexT) = fVal.isBottomResolutionDependent;
                        g_re_isNextBottom(indexT) = fVal.isNextBottom;
                        g_re_interpolated(indexT) = fVal.interpolated;
                        g_re_bottom_lg(indexT) = fVal.bottom_lg;
                        g_re_next_lg(indexT) = fVal.next_lg;
                        g_re_dif2Bottom_lg(indexT) = fVal.dif2Bottom_lg;
                        g_re_dif2Last_lg(indexT) = fVal.dif2Last_lg;
                        g_re_bottom(indexT) = fVal.bottom;
                        g_re_next(indexT) = fVal.next;
                        g_re_bottom_orig(indexT) = fVal.bottom_orig;
                        g_re_next_orig(indexT) = fVal.next_orig;
                        g_re_bottom_prob(indexT) = fVal.bottom_prob;
                        g_re_next_prob(indexT) = fVal.next_prob;
                        g_re_dif2Bottom(indexT) = fVal.dif2Bottom;
                        g_re_dif2Last(indexT) = fVal.dif2Last; 
                    end

                end
            else
                correctedVal = true;
            end
        end
    end
end

function op = interpolate_missing_values()
    global g_re_bottom 
    %%% this function interpolates missing values in g_re_bottom (a subset of g_bottom)

    allCount = size(g_re_bottom,1);

    if allCount < 2
        return;
    end
    
    aX_all = 1:allCount;                %%% index of all photons (incl. nan)
    aX_or = aX_all;                     %%% create copy of aX_all
    aX_or(isnan(g_re_bottom)) = [];     %%% remove nan values in copy of aX_all
    or_Bottom = g_re_bottom;            %%% create opy of bottom estimates
    or_Bottom(isnan(or_Bottom)) = [];   %%% remove nan values in bottom estimates

    allCount = size(or_Bottom,1);

    if allCount < 2                     %%% less than two valid points available remaining, interpolation not possible
        return;
    end

    g_re_bottom = interp1(aX_or, or_Bottom, aX_all);
    g_re_bottom = g_re_bottom';
end

function op = final_smoothing()
    global g_re_top g_re_bottom  g_subsetOfHeights g_photon_count g_re_TopSmooth g_re_BottomSmooth ...
    % This function smooths the bottom using the rloess method with
    % a span of g_photon_count*5/(total number of photons in subset)        
    
    sZ = size(g_subsetOfHeights,1);
    xV = 1:sZ;
    
    smoothingValue = (g_photon_count*5)/sZ;
    
    g_re_BottomSmooth = smooth(xV,g_re_bottom,smoothingValue,'rloess');
end
