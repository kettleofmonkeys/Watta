function op = post_replace_top_values()
%%% Given weights as determined in post_evaluate_pond_top_density, rearrange top,  bottom and potential ice layers. Interpolate when necessary
    global  g_top_new g_top_new_interp g_top g_heights_orig ...
       g_unreliablePhDens_surf g_weight_dist_top g_weight_dens_top ...
       g_bottom g_next g_bottom_new g_next_new g_altered_top ...
       g_subsurf g_subsurf_new g_impossibleHeight g_acceptable_heights

    lHeights = size(g_top,1);
    g_subsurf_new = g_subsurf;
    tBestHeights = g_top;                               %%% copy g_top to tBestHeights
    tBestHeights(g_weight_dens_top <= 3) = NaN;         %%% remove photons with low photon density at top
    tBestHeights(g_weight_dist_top <= 3) = NaN;         %%% remove photons where top is far away from peak in photon density histogram
    tBestHeights(g_impossibleHeight == 1) = NaN;        %%% idem
    
    g_acceptable_heights = zeros(size(g_top));          %%% define which photons have a reliable top estimate
    g_acceptable_heights(~isnan(tBestHeights)) = 1;
    
    %%% copy arrays
    g_top_new = g_top;
    g_bottom_new = g_bottom;
    g_next_new = g_next;
    g_altered_top = zeros(size(g_top));
    
    for iT = 1:lHeights
        if(g_unreliablePhDens_surf(iT))     % low photon density around estimated top -> signal for unusability, nothing around here is usable
            continue;
        end
        
        if( isnan(tBestHeights(iT)))
            % due for replacement, using 75 photons left and right
            wideBerth = 75;           
            startOrig = max(iT - wideBerth,1);
            endOrig = min(iT + wideBerth, lHeights);

            g_top_new(iT) = NaN;
            botVal = g_bottom(iT);
            nextVal = g_next(iT);
            meanHtHighDens = nanmean(tBestHeights(startOrig:endOrig));
            stdHtHighDens = nanstd(tBestHeights(startOrig:endOrig));  
            
            if(botVal > nextVal)
                 %%% if bottom value is within 2 standard devations + 1 m
                 %%% of mean top surface height in bin, we assume that the
                 %%% bottom is actually the top surface
                if( botVal < meanHtHighDens + 2*stdHtHighDens+1 && ...
                        botVal > meanHtHighDens - 2*stdHtHighDens-1 )  % within range
                    g_top_new(iT) = botVal;
                    g_bottom_new(iT) = nextVal;
                    g_next_new(iT) = NaN;
                    g_altered_top(iT) = 1;
                end                
            else
                 %%% if next value is within 2 standard devations + 1 m
                 %%% of mean top surface height in bin, we assume that the
                 %%% next value is actually the top surface
                if( nextVal < meanHtHighDens + 2*stdHtHighDens && ...
                        nextVal > meanHtHighDens - 2*stdHtHighDens )  % within range
                    g_top_new(iT) = nextVal;  %this was orginally misconstrued as an ice layer
                    g_bottom_new(iT) = NaN;
                    g_next_new(iT) = NaN;
                    g_altered_top(iT) = 2;
                end
            end
        end
    end
    
    %%% Remove outliers using moving mean
    WindowScale = 75*3;ThresholdFactor = 2;
    outlierList = isoutlier(g_top_new,'movmean',WindowScale,'ThresholdFactor',ThresholdFactor);
    g_top_new(outlierList == 1) = NaN;  %kill outliers
    g_acceptable_heights(outlierList == 1) = 0;
    g_subsurf_new(g_acceptable_heights == 0) = NaN;
    
    allCount = size(g_heights_orig,1);
    if sum(~isnan(g_heights_orig)) < 1
        return;         %%% no photons remaining
    end
    
    
    %%% interpolate missing values
    aX_all = 1:allCount;
    aX_or = aX_all;
    aX_or(isnan(g_top_new)) = [];
    or_Top = g_top_new;
    or_Top(isnan(or_Top)) = [];    
    InterpolatedTops = isnan(g_top_new);
    g_altered_top(InterpolatedTops == 1) = 3;
    g_altered_top(outlierList) = 4;
    
    allCount = size(or_Top,1);
    if sum(~isnan(or_Top)) < 2
        return;     %%% not enough photons remaining for interpolation
    end
    g_top_new_interp= interp1(aX_or, or_Top, aX_all);
    g_top_new_interp = g_top_new_interp';

    %%% repeat interpolation for subsurface layer
    allCount = size(g_heights_orig,1);
    aX_all = 1:allCount;
    aX_or = aX_all;
    aX_or(isnan(g_subsurf_new)) = [];
    or_Ssurf = g_subsurf_new;
    or_Ssurf(isnan(or_Ssurf)) = [];        
    allCount = size(or_Ssurf,1);
    if sum(~isnan(or_Ssurf)) < 2
        return;
    end
    g_subsurf_new= interp1(aX_or, or_Ssurf, aX_all);
    g_subsurf_new = g_subsurf_new';
    
 
end