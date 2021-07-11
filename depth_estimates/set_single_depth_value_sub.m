function op = set_single_depth_value_sub( iM, photon_count)
    % primary function where we determine the elevation of the water/ice surface and bottom 
    
    % Look for peaks in the proability density function and return
    % elevation of peak (hts), second derivative of pdf at peak (df2) and probability (probs).
    % Also identify the lower elevation of the lake top (subSurf)
    % of 
    bandwidth = 0.1;
    [hts_01, df2_01, probs_01, subSurf_01, isSlope_01] = subfn_pull_subset_for_estimate_il(iM,photon_count, bandwidth);         

    bandwidth = 0.3;
    [hts_03, df2_03, probs_03, subSurf_03, isSlope_03] = subfn_pull_subset_for_estimate_il(iM,photon_count, bandwidth);         
    %Define a possible ice layer
    %If the bottom and the top are within 1.2m of each other
    %and if it doesn't show up in the wider bandwidth attempt
    %and if the height of the next probability is at the same height
    %then it's probably not an ice layer
    op.notEnoughInformation = 0;
    %illegitimate information 
    if (isSlope_01 == false)
        %means insufficient variability in the photon return (peak identification in kernel density not possible)
        op.notEnoughInformation = 1;
        return;
    end
        
    %is narrow layer
    isNarrowLayer = (abs(hts_01(1) - hts_01(2))) < 1.2;                 % Is the second layer too close to the surface? 
    isBottomResolutionDependent = (abs(hts_01(2) - hts_03(2))) > 0.3;   % Does the second layer change based on resolution?
    isNextBottom = (abs(hts_01(3) - hts_03(2))) < 0.3;                  % Third peak from narrow bandwith close to the second peak of broad bandwith -> next is actually the bottom
    %%% store all information
    op.isNarrowLayer = isNarrowLayer;
    op.isBottomResolutionDependent = isBottomResolutionDependent;
    op.isNextBottom = isNextBottom;
    op.interpolated = false;
    op.top_lg = hts_03(1);
    op.bottom_lg = hts_03(2);
    op.next_lg = hts_03(3);
    op.dif2Top_lg = df2_03(1);
    op.dif2Bottom_lg = df2_03(2);
    op.dif2Last_lg = df2_03(3);
    op.subsurf = subSurf_01.height;
   
    if( isNarrowLayer && isBottomResolutionDependent)
        if( isNextBottom )  
            % Third peak from narrow bandwith close to the second peak of broad bandwith -> next is actually the bottom. 
            % Just switch is so that next (lowest elevation out of the three) is the ice layer..
            op.top = hts_01(1);
            op.bottom = hts_01(3);
            op.next = hts_01(2);
            op.top_orig = hts_01(1);
            op.bottom_orig = hts_01(3);
            op.next_orig = hts_01(2);
            op.top_prob = probs_01(1);
            op.bottom_prob = probs_01(3);
            op.next_prob = probs_01(2);
            op.dif2Top = df2_01(1);
            op.dif2Bottom = df2_01(3);
            op.dif2Last = df2_01(2); 
        else
            %%% detected second peak too close to surface, needs to be
            %%% interpolated from neighbouring values
            op.interpolated = true;
            op.top = hts_01(1);
            op.bottom = NaN;
            op.next = NaN;
            op.top_orig = hts_01(1);
            op.bottom_orig = hts_01(2);
            op.next_orig = hts_01(3);
            op.top_prob = probs_01(1);
            op.bottom_prob = probs_01(2);
            op.next_prob = probs_01(3);
            op.dif2Top = df2_01(1);
            op.dif2Bottom = df2_01(2);
            op.dif2Last = df2_01(3);          
        end
    else        
        op.top = hts_01(1);
        op.bottom = hts_01(2);
        op.next = hts_01(3);
        op.top_orig = hts_01(1);
        op.bottom_orig = hts_01(2);
        op.next_orig = hts_01(3);
        op.top_prob = probs_01(1);
        op.bottom_prob = probs_01(2);
        op.next_prob = probs_01(3);
        op.dif2Top = df2_01(1);
        op.dif2Bottom = df2_01(2);
        op.dif2Last = df2_01(3);          
    end        

end