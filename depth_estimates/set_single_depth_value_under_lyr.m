function op = set_single_depth_value_under_lyr( iM, photon_count, origTop, s_heights)
    % primary function where we determine the elevation of the water/ice surface and bottom 
    
    bandwidth = 0.1;
    [hts_01, df2_01, probs_01, subSurf_01, isSlope_01] = subfn_pull_subset_for_estimate_il_sub(iM,photon_count, bandwidth, s_heights);  
    %%% redo, but with wider bandwidth
    bandwidth = 0.3;
    [hts_03, df2_03, probs_03, subSurf_03, isSlope_03] = subfn_pull_subset_for_estimate_il_sub(iM,photon_count, bandwidth, s_heights);         
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
    
 
    %%% copy of information returned by subfn_pull_subset_for_estimate_il_sub
    htsOrig_01 = hts_01;
    htsOrig_03 = hts_03;
    df2Orig_01 = df2_01;
    df2Orig_03 = df2_03;
    %%% keep information about the the original top layer, ad information
    %%% from the first two peaks identified in subfn_pull_subset_for_estimate_il_sub
    hts_01(1) = origTop;
    hts_01(2) = htsOrig_01(1);
    hts_01(3) = htsOrig_01(2);
    df2_01(1) = NaN;
    df2_01(2) = df2Orig_01(1);
    df2_01(3) = df2Orig_01(2);

    hts_03(1) = origTop;
    hts_03(2) = htsOrig_03(1);
    hts_03(3) = htsOrig_03(2);
    df2_03(1) = NaN;
    df2_03(2) = df2Orig_03(1);
    df2_03(3) = df2Orig_03(2);
    
    %is narrow layer
    isNarrowLayer = (abs(hts_01(1) - hts_01(2))) < 1.2;                 % Is the second peak too close to the surface? 
    isBottomResolutionDependent = (abs(hts_01(2) - hts_03(2))) > 0.3;   % Does the second peak change based on resolution?
    isNextBottom = (abs(hts_01(3) - hts_03(2))) < 0.3;                  % Third peak from narrow bandwith close to the second peak of broad bandwith 
    
    %%% store all information
    op.isNarrowLayer = isNarrowLayer;
    op.isBottomResolutionDependent = isBottomResolutionDependent;
    op.isNextBottom = isNextBottom;
    op.interpolated = false;
    op.bottom_lg = hts_03(2);
    op.next_lg = hts_03(3);
    op.dif2Bottom_lg = df2_03(2);
    op.dif2Last_lg = df2_03(3);
    op.subsurf = subSurf_01.height;
    
    if( isNarrowLayer || isNextBottom )
        % if second peak is close to surface, or third peak from narrow bandwith close to the second peak of broad bandwith
        % -> then it's definitely an ice layer or they match depending on resolution
        op.bottom = hts_01(3);
        op.next = hts_01(2);
        op.bottom_orig = hts_01(3);
        op.next_orig = hts_01(2);
        op.bottom_prob = probs_01(3);
        op.next_prob = probs_01(2);
        op.dif2Bottom = df2_01(3);
        op.dif2Last = df2_01(2); 
    else
        if( isBottomResolutionDependent )
            op.interpolated = true;
            op.bottom = NaN;
            op.next = NaN;
            op.bottom_orig = hts_01(2);
            op.next_orig = hts_01(3);
            op.bottom_prob = probs_01(2);
            op.next_prob = probs_01(3);
            op.dif2Bottom = df2_01(2);
            op.dif2Last = df2_01(3);  
            %location of the second peak changes based on resolution
            %no matches and it's probably nonsense
        else
            op.bottom = hts_01(2);
            op.next = hts_01(3);
            op.top_orig = hts_01(1);
            op.bottom_orig = hts_01(2);
            op.next_orig = hts_01(3);
            op.bottom_prob = probs_01(2);
            op.next_prob = probs_01(3);
            op.dif2Bottom = df2_01(2);
            op.dif2Last = df2_01(3);          
        end
    end        

end