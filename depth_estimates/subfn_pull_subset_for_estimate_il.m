function [heights, diff2_line, probabilities, isLake, isFound] = subfn_pull_subset_for_estimate_il( midPoint_index, photonCount, bandwidth)
    global g_bandwidth
    % this function is the meat of the operation
    % it selects the number of photons within a width and then fits
    % it to a kernel density estimate with a bandwidth (think bin width)
    % picks off the first three peaks of the curve (high 2nd deriv)
    % and then organizes them by height to determine a top and a bottom
    
    [startOrig, endOrig, heights_list_sub] = subfn_select_subset_photons(midPoint_index, photonCount);  %%% select a number of photons (photonCount) left and right of the midpoint photon
    
    if (endOrig<1)          % no valid photons found -> return
         heights=[];
         diff2_line=[];
         probabilities = [];
         isLake=[];
         isFound=false;
         return;
    end
    
    %%% range of heights and number of steps
    minVal = min(heights_list_sub);
    maxVal = max(heights_list_sub);
    xHeight = minVal:0.1:maxVal;
    %%% create probability density function
    pdSix = fitdist(heights_list_sub,'Kernel','BandWidth',bandwidth);
    ySix_orig = pdf(pdSix,xHeight);
    nSum = nansum(ySix_orig);
    ySix = ySix_orig ./ nSum;                           %%% divide by number of valid points
    diff1 = gradient(ySix) ./ gradient(xHeight);        %%% first derivative
    diff2 = gradient(diff1) ./ gradient(xHeight);       %%% second derivative
    [heights, diff2_line, probabilities, isLake, isFound] = find_peaks_of_curve_simple_Tr(diff1, diff2, xHeight, ySix); 
    %This is returning heights (in order of elevation) and the corresponding
    %second derivative of the pdf (diff2_line) and probabilities, which will indicate
    %the strength of the peak.   
    %isLake contains the elevation (isLake.height), probability
    %(isLake.prob) and second derivative of the pdf (isLake.diff2)
    %if isFound is true, detection of at least one peak was sucessful
end


