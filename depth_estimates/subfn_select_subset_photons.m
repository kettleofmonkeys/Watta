function [startOrig, endOrig, subList] = subfn_select_subset_photons(midPoint_index, along_track_width )
    global g_heights  g_recalc g_St g_Ed g_subsurf_new
    % this function is designed to pick off a number of photons to the left
    % and right of the midPoint_index of width along_track_width (measured
    % as a number of photons, not a distance!), but excluding NaN values 
    % (those which were likely excluded due to being too high in the
    % atmosphere).
    % This returns a subList of photons and the startOrig, endOrig
    % which are indices within the original array of photons (which
    % includes NaNs) 

    subheights = g_heights;

    if( g_recalc == true)
        midPoint_index = midPoint_index+g_St-1;
        subheights = g_heights;
        subheights(subheights >= g_subsurf_new) = NaN;        
    end
    
    
    %choose a larger window (where there will inherently be NaN values
    wideBerth = 3000;           
    maxStart = max(midPoint_index - wideBerth,1);
    maxEnd = min(midPoint_index + wideBerth, size(subheights,1));
    s_heights = subheights(maxStart:maxEnd); %create a subarray 

    %defining which are valid and which are not
    s_mid = midPoint_index - maxStart+1;        %reset the midpoint within the new sub-array
    validValues = ~isnan(s_heights);            %which values are valid (notNaN)
    nonValidValues = isnan(s_heights);          %which values are not valid (are NaN)
    nonValid_indices = find(isnan(s_heights));  %find indices which are NaN values
    readj = nansum(nonValid_indices < s_mid);   % how many indices previous to the midpoint are nan values
    
    s_mid_new = s_mid - readj;                  %the new midpoint accounting for the absence of nan values
    s_heights(nonValidValues) = [];             %remove nan elements
	
    if size(s_heights,1) == 0                   %this basically implies that there were no valid photons  (all NaNs)
        endOrig = 0;
        startOrig= 0;
        subList = [];
        return;
    end	
    
    % Capture the window of width 2*along_track_width (unless there's not enough room on either side)
    % this is within the larger window, and limited to actual valid height
    % values (notNaN)
    sStart = max(s_mid_new - along_track_width,1);
    sEnd = min(s_mid_new + along_track_width, size(s_heights,1));
    subList = s_heights(sStart:sEnd); 
    
    % now readjust to fit into the original g_heights list
    % find the start point accounting for the number of NaN values
    validValues_tStart = find((validValues == 1), sStart);
    nonValidValues_tStart = nonValidValues(1:max(validValues_tStart));
    nonValidValues_tStart = nansum(nonValidValues_tStart);
    startOrig = nonValidValues_tStart + sStart + maxStart-1;
   	
    validValues_tEnd = find((validValues == 1), sEnd);
    nonValidValues_End = nonValidValues(1:max(validValues_tEnd));
    nonValidValues_End = nansum(nonValidValues_End);
    endOrig = sEnd + nonValidValues_End + maxStart-1;
    
end
