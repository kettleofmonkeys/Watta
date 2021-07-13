function [startOrig, endOrig, subList] = subfn_select_subset_photons_redo(midPoint_index, along_track_width, s_heights )
    % this function is designed to pick off a number of photons to the left
    % and right of the midPoint_index of width along_track_width within the
    % the subset of photons stored in s_heigths
    % but excluding NaN values (those which were likely excluded due to
    % being too high in the atmosphere
    % this returns a subList of photons and the startOrig, endOrig
    % which are indices within the original array of photons (which
    % includes NaNs) 
    %
    % The code is similar to subfn_select_subset_photons, but is altered to 
    % account for the fact that maxStart is now always in s_heights
    s_mid = midPoint_index;
    validValues = ~isnan(s_heights);            %which values are valid (notNaN)
    nonValidValues = isnan(s_heights);          %which values are not valid (are NaN)
    nonValid_indices = find(isnan(s_heights));  %find indices which are NaN values
    readj = nansum(nonValid_indices < s_mid);   % how many indices previous to the midpoint are nan values
    
    s_mid_new = s_mid - readj;                  %the new midpoint accounting for hte absence of nan values
    s_heights(isnan(s_heights)) = [];           %remove nan elements
	
    if size(s_heights,1) == 0                   %this basically implies that there were no valid photons  (all NaNs)
        endOrig = 0;
        startOrig= 0;
        subList = [];
        return;
    end	
    
    % Capture the window of width 2*along_track_width (unless there's not enough room on eihter side)
    % this is within the larger window, and limited to actual height
    % values (notNaN)
    sStart = max(s_mid_new - along_track_width,1);
    sEnd = min(s_mid_new + along_track_width, size(s_heights,1));
    subList = s_heights(sStart:sEnd); 
    
    % now readjust to fit into the original g_heights list
    % find the start point accounting for the number of NaN values
    validValues_tStart = find((validValues == 1), sStart);
    nonValidValues_tStart = nonValidValues(1:max(validValues_tStart));
    nonValidValues_tStart = nansum(nonValidValues_tStart);
    startOrig = nonValidValues_tStart + sStart;
   	
    validValues_tEnd = find((validValues == 1), sEnd);
    nonValidValues_End = nonValidValues(1:max(validValues_tEnd));
    nonValidValues_End = nansum(nonValidValues_End);
    endOrig = sEnd + nonValidValues_End;  
    
end
