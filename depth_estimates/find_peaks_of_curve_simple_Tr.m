function [hts, df2s, probabilities, subSurf, foundN] = find_peaks_of_curve_simple_Tr(df1, df2, xHeight, ySix)
    % finds the first three peaks with the most negative second derivative
    % and then organizes them by height to determine a top and a bottom

    
    foundN = true;
    %find the points where the first derivative shifts from positive to
    %negative
    %and then find the min 2nd derivative at that point
    df1_sub = df1(1:size(df1,2)-1) >= 0;        % points where first derivative is positive
    df1_sub2 = df1(2:end) < 0;                  % points where first derivative is strictly negative
    dfCheck = df1_sub2 - df1_sub;               
    changePoint = dfCheck == 0;                 %the point where a shift occurs
   
%     %%% Note that above method does not always identifies the peaks correctly
%     %%% find peaks by looking at change in slope between consecutive points. If product of slope < 0 then a peak or through is found.  
%     changePoint=zeros(1,length(df1));
%     d1=diff(ySix);
%     d1_2=d1(1:(end-1)).*d1(2:end);
%     idx_peaks=find(d1_2<0&d1(1:(end-1))>0)+1;
%     changePoint(idx_peaks)=1;
    
    diff2ChPoint = df2(1:size(df1,2)-1);        %second derivative at the change points
    diff2ChPoint(~changePoint) = NaN;
    
    tSlope = sum(isnan(diff2ChPoint));          % if less than 10 non-peak points remain we just don't have enough data to make a determination
    if(tSlope < 10)
	hts=[]; df2s=[]; probabilities=[]; subSurf=[];
        foundN = false;
        return;
    end
        
    %find the 10 lower 2nd derivatives
    [vals, inds] = sort(diff2ChPoint);
    peak2ndDifs = vals(1:10); 
    peak2ndDifs_ind = inds(1:10);
    %remove NaN values
    peak2ndDifs_ind=peak2ndDifs_ind(isfinite(peak2ndDifs));
    peak2ndDifs=peak2ndDifs(isfinite(peak2ndDifs));
    
    %find the corresponding probabilities    
    probInd = nan(1,length(peak2ndDifs));
    probInd = ySix(peak2ndDifs_ind);
    %and sort them
    [valsP, indsProb] = sort(probInd);
 
    %find the indices in diff2ChPoint which correspond
    %to the three highest-probability cases 
    tIndex_A=nan(1,3);
    tIndex_A(1) = peak2ndDifs_ind(indsProb(end));
    if (length(indsProb)>1); tIndex_A(2) = peak2ndDifs_ind(indsProb(end-1)); end
    if (length(indsProb)>2); tIndex_A(3) = peak2ndDifs_ind(indsProb(end-2)); end
    %%% create arrays
    diff2_line=nan(1,3);
    heights=nan(1,3);
    probs=nan(1,3);
    %second derivatives
    diff2_line(1) = diff2ChPoint(tIndex_A(1));
    if (~isnan(tIndex_A(2))); diff2_line(2) = diff2ChPoint(tIndex_A(2));end
    if (~isnan(tIndex_A(3))); diff2_line(3) = diff2ChPoint(tIndex_A(3)); end
    %heights
    heights(1) = xHeight(tIndex_A(1));
    if (~isnan(tIndex_A(2))); heights(2) = xHeight(tIndex_A(2));end
    if (~isnan(tIndex_A(3))); heights(3) = xHeight(tIndex_A(3));end
    %probability
    probs(1) = ySix(tIndex_A(1));
    if (~isnan(tIndex_A(2))); probs(2) = ySix(tIndex_A(2));end
    if (~isnan(tIndex_A(3))); probs(3) = ySix(tIndex_A(3));end

    %now just organize by height
    [hts,df2s, probabilities, newInd] = select_Top_Bottom_spread_rev(heights, diff2_line, probs);
    
    %new Ind has the new order of heights
    indexTop = tIndex_A(newInd(3));
    indexBotMaybe = tIndex_A(newInd(2));
    
    %find the local minima in between to indicate when surface-esque
    %photons drop off, everything below this can be considered as
    %non-surface
    aSt = min(indexTop,indexBotMaybe);
    aEd = max(indexTop,indexBotMaybe);
    df2_inBetween = df2(1:size(df1,2)-1);
    df2_inBetween(1:aSt) = NaN;     %surface
    df2_inBetween(aEd:end) = NaN;     %surface
    [B, I] = sort(df2_inBetween);
    I(isnan(B)) = [];
    if(isempty(I))
        subSurf.diff2 = NaN;
        subSurf.height = NaN;
        subSurf.prob = NaN;
    else
        indexOf_SubSurface = I(end);
        subSurf.diff2 = df2(indexOf_SubSurface);
        subSurf.height = xHeight(indexOf_SubSurface);
        subSurf.prob = ySix(indexOf_SubSurface);
    end
    
    
end

function [hts,df2s, probabilities, newInd] = select_Top_Bottom_spread_rev(heightsVal, diff2_lineV, probs)
    %this just picks off the top and the bottom based on height
    
    [B,I] = sort(heightsVal,'MissingPlacement','first');
    hts = [heightsVal(I(3)) heightsVal(I(2)) heightsVal(I(1))];
    df2s = [diff2_lineV(I(3)) diff2_lineV(I(2)) diff2_lineV(I(1))];
    probabilities = [probs(I(3)) probs(I(2)) probs(I(1))];
    newInd = I;
    tCheck = nansum(~isnan(heightsVal));
    isLake = (tCheck > 1);                  %%% if more than one valid peak set isLake to TRUE 
    
end
