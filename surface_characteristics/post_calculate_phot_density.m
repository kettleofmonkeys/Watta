function densityTot = post_calculate_phot_density(lowerBound_m, upperBound_m, startOrig, endOrig)
    global  g_heights_orig g_dist g_re_dist g_recalc
    
    subSet = g_heights_orig(startOrig:endOrig); %%% collect photons in subset
   
    subSet(subSet < lowerBound_m) = NaN;    %%% remove photons below lower bound
    subSet(subSet > upperBound_m) = NaN;    %%% remove photons above lower bound
    lSz = size(subSet,1);                   %%% number of photons in subset
    if( g_recalc == true)
        tGr = abs(g_re_dist(endOrig) - g_re_dist(startOrig)) / lSz;        %%% distance covered in subset divided by number of photons [units: m/photon]
    else
        tGr = abs(g_dist(endOrig) - g_dist(startOrig)) / lSz;              %%% distance covered in subset divided by number of photons [units: m/photon]
    end
    tGr = repmat(tGr,(size(subSet,1)),1);             
    tHeight = max(upperBound_m -lowerBound_m,0);    %%% verticale range of min/max bound [units: m]
    AreaV = nansum(abs(tGr) .* tHeight);            %%% distance * height 
    densityTot = nansum(~isnan(subSet)) / AreaV;    %%% number of valid photons in subset / area
end
    