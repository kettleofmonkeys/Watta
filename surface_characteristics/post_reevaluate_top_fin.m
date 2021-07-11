function op = post_reevaluate_top_fin()
    global g_top_new_interp g_top

    g_top_new_interp = nan(size(g_top));
    
    develop_signal_spread();        %%% estimate local density of photons for a surface density and for background
    post_evaluate_pond_top_density();  %%% Compare the difference between the local photon density at any given pond top estimate and that of the background rate/density around
    post_replace_top_values();          %%% rearrange top,  bottom and potential ice layers.
    post_evaluate_pond_bottom_density_bob();  %%% produces weighted likelihood that a subsurface is a lake bottom 
end

function op = develop_signal_spread()
    global g_top_density_ph g_heights_orig g_background_rate g_surf_dens  g_dist ...
        g_binSurf g_binBack g_dist
    %%% for each point along the track, assess the density of photons for background vs near the surface

    szAll = size(g_heights_orig,1);                 %%% size of array
    g_background_rate = nan(size(g_heights_orig));  %%% create output arrays
    g_surf_dens = nan(size(g_heights_orig));
    g_binSurf = nan(size(g_heights_orig));
    g_binBack = nan(size(g_heights_orig));
    
    %%% cut off 5 m above and below lowest and highest photon
    maxM = max(g_heights_orig);
    minM = min(g_heights_orig);
    minM_base = floor(minM/10)*10 + 5;
    maxM_base = maxM - 5;
    %%% define edges of histograms in steps of 5 m; histograms is computed
    %%% twice, with a 2.5 vertical shift in the second round
    nedges = minM_base:5:maxM_base;
    nedges2 = minM_base+2.5:5:maxM_base;    
    
    tSep = 50;                                  %%% steps of 50 photons
    for iT = 1:tSep:szAll
        mSep = min(iT+tSep-1,szAll);
        wideBerth = 750;                        %%% collect 750 photons left and right        
        startOrig = max(iT - wideBerth,1);      %%% index of start photon
        endOrig = min(iT + wideBerth, szAll);   %%% index of end photon   <-------- why not iT+tsep+wideBerth ?

        dist1 = g_dist(startOrig);              %%% coordinate of start photon
        dist2 = g_dist(endOrig);                %%% coordinate of end photon
        areaV = abs(dist2 - dist1) .* 5;        %<-------- why divide by 2

        allHeights = g_heights_orig(startOrig: endOrig);   %%% photons in window 
        %basic photon density
        [N1,edges] = histcounts(allHeights,nedges);        %%% first histogram
        [N2,edges2] = histcounts(allHeights,nedges2);      %%% shifted histogram
        [surfBinDens,surfBinEd, backgroundBinDens, backgroundBinEdge] = find_max( N1,edges);        %%% find bin with maximum density (surface), return count in that bin and its edge. Estimate median background density and find the corresponding edge 
        [surfBinDens2,surfBinEd2, backgroundBinDens2, backgroundBinEdge2] = find_max( N2,edges2);   %%% repeat for shifted histogram
        %%% use histrogram with strongest peak; convert peak (surface) and
        %%% background count to density; assign to photons iT:mSep in
        %%% outpur arrays
        if( surfBinDens2 > surfBinDens)                     %%% shifted histogram has highest peak                
            g_surf_dens(iT:mSep) = surfBinDens2 / areaV;     
            g_background_rate(iT:mSep) = backgroundBinDens2 / areaV;            
            g_binSurf(iT:mSep) = surfBinEd2;
            g_binBack(iT:mSep) = backgroundBinEdge2;
        else                                                %%% original histogram has highest peak
            g_surf_dens(iT:mSep) = surfBinDens / areaV;
            g_background_rate(iT:mSep) = backgroundBinDens / areaV;            
            g_binSurf(iT:mSep) = surfBinEd;
            g_binBack(iT:mSep) = backgroundBinEdge;
        end
    end

    
end


function [surfBinDens,surfBinEd, backgroundBinDens, backgroundBinEdge] = find_max( N,edges)
    
    %the maximum case will be where the surface is
    [vals,inds] = sort(N);
    surfBinDens = vals(end);        %%% count in maximum bin 
    surfBinEd = edges(inds(end));   %%% edge of maximum bin
    tVals = vals;                   %%% temporary copy 
    vals(tVals == 0) = [];          %%% remove bins with no photons
    inds(tVals == 0) = [];
    vals(end) = [];                 %%% remove maximum bin
    inds(end) = [];
    if(isempty(inds))               %%% if no valid bins remain after maximum is removed, assign zero to background density
        backgroundBinDens = 0;
        backgroundBinEdge = edges(end-2);
    else        
        %use the median of the remaining bins for the background density and find the corresponding edge
        szT = size(vals,2);
        indV = floor(szT/2);
        backgroundBinDens = vals(indV);
        backgroundBinEdge = edges(inds(indV));
    end
end

