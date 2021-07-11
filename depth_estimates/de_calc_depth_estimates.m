function op = de_calc_depth_estimates( dirIN, dirPondOut, shortName, startIndex, endIndex, indexName )
    global g_photon_count g_fileName g_dirIn ...
        g_dirPondOut g_startIndex g_endIndex g_shortName g_recalc
    
    %%% define global variables
    g_recalc = false;
    g_dirIn = dirIN;
    g_dirPondOut = dirPondOut; 
    g_shortName = shortName;
    g_startIndex = startIndex;
    g_endIndex = endIndex;
    g_photon_count = 75;        %this*2 defines the number of VALID (non atmosphere) photons that will be considered for calculation

    fileNameOut = sprintf('%s/SubDepthEst_%s_%s.csv',g_dirPondOut, shortName, indexName);
    g_fileName = fileNameOut;

    extract_initial_estimates( startIndex, endIndex);       %calculates the first pond tops/bottoms as well as the derivatives
    create_smoothed_pond_mn();                              %this chops off erroneous values
    save_file(fileNameOut);                                 %save in .mat file
end

function op = calc_surface_strength()
    global g_surfaceStrength g_heights
    
    g_surfaceStrength = nan(size(g_heights));
    szV = size(g_heights,1);
    for iT = 300:szV-300
        g_surfaceStrength(iT) = nansum(~isnan(g_heights(iT-299:iT+299)));        
    end
    
end

function op = extract_initial_estimates(startIndex, endIndex)
    global g_heights g_heights_orig g_photon_count  ...
        g_top g_bottom g_next g_subsurf g_dif2Top g_dif2Bottom g_dif2Last  ...
        g_lats g_lons g_ph_conf  ...
        g_top_prob g_bottom_prob g_next_prob ...
        g_isNarrowLayer g_isBottomResolutionDependent g_isNextBottom ...
        g_top_orig g_bottom_orig g_next_orig g_interpolated ...
        g_top_lg g_bottom_lg g_next_lg g_dif2Top_lg g_dif2Bottom_lg g_dif2Last_lg ...
        g_dirIn g_shortName  ...
        g_notEnoughInformation 

    %%% load photon height data 
    heightFl = sprintf('%s/HeightFn/HeightFn_%s.mat', g_dirIn, g_shortName);
    load(heightFl);
    g_heights_orig = heights_val.heights(startIndex:endIndex);
    g_heights = g_heights_orig;
    clear('heights_val');
    %%% load photon location data
    locFl = sprintf('%s/LocFn/LocFn_%s.mat', g_dirIn, g_shortName);
    load(locFl);
    g_lats = loc_val.latitide(startIndex:endIndex);
    g_lons = loc_val.longitude(startIndex:endIndex);
    clear('loc_val');
    %%% load photon confidence data
    phConfFl = sprintf('%s/PhConfFn/PhConfFn_%s.mat', g_dirIn, g_shortName);
    load(phConfFl);
    g_ph_conf = ph_conf_val.ph_conf(startIndex:endIndex);
    clear('ph_conf_val');
    %%% load atl06 data
    atl06lfl = sprintf('%s/Heights_ATL06/ATL06hts_%s.mat', g_dirIn, g_shortName);
    load(atl06lfl);
    atl6heights = atl06heights.atl06_height_interp(startIndex:endIndex);
    atl6heights( atl6heights > 15000) = NaN;    %this is higher than any point on earth -> remove
    
    %%% remove photons outside of +/- 50 m of ATL06 
    g_heights(g_heights > (atl6heights+50)) = NaN;
    g_heights(g_heights < (atl6heights-50)) = NaN;
    calc_surface_strength();

    %%% create arrays for later use
    g_notEnoughInformation = zeros(size(g_heights_orig));
    g_top_orig = nan(size(g_heights_orig));
    g_bottom_orig = nan(size(g_heights_orig));
    g_next_orig = nan(size(g_heights_orig));

    g_top = nan(size(g_heights_orig));
    g_subsurf = nan(size(g_heights_orig));

    g_bottom = nan(size(g_heights_orig));
    g_next = nan(size(g_heights_orig));
    g_isNarrowLayer = nan(size(g_heights_orig));
    g_isBottomResolutionDependent = nan(size(g_heights_orig));
    g_isNextBottom = nan(size(g_heights_orig));
    g_interpolated = nan(size(g_heights_orig));
    
    g_top_prob = nan(size(g_heights_orig));
    g_bottom_prob = nan(size(g_heights_orig));
    g_next_prob = nan(size(g_heights_orig));

    g_dif2Top = nan(size(g_heights_orig));
    g_dif2Bottom = nan(size(g_heights_orig));
    g_dif2Last = nan(size(g_heights_orig));
    
    g_top_lg = nan(size(g_heights_orig));
    g_bottom_lg = nan(size(g_heights_orig));
    g_next_lg = nan(size(g_heights_orig));
    g_dif2Top_lg = nan(size(g_heights_orig));
    g_dif2Bottom_lg = nan(size(g_heights_orig));
    g_dif2Last_lg = nan(size(g_heights_orig));

    % go through each photon and calculate the 
    % the top, the bottom and the next peak (alternative bottom)
    % as well as the 2nd differential at that location
    allCount = size(g_heights,1);
    for iM = 1:allCount-1 
        set_single_depth_value( iM, g_photon_count );
    end    
    g_isNarrowLayer(isnan(g_isNarrowLayer)) = false;
    g_isBottomResolutionDependent(isnan(g_isBottomResolutionDependent)) = false;
    g_isNextBottom(isnan(g_isNextBottom)) = false;

end

function op = save_file(flNm )
    global g_dif2Top g_dif2Bottom g_dif2Last   ...
        g_top g_bottom g_next g_subsurf  g_IsIce ...
        g_csvList g_TopSmooth g_BottomSmooth ...
        g_top_prob g_bottom_prob g_next_prob ...
        g_isNarrowLayer g_isBottomResolutionDependent ...
        g_interpolated  g_isNextBottom ...
        g_top_lg g_bottom_lg g_next_lg g_dif2Top_lg g_dif2Bottom_lg g_dif2Last_lg ...
        g_startIndex g_endIndex ...
        g_notEnoughInformation g_surfaceStrength ...
        g_heights_orig g_heights g_lats g_lons g_ph_conf

    g_csvList = nan(g_endIndex-g_startIndex+1, 29);

    g_csvList(:,1) = g_heights_orig;
    g_csvList(:,2) = g_heights;
    g_csvList(:,3) = g_lats;
    g_csvList(:,4) = g_lons;
    g_csvList(:,5) = g_ph_conf;
    g_csvList(:,6) = g_top;
    g_csvList(:,7) = g_bottom;
    g_csvList(:,8) = g_next; 
    g_csvList(:,9) = g_TopSmooth;
    g_csvList(:,10) = g_BottomSmooth;     
    g_csvList(:,11) = g_dif2Top;
    g_csvList(:,12) = g_dif2Bottom;
    g_csvList(:,13) = g_dif2Last;
    g_csvList(:,14) = g_top_prob;
    g_csvList(:,15) = g_bottom_prob;
    g_csvList(:,16) = g_next_prob;
    g_csvList(:,17) = g_isNarrowLayer;
    g_csvList(:,18) = g_isBottomResolutionDependent;
    g_csvList(:,19) = g_isNextBottom;
    g_csvList(:,20) = g_interpolated;    
    g_csvList(:,21) = g_top_lg;
    g_csvList(:,22) = g_bottom_lg;
    g_csvList(:,23) = g_next_lg;
    g_csvList(:,24) = g_dif2Top_lg;
    g_csvList(:,25) = g_dif2Bottom_lg;
    g_csvList(:,26) = g_dif2Last_lg;
    g_csvList(:,27) = g_subsurf;
    g_csvList(:,28) = g_notEnoughInformation;
    g_csvList(:,29) = g_surfaceStrength;
    
    writematrix(g_csvList,flNm);
end
