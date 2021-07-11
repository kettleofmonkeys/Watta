    function op = post_pull_in_pond_file( fileName )

    %pulls in contents of csv file which generated surface tops and bottoms
    
    global g_dif2Top g_dif2Bottom g_dif2Last   ...
        g_top g_bottom g_next g_heights g_IsIce  ...
        g_fileName g_dirN g_TopSmooth g_BottomSmooth ...
        g_heights_orig g_lats g_lons g_time g_ph_conf ...
        g_top_prob g_bottom_prob g_next_prob ...
        g_isNarrowLayer g_isBottomResolutionDependent ...
        g_interpolated  g_isNextBottom ...
        g_top_lg g_bottom_lg g_next_lg g_dif2Top_lg g_dif2Bottom_lg g_dif2Last_lg ...
        g_subsurf g_notEnoughInformation g_surfaceStrength


    g_fileName = fileName;
    aFl = sprintf('%s/%s', g_dirN, g_fileName);
    g_csvList = load(aFl);
    g_heights_orig = g_csvList(:,1);
    g_heights = g_csvList(:,2);
    g_lats = g_csvList(:,3);
    g_lons = g_csvList(:,4);
    g_ph_conf = g_csvList(:,5);

    sz = size(g_heights,1);
    g_top = g_csvList(:,6);
    g_top(1) = g_top(2);
    g_top(sz) = g_top(sz-1);
    g_bottom = g_csvList(:,7);
    g_bottom(1) = g_bottom(2);
    g_bottom(sz) = g_bottom(sz-1);
    g_next = g_csvList(:,8); 
    g_next(1) = g_next(2);
    g_next(sz) = g_next(sz-1);
    g_TopSmooth = g_csvList(:,9);
    g_BottomSmooth= g_csvList(:,10); 
    
    g_dif2Top = g_csvList(:,11);
    g_dif2Bottom = g_csvList(:,12);
    g_dif2Last = g_csvList(:,13);
    g_top_prob = g_csvList(:,14);
    g_bottom_prob = g_csvList(:,15);
    g_next_prob = g_csvList(:,16);
    
    g_isNarrowLayer = g_csvList(:,17);
    g_isBottomResolutionDependent = g_csvList(:,18);
    g_isNextBottom = g_csvList(:,19);
    g_interpolated = g_csvList(:,20);      
    g_top_lg = g_csvList(:,21);
    g_bottom_lg = g_csvList(:,22);
    g_next_lg = g_csvList(:,23);
    g_dif2Top_lg = g_csvList(:,24);
    g_dif2Bottom_lg = g_csvList(:,25);
    g_dif2Last_lg = g_csvList(:,26);
    g_subsurf = g_csvList(:,27);
    g_notEnoughInformation = g_csvList(:,28);
    g_surfaceStrength = g_csvList(:,29);

end
