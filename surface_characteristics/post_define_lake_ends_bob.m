function op = post_define_lake_ends_bob(lake_edges )
    global  g_lakes  g_Flat_Change

    %Assigns a number to the lakes, stored in g_lakes
    g_lakes = nan(size(g_Flat_Change));
    whichVal = 1;
    for i = 1:length(lake_edges.idx_breaks_left)
        leftT = lake_edges.idx_breaks_left(i);
        rightT = lake_edges.idx_breaks_right(i);
        g_lakes(leftT:rightT) = whichVal;
        whichVal = whichVal+1;

    end
     g_lakes(isnan(g_lakes)) = 0;

end