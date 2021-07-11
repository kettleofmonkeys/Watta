function op = create_dist()
    %%% convert longitude/latitude to polar stereographic cooordinates
    %%% (referenced to first photon)
    global g_dist g_lats g_lons
    if (nanmean(g_lats)>0)
     [g_x,g_y]=polarstereo_fwd(g_lats,g_lons,6378137.0,0.08181919,70);
    else
     [g_x,g_y]=polarstereo_fwd(g_lats,g_lons);
    end
    
    g_dist=sqrt(g_x.^2+g_y.^2);
%    g_dist = g_dist - min(g_dist);
    g_dist=abs(g_dist-min(g_dist(1)));
end
