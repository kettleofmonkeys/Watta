function op = post_find_lake_surf_newest_bob()
    slope_lake_top_newest();
end


function op = slope_lake_top_newest()
    % define the lake surface
    global g_top g_top_new_interp g_classes g_top_filt g_lake_score g_heights g_dist g_Flat_Change  ...
       g_Slope_a g_Curve_a  g_interc_a g_index_change g_peak_ratio g_classes   g_index_change_slope ...

    %%% clean extreme outliers in g_heigths
    length_track=range(g_dist);
    length_seg=2e3;
    n_seg_track=ceil(length_track/length_seg);
    perc_lim_low=1;     %%% used to define upper part of lake surface elevation
    perc_lim_hi=99;     %%% futher outlier removal for lake surface
    
    g_heights_copy = g_heights;
    for i = 1:n_seg_track
        idx_ph_seg=find(g_dist>=(i-1)*length_seg&g_dist<i*length_seg);
        %%% find percentiles in segment
        lim_hi=prctile(g_heights_copy(idx_ph_seg),perc_lim_hi);
        lim_low=prctile(g_heights_copy(idx_ph_seg),perc_lim_low);
        idx_hi=find(g_heights_copy(idx_ph_seg)>lim_hi);
        idx_lo=find(g_heights_copy(idx_ph_seg)<lim_low);
        g_heights_copy(idx_ph_seg(idx_hi))=nan;     % <------------------
        g_heights_copy(idx_ph_seg(idx_lo))=nan;     % <------------------
    end

    %%% find breaks in surface layer
    [a,~]=ischange(g_top_new_interp,'mean');      %%% array a has the same size as g_top_new_interp, a value of 1 indicates an abrupt change in the mean of g_top_new_interp      
        g_index_change= [ 1 ; find(a==1) ;length(a)];  %%% g_index_change contains the locations of the breakpoints, add first and last photon
        g_top_filt=g_top_new_interp; %g_top
    %%% outlier removal based on cumulative density function 
        min_seg = 100;                    %%% minimum segment length
        cdf_low = 0.30;                   %%% lower limit for cdf
        cdf_high = 0.70;                  %%% upper limit for cdf; photons outside these limits are considered to be atmospher/noise
        min_photons = 200;                %%% minimun number of photons required to compute alternative top elevation  
        exp_perc=0.02;
        
        for i_a = 1:(size(g_index_change,1)-1)                              %%% Loop through breakpoints   
            idx_fit=g_index_change(i_a):(g_index_change(i_a+1)-1);              %%% index of points within bin 
            idx_fit_orig=idx_fit;                                           %%% save for later
            if ( max(g_dist(idx_fit))-min(g_dist(idx_fit))<min_seg)         %%% expand if necessary
                mid_seg=(min(g_dist(idx_fit))+max(g_dist(idx_fit)))/2;
                idx_fit=find(g_dist>=mid_seg-min_seg/2&g_dist<=mid_seg+min_seg/2);
            end
            if (sum(isfinite(g_heights_copy(idx_fit)))>0)
             [ cdf_seg,x_cdf_seg]=ecdf(g_heights_copy(idx_fit));                  %%% compute cdf
                idx_cdf_seg0=find(cdf_seg>=cdf_low);idx_cdf_seg0=idx_cdf_seg0(1);   %%% elevation corresponding to lower limit
                idx_cdf_seg1=find(cdf_seg<=cdf_high);idx_cdf_seg1=idx_cdf_seg1(end);%%% elevation corresponding to upper limit
                idx_cdf_seg1_exp=find(x_cdf_seg<=x_cdf_seg(idx_cdf_seg1)+range(x_cdf_seg)*exp_perc);    %%% expand top limit for very narrow peaks in cdf (multipe surfaces)    
                idx_cdf_seg1=idx_cdf_seg1_exp(end);        
                idx_bad=find(g_top_new_interp(idx_fit_orig)<x_cdf_seg(idx_cdf_seg0)|g_top_new_interp(idx_fit_orig)>x_cdf_seg(idx_cdf_seg1));  %%% find g_top photons outside these limits
                idx_bad=idx_fit_orig(idx_bad);
                if ( ~isempty(idx_bad) )                                          %%% if top photons found outside limits, replace with median of elevations within cdf limits of ORIGINAL bin
                    while (sum(isfinite(g_heights_copy(idx_fit_orig)))<min_photons)      %%% check if sufficient photons are availabe
                        add=round(min_photons-sum(isfinite(g_heights_copy(idx_fit_orig)))/2);    %%% if not, expand bin
                        idx_fit_orig=(max(idx_fit_orig(1)-add,1):min(idx_fit_orig(end)+add,length(g_top_new_interp)));
                    end
                    [ cdf_seg_orig,x_cdf_seg_orig]=ecdf(g_heights_copy(idx_fit_orig));      %%% cdf of (expanded) original bin photons
                    idx_cdf_seg0_orig=find(cdf_seg_orig>=cdf_low);idx_cdf_seg0_orig=idx_cdf_seg0_orig(1);   %%% elevation of first photon above lower threshold
                    idx_cdf_seg1_orig=find(cdf_seg_orig<=cdf_high);idx_cdf_seg1_orig=idx_cdf_seg1_orig(end); %%% elevation of last photon below upper threshold
                    idx_cdf_seg1_exp_orig=find(x_cdf_seg_orig<=x_cdf_seg_orig(idx_cdf_seg1_orig)+range(x_cdf_seg_orig)*exp_perc);    %%% expand top limit for very narrow peaks in cdf (multipe surfaces)    
                    idx_cdf_seg1_orig=idx_cdf_seg1_exp_orig(end); 
                    g_top_filt(idx_bad)=median(x_cdf_seg_orig(idx_cdf_seg0_orig:idx_cdf_seg1_orig)); %%% replace by median of photons within cdf limits  
                end
           end
        end

%%% new breakpoint analysis, based on slope 
    [a_slope,~]=ischange(g_top_filt,'linear','Threshold',100); 
    g_index_change_slope= [ 1 ; find(a_slope==1) ;length(a_slope)];  %%% g_index_change_slope contains the locations of the breakpoints (based on slope now), add first and last photon

%%% lake classification. This part computes a histogram for all bins, excluding the photons within +/- 0.5 m of the g_top_filt surface. Lakes are expected to have a prominent double peak. Double peaks may also occur in sloping areas, but are less prominent and gradually taper off. Detection is based on the ratio of first peak and second peak. 
%%% Additionally slope is also calculated to identify lakes with a weaker
%%% double peak, but flat surface
    min_depth = 0.5; % only consider photons 0.5 m below g_top
    max_depth = 30; % only consider photons no more than 30 m below g_top
    hist_step = 2 ; % bin width of histogram
    %%% construct arrays
    g_Flat_Change=nan(size(g_top_new_interp));
    g_peak_ratio=nan(size(g_top_filt));
    g_Slope_a=nan(size(g_top_filt));
    g_Curve_a=nan(size(g_top_filt));    %%% currently not used for classification
    g_interc_a=nan(size(g_top_filt));
    g_lake_score=nan(size(g_top_filt));
    %%% loop through segments between breakpoints and detect strength of
    %%% second peak
    for i_a = 1:(size(g_index_change_slope,1)-1) 
        idx_fit=g_index_change_slope(i_a):(g_index_change_slope(i_a+1)-1); 
        idx_fit_orig=idx_fit;                                           %%% save for later
        if ( max(g_dist(idx_fit))-min(g_dist(idx_fit))<min_seg)         %%% expand if necessary
            mid_seg=(min(g_dist(idx_fit))+max(g_dist(idx_fit)))/2;
            idx_fit=find(g_dist>=mid_seg-min_seg/2&g_dist<=mid_seg+min_seg/2);
        end
        while (sum(isfinite(g_heights_copy(idx_fit_orig)))<min_photons)      %%% check if sufficient photons are availabe
                add=round(min_photons-sum(isfinite(g_heights_copy(idx_fit_orig)))/2);    %%% if not, expand bin
                idx_fit_orig=(max(idx_fit_orig(1)-add,1):min(idx_fit_orig(end)+add,length(g_top_new_interp)));
        end
        diff_ph=(g_heights_copy(idx_fit)-g_top_filt(idx_fit));
        idx_ph_nosurf=find(abs(diff_ph)>min_depth&abs(diff_ph)<max_depth) ;
        hist_ph=hist(abs(diff_ph(idx_ph_nosurf)),(0:hist_step:max_depth));
        [~, a]=peakseek(hist_ph,3);a=sort(a,"descend"); % detect peaks in histogram; with minimum peak separation

        if (length(a)>1)
            g_peak_ratio(idx_fit_orig)=a(1)/a(2);       % compute ratio of strongest and second-strongest peak
        end
        fit=polyfit(g_dist(idx_fit)-mean(g_dist(idx_fit)),g_top_filt(idx_fit)-mean(g_top_filt(idx_fit)),2);    %%% get slope and curvature
        g_Slope_a(idx_fit_orig)=fit(2);
        g_Curve_a(idx_fit_orig)=fit(1);
        g_interc_a(idx_fit_orig)=fit(3);
    end
%%% classify
    peak_high=3.5;      %%% if ratio of strongest and second-strongest peaks exceeds this ratio, the bottom is well-defined
    peak_low=2.;        %%% if ratio of strongest and second-strongest peaks exceeds this ratio, the bottom is reasonably well-defined
    slope_thresh=3e-04; %%% if surface slope of segment is below this value, the surface is considered to be flat 
    curve_thresh=2e-05; %%% currently not used  
    %%% 1,2 and 4 most likely candidates for lakes, 3 could be streams on
    %%% inclined surface but also picks up lots of noise
    idx1=find(g_peak_ratio>=peak_high&abs(g_Slope_a)<=slope_thresh);g_lake_score(idx1)=1;  %%% double peak, flat 
    idx2=find(g_peak_ratio>=peak_high&abs(g_Slope_a)>slope_thresh&abs(g_Slope_a)<=10*slope_thresh);g_lake_score(idx2)=2;  %%% double peak, relatively flat
    idx3=find(g_peak_ratio>=peak_high&abs(g_Slope_a)>10*slope_thresh);g_lake_score(idx3)=3; %%% double peak, not flat
    idx4=find(g_peak_ratio>=peak_low&g_peak_ratio<peak_high&abs(g_Slope_a)<=slope_thresh);g_lake_score(idx4)=4; % weaker double peak, flat
    idx5=find(g_peak_ratio>=peak_low&g_peak_ratio<peak_high&abs(g_Slope_a)>slope_thresh&abs(g_Slope_a)<=10*slope_thresh);g_lake_score(idx5)=5; % weaker double peak, relatively flat
    idx6=find(g_peak_ratio>=peak_low&g_peak_ratio<peak_high&abs(g_Slope_a)>10*slope_thresh);g_lake_score(idx6)=6; % weaker double peak, not flat
    idx7=find(g_peak_ratio<peak_low&abs(g_Slope_a)<slope_thresh);g_lake_score(idx7)=7; % flat, but no obvious double peak 

    
    [idx_seg,~]=find(g_classes==g_lake_score);
    g_Flat_Change(idx_seg) = 1;
 

end