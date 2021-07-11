function op = surface_characteristics_mn_tri(dirIn, dirOut, fileNm, shName )
    global gName g_photon_count gDirOut g_heights g_recalc g_subsurf g_subsurf_new g_classes

    g_recalc = false;                               %%% switch for second round of calculating the bottom
    warning('off','all');
    g_photon_count_normal=75;                       %%% photon window width for detection of top/bottom
    g_photon_count_recalc=30;                       %%% photon window width for second round of calculating the bottom
    g_classes = [1,2,4,7];                              %%% only consider these classes of lakes (see post_find_lake_surf_newest for definition of classes)

    
    g_photon_count = g_photon_count_normal;

    gDirOut = dirOut;                               %%% output directory
    gName = sprintf('SType_%s', shName);            %%% input file
    outFl1 = sprintf('%s/%s_gph%d_%d_surf_wt.csv', dirOut, gName,g_photon_count,g_photon_count_recalc); %%% output file 1
    outFl2 = sprintf('%s/%s_gph%d_%d_surf_lt.csv', dirOut, gName,g_photon_count,g_photon_count_recalc); %%% output file 2

    disp('pulling file');                           
    post_pull_in_pond_file( fileNm );%, sName );    %%% Read in file     
    disp('create_dist');
    create_dist();                                  %%% Convert longitude/latitude to polar stereographic cooordinate
    disp('reeva luate top');
    post_reevaluate_top_fin();                      %%% reduce outliers substantially 
    disp('lake Surf');    
    post_find_lake_surf_newest_bob();  
    disp('lake edges');
    lake_edges = post_find_edges_of_lake_newest_bob();

    g_photon_count = g_photon_count_recalc
    g_recalc = true;
    
    n_nodes = 12  
%   parpool(n_nodes);       %uncomment for parallel computing and change 'for' to 'parfor' in post_reset_kerndensity_lakes.m
    for i  = 1 : length(lake_edges.idx_breaks_left)    
       msgDsp = sprintf('Running %d / %d', i, length(lake_edges.idx_breaks_left));
       disp(msgDsp);
       sub_recalculation(lake_edges.idx_breaks_left(i),lake_edges.idx_breaks_right(i));
       refit_sub_section(lake_edges.idx_breaks_left(i),lake_edges.idx_breaks_right(i))
    end 
%    delete(gcp('nocreate')); %uncomment for parallel computing

     g_recalc = false;
     g_photon_count = g_photon_count_normal;
     
     disp('reevaluate top');
     post_reevaluate_top_fin();
     disp('saturation');
     pull_saturation(dirIn, shName);
     disp('lake edges');
     lake_edges = post_find_edges_of_lake_newest_bob();
     disp('lake ends');
     post_define_lake_ends_bob(lake_edges );
     disp('ice layer');
     post_detect_ice_layer_bob(lake_edges );
     disp('reset edges');
     post_reset_edges_for_ice_bob(lake_edges );
     disp('final smoothing');
     for i  = 1 : length(lake_edges.idx_breaks_left)  
        post_final_smoothing(lake_edges.idx_breaks_left(i),lake_edges.idx_breaks_right(i))
     end      
    disp('final estimates');
    post_produce_final_estimates_bob();
    disp('calculate lake probabilities');
    post_calculate_lake_prob();
    disp('write final files');
    write_final_files_weighted(outFl1);
    write_final_files_lakes(outFl2);   
end

function op = pull_saturation(dirIn, shortName)
    global g_full_sat_fract_interp g_near_sat_fract_interp
    %%% load saturation data contained in ATL03
    heightFl = sprintf('%s/SaturFn/SaturFn_%s.mat', dirIn, shortName);
    load(heightFl);
    g_full_sat_fract_interp = sat_val.full_sat_fract_interp;
    g_near_sat_fract_interp = sat_val.near_sat_fract_interp;
end

function op = refit_sub_section(st,ed)
%%% put the data produced in the recalculation in the right place
    global g_heights g_St g_Ed ... 
        g_isNarrowLayer g_re_isNarrowLayer ...
        g_isBottomResolutionDependent g_re_isBottomResolutionDependent ...
        g_isNextBottom g_re_isNextBottom ...
        g_interpolated g_re_interpolated ...
        g_top g_re_top ...
        g_bottom g_re_bottom ...
        g_next g_re_next ...
        g_subsurf g_re_subsurf...
        g_top_orig g_re_top_orig ...
        g_bottom_orig g_re_bottom_orig ...
        g_next_orig g_re_next_orig ...
        g_top_prob g_re_top_prob ...
        g_bottom_prob g_re_bottom_prob ...
        g_next_prob g_re_next_prob ...
        g_dif2Top g_re_dif2Top ...
        g_dif2Bottom g_re_dif2Bottom ...
        g_dif2Last g_re_dif2Last ...
        g_top_lg g_re_top_lg ...
        g_bottom_lg g_re_bottom_lg ... 
        g_next_lg g_re_next_lg  ....
        g_dif2Top_lg g_re_dif2Top_lg ...
        g_dif2Bottom_lg g_re_dif2Bottom_lg ...
        g_dif2Last_lg g_re_dif2Last_lg ...
        g_notEnoughInformation g_re_notEnoughInformation ...
        g_surfaceStrength g_re_surfaceStrength ...        
        g_TopSmooth g_re_TopSmooth ...
        g_BottomSmooth g_re_BottomSmooth ...
        g_top_density_ph g_re_top_density_ph...
        g_unreliablePhDens_surf g_re_unreliablePhDens_surf ...
        g_impossibleHeight g_re_impossibleHeight ...
        g_weight_dist_top g_re_weight_dist_top ...
        g_threshold_spreadDens g_re_threshold_spreadDens ...
        g_top_new  g_re_top_new ...
        g_top_new_interp g_re_top_new_interp  ...
        g_bottom_new g_re_bottom_new ...
        g_next_new g_re_next_new ...
        g_altered_top g_re_altered_top ...
        g_subsurf_new g_re_subsurf_new ...
        g_acceptable_heights g_re_acceptable_heights

        g_isNarrowLayer(st:ed) =  g_re_isNarrowLayer;
        g_isBottomResolutionDependent(st:ed) =   g_re_isBottomResolutionDependent;
        g_isNextBottom(st:ed) =   g_re_isNextBottom;
        g_interpolated(st:ed) =   g_re_interpolated;
        g_bottom(st:ed) =   g_re_bottom;
        g_next(st:ed) =   g_re_next;
        g_bottom_orig(st:ed) =   g_re_bottom_orig;
        g_next_orig(st:ed) =   g_re_next_orig;
        g_bottom_prob(st:ed) =   g_re_bottom_prob;
        g_next_prob(st:ed) =   g_re_next_prob;
        g_dif2Bottom(st:ed) =   g_re_dif2Bottom;
        g_dif2Last(st:ed) =   g_re_dif2Last;
        g_bottom_lg(st:ed) =   g_re_bottom_lg; 
        g_next_lg(st:ed) =   g_re_next_lg;
        g_dif2Bottom_lg(st:ed) =   g_re_dif2Bottom_lg;
        g_dif2Last_lg(st:ed) =   g_re_dif2Last_lg;
        g_notEnoughInformation(st:ed) =   g_re_notEnoughInformation;

end





function op = write_final_files_lakes( saveNm )
% write output le including lake surface, bottom, depth, surface ice, subsurface ice and saturation fraction (from the original le)
    global g_lakesurface_final g_subsurf_ice g_surface_ice ...
        g_depth_all g_depth_real g_depth_corr_all g_depth_corr_real ...
        g_dist g_bottom_final g_bottom_final_real  ...
        g_iceLayer g_lats g_lons g_ph_conf g_heights_orig ...
        g_bottom_refit_edges  g_lakes ...
        g_next_new g_bottom_new  g_top_new_interp ...
        g_subsurf_new  ...
        g_full_sat_fract_interp g_near_sat_fract_interp ...
        g_Slope_a g_Curve_a g_interc_a g_BottomSmooth g_BottomSmooth_low g_top_filt
    
    sz = size(g_dist,1);
    sV = nan(sz,18);
    sV(:,1) = g_lats;
    sV(:,2) = g_lons;
    sV(:,3) = g_heights_orig;
    sV(:,4) = g_ph_conf;
    sV(:,5) = g_dist;
    
    %top/bot metrics
    sV(:,6) = g_top_new_interp;
    sV(:,7) = g_next_new;
    sV(:,8) = g_bottom_new;
    sV(:,9) = g_subsurf_new;        %just under the surface (the trough in the distribution)
    sV(:,10) = g_iceLayer;          %used to calculate subsurface ice
    sV(:,11) = g_subsurf_ice;       %not to be confused with g_subsurf_new
    sV(:,12) = g_surface_ice;    
    sV(:,13) = g_full_sat_fract_interp;
    sV(:,14) = g_near_sat_fract_interp;

    
    sV(:,15) = g_lakesurface_final;
    sV(:,16) = g_lakes;    
    sV(:,17) = g_bottom_refit_edges;        
    sV(:,18) = g_bottom_final;
    sV(:,19) = g_bottom_final_real;         %ignoring ice layer areas
    sV(:,20) = g_BottomSmooth;
    sV(:,21) = g_depth_all;
    sV(:,22) = g_depth_corr_all;
    sV(:,23) = g_depth_real;                %this excludes interpolated elements
    sV(:,24) = g_depth_corr_real;
    
    sV(:,25) = g_Slope_a;
    sV(:,26) = g_Curve_a;
    sV(:,27) = g_interc_a;
    sV(:,28) = g_top_filt;
    if (length(g_BottomSmooth_low)==0)
        g_BottomSmooth_low=g_BottomSmooth;
    end
    sV(:,29) = g_BottomSmooth_low; 
    writematrix(sV,saveNm);
    
end

function op = write_final_files_weighted( saveNm )
% write output including top/bottom values as well as metrics used to calculate the lake score
    global g_lakesurface_final g_subsurf_ice g_surface_ice ...
        g_depth_all g_depth_real g_depth_corr_all g_depth_corr_real ...
        g_dist g_x_dist g_BottomSmooth_bin_fix g_Flat_Change ...
        g_bottom_final g_bottom_final_real  ...
        g_iceLayer g_interpolated ...
        g_next g_bottom g_lats g_lons g_ph_conf g_heights_orig g_csvList gName ...
        g_top g_tot_bottom_refit_edges g_lakebottom_final_splined g_lakes ...
        g_bot_density_interpolated_ph g_top_density_ph g_bot_density_ph ...
        g_binSurf g_surf_dens g_background_rate g_impossibleHeight ...
        g_weight_dist_top g_weight_dens_top ...
        g_top_new g_next_new g_bottom_new g_altered_top g_top_new_interp ...
        g_subsurf_new g_impossibleDepth g_weight_dens_bottom ...
        g_bottom_density_ph g_unreliablePhDens_bottom ...
        g_topTest1 g_topTest2 g_botTest1 g_botTest2 g_botTest3 g_botTest4 ...
        g_stdDevBot g_lake_score
    
    sz = size(g_dist,1);
    sV = nan(sz,34);
    sV(:,1) = g_lats;
    sV(:,2) = g_lons;
    sV(:,3) = g_heights_orig;
    sV(:,4) = g_ph_conf;
    sV(:,5) = g_top;
    sV(:,6) = g_bottom;
    sV(:,7) = g_next;  
    
    %calculation of hist heights
    sV(:,8) = g_binSurf;
    sV(:,9) = g_surf_dens;
    sV(:,10) = g_background_rate;
    
    %top metrics
    sV(:,11) = g_impossibleHeight;
    sV(:,12) = g_weight_dist_top;
    sV(:,13) = g_weight_dens_top;
    sV(:,14) = g_top_density_ph;
    sV(:,15) = g_top_new;
    sV(:,16) = g_next_new;
    sV(:,17) = g_bottom_new;
    sV(:,18) = g_altered_top;
    sV(:,19) = g_top_new_interp;
    
    %bottom values
    sV(:,20) = g_subsurf_new;
    sV(:,21) = g_subsurf_ice;
    sV(:,22) = g_impossibleDepth;
    sV(:,23) = g_weight_dens_bottom;
    sV(:,24) = g_bottom_density_ph;
    sV(:,25) = g_unreliablePhDens_bottom;
    sV(:,26) = g_bottom_new; 
    % test lake or not a lake
    sV(:,27) = g_topTest1;
    sV(:,28) = g_topTest2;
    sV(:,29) = g_botTest1;
    sV(:,30) = g_botTest2;
    sV(:,31) = g_botTest3;
    sV(:,32) = g_botTest4;
    sV(:,33) = g_stdDevBot;
    sV(:,34) = g_lake_score;
  
    writematrix(sV,saveNm);
    
end
