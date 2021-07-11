function op = read_ATL03_ATL06_files()
    global gDir

    % create subdirectories and matfiles within from ATL03 and ATL06 files
    % which can then be used for the Watta algorithm

    gDir =  '';                             %set location for output files
    ATL0306Directory = '';                  %set location for ATL03/ATL06 files
    ATL03_filename = 'ATL03_20190729171304_04830405_003_01.h5'; %name of ATL03 file
    ATL06_filename = 'ATL06_20190729171304_04830405_003_01.h5'; %name of ATL06 file
    designation = '20190729_04830405_003';                      %a unique identifier for the IS2 file
    startLatitude = 62;                                         %starting latidue
    endLatitude = 65;                                           %ending latitude
    

    
    FILE_NAME = sprintf('%s/%s', ATL0306Directory, ATL03_filename);
    FILE_NAME_2 = sprintf('%s/%s', ATL0306Directory,ATL06_filename );
        
    run_all_beams_to_mat_ATL03(FILE_NAME, designation, startLatitude, endLatitude );
    run_all_beams_to_mat_ATL06(FILE_NAME,FILE_NAME_2, designation,startLatitude, endLatitude);
  
end

function op = run_all_beams_to_mat_ATL03( FILE_NAME, subT, stLat, endLat)
    global gStartLat gEndLat
    
    gStartLat = stLat;
    gEndLat = endLat;
    
    beams = {'gt1l'; 'gt1r'; 'gt2l'; 'gt2r'; 'gt3l'; 'gt3r'};
    for i=1:6
        beamN = beams{i};
        msgR = sprintf('Running %s %s', subT, beamN);
        disp(msgR);
        [stLatInd, endLatInd] = read_file_all_loc( FILE_NAME, beamN, subT );
        read_file_all_heights( FILE_NAME, beamN, subT, stLatInd, endLatInd );
        read_file_sat( FILE_NAME, beamN, subT, stLatInd, endLatInd );
        read_file_all_phconf(FILE_NAME, beamN, subT, stLatInd, endLatInd );
        read_file_surf_type( FILE_NAME, beamN, subT, stLatInd, endLatInd );
%        read_file_segment(FILE_NAME, beamN, subT );
    end
end



function op = run_all_beams_to_mat_ATL06( FILE_NAME,FILE_NAME_ATL06, subT, stLat, endLat)
    global gStartLat gEndLat
    
    gStartLat = stLat;
    gEndLat = endLat;
    
    beams = {'gt1l'; 'gt1r'; 'gt2l'; 'gt2r'; 'gt3l'; 'gt3r'};
    for i=1:6
        beamN = beams{i};
        msgR = sprintf('Running %s %s', subT, beamN);
        disp(msgR);
        [stLatInd, endLatInd] = read_file_all_loc( FILE_NAME, beamN, subT );
        read_file_segment_ATL06(FILE_NAME, FILE_NAME_ATL06, beamN, subT, stLatInd, endLatInd );
    end
end
function op = read_file_segment_ATL06( FILE_NAME, FILE_NAME_2, gtl_name, outStr, stLatInd, endLatInd )
    global gDir
    fileNameOut = sprintf('%s/Heights_ATL06/ATL06hts_%s_%s', gDir, outStr, gtl_name);
    fileNameOutQu = sprintf('%s/Heights_ATL06/ATL06Qu_%s_%s',gDir, outStr, gtl_name);
    
    file_id = H5F.open (FILE_NAME, 'H5F_ACC_RDONLY', 'H5P_DEFAULT');
    % Open the datasets.
    HEIGHTS_NAME=sprintf('%s/heights/h_ph', gtl_name);
    heights_id=H5D.open(file_id, HEIGHTS_NAME);

    PH_INDEX_BEG_NAME=sprintf('%s/geolocation/ph_index_beg', gtl_name);
    ph_index_beg_id=H5D.open(file_id, PH_INDEX_BEG_NAME);

    SEGMENT_PH_CNT_NAME=sprintf('%s/geolocation/segment_ph_cnt', gtl_name);
    seg_ph_cnt_id=H5D.open(file_id, SEGMENT_PH_CNT_NAME);

    SEGMENT_ID_NAME=sprintf('%s/geolocation/segment_id', gtl_name);
    segment_id_id=H5D.open(file_id, SEGMENT_ID_NAME);


    % Read the datasets.
    heights=H5D.read(heights_id,'H5T_NATIVE_DOUBLE', 'H5S_ALL', 'H5S_ALL', ...
                 'H5P_DEFAULT');
    ph_index_beg=H5D.read(ph_index_beg_id,'H5T_NATIVE_DOUBLE', 'H5S_ALL', 'H5S_ALL', ...
                 'H5P_DEFAULT');
    seg_ph_cnt=H5D.read(seg_ph_cnt_id,'H5T_NATIVE_DOUBLE', 'H5S_ALL', 'H5S_ALL', ...
                 'H5P_DEFAULT');
    seg_id=H5D.read(segment_id_id,'H5T_NATIVE_DOUBLE', 'H5S_ALL', 'H5S_ALL', ...
                 'H5P_DEFAULT');

    % Close and release resources.
    H5D.close (segment_id_id);
    H5D.close (seg_ph_cnt_id);
    H5D.close (ph_index_beg_id);

    H5F.close (file_id);

    %now the ATL06
    file_id = H5F.open (FILE_NAME_2, 'H5F_ACC_RDONLY', 'H5P_DEFAULT');

    % Open the datasets.

    HEIGHTS_ATL06_NAME=sprintf('%s/land_ice_segments/h_li', gtl_name);
    heights_atl06_id=H5D.open(file_id, HEIGHTS_ATL06_NAME);
    
    HEIGHTS_ATL06_QNAME_NAME=sprintf('%s/land_ice_segments/atl06_quality_summary', gtl_name);
    quality_atl06_id=H5D.open(file_id, HEIGHTS_ATL06_QNAME_NAME);
    

    SEG_ID_ATL06_NAME=sprintf('%s/land_ice_segments/segment_id', gtl_name);
    segid_atl06_id=H5D.open(file_id, SEG_ID_ATL06_NAME);

    % Read the datasets.
     atl06_heights=H5D.read(heights_atl06_id,'H5T_NATIVE_DOUBLE', 'H5S_ALL', 'H5S_ALL', ...
                 'H5P_DEFAULT');

     atl06_quality=H5D.read(quality_atl06_id,'H5T_NATIVE_DOUBLE', 'H5S_ALL', 'H5S_ALL', ...
                 'H5P_DEFAULT');

     atl06_segids=H5D.read(segid_atl06_id,'H5T_NATIVE_DOUBLE', 'H5S_ALL', 'H5S_ALL', ...
                 'H5P_DEFAULT');

     % Close and release resources.
    H5D.close (heights_atl06_id);
    H5D.close (quality_atl06_id);
    H5D.close (segid_atl06_id);
    H5F.close (file_id);

    begin_segment_interp = nan(size(heights));

    sZ = size(ph_index_beg,1);
    sz06 = size(atl06_segids,1);
    atl06_height_interp = nan(size(heights));
    atl06_quality_interp = nan(size(heights));
    
    for iT=1:sz06
%        msg = sprintf('Does: %d', iT);
%        disp(msg);
        %find atl03 segs
        atl06seg = atl06_segids(iT);
        ind2 = find(seg_id == atl06seg);
        ind1 = find(seg_id == (atl06seg-1));

        tBeg = ph_index_beg(ind1);
        if(tBeg == 0) 
            tBeg = ph_index_beg(ind2);
            if( tBeg == 0)
                continue;
            end
        end
        tEnd = ph_index_beg(ind2) + seg_ph_cnt(ind2);
        atl06_height_interp(tBeg:tEnd-1) = repmat(atl06_heights(iT),1,tEnd-tBeg);
        atl06_quality_interp(tBeg:tEnd-1) = repmat(atl06_quality(iT),1,tEnd-tBeg);
    end
    atl06heights.atl06_height_interp = atl06_height_interp(stLatInd:endLatInd);
    atl06quality.atl06_quality_interp = atl06_quality_interp(stLatInd:endLatInd);

    save(fileNameOut, 'atl06heights');
    save(fileNameOutQu, 'atl06quality');

end


function op = read_file_sat( FILE_NAME, gtl_name, outStr,stLatInd, endLatInd)
    global gDir


    fileNameOut_65 = sprintf('%s/SaturFn/SaturFn_%s_%s',gDir, outStr, gtl_name);
    
    file_id = H5F.open (FILE_NAME, 'H5F_ACC_RDONLY', 'H5P_DEFAULT');

    % Open the datasets.
    HEIGHTS_NAME=sprintf('%s/heights/h_ph', gtl_name);
    heights_id=H5D.open(file_id, HEIGHTS_NAME);
    
    PH_INDEX_BEG_NAME=sprintf('%s/geolocation/ph_index_beg', gtl_name);
    ph_index_beg_id=H5D.open(file_id, PH_INDEX_BEG_NAME);

    SEGMENT_PH_CNT_NAME=sprintf('%s/geolocation/segment_ph_cnt', gtl_name);
    seg_ph_cnt_id=H5D.open(file_id, SEGMENT_PH_CNT_NAME);

    FULL_SAT_FRACT_NAME=sprintf('%s/geolocation/full_sat_fract', gtl_name);
    full_sat_fract_id=H5D.open(file_id, FULL_SAT_FRACT_NAME);

    NEAR_SAT_FRACT_NAME=sprintf('%s/geolocation/near_sat_fract', gtl_name);
    near_sat_fract_id=H5D.open(file_id, NEAR_SAT_FRACT_NAME);

    % Read the datasets.
    heights=H5D.read(heights_id,'H5T_NATIVE_DOUBLE', 'H5S_ALL', 'H5S_ALL', ...
                 'H5P_DEFAULT');
    ph_index_beg=H5D.read(ph_index_beg_id,'H5T_NATIVE_DOUBLE', 'H5S_ALL', 'H5S_ALL', ...
                 'H5P_DEFAULT');
    seg_ph_cnt=H5D.read(seg_ph_cnt_id,'H5T_NATIVE_DOUBLE', 'H5S_ALL', 'H5S_ALL', ...
                 'H5P_DEFAULT');
    full_sat_fract=H5D.read(full_sat_fract_id,'H5T_NATIVE_DOUBLE', 'H5S_ALL', 'H5S_ALL', ...
                 'H5P_DEFAULT');
    near_sat_fract=H5D.read(near_sat_fract_id,'H5T_NATIVE_DOUBLE', 'H5S_ALL', 'H5S_ALL', ...
                 'H5P_DEFAULT');

    % Close and release resources.
    H5D.close (seg_ph_cnt_id);
    H5D.close (near_sat_fract_id);
    H5D.close (full_sat_fract_id);
    H5D.close (ph_index_beg_id);

    H5F.close (file_id);

    full_sat_fract_interp = nan(size(heights));
    near_sat_fract_interp = nan(size(heights));
    begin_segment = nan(size(heights));

    sZ = size(ph_index_beg,1);
    
    for iT=1:sZ
        tBeg = ph_index_beg(iT);
        begin_segment(iT) = 1;
        tAdd = seg_ph_cnt(iT);
        full_sat_fract_interp(tBeg:tBeg+tAdd-1) = repmat(full_sat_fract(iT),1,tAdd);
        near_sat_fract_interp(tBeg:tBeg+tAdd-1) = repmat(near_sat_fract(iT),1,tAdd);
    end
    
    sat_val.full_sat_fract_interp = full_sat_fract_interp(stLatInd:endLatInd);
    sat_val.near_sat_fract_interp = near_sat_fract_interp(stLatInd:endLatInd);
    save(fileNameOut_65, 'sat_val');
end


function op = read_file_surf_type( FILE_NAME, gtl_name, outStr, stLatInd, endLatInd)
    global gDir


    %land, ocean, sea ice, land
    %ice and inland water
    
    fileNameOut_65 = sprintf('%s/SurfTypeFn/SurfTypeFn_%s_%s', gDir, outStr, gtl_name);

    file_id = H5F.open (FILE_NAME, 'H5F_ACC_RDONLY', 'H5P_DEFAULT');

    % Open the datasets.
    HEIGHTS_NAME=sprintf('%s/heights/h_ph', gtl_name);
    heights_id=H5D.open(file_id, HEIGHTS_NAME);

    PH_INDEX_BEG_NAME=sprintf('%s/geolocation/ph_index_beg', gtl_name);
    ph_index_beg_id=H5D.open(file_id, PH_INDEX_BEG_NAME);

    SEGMENT_PH_CNT_NAME=sprintf('%s/geolocation/segment_ph_cnt', gtl_name);
    seg_ph_cnt_id=H5D.open(file_id, SEGMENT_PH_CNT_NAME);

    SURF_TYPE_NAME=sprintf('%s/geolocation/surf_type', gtl_name);
    surf_type_id=H5D.open(file_id, SURF_TYPE_NAME);

    % Read the datasets.
    heights=H5D.read(heights_id,'H5T_NATIVE_DOUBLE', 'H5S_ALL', 'H5S_ALL', ...
                 'H5P_DEFAULT');
    ph_index_beg=H5D.read(ph_index_beg_id,'H5T_NATIVE_DOUBLE', 'H5S_ALL', 'H5S_ALL', ...
                 'H5P_DEFAULT');
    seg_ph_cnt=H5D.read(seg_ph_cnt_id,'H5T_NATIVE_DOUBLE', 'H5S_ALL', 'H5S_ALL', ...
                 'H5P_DEFAULT');
    surf_type=H5D.read(surf_type_id,'H5T_NATIVE_DOUBLE', 'H5S_ALL', 'H5S_ALL', ...
                 'H5P_DEFAULT');

    % Close and release resources.
    H5D.close (seg_ph_cnt_id);
    H5D.close (surf_type_id);
    H5D.close (ph_index_beg_id);

    H5F.close (file_id);

    surf_type_interp = nan(5,size(heights,1));
    sZ = size(ph_index_beg,1);
    for iT=1:sZ
        tBeg = ph_index_beg(iT);
        tAdd = seg_ph_cnt(iT);
        surf_type_interp(:,tBeg:tBeg+tAdd-1) = repmat(surf_type(:,iT),1,tAdd);
    end
    
    surftype_val.surf_type_interp = surf_type_interp(stLatInd:endLatInd);
    save(fileNameOut_65, 'surftype_val');
end


function op = read_file_all_heights( FILE_NAME, gtl_name, outStr,stLatInd, endLatInd)
    global gDir

    fileNameOut_65 = sprintf('%s/HeightFn/HeightFn_%s_%s',gDir, outStr, gtl_name);

    file_id = H5F.open (FILE_NAME, 'H5F_ACC_RDONLY', 'H5P_DEFAULT');

    % Open the datasets.

    HEIGHTS_NAME=sprintf('%s/heights/h_ph', gtl_name);
    heights_id=H5D.open(file_id, HEIGHTS_NAME);

    % Read the datasets.
    heights=H5D.read(heights_id,'H5T_NATIVE_DOUBLE', 'H5S_ALL', 'H5S_ALL', ...
                 'H5P_DEFAULT');

    % Close and release resources.
    H5D.close (heights_id);
    H5F.close (file_id);

    
    heights_val.heights = heights(stLatInd:endLatInd);
    save(fileNameOut_65, 'heights_val');
        


end

function op = read_file_all_phconf( FILE_NAME, gtl_name, outStr, stLatInd, endLatInd)
    global gDir


    fileNameOut_65 = sprintf('%s/PhConfFn/PhConfFn_%s_%s', gDir, outStr, gtl_name);
    file_id = H5F.open (FILE_NAME, 'H5F_ACC_RDONLY', 'H5P_DEFAULT');

    % Open the datasets.
    PH_CONF_NAME=sprintf('%s/heights/signal_conf_ph', gtl_name);
    ph_conf_id=H5D.open(file_id, PH_CONF_NAME);
    % Read the datasets.
    ph_conf=H5D.read(ph_conf_id,'H5T_NATIVE_DOUBLE', 'H5S_ALL', 'H5S_ALL', ...
                 'H5P_DEFAULT');
    % Close and release resources.
    H5D.close (ph_conf_id);

    H5F.close (file_id);

    ph_conf_val.ph_conf = ph_conf(stLatInd:endLatInd);
    save(fileNameOut_65, 'ph_conf_val');

end

function [stLatInd, endLatInd ] = read_file_all_loc( FILE_NAME, gtl_name, outStr)
    global gStartLat gEndLat gDir



%    fileNameOut = sprintf('Loc_%s_%s', outStr, gtl_name);
    fileNameOut_65 = sprintf('%s/LocFn/LocFn_%s_%s', gDir, outStr, gtl_name);
    
    file_id = H5F.open (FILE_NAME, 'H5F_ACC_RDONLY', 'H5P_DEFAULT');

    % Open the datasets.
    LATFIELD_NAME= sprintf('%s/heights/lat_ph', gtl_name);
    lat_id=H5D.open(file_id, LATFIELD_NAME);

    LONFIELD_NAME=sprintf('%s/heights/lon_ph', gtl_name);
    lon_id=H5D.open(file_id, LONFIELD_NAME);


    % Read the datasets.
    latitude=H5D.read(lat_id,'H5T_NATIVE_DOUBLE', 'H5S_ALL', 'H5S_ALL',...
                    'H5P_DEFAULT');
    longitude=H5D.read(lon_id,'H5T_NATIVE_DOUBLE', 'H5S_ALL', 'H5S_ALL', ...
                 'H5P_DEFAULT');
             
    ascending = latitude(100) > latitude(2);
    if( ascending)   %ascending
        stLatInd =  max(find(latitude < gStartLat));
        endLatInd =  min(find(latitude >=  gEndLat));
        loc_val.latitide = latitude(stLatInd:endLatInd);
        loc_val.longitude = longitude(stLatInd:endLatInd);
        save(fileNameOut_65, 'loc_val');
    else        %descending
        stLatInd =  max(find(latitude > gEndLat));        
        endLatInd =  min(find(latitude <=  gStartLat));
        loc_val.latitide = latitude(stLatInd:endLatInd);
        loc_val.longitude = longitude(stLatInd:endLatInd);        
        save(fileNameOut_65, 'loc_val');

    end
    
    % Close and release resources.
    H5D.close (lon_id);
    H5D.close (lat_id);


    H5F.close (file_id);    
end

