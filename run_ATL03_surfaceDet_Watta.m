function op = run_ATL03_surfaceDet_Watta()


%%% define necessary paths
    dirATL03='/YOUR_DATA_DIR_CONTAINING_ATL03_H5_files/'             % directory where the ATL03 .h5 files are stored
    dirATL06='/YOUR_DATA_DIR_CONTAINING_ATL06_H5_files/'             % directory where the ATL06 .h5 files are stored
    dirMain = '/YOUR_WATTA_OUTPUT_DIRECTORY/';                                 % this should be the main output folder, where all the data output will be stored in subdirectories
    dirOut = '/YOUR_DIRECTORY_WHERE_YOU_STORE_THE_FINAL_OUTPUT/';                     % the directory where the products of Step 1: surface detection and Step 2: interpretation are stored (can be identical to dirMain)
    
%%% specify the data to be processed    
    ProcessBeams = { ...                        % array containing tracks and beam in input directory to be processed. Naming convention yyyymmdd_ttttcc_rrr where tttt is the reference ground track (0001 to 1387), cc the cycle number and rrr the data release  
        '20190729_04830405_003', '1l'; ...      % These are just examples
        '20190729_04830405_003', '3l'; ...
    };

%%% STEP 0: preprocess the ATL03/ATL06 to .h5 files to .mat format for the region of interest (saved in subfolrder HeightFn, LocFn, PhConfFn, SaturFn, SurfTypeFn in the main     revision = '01';; % revision number, e.g. '01' for files ATL06_/ATL03_20190729171304_04830405_003_01.h5 
   startLatitude = 62; % starting latitude [-90:90]
   endLatitude = 65;  % endlatitude [-90:90] 
   revision = '01' % revision number of file, e.g. '01' for file ATL03_20190729171304_04830405_003_01.h5
   for iT = 1:size(ProcessBeams,1)
        read_ATL03_ATL06_files(dirATL03,dirATL06,ProcessBeams{iT,1},revision,ProcessBeams{iT,2},62,65,dirMain)
   end    

%%% loop through the data to be processed    
   for iT = 1:size(ProcessBeams,1)
        shortName = ProcessBeams{iT,1};
        beamN = sprintf('gt%s',ProcessBeams{iT,2});
        shortName_beam = sprintf('%s_%s',shortName, beamN);
       
        dirMain_file = sprintf('%s/%s', dirMain, shortName);
        
        %%% part I: surface detection
        dirOut_depth = sprintf('%s/%s/%s',dirOut,'outputDepths',shortName);
        mkdir(dirOut_depth);
        arraySpread = extract_tot_size(dirMain, shortName_beam, 30000);
        run_make_depthFiles(dirMain, dirOut_depth,shortName_beam, arraySpread);
        fileNameFin = sprintf('%s/DepthEst_%s.csv',dirOut_depth, shortName_beam);
        realign_overlaps(dirOut_depth, shortName_beam, arraySpread)
        writeFinalFiles(dirOut_depth, shortName_beam, fileNameFin, arraySpread); 
        
        %%% part II: interpretative part
        surface_characteristics_mn_tri(dirMain, dirOut, fileNameFin, shortName_beam)
   end          
end


function op = run_make_depthFiles( dirMain, dirOut, shortName, arraySpread )
   totCases = size(arraySpread,1);
   n_nodes = 12
   parpool(n_nodes)		   %%% run loop in parallel. Comment out if parellel computing is not desired
   parfor iT = 1: totCases         %%% create depth files for chunks of photons within start and end indices as defined in arraySpread
        startIndex = arraySpread(iT,1);
        endIndex = arraySpread(iT,2);
        msgDisp = sprintf('%d - %d', startIndex, endIndex);
        disp(msgDisp);
        stIndex = sprintf('Ind_%d', iT);
        de_calc_depth_estimates(dirMain,dirOut, shortName, startIndex, endIndex, stIndex);
                
        if( iT < totCases)      %%% make depth files for overlapping sections
            startIndex = endIndex-3000;
            endIndex = endIndex+3000;
            stInter = sprintf('Inter_it%d', iT);
            de_calc_depth_estimates(dirMain,dirOut, shortName, startIndex, endIndex, stInter);           
        end        
    end    
    delete(gcp('nocreate'));

end

function  op = realign_overlaps( dirOut, shortName, arraySpread)   
    %%% put segments back together
    numEntries_Per_File = 30000;
    totCases = size(arraySpread,1);
    midIndex = 3001;
    
    Iteration = 1;
    fl1 = sprintf('%s/SubDepthEst_%s_Ind_%d.csv', dirOut, shortName, Iteration);
    flInter2 = sprintf('%s/SubDepthEst_%s_Inter_it%d.csv',dirOut, shortName, Iteration);
    Index_1 = load(fl1);
    Inter_2 = load(flInter2);
    Index_1(end-750:end,:) = Inter_2(midIndex-750:midIndex,:);
    flFin = sprintf('%s/DepthEstSection_%s_%d.csv', dirOut, shortName, Iteration);
    writematrix(Index_1,flFin);
    clear('Inter_2', 'Inter_1', 'Index_1');
    
    for Iteration = 2: totCases-1       
        flInter1 = sprintf('%s/SubDepthEst_%s_Inter_it%d.csv',dirOut, shortName, Iteration-1);
        fl1 = sprintf('%s/SubDepthEst_%s_Ind_%d.csv', dirOut, shortName, Iteration);
        flInter2 = sprintf('%s/SubDepthEst_%s_Inter_it%d.csv',dirOut, shortName, Iteration);
        Index_1 = load(fl1);
        Inter_1 = load(flInter1);
        Index_1(1:750,:) = Inter_1(midIndex+1:midIndex+1+750-1,:);
        Inter_2 = load(flInter2);
        Index_1(end-750:end,:) = Inter_2(midIndex-750:midIndex,:);
        clear('Inter_2', 'Inter_1');
        flFin = sprintf('%s/DepthEstSection_%s_%d.csv', dirOut, shortName, Iteration);
        writematrix(Index_1,flFin);
    end
    Iteration = totCases;
    fl1 = sprintf('%s/SubDepthEst_%s_Ind_%d.csv', dirOut, shortName, Iteration);
    flInter1 = sprintf('%s/SubDepthEst_%s_Inter_it%d.csv',dirOut, shortName, Iteration-1);
    Index_1 = load(fl1);
    Inter_1 = load(flInter1);
    Index_1(1:750,:) = Inter_1(midIndex+1:midIndex+1+750-1,:);
    flFin = sprintf('%s/DepthEstSection_%s_%d.csv', dirOut, shortName, Iteration);
    writematrix(Index_1,flFin);
    clear('Inter_2', 'Inter_1', 'Index_1'); 
end

function op = writeFinalFiles(dirOut, shortName, finalFileName, arraySpread)  
    totCases = size(arraySpread,1);
    %%% write final results of Step I (surface detection) to file
    aFin = nan(1,29);
    for iT = 1: totCases                
        flFin = sprintf('%s/DepthEstSection_%s_%d.csv', dirOut, shortName, iT);
        cvLoad = load(flFin);
        aFin = [aFin; cvLoad];
        clear('cvLoad');
    end
    aFin(1,:) = [];
    writematrix(aFin,finalFileName);
end

function arraySpread = extract_tot_size(dirMain, shortName, numEntries_Per_File)
    %%% divide tracks in segments
    heightFl = sprintf('%s/HeightFn/HeightFn_%s.mat', dirMain, shortName);
    load(heightFl);
    heights_val = heights_val.heights;
    totSize = size(heights_val, 1);
    totCases = floor(totSize / numEntries_Per_File);
    aTot = 1:totCases;
    aStart = (aTot-1)*numEntries_Per_File +1;
    aEnd = aStart + numEntries_Per_File-1;
    aEnd(end) = totSize;
    
    arraySpread = nan(totCases,2);
    arraySpread(:,1) = aStart;
    arraySpread(:,2) = aEnd;
end


