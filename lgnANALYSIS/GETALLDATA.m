function [ROISCAT, DISTANCECAT]=  GETALLDATA

%%  COMBINE DATA TABLES FROM ALL IMAGING FIELDS


%% DATA FILES 2 PROCESS: HARD CODED
mouseNameList           = {'lgn02', 'lgn03', 'lgn04', 'lgn05'};
imfield_expNames_All    = cell(numel(mouseNameList),1);
 
% Lgn02 mouse
imfield_expNames_All{1}     = {'_000_004', '_000_003';... % 1
                              '_001_000', '_001_001';... %2
                              '_002_000', '_002_002';... % this 3rd field has a really bad match
                              '_003_000', '_003_001'};   %4
% Lgn03 mouse                        
imfield_expNames_All{2}     = {'_001_001', '_001_000';...%5
                               '_001_003', '_001_002';...%6
                               '_002_001', '_002_000'}; %7
% Lgn04 mouse                                               
imfield_expNames_All{3}     = {'_002_000', '_002_003'};  %8  

% Lgn05 mouse  
imfield_expNames_All{4}     = {'_001_001', '_001_000';... %9
                                '_008_001', '_008_000';... %10
                                '_011_000', '_011_001'};  %11       

%% CREATE TABLE  
ROISCAT     = [];
DISTANCECAT = [];
unique_Im_field_number = 0;
for i = 1:numel(mouseNameList)
    for j = 1:size(imfield_expNames_All{i},1)
        
        mouse = mouseNameList{i};
        houghExpId = imfield_expNames_All{i}{j,1};
        orisfExpId = imfield_expNames_All{i}{j,2};
        unique_Im_field_number = unique_Im_field_number + 1;
       
        [rois, uifInfo] = getRoiTable(mouse, houghExpId, orisfExpId, unique_Im_field_number);
        ROISCAT     = [ROISCAT;rois];   
    end
end

ROISCAT = [table((1:height(ROISCAT))', 'VariableNames', {'masterEntry'} ), ROISCAT];




for uif = 1:numel(unique(ROISCAT.uniqImFieldNum))
    ROIS = ROISCAT(ROISCAT.uniqImFieldNum == uif,:);
    [distances, distSelection] = getDistanceTable(ROIS);  
    DISTANCECAT = [DISTANCECAT; distances];        
end

close all
DISTANCECAT = [table((1:height(DISTANCECAT))', 'VariableNames', {'masterEntry'} ), DISTANCECAT];   



end

%% CREATE ROI TABLE FOR 1 IMAGING FIELD
function [ROIS, UIFinfo] = getRoiTable(mouse, houghExpId, orisfExpId, uniqueImagFieldNumber)

% mouse       = 'lgn02';
% houghExpId  = '_000_004';
% orisfExpId  = '_000_003';
% close all;
      
    % Format File Names    
    RET_fileNameBase    = [lojGetDataPath, mouse,'/', mouse, houghExpId];   
    kernRET_fName    = [RET_fileNameBase , '_rigid.houghkernstats'];
    quadRET_fName    = [RET_fileNameBase , '_quadrature.mat'];
    segmentRET_fName = [RET_fileNameBase , '_rigid.segment'];
    alignRET_fName   = [RET_fileNameBase,  '.align'];
    spksRET_fName    = [RET_fileNameBase , '_rigid.signals'];
    
    TUN_fileNameBase    = [lojGetDataPath, mouse,'/', mouse, orisfExpId];
    kernTUN_fName    = [TUN_fileNameBase ,'_rigid.orisf'];
    quadTUN_fName    = [TUN_fileNameBase , '_quadrature.mat'];
    segmentTUN_fName = [TUN_fileNameBase , '_rigid.segment'];
    alignTUN_fName   = [TUN_fileNameBase,  '.align'];
    spksTUN_fName    = [TUN_fileNameBase , '_rigid.signals'];
    

    % Get Kernels from this population
    RET_Kern     = load(kernRET_fName, '-mat');   
    RET_Kern     = RET_Kern.S;
    TUN_Kern     = load(kernTUN_fName, '-mat');   
    TUN_Kern     = TUN_Kern.stat;
    
    % Get Quad Data From This Population
    if exist(quadRET_fName, 'file')
        RET_quad = load(quadRET_fName, '-mat');   
        RET_quad = {double(RET_quad.quad_data)}; 
        hasQuadData_Ret = true;
    else
        RET_quad = {nan};
        hasQuadData_Ret = false;               
    end
    
    if exist(quadTUN_fName, 'file')
        TUN_quad = load(quadTUN_fName, '-mat');   
        TUN_quad = {double(TUN_quad.quad_data)};   
        hasQuadData_Tun = true;                       
    else
        TUN_quad = {nan};
        hasQuadData_Tun = false;                       
    end
        
    % Load Segment Images
    RET_seg = load(segmentRET_fName, '-mat', 'mask'); 
    RET_seg = RET_seg.mask;
    TUN_seg = load(segmentTUN_fName, '-mat', 'mask'); 
    TUN_seg = TUN_seg.mask;
    
    % Load Alignment Images
    % Get Quad Data From This Population
    if exist(alignRET_fName, 'file')
        RET_align = load(alignRET_fName, '-mat', 'm'); 
        RET_align = RET_align.m;
    else
        RET_align = {nan};
    end
    
    if exist(alignTUN_fName, 'file')
        TUN_align = load(alignTUN_fName, '-mat', 'm'); 
        TUN_align = TUN_align.m;    
    else
        TUN_align = {nan};
    end 


    
    
    % Match rois in this population
    MTC         = sbxmatchfields([RET_fileNameBase, '_rigid'], [TUN_fileNameBase,'_rigid'],.2);
    
    % Create Table for Population
    S         = struct2table(RET_Kern);
    S         = S(:,{'kern','kurtmax', 'sig', 'xy'});
    S.Properties.VariableNames = {'kern_ret', 'kurt_ret', 'sig_ret', 'xy_ret'};
    S.sig_ret = num2cell(S.sig_ret,2);
    
    T         =  struct2table(TUN_Kern);
    T         = T(:,{'kern_tune', 'sig_tune', 'signal_orisf', 'oriresp_tune', 'oriest_tune', 'sfresp_tune', 'osi_tune', 'sfest_tune', 'xy_soma', 'area_soma'});
    T.Properties.VariableNames = {'kern_tun', 'signif_tun', 'sig_tun', 'oricurve', 'oriest', 'sfcurve', 'osi', 'sfPeak',  'segment_xy_tun', 'segment_area_tun'};
   
    
    % Create No Match Table For Retinotpoy Kernels  
    rois_haveRetData    = array2table(MTC.AnotB', 'variablenames', {'roi_ID_Ret'});
    num_haveRetData     = height(rois_haveRetData);   
    hasRetData         = table(true(num_haveRetData ,1), 'variablenames', {'hasRetData'});

    rois_haveTunData    = array2table(nan(num_haveRetData, 1), 'variablenames', {'roi_ID_Tun'});    
    num_haveTunData     = height(rois_haveTunData);  
    hasTunData         = table(false(num_haveTunData ,1), 'variablenames', {'hasTunData'});
    
    K_onlyRetData       = [hasRetData, hasTunData, rois_haveRetData, rois_haveTunData,...
                                [S(MTC.AnotB,:), makeEmptyTable(T(1,:),num_haveRetData)] ]; 
 
    % Create No Match Table For Tuning Kernels                              
    rois_haveTunData    = array2table(MTC.BnotA', 'variablenames', {'roi_ID_Tun'}); 
    num_haveTunData     = height(rois_haveTunData);         
    hasTunData         = table(true(num_haveTunData ,1), 'variablenames', {'hasTunData'});
    
    rois_haveRetData    = array2table(nan(num_haveTunData, 1), 'variablenames', {'roi_ID_Ret'});    
    num_haveRetData     = height(rois_haveRetData);         
    hasRetData         = table(false(num_haveRetData ,1), 'variablenames', {'hasRetData'});
    
    K_onlyTunData       = [hasRetData, hasTunData, rois_haveRetData, rois_haveTunData,...
                                [makeEmptyTable(S(1,:), num_haveTunData), T(MTC.BnotA,:)] ]; 
     
                            
    % Create Match Table                    
    rois_haveRetData    = array2table(MTC.match(:,1), 'variablenames', {'roi_ID_Ret'});
    rois_haveTunData    = array2table(MTC.match(:,2), 'variablenames', {'roi_ID_Tun'});
    num_haveRetData     = height(rois_haveRetData);         
    num_haveTunData     = height(rois_haveTunData);         
    hasRetData         = table(true(num_haveRetData ,1), 'variablenames', {'hasRetData'});
    hasTunData         = table(true(num_haveTunData ,1), 'variablenames', {'hasTunData'});
    K_haveRetAndTuneData       = [hasRetData, hasTunData, rois_haveRetData, rois_haveTunData,...
                                [[S(MTC.match(:,1),:) T(MTC.match(:,2),:)]]];   

    % Combine All & assign Mouse Id, Unique Population ID, and Append
    ROIS     = [K_onlyRetData; K_onlyTunData; K_haveRetAndTuneData];
    mse         = table(repmat({mouse}, height(ROIS) ,1), 'variablenames', {'mouse'}); 
    uniqImFieldNum  = table(repmat(uniqueImagFieldNumber, height(ROIS) ,1), 'variablenames', {'uniqImFieldNum'});
    entryNumber   = table((1:height(ROIS))', 'variablenames', {'entryNumber'});
    
    ROIS     = [entryNumber, mse, uniqImFieldNum, ROIS] ;
    
%     Create Imaging Field
    UIFinfo = table({mouse},  uniqueImagFieldNumber,  {houghExpId},    {orisfExpId},...
                            hasQuadData_Ret,  RET_quad,  hasQuadData_Tun, TUN_quad,...
                            {RET_seg}, {TUN_seg}, {RET_align},{TUN_align}, ...
        'variablenames',{'mouseName', 'uniqImFieldNum', 'unit_exp_hough','unit_exp_orisf',...
                        'hasQuadData_Ret', 'hough_quad',   'hasQuadData_Orisf', 'orisf_quad',...
                        'hough_segMask', 'orisf_segMask', 'hough_alignImage', 'orisf_alignImage'});
                    
                    
    function T_empty = makeEmptyTable(table2rep, heightOfTable)

        for i = 1:size(table2rep,2)
            switch class(table2rep{1,i})
                case 'cell'
                    table2rep{1,i} = {nan};
                case {'single', 'double'}
                    table2rep{1,i} = nan;
                case 'logical'
                    table2rep{1,i} = false;
                otherwise 
                    error('');
            end    
        end

        T_empty = repmat(table2rep, heightOfTable,1);

    end    

    
    function fixMatchTable(MTC)
    end      

                    
end

%% CREATE DISTANCE TABLE FOR AN IMAGING FIELD
function [D, roiConstraints] = getDistanceTable(ROIS)

isInUpperQuart = @(x) x(:)>prctile(x(:), 75);
isInLowerQuart = @(x) x(:)<prctile(x(:), 25);
isInMidFifty  = @(x) x(:)>prctile(x(:), 25) & x(:)<prctile(x(:), 75);

% SELECT ROIS FOR DISTANCE TABLE 
% Must have tuning and ret data, and tuning and ret kernels must be significant

roiConstraints = ROIS.hasRetData &...
    ROIS.hasTunData &...
    ROIS.kurt_ret>5 &...
    ROIS.signif_tun;

R = ROIS(roiConstraints,:);

% Compute Distance
D    = [];
    
        % Compute Retinotopic Overlap Measures
        dRetKern_c2cDist    = table(computeDistance(R.kern_ret, 'norm_dRfCenters'), 'variablenames', {'dRetKern_c2cDist'});   % Find Retinotopy Kernel Similirity
    
        dRetKern_corr       = table(computeDistance(R.kern_ret, 'kerncorr'), 'variablenames', {'dRetKern_corr'});   % Find Retinotopy Kernel Similirity        
        dRetKern_corrStrong = table(dRetKern_corr.dRetKern_corr>.55, 'variablenames', {'dRetKern_corrStrong'});   % Find Retinotopy Kernel Similirity        
        dRetKern_corrWeak   = table(dRetKern_corr.dRetKern_corr<.15, 'variablenames', {'dRetKern_corrWeak'});   % Find Retinotopy Kernel Similirity        
        
        dRetKern_Overlap      = table(computeDistance(R.kern_ret, 'overlap'), 'variablenames', {'dRetKern_Overlap'});   % Find Retinotopy Kernel Similirity
        dRetKern_OverlapFiftyPrc = table(dRetKern_Overlap.dRetKern_Overlap > .5, 'variablenames', {'dRetKern_OverlapFiftyPrc'}); % pairs that have high ret overlap         
        dRetKern_OverlapYes   = table(dRetKern_Overlap.dRetKern_Overlap > 0, 'variablenames', {'dRetKern_OverlapYes'}); % pairs that have high ret overlap 
        dRetKern_OverlapNo    = table(dRetKern_Overlap.dRetKern_Overlap == 0, 'variablenames', {'dRetKern_OverlapNo'}); % pairs that do not overlap
             
        % Compute Tuning Similarity Measures
        dTunKern_cos  = table(computeDistance(R.kern_tun, 'cosangle'), 'variablenames', {'dTunKern_cos'}); % Tuning Kernel  Similirity
        dTunKern_corr = table(computeDistance(R.kern_tun, 'kerncorr'), 'variablenames', {'dTunKern_corr'}); % Tuning Kernel  Similirity
       
        dSfResp            = table(computeDistance(R.sfcurve,'respcorr')', 'variablenames', {'dSfResp'}); % Diff in Log of SF        
        dLogSfEst           = table(computeDistance(log2(R.sfPeak),'absdiff'), 'variablenames', {'dLogSfEst'}); % Diff in Log of SF
        
        dLogSfEst_lessThanHalfOct          = table(dLogSfEst.dLogSfEst <= .5, 'variablenames', {'dLogSfEst_lessThanHalfOct'}); % Pairs with Same SF (diff < 1 octave)        
        dLogSfEst_lessThanOneOct          = table(dLogSfEst.dLogSfEst <= 1, 'variablenames', {'dLogSfEst_lessThanOneOct'}); % Pairs with Same SF (diff < 1 octave)
        dLogSfEst_lessThanTwoOct  = table(dLogSfEst.dLogSfEst <= 2, 'variablenames', {'dLogSfEst_lessThanTwoOct'}); % Pairs with Same SF (diff < 1 octave)
        dLogSfEst_lessThanThreeOct      = table(dLogSfEst.dLogSfEst <= 3, 'variablenames', {'dLogSfEst_lessThanThreeOct'}); % Pairs with Same SF (diff < 1 octave)
        dLogSfEst_AnyDiff           = table(dLogSfEst.dLogSfEst <= max(dLogSfEst.dLogSfEst), 'variablenames', {'dLogSfEst_AnyDiff'}); % Pairs with Different SF (diff > 1 octave)
        dLogSfEst_moreThanHalfOct          = table(dLogSfEst.dLogSfEst > .5, 'variablenames', {'dLogSfEst_moreThanHalfOct'}); % Pairs with Same SF (diff < 1 octave)                
        dLogSfEst_moreThanOneOct      = table(dLogSfEst.dLogSfEst > 1, 'variablenames', {'dLogSfEst_moreThanOneOct'}); % Pairs with Same SF (diff < 1 octave)
        dLogSfEst_moreThanTwoOct      = table(dLogSfEst.dLogSfEst > 2, 'variablenames', {'dLogSfEst_moreThanTwoOct'}); % Pairs with Same SF (diff < 1 octave)
        dLogSfEst_moreThanThreeOct      = table(dLogSfEst.dLogSfEst > 3, 'variablenames', {'dLogSfEst_moreThanThreeOct'}); % Pairs with Same SF (diff < 1 octave)
       
        dOriResp_corr       = table(computeDistance(R.oricurve,'respcorr')', 'variablenames', {'dOriResp_corr'}); % Diff in Log of SF
        dOriResp_corrStrong = table(dOriResp_corr.dOriResp_corr>.55, 'variablenames', {'dOriResp_corrStrong'});   % Find Retinotopy Kernel Similirity        
        dOriResp_corrWeak   = table(dOriResp_corr.dOriResp_corr<.15, 'variablenames', {'dOriResp_corrWeak'});   % Find Retinotopy Kernel Similirity        
                
        
        dOriEst             = table(computeDistance(R.oriest,'distcirc'), 'variablenames', {'dOriEst'}); % Diff in Log of SF
        dOriEst_Same        = table(dOriEst.dOriEst<=15, 'variablenames', {'dOriEst_Same'}); % Diff in Log of SF
        dOriEst_Different   = table(dOriEst.dOriEst>=45, 'variablenames', {'dOriEst_Different'}); % Diff in Log of SF
       
        
        % Compute Roi Segment Position Distance
        dSegRoi = computeDistance({[R.segment_xy_tun, R.segment_area_tun]}, 'distance2dNorm');
        dSegRoi = table(dSegRoi, 'variablenames', {'dSegRoi'});
        
        % Get Roi Pairing Info
        [~, ii, jj]   = computeDistance(R.sfPeak,'absdiff'); % RetrunedsPair Index
        pair_roiMasterEntry = table([R.masterEntry(ii), R.masterEntry(jj)], 'VariableNames', {'pair_roiMasterEntry'});
        pair_roiUifEntry    = table([R.entryNumber(ii), R.entryNumber(jj)], 'VariableNames', {'pair_roiUifEntry'});
        pair_roiDistMatIdx  = table([ii,jj], 'VariableNames', {'pair_roiDistMatIdx'});
        pair_roiOrisfId     = table([R.roi_ID_Tun(ii), R.roi_ID_Tun(jj)], 'VariableNames', {'pair_roiOrisfId'});
        pair_roiHoughId     = table([R.roi_ID_Ret(ii), R.roi_ID_Ret(jj)], 'VariableNames', {'pair_roiHoughId'});        
        numPairs            =  height(pair_roiDistMatIdx);
        
        % Create Table
        mouse       = table(repmat(ROIS.mouse(1), numPairs ,1), 'variablenames', {'mouse'}); 
        uif         = table(repmat(ROIS.uniqImFieldNum(1), numPairs ,1), 'variablenames', {'uif'});
        uifEntry    = table((1:numPairs)', 'variablenames', {'uifEntry'});
        
        D = [uifEntry, mouse, pair_roiMasterEntry, uif, pair_roiUifEntry, pair_roiDistMatIdx, pair_roiOrisfId, pair_roiHoughId,...
                        dRetKern_c2cDist,...
                        dRetKern_corr, dRetKern_corrStrong, dRetKern_corrWeak,...
                        dRetKern_Overlap, dRetKern_OverlapFiftyPrc, dRetKern_OverlapYes, dRetKern_OverlapNo,...
                        dTunKern_cos, dTunKern_corr,...
                        dSfResp,...
                        dLogSfEst, dLogSfEst_lessThanHalfOct, dLogSfEst_lessThanOneOct , dLogSfEst_lessThanTwoOct ,dLogSfEst_lessThanThreeOct,dLogSfEst_AnyDiff,...
                                   dLogSfEst_moreThanHalfOct, dLogSfEst_moreThanOneOct, dLogSfEst_moreThanTwoOct, dLogSfEst_moreThanThreeOct,...
                        dOriResp_corr, dOriResp_corrStrong, dOriResp_corrWeak,...
                        dOriEst,    dOriEst_Same ,  dOriEst_Different,...
                        dSegRoi];    


end


%% DISTANCE MEASURMENTS 
function [D, ii, jj] = computeDistance(X, distanceMeasure)

switch distanceMeasure
    
    case 'overlap'
        D = pixOverlap(X);
    case 'cosangle'
        D = cosangle(X);
    case 'kerncorr'
        stretchKernel   = @(kernelCellArray) cellfun(@(k) k(:), kernelCellArray, 'UniformOutput' , false);
        stretchedRetKernels = stretchKernel(X);
        D = 1- pdist(cat(2,stretchedRetKernels{:})' , 'correlation')';  

    case 'respcorr'
        D = 1- pdist(X , 'correlation'); 
        
    case 'absdiff'
        D  = abs(takeDiff(X));
        [~, ii, jj] = takeDiff(X); 
    case 'distance2dNorm' % used for 
        xy      = X{:}(:, 1:2);
        area    = X{:}(:, 3);
        radius = (2*sqrt(mean(area)/pi));
        D = pdist(xy, 'euclidean')'/radius;   
        
    case 'norm_dRfCenters'
        D = norm_dRfCenters(X);
        
    case 'distcirc'
        D  = abs(takeDiff(X));
        D(D>90)   = 180 - D(D>90);        
        [~, ii, jj] = takeDiff(X); 
        
    otherwise
        error('Not a defined distance measure')
end

%% Cosinge Angle
function [d, ii, jj] = cosangle(KcellArray)
        Kstretched = cellfun(@(k) k(:), KcellArray, 'UniformOutput' , false);
        Kstretched = cat(2,Kstretched{:});


        numUnits = numel(KcellArray);
        numPairs = (numUnits *(numUnits-1))/2;
        [ii, jj] = find(tril(ones(numUnits),-1)) ;

        d = nan(numPairs,1);
        for p = 1:numPairs

            Ki = Kstretched(:, ii(p));
            Kj = Kstretched(:, jj(p));

            d(p) = sum(Ki.*Kj)/(norm(Ki)*norm(Kj));        
            sprintf('%d/%d', p, numPairs)
        end   
        
        
end

%% Kernel Pixel Overlap
function [D,  ii, jj] = pixOverlap(KcellArray)
    Kstretched = cellfun(@(k) k(:), KcellArray, 'UniformOutput' , false);
    Kstretched = cat(2,Kstretched{:});
    
    
    numUnits = numel(KcellArray);
    numPairs = (numUnits *(numUnits-1))/2;
    [ii, jj] = find(tril(ones(numUnits),-1)) ;
    
    nrm     = @(s) (s - min(s(:)))/ (max(s(:)) - min(s(:)) );

    D = nan(numPairs,1);
    for p = 1:numPairs
        
        Ki = Kstretched(:, ii(p));
        Kj = Kstretched(:, jj(p));
        
        Knrm_i  = nrm(Ki);
        Kidx_i  = find(Knrm_i >.75);

        Knrm_j    = nrm(Kj);
        Kidx_j  = find(Knrm_j >.75);

        D(p) = length(intersect(Kidx_i,Kidx_j))/length(union(Kidx_i,Kidx_j)) ;
        
        sprintf('%d/%d', p, numPairs)
    end
    
end

%% Differences
function [D, ii, jj] = takeDiff(array)

    
    numUnits = numel(array);
    numPairs = (numUnits *(numUnits-1))/2;
    [ii, jj] = find(tril(ones(numUnits),-1)) ;
    
    D = nan(numPairs,1);
    for p = 1:numPairs
        Ki = array(ii(p));
        Kj = array(jj(p));
        D(p) = (Ki-Kj);        
    end
    
end

%% Normalized Distance Between RF Centers
function [D, ii, jj] = norm_dRfCenters(X)
        numUnits = numel(X);
        numPairs = (numUnits *(numUnits-1))/2;
        [ii, jj] = find(tril(ones(numUnits),-1)) ;

        nrm     = @(s) (s - min(s(:)))/ (max(s(:)) - min(s(:)) );
        D = nan(numPairs,1);
        
        for p = 1:numPairs

            rpi = regionprops(nrm(X{ii(p)})>.85, 'Area', 'centroid');
            rpj = regionprops(nrm(X{jj(p)})>.85, 'Area', 'centroid');
            Ki = rpi( [rpi.Area] == max([rpi.Area]) );
            Kj = rpj( [rpj.Area] == max([rpj.Area]) );

            xi = Ki.Centroid(1);
            yi = Ki.Centroid(2);
            ri = (2*sqrt(mean(Ki.Area)/pi)) ;            

            xj = Kj.Centroid(1);
            yj = Kj.Centroid(2);
            rj = (2*sqrt(mean(Kj.Area)/pi)) ;            
            
            
            D(p) = sqrt((xi-xj)^2 + (yi-yj)^2) / mean([ri, rj]);
            
            sprintf('%d/%d', p, numPairs)
        end 

end






end