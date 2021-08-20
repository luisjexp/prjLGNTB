function DISTCAT = onLeave(ROISCAT)

ROISCAT(448,:) = [];
DISTCAT = [];
roisSelected = [];
for uif = 1:numel(unique(ROISCAT.uif))
    ROIS = ROISCAT(ROISCAT.uif == uif,:);
    [distances, rselection] = getDistanceTable(ROIS); 
    fprintf('\n%d rois (%d Pairs) for uif %d\n', height(rselection), height(distances), uif)
    DISTCAT = [DISTCAT; distances]; 
    roisSelected = [roisSelected;rselection];
end

DISTCAT = [table((1:height(DISTCAT))', 'VariableNames', {'masterEntry'} ), DISTCAT];   
end


%% CREATE DISTANCE TABLE FOR AN IMAGING FIELD

function [D, R] = getDistanceTable(ROIS)
% Categorization
isInUpperQuart = @(x) x(:)>prctile(x(:), 75);
isNotUpperTenth = @(x) x(:)<prctile(x(:), 90);
isInMidFifty  = @(x)  x(:)>prctile(x(:), 25) & x(:)<prctile(x(:), 75);
isInLowerQuart = @(x) x(:)<prctile(x(:), 25);



% SELECT ROIS FOR DISTANCE TABLE 
%  tuning and ret kernels must be significant
% roiConstraints = ROIS.kurt_ret>7 & ROIS.signif_tun;

roiConstraints =    ROIS.kurt_ret>13 & ROIS.signif_tun &...
    ROIS.kern_ret_binCentroid(:,1) > size(ROIS.kern_ret{2},2)/2-30;

% roiConstraints =    ROIS.kurt_ret>13 &...
%     ROIS.kern_ret_binCentroid(:,1) > size(ROIS.kern_ret{2},2)/2-30;
R = ROIS(roiConstraints,:);
    

% Compute Retinotopic Overlap 
dRetKern_c2cDist        = table(getDist(R.kern_ret, 'norm_dRfCenters'), 'variablenames', {'dRetKern_c2cDist'});   % Find Retinotopy Kernel Similirity

dRetKern_c2cDistFar     = table(isInUpperQuart(dRetKern_c2cDist.dRetKern_c2cDist), 'variablenames', {'dRetKern_c2cDistFar'});   % Find Retinotopy Kernel Similirity        
dRetKern_c2cDistClose   = table(isInLowerQuart(dRetKern_c2cDist.dRetKern_c2cDist), 'variablenames', {'dRetKern_c2cDistClose'});   % Find Retinotopy Kernel Similirity        


dRetKern_corr       = table(getDist(R.kern_ret, 'kerncorr'), 'variablenames', {'dRetKern_corr'});   % Find Retinotopy Kernel Similirity        
dRetKern_corrStrong = table(dRetKern_corr.dRetKern_corr>.55, 'variablenames', {'dRetKern_corrStrong'});   % Find Retinotopy Kernel Similirity        
dRetKern_corrWeak   = table(dRetKern_corr.dRetKern_corr<.15, 'variablenames', {'dRetKern_corrWeak'});   % Find Retinotopy Kernel Similirity        

% dRetKern_ThreshCorr       = table(getDist([R.kern_ret,R.kern_ret_binPixelList], 'threshKernCorr'), 'variablenames', {'dRetKern_ThreshCorr'});   % Find Retinotopy Kernel Similirity        

dRetKern_cos  = table(getDist(R.kern_ret, 'cosangle'), 'variablenames', {'dRetKern_cos'}); % Tuning Kernel  Similirity

dRetKern_Overlap            = table(getDist(R.kern_ret, 'overlap'), 'variablenames', {'dRetKern_Overlap'});   % Find Retinotopy Kernel Similirity
dRetKern_OverlapFiftyPrc    = table(dRetKern_Overlap.dRetKern_Overlap > .5, 'variablenames', {'dRetKern_OverlapFiftyPrc'}); % pairs that have high ret overlap         
dRetKern_OverlapYes         = table(dRetKern_Overlap.dRetKern_Overlap > 0, 'variablenames', {'dRetKern_OverlapYes'}); % pairs that have high ret overlap 
dRetKern_OverlapNo          = table(dRetKern_Overlap.dRetKern_Overlap == 0, 'variablenames', {'dRetKern_OverlapNo'}); % pairs that do not overlap


% Compute Joint Tuning Similarity 
dTunKern_cos  = table(getDist(R.kern_tun, 'cosangle'), 'variablenames', {'dTunKern_cos'}); % Tuning Kernel  Similirity
dTunKern_corr = table(getDist(R.kern_tun, 'kerncorr'), 'variablenames', {'dTunKern_corr'}); % Tuning Kernel  Similirity
dTunKernSmth_cos  = table(getDist(R.kern_tun_smth, 'cosangle'), 'variablenames', {'dTunKernSmth_cos'}); % Tuning Kernel  Similirity
dTunKernSmth_corr  = table(getDist(R.kern_tun_smth, 'kerncorr'), 'variablenames', {'dTunKernSmth_corr'}); % Tuning Kernel  Similirity

% Compute Spat Freq Similarity
dSfResp            = table(getDist(R.sfcurve,'respcorr')', 'variablenames', {'dSfResp'}); % Diff in Log of SF        
dLogSfEst           = table(getDist(log2(R.sfPeak),'absdiff'), 'variablenames', {'dLogSfEst'}); % Diff in Log of SF

dLogSfEst_lessThanHalfOct   = table(dLogSfEst.dLogSfEst <= .5, 'variablenames', {'dLogSfEst_lessThanHalfOct'}); % Pairs with Same SF (diff < 1 octave)        
dLogSfEst_lessThanOneOct    = table(dLogSfEst.dLogSfEst <= 1, 'variablenames', {'dLogSfEst_lessThanOneOct'}); % Pairs with Same SF (diff < 1 octave)
dLogSfEst_lessThanOneAndHalfOct    = table(dLogSfEst.dLogSfEst <= 1.5, 'variablenames', {'dLogSfEst_lessThanOneAndHalfOct'}); % Pairs with Same SF (diff < 1 octave)
dLogSfEst_lessThanTwoOct    = table(dLogSfEst.dLogSfEst <= 2, 'variablenames', {'dLogSfEst_lessThanTwoOct'}); % Pairs with Same SF (diff < 1 octave)
dLogSfEst_lessThanThreeOct  = table(dLogSfEst.dLogSfEst <= 3, 'variablenames', {'dLogSfEst_lessThanThreeOct'}); % Pairs with Same SF (diff < 1 octave)
dLogSfEst_moreThanHalfOct   = table(dLogSfEst.dLogSfEst > .5, 'variablenames', {'dLogSfEst_moreThanHalfOct'}); % Pairs with Same SF (diff < 1 octave)                
dLogSfEst_moreThanOneOct    = table(dLogSfEst.dLogSfEst > 1, 'variablenames', {'dLogSfEst_moreThanOneOct'}); % Pairs with Same SF (diff < 1 octave)
dLogSfEst_moreThanTwoOct    = table(dLogSfEst.dLogSfEst > 2, 'variablenames', {'dLogSfEst_moreThanTwoOct'}); % Pairs with Same SF (diff < 1 octave)
dLogSfEst_moreThanThreeOct  = table(dLogSfEst.dLogSfEst > 3, 'variablenames', {'dLogSfEst_moreThanThreeOct'}); % Pairs with Same SF (diff < 1 octave)

% Compute Ori Pref Similarity
dOriResp_corr       = table(getDist(R.oricurve,'respcorr')', 'variablenames', {'dOriResp_corr'}); % Diff in Log of SF
dOriResp_corrStrong = table(dOriResp_corr.dOriResp_corr>.55, 'variablenames', {'dOriResp_corrStrong'});   % Find Retinotopy Kernel Similirity        
dOriResp_corrWeak   = table(dOriResp_corr.dOriResp_corr<.15, 'variablenames', {'dOriResp_corrWeak'});   % Find Retinotopy Kernel Similirity        

dOriEst             = table(getDist(R.oriest,'distcirc'), 'variablenames', {'dOriEst'}); % 
dOriEst_lessThan10  = table(dOriEst.dOriEst<=10, 'variablenames', {'dOriEst_lessThan10'}); % 
dOriEst_lessThan15  = table(dOriEst.dOriEst<=15, 'variablenames', {'dOriEst_lessThan15'}); % 
dOriEst_lessThan30  = table(dOriEst.dOriEst<=30, 'variablenames', {'dOriEst_lessThan30'}); % 
dOriEst_lessThan45  = table(dOriEst.dOriEst<=45, 'variablenames', {'dOriEst_lessThan45'}); % 
dOriEst_moreThan15  = table(dOriEst.dOriEst>15, 'variablenames', {'dOriEst_moreThan15'}); % 
dOriEst_moreThan30  = table(dOriEst.dOriEst>30, 'variablenames', {'dOriEst_moreThan30'}); % 
dOriEst_moreThan45  = table(dOriEst.dOriEst>45, 'variablenames', {'dOriEst_moreThan45'}); % 
dOriEst_moreThan60  = table(dOriEst.dOriEst>60, 'variablenames', {'dOriEst_moreThan60'}); % 
dOriEst_moreThan75  = table(dOriEst.dOriEst>75, 'variablenames', {'dOriEst_moreThan75'}); % 

% Compute Roi Segment Position Distance
% dSegRoi = computeDistance({[R.segment_xy_tun, R.segment_area_tun]}, 'distance2dNorm');
% dSegRoi = table(dSegRoi, 'variablenames', {'dSegRoi'});

% Get Roi Id Pairing Info
[~, ii, jj]   = getDist(R.sfPeak,'absdiff'); % RetrunsPair Index
pair_roiMasterEntry = table([R.masterEntry(ii), R.masterEntry(jj)], 'VariableNames', {'pair_roiMasterEntry'});
pair_roiUifEntry    = table([R.entryNumber(ii), R.entryNumber(jj)], 'VariableNames', {'pair_roiUifEntry'});
pair_roiDistMatIdx  = table([ii,jj], 'VariableNames', {'pair_roiDistMatIdx'});
numPairs            =  height(pair_roiDistMatIdx);

% Create Table
mouse       = table(repmat(ROIS.mouse(1), numPairs ,1), 'variablenames', {'mouse'}); 
uif         = table(repmat(ROIS.uif(1), numPairs ,1), 'variablenames', {'uif'});
uifEntry    = table((1:numPairs)', 'variablenames', {'uifEntry'});

D = [uifEntry, mouse, pair_roiMasterEntry, uif, pair_roiUifEntry, pair_roiDistMatIdx,...
    dRetKern_cos,...
    dRetKern_c2cDist, dRetKern_c2cDistFar,dRetKern_c2cDistClose...
    dRetKern_corr, dRetKern_corrStrong, dRetKern_corrWeak,...
    dRetKern_Overlap, dRetKern_OverlapFiftyPrc, dRetKern_OverlapYes, dRetKern_OverlapNo,...
    dTunKern_cos, dTunKern_corr,...
    dTunKernSmth_cos, dTunKernSmth_corr,...                
    dSfResp,...
    dLogSfEst, dLogSfEst_lessThanHalfOct, dLogSfEst_lessThanOneOct , dLogSfEst_lessThanOneAndHalfOct, dLogSfEst_lessThanTwoOct ,dLogSfEst_lessThanThreeOct,...
    dLogSfEst_moreThanHalfOct, dLogSfEst_moreThanOneOct, dLogSfEst_moreThanTwoOct, dLogSfEst_moreThanThreeOct,...
    dOriResp_corr, dOriResp_corrStrong, dOriResp_corrWeak,...
    dOriEst,...
    dOriEst_lessThan10, dOriEst_lessThan15, dOriEst_lessThan30, dOriEst_lessThan45,...
    dOriEst_moreThan15, dOriEst_moreThan30, dOriEst_moreThan45, dOriEst_moreThan60, dOriEst_moreThan75,];
%                 dSegRoi];    


end


%% DISTANCE MEASURMENTS 
function [D, ii, jj] = getDist(X, distanceMeasure)

switch distanceMeasure
    
    case 'overlap'
        D = pixOverlap(X);
        
    case 'cosangle'
        D = cosangle(X);
        
    case 'kerncorr'
        stretchKernel   = @(kernelCellArray) cellfun(@(k) k(:), kernelCellArray, 'UniformOutput' , false);
        stretchedRetKernels = stretchKernel(X);
        D = 1- pdist(cat(2,stretchedRetKernels{:})' , 'correlation')';
        
    case 'threshKernCorr'
        D = [1];
        
        

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

%%

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
            sprintf('%d/%d', p, numPairs);
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
        
        sprintf('%d/%d', p, numPairs);
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
            
            sprintf('%d/%d', p, numPairs);
        end 

end



end
