classdef lgnanz < handle
    % LGN Local Tuning Biases Analysis
    %   This class performs the secondary and final component of the
    %   processing pipeline (after imaging, alignment, segmentation,
    %   merging and signal pulling) for the LGN tuning bias project
    
    
    properties (Constant)
        path_raw_data       = lgnanz.get_path_raw_data
        uif_tbl_base        = lgnanz.get_uif_table_base        
    end
    
    properties 
        roi_stack      
        roi_stack_info  
        dist_stack
    end
    
    methods (Static)
        
        function uif_tbl_base = get_uif_table_base()
            %% generate a table of base file names for each imaging sessions 
            % rows correspond to unique imaging field
            % variables...
            %   mouse_id: id of the mouse
            %   fname_base_orisf: the base name of data files that
            %       correspond to imaging data under grating stimulation
            %       used to measure orientation * spatfreq tuning
            %   fname_base_orisf: the base name of data files that
            %       correspond to imaging data under bar/hough stimulation
            %       used to measure retinotopy
            
            uif_1 = table(1, {'lgn02'}, {'lgn02_000_004'}, {'lgn02_000_003'});
            uif_2 = table(2, {'lgn02'}, {'lgn02_001_000'}, {'lgn02_001_001'});
            uif_3 = table(3, {'lgn02'}, {'lgn02_002_000'}, {'lgn02_002_002'});
            uif_4 = table(4, {'lgn02'}, {'lgn02_003_000'}, {'lgn02_003_001'});


            uif_5 = table(5, {'lgn03'}, {'lgn03_001_001'}, {'lgn03_001_000'});
            uif_6 = table(6, {'lgn03'}, {'lgn03_001_003'}, {'lgn03_001_002'});
            uif_7 = table(7, {'lgn03'}, {'lgn03_002_001'}, {'lgn03_002_000'});


            uif_8 = table(8, {'lgn04'}, {'lgn04_002_000'}, {'lgn04_002_003'});
            
            uif_9 = table(9, {'lgn05'}, {'lgn05_001_001'}, {'lgn05_001_000'});
            uif_10 = table(10, {'lgn05'}, {'lgn05_008_001'}, {'lgn05_008_000'});
            uif_11 = table(11, {'lgn05'}, {'lgn05_011_000'}, {'lgn05_011_001'});


            uif_tbl_base = [...
                uif_1;...
                uif_2;...
                uif_3;...
                uif_4;...
                uif_5;...
                uif_6;...
                uif_7;...
                uif_8;...
                uif_9;...
                uif_10;...
                uif_11];          
            
            uif_tbl_base.Properties.VariableNames = {'uif', 'mouse_id', 'fname_base_ret', 'fname_base_orisf'};
            
            
        end
            
        
        function path_raw_data = get_path_raw_data() 
            %% Locate Data Storage Folder
            % NOTE: replaces ''lojGetDataPath'' function in retired folder
            
            if ismac && strcmp(getenv('LOGNAME'), 'luis')
                path_raw_data = '/Users/luis/Box/boxPHD/Toolbox/2pdata/';
            elseif ispc && strcmp(getenv('username'), 'Luis')
                path_raw_data = 'C:\Users\Luis\Box Sync\boxEXPER\Toolbox\2pdata\';
            elseif ispc && strcmp(getenv('username'), 'dario')
                path_raw_data = 'C:\2pdata\';
            else
                error('unknown computer')
            end
            
        end 
        
    end
    
    methods 
        function [roi_stack, roi_stack_info] = get_stacked_roi_tables(obj, uif_list)
            %% Concatenate ROI tables from multiple imaging fields            
            arguments
                obj
                uif_list = [6 7]
            end
            
            roi_stack       = [];
            roi_stack_info        = [];
            for uif = uif_list
                [roi_tbl_ith, uif_info_ith] = get_single_imfield_roi_table(uif) ;              
                roi_stack                   = [roi_stack; roi_tbl_ith];   
                roi_stack_info              = [roi_stack_info;  uif_info_ith];
            end
            
            obj.roi_stack       = [table((1:height(roi_stack))', 'VariableNames', {'masterEntry'} ), roi_stack];
            obj.roi_stack_info  = roi_stack_info;
            
            function [roi_tbl, imaging_field_info] = get_single_imfield_roi_table(uif)
                uif_idx             = find(lgnanz.uif_tbl_base.uif == uif);
                mouse               = lgnanz.uif_tbl_base.mouse_id{uif_idx};
                fname_base_ret      = lgnanz.uif_tbl_base.fname_base_ret{uif_idx};
                fname_base_orisf    = lgnanz.uif_tbl_base.fname_base_orisf{uif_idx}   ;         

                % Build data File names created from retinotopic imaging sessions    
                ret_fname_base    = [lgnanz.path_raw_data, mouse,'/' fname_base_ret];   
                ret_fname_kernels = [ret_fname_base , '_rigid.houghkernstats'];
                ret_fname_quad    = [ret_fname_base , '_quadrature.mat'];
                ret_fname_segment = [ret_fname_base , '_rigid.segment'];
                ret_fname_align   = [ret_fname_base,  '.align'];
                ret_fname_spikes  = [ret_fname_base , '_rigid.signals'];

                % Build File Names from orisf imaging sessions                
                tun_fname_base    = [lgnanz.path_raw_data, mouse,'/', fname_base_orisf];
                tun_fname_kernels = [tun_fname_base ,'_rigid.orisf'];
                tun_fname_quad    = [tun_fname_base , '_quadrature.mat'];
                tun_fname_segment = [tun_fname_base , '_rigid.segment'];
                tun_fname_align   = [tun_fname_base,  '.align'];
                tun_fname_spikes  = [tun_fname_base , '_rigid.signals'];


                % Get Kernels from this population
                ret_kernels     = load(ret_fname_kernels, '-mat');   
                ret_kernels     = ret_kernels.S;
                tun_kernels     = load(tun_fname_kernels, '-mat');   
                tun_kernels     = tun_kernels.stat;

                % Get Quad Data From This Population
                if exist(ret_fname_quad, 'file')
                    ret_quad = load(ret_fname_quad, '-mat');   
                    ret_quad = {double(ret_quad.quad_data)}; 
                    ret_has_quad_data = true;
                else
                    ret_quad = {nan};
                    ret_has_quad_data = false;               
                end

                % Get Quad Data From This Population
                if exist(tun_fname_quad, 'file')
                    tun_quad = load(tun_fname_quad, '-mat');   
                    tun_quad = {double(tun_quad.quad_data)};   
                    tun_has_quad_data = true;                       
                else
                    tun_quad = {nan};
                    tun_has_quad_data = false;                       
                end

                % Load Segmented Images from this imaging session
                ret_segment_mask = load(ret_fname_segment, '-mat', 'mask'); 
                ret_segment_mask = ret_segment_mask.mask;
                tun_segment_mask = load(tun_fname_segment, '-mat', 'mask'); 
                tun_segment_mask = tun_segment_mask.mask;

                % Load Alignment Images from this population
                if exist(ret_fname_align, 'file')
                    ret_align_image = load(ret_fname_align, '-mat', 'm'); 
                    ret_align_image = ret_align_image.m;
                else
                    ret_align_image = {nan};
                end

                if exist(tun_fname_align, 'file')
                    tun_align_image = load(tun_fname_align, '-mat', 'm'); 
                    tun_align_image = tun_align_image.m;    
                else
                    tun_align_image = {nan};
                end 



                % Match rois in this population
                MTC         = sbxmatchfields([ret_fname_base, '_rigid'], [tun_fname_base,'_rigid'],.2);

                % Create data table for the population of rois
                S         = struct2table(ret_kernels);
                S         = S(:,{'kern','kurtmax', 'sig', 'xy'});
                S.Properties.VariableNames = {'kern_ret', 'kurt_ret', 'sig_ret', 'xy_ret'};
                S.sig_ret = num2cell(S.sig_ret,2);

                T         =  struct2table(tun_kernels);
                T         = T(:,{'kern_tune', 'sig_tune', 'signal_orisf', 'oriresp_tune', 'oriest_tune', 'sfresp_tune', 'osi_tune', 'sfest_tune', 'xy_soma', 'area_soma'});
                T.Properties.VariableNames = {'kern_tun', 'signif_tun', 'sig_tun', 'oricurve', 'oriest', 'sfcurve', 'osi', 'sfPeak',  'segment_xy_tun', 'segment_area_tun'};


                % Create Table of ROIS that ONLYL have Retinotpic Kernels
                % (are missing tuning kernels)
                rois_haveRetData    = array2table(MTC.AnotB', 'variablenames', {'roi_ID_Ret'});
                num_haveRetData     = height(rois_haveRetData);   
                hasRetData          = table(true(num_haveRetData ,1), 'variablenames', {'hasRetData'});

                rois_haveTunData    = array2table(nan(num_haveRetData, 1), 'variablenames', {'roi_ID_Tun'});    
                num_haveTunData     = height(rois_haveTunData);  
                hasTunData          = table(false(num_haveTunData ,1), 'variablenames', {'hasTunData'});

                K_onlyRetData       = [hasRetData, hasTunData, rois_haveRetData, rois_haveTunData,...
                                            [S(MTC.AnotB,:), makeEmptyTable(T(1,:),num_haveRetData)] ]; 

                % Create Table of ROIS that ONLY have tuning Kernels
                % (are missing retinotopic kernels)                             
                rois_haveTunData    = array2table(MTC.BnotA', 'variablenames', {'roi_ID_Tun'}); 
                num_haveTunData     = height(rois_haveTunData);         
                hasTunData         = table(true(num_haveTunData ,1), 'variablenames', {'hasTunData'});

                rois_haveRetData    = array2table(nan(num_haveTunData, 1), 'variablenames', {'roi_ID_Ret'});    
                num_haveRetData     = height(rois_haveRetData);         
                hasRetData         = table(false(num_haveRetData ,1), 'variablenames', {'hasRetData'});

                K_onlyTunData       = [hasRetData, hasTunData, rois_haveRetData, rois_haveTunData,...
                                            [makeEmptyTable(S(1,:), num_haveTunData), T(MTC.BnotA,:)] ]; 


                % Create Table of ROIS that have BOTH tuning and retinotopic kernels
                rois_haveRetData    = array2table(MTC.match(:,1), 'variablenames', {'roi_ID_Ret'});
                rois_haveTunData    = array2table(MTC.match(:,2), 'variablenames', {'roi_ID_Tun'});
                num_haveRetData     = height(rois_haveRetData);         
                num_haveTunData     = height(rois_haveTunData);         
                hasRetData         = table(true(num_haveRetData ,1), 'variablenames', {'hasRetData'});
                hasTunData         = table(true(num_haveTunData ,1), 'variablenames', {'hasTunData'});
                K_haveRetAndTuneData       = [hasRetData, hasTunData, rois_haveRetData, rois_haveTunData,...
                                            [[S(MTC.match(:,1),:) T(MTC.match(:,2),:)]]];   

                % Stack tables, assign Mouse Id, Unique Population ID
                roi_tbl         = [K_onlyRetData; K_onlyTunData; K_haveRetAndTuneData];
                mse             = table(repmat({mouse}, height(roi_tbl) ,1), 'variablenames', {'mouse'}); 
                uniqImFieldNum  = table(repmat(uif, height(roi_tbl) ,1), 'variablenames', {'uniqImFieldNum'});
                entryNumber     = table((1:height(roi_tbl))', 'variablenames', {'entryNumber'});


                roi_tbl     = [entryNumber, mse, uniqImFieldNum, roi_tbl] ;

            %    Create Imaging Field
                imaging_field_info = table({mouse},  uif,  {fname_base_ret},    {fname_base_orisf},...
                                        ret_has_quad_data,  ret_quad,  tun_has_quad_data, tun_quad,...
                                        {ret_segment_mask}, {tun_segment_mask}, {ret_align_image},{tun_align_image}, ...
                    'variablenames',{'mouseName', 'uniqImFieldNum', 'unit_exp_hough','unit_exp_orisf',...
                                    'hasQuadData_Ret', 'hough_quad',   'hasQuadData_Orisf', 'orisf_quad',...
                                    'hough_segMask', 'orisf_segMask', 'hough_alignImage', 'orisf_alignImage'});


                % For making empty table                
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
            end 
            
         end       
 
        function dist_stack = get_stacked_distance_table(obj)
            %% Compute Distance measures for rois within each imaging field                      

            dist_stack = [];
            uif_list = unique(obj.roi_stack.uniqImFieldNum)';
            
            for uif_id = uif_list
                r = obj.roi_stack(obj.roi_stack.uniqImFieldNum == uif_id,:);
                [distances, dist_selection] = get_single_imfield_dist_table(r);  
                dist_stack = [dist_stack; distances];        
            end

            close all
            obj.dist_stack = [table((1:height(dist_stack))', 'VariableNames', {'masterEntry'} ), dist_stack];   
            
            
            function [D, roiConstraints] = get_single_imfield_dist_table(roi_pop)
                %% create distance table for a imaging field of rois

                isInUpperQuart = @(x) x(:)>prctile(x(:), 75);
                isInLowerQuart = @(x) x(:)<prctile(x(:), 25);
                isInMidFifty  = @(x) x(:)>prctile(x(:), 25) & x(:)<prctile(x(:), 75);

                % SELECT ROIS FOR DISTANCE TABLE 
                % Must have tuning and ret data, and tuning and ret kernels must be significant

                roiConstraints = roi_pop.hasRetData &...
                    roi_pop.hasTunData &...
                    roi_pop.kurt_ret>2 &...
                    roi_pop.signif_tun;

                R = roi_pop(roiConstraints,:);

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
                mouse       = table(repmat(roi_pop.mouse(1), numPairs ,1), 'variablenames', {'mouse'}); 
                uif         = table(repmat(roi_pop.uniqImFieldNum(1), numPairs ,1), 'variablenames', {'uif'});
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
                            
            function [D, ii, jj] = computeDistance(X, distanceMeasure)
            %% Compute distances between properties of rois
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

            function [d, ii, jj] = cosangle(KcellArray)
                %% Cosinge Angle
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

            end

        end
        



        
    end



end
