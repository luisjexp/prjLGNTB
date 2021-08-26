classdef lgnanz < tbanz
    % LGN Local Tuning Biases Analysis
    %   This class performs the secondary and final component of the
    %   processing pipeline (after imaging, alignment, segmentation,
    %   merging and signal pulling) for the LGN tuning bias project
    

    %% PREPARE DATA FILES
    properties (Constant)
%         path_raw_data       = lgnanz.get_path_raw_data
        uif_tbl_base        = lgnanz.get_uif_table_base   
        all_available_uif   = lgnanz.uif_tbl_base.uif
    end
    
    methods (Static)
       %% - - - Locate Data Storage Folder
        function path_raw_data = get_path_raw_data() 
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
        
        %% - - - Generate a table of base file names for each imaging sessions 
        function uif_tbl_base = get_uif_table_base()
            % rows correspond to unique imaging field
            % variables...
            %   mouse_id: id of the mouse
            %   fname_base_orisf: the base name of data files that
            %       correspond to imaging data under grating stimulation
            %       used to measure orientation * spatfreq tuning
            %   fname_base_orisf: the base name of data files that
            %       correspond to imaging data under bar/hough stimulation
            %       used to measure retinotopy

%             uif_1 = table(1, {'lgn02'}, {'lgn02_000_004'}, {'lgn02_000_003'});
            uif_1 = table(1, {'lgn02'}, {'lgn02_001_000'}, {'lgn02_001_001'});
            uif_2 = table(2, {'lgn02'}, {'lgn02_002_000'}, {'lgn02_002_002'});
            uif_3 = table(3, {'lgn02'}, {'lgn02_003_000'}, {'lgn02_003_001'});


            uif_4 = table(4, {'lgn03'}, {'lgn03_001_001'}, {'lgn03_001_000'});
            uif_5 = table(5, {'lgn03'}, {'lgn03_001_003'}, {'lgn03_001_002'});
            uif_6 = table(6, {'lgn03'}, {'lgn03_002_001'}, {'lgn03_002_000'});


            uif_7 = table(7, {'lgn04'}, {'lgn04_002_000'}, {'lgn04_002_003'});
            
            uif_8 = table(8, {'lgn05'}, {'lgn05_001_001'}, {'lgn05_001_000'});
            uif_9 = table(9, {'lgn05'}, {'lgn05_008_001'}, {'lgn05_008_000'});
            uif_10 = table(10, {'lgn05'}, {'lgn05_011_000'}, {'lgn05_011_001'});

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
                uif_10]   ;      
            
            uif_tbl_base.Properties.VariableNames = {'uif', 'mouse_id', 'fname_base_ret', 'fname_base_orisf'};        
        end
           
    end
    
    %% GENERATE DATA TABLES
    properties     
        roi_stack      
        roi_stack_info 
        roi_stack_uif_list 
        
        dist_stack
        dist_stack_uif_list
    end
    
    methods 
        %% - - - For each imaging field, generate an roi table. and stack them together. 
        function [roi_stack_new, roi_stack_info_new] = get_roi_tables(obj, uif_list)
            arguments
                obj
                uif_list = obj.all_available_uif % default run through all available imaging fields
            end
            
            % clear the previous data 
            obj.roi_stack           = [];
            obj.roi_stack_info      = [];
            obj.roi_stack_uif_list  = [];
            
            obj.dist_stack          = [];
            obj.ranksum_stack       = [];
            obj.lmfit_stack         = [];
            
            roi_stack_new = [];
            roi_stack_info_new = [];
            for uif = uif_list(:)'
                [roi_tbl_ith, uif_info_ith] = get_one_roi_table(uif) ;              
                roi_stack_new                   = [roi_stack_new; roi_tbl_ith];   
                roi_stack_info_new              = [roi_stack_info_new;  uif_info_ith];
            end
            
            obj.roi_stack           = [table((1:height(roi_stack_new))', 'VariableNames', {'masterEntry'} ), roi_stack_new];
            obj.roi_stack_info      = roi_stack_info_new;
            obj.roi_stack_uif_list  = uif_list;
            
            function [roi_tbl, imaging_field_info] = get_one_roi_table(uif)
                % - - - run
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



                % MATCH rois in this population
                MTC         = sbxmatchfields([ret_fname_base, '_rigid'], [tun_fname_base,'_rigid'],.2);
                
                % Create data table for the population of rois
                S         = struct2table(ret_kernels);
                S         = S(:,{'kern','kurtmax', 'sig', 'xy'});
                S.Properties.VariableNames = {'kern_ret', 'kurt_ret', 'sig_ret', 'xy_ret'};
                S.sig_ret = num2cell(S.sig_ret,2);

                T         =  struct2table(tun_kernels);
                T         = T(:,{'kern_tune', 'sig_tune', 'signal_orisf', 'oriresp_tune', 'oriest_tune', 'sfresp_tune', 'osi_tune', 'sfest_tune', 'xy_soma', 'area_soma'});
                T.Properties.VariableNames = {'kern_tun', 'signif_tun', 'sig_tun', 'oricurve', 'oriest', 'sfcurve', 'osi', 'sfPeak',  'segment_xy_tun', 'segment_area_tun'};


                % Create Table of ROIS that ONLY have Retinotpic Kernels
                % (are missing tuning kernels)
                rois_haveRetData    = array2table(MTC.AnotB', 'variablenames', {'roi_ID_Ret'});
                num_haveRetData     = height(rois_haveRetData);   
                hasRetData          = table(true(num_haveRetData ,1), 'variablenames', {'hasRetData'});

                rois_haveTunData    = array2table(nan(num_haveRetData, 1), 'variablenames', {'roi_ID_Tun'});    
                num_haveTunData     = height(rois_haveTunData);  
                hasTunData          = table(false(num_haveTunData ,1), 'variablenames', {'hasTunData'});

                K_onlyRetData       = [hasRetData, hasTunData, rois_haveRetData, rois_haveTunData,...
                                            [S(MTC.AnotB,:), make_empty_table(T(1,:),num_haveRetData)] ]; 

                % Create Table of ROIS that ONLY have tuning Kernels
                % (are missing retinotopic kernels)                             
                rois_haveTunData    = array2table(MTC.BnotA', 'variablenames', {'roi_ID_Tun'}); 
                num_haveTunData     = height(rois_haveTunData);         
                hasTunData         = table(true(num_haveTunData ,1), 'variablenames', {'hasTunData'});

                rois_haveRetData    = array2table(nan(num_haveTunData, 1), 'variablenames', {'roi_ID_Ret'});    
                num_haveRetData     = height(rois_haveRetData);         
                hasRetData         = table(false(num_haveRetData ,1), 'variablenames', {'hasRetData'});

                K_onlyTunData       = [hasRetData, hasTunData, rois_haveRetData, rois_haveTunData,...
                                            [make_empty_table(S(1,:), num_haveTunData), T(MTC.BnotA,:)] ]; 


                % Create Table of ROIS that have BOTH tuning and retinotopic kernels
                rois_haveRetData    = array2table(MTC.match(:,1), 'variablenames', {'roi_ID_Ret'});
                rois_haveTunData    = array2table(MTC.match(:,2), 'variablenames', {'roi_ID_Tun'});
                num_haveRetData     = height(rois_haveRetData);         
                num_haveTunData     = height(rois_haveTunData);         
                hasRetData          = table(true(num_haveRetData ,1), 'variablenames', {'hasRetData'});
                hasTunData          = table(true(num_haveTunData ,1), 'variablenames', {'hasTunData'});
                K_haveRetAndTuneData       = [hasRetData, hasTunData, rois_haveRetData, rois_haveTunData,...
                                            [[S(MTC.match(:,1),:) T(MTC.match(:,2),:)]]];   

                % Stack tables, assign Mouse Id, Unique Population ID
                roi_tbl         = [K_onlyRetData; K_onlyTunData; K_haveRetAndTuneData];
                mse             = table(repmat({mouse}, height(roi_tbl) ,1), 'variablenames', {'mouse'}); 
                uniqImFieldNum  = table(repmat(uif, height(roi_tbl) ,1), 'variablenames', {'uniqImFieldNum'});
                entryNumber     = table((1:height(roi_tbl))', 'variablenames', {'entryNumber'});


                roi_tbl     = [entryNumber, mse, uniqImFieldNum, roi_tbl] ;
                
                % Create Table of Rois From this imaging field 
                imaging_field_info = table({mouse},  uif,  {fname_base_ret},    {fname_base_orisf},...
                    ret_has_quad_data,  ret_quad,  tun_has_quad_data, tun_quad,...
                    {ret_segment_mask}, {tun_segment_mask}, {ret_align_image},{tun_align_image}, ...
                    'variablenames',{'mouseName', 'uniqImFieldNum', 'unit_exp_hough','unit_exp_orisf',...
                    'hasQuadData_Ret', 'hough_quad',   'hasQuadData_Orisf', 'orisf_quad',...
                    'hough_segMask', 'orisf_segMask', 'hough_alignImage', 'orisf_alignImage'});
                
                                  
                % For making empty table                
                function T_empty = make_empty_table(table2rep, heightOfTable)
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
             
        %% - - - For each imaging field, generate a distance table. and stack together. 
        function dist_stack = get_distance_table(obj,options)
            arguments
                obj
                options.uif_list = obj.roi_stack_uif_list % default run all imaging fields with an roi table 
            end
            dist_stack  = [];
            
            for uif_val = options.uif_list(:)'
                r = obj.roi_stack(obj.roi_stack.uniqImFieldNum == uif_val,:);
                mouse_id = r.mouse{1};
                [distances] = get_one_dist_table(r,mouse_id,uif_val);  
                dist_stack = [dist_stack; distances];        
            end

            
            obj.dist_stack          = dist_stack;
            
            % bad rois from imaging fields are tossed when creating distance table,
            % so tables cannot be created for some imaging fields. so only
            % note those that are left
            obj.dist_stack_uif_list = unique(dist_stack.uniqImFieldNum); 
            
            function D = get_one_dist_table(roi_table_from_one_uif, mouse_id, uif_id)

                isInUpperQuart = @(x) x(:)>prctile(x(:), 75);
                isInLowerQuart = @(x) x(:)<prctile(x(:), 25);
                isInMidFifty  = @(x) x(:)>prctile(x(:), 25) & x(:)<prctile(x(:), 75);
                
                % Lets Filter Out Bad Rois
                roi_requirements = roi_table_from_one_uif.kurt_ret>13 & roi_table_from_one_uif.signif_tun;
                R                = roi_table_from_one_uif(roi_requirements,:);

                % Compute retinotopic Overlap/similarity
                dRetKern_c2cDist        = droi(R.kern_ret, 'norm_dRfCenters');   % Find Retinotopy Kernel Similirity
                dRetKern_corr           = droi(R.kern_ret, 'kerncorr');   % Find Retinotopy Kernel Similirity        
                dRetKern_corrStrong     = dRetKern_corr>.55;   % Find Retinotopy Kernel Similirity        
                dRetKern_corrWeak       = dRetKern_corr<.15;   % Find Retinotopy Kernel Similirity        

                dRetKern_Overlap            = droi(R.kern_ret, 'overlap');   % Find Retinotopy Kernel Similirity
                dRetKern_OverlapFiftyPrc    = dRetKern_Overlap > .5;         % pairs that have high ret overlap         
                dRetKern_OverlapYes         = dRetKern_Overlap > 0;         % pairs that have high ret overlap 
                dRetKern_OverlapNo          = dRetKern_Overlap == 0;        % pairs that do not overlap

                % Compute Tuning Similarity Measures
                dTunKern_cos    = droi(R.kern_tun, 'cosangle'); % Tuning Kernel  Similirity
                dTunKern_corr   = droi(R.kern_tun, 'kerncorr'); % Tuning Kernel  Similirity

                dSfResp            = droi(R.sfcurve,'respcorr'); % Diff in Log of SF        
                dLogSfEst          = droi(log2(R.sfPeak),'absdiff'); % Diff in Log of SF

                dLogSfEst_lessThanHalfOct   = dLogSfEst <= .5; % Pairs with Same SF (diff < 1 octave)        
                dLogSfEst_lessThanOneOct    = dLogSfEst <= 1; % Pairs with Same SF (diff < 1 octave)
                dLogSfEst_lessThanTwoOct    = dLogSfEst <= 2; % Pairs with Same SF (diff < 1 octave)
                dLogSfEst_lessThanThreeOct  = dLogSfEst <= 3; % Pairs with Same SF (diff < 1 octave)
                dLogSfEst_AnyDiff           = dLogSfEst <= max(dLogSfEst); % Pairs with Different SF (diff > 1 octave)
                dLogSfEst_moreThanHalfOct   = dLogSfEst > .5; % Pairs with Same SF (diff < 1 octave)                
                dLogSfEst_moreThanOneOct    = dLogSfEst > 1; % Pairs with Same SF (diff < 1 octave)
                dLogSfEst_moreThanTwoOct    = dLogSfEst > 2; % Pairs with Same SF (diff < 1 octave)
                dLogSfEst_moreThanThreeOct  = dLogSfEst > 3; % Pairs with Same SF (diff < 1 octave)

                dOriResp_corr       = droi(R.oricurve,'respcorr'); % Diff in Log of SF
                dOriResp_corrStrong = dOriResp_corr>.55;   % Find Retinotopy Kernel Similirity        
                dOriResp_corrWeak   = dOriResp_corr<.15;  % Find Retinotopy Kernel Similirity        


                dOriEst             = droi(R.oriest,'distcirc'); % Diff in Log of SF
                dOriEst_Same        = dOriEst<=15; % Diff in Log of SF
                dOriEst_Different   = dOriEst>=45; % Diff in Log of SF                
                
                % Compute Roi Segment Position Distance
                dSegRoi = droi({[R.segment_xy_tun, R.segment_area_tun]}, 'distance2dNorm');

                % GENERATE TABLE
                % Get Roi Pairing Info
                [~, ii, jj]         = droi(R.sfPeak,'absdiff'); % RetrunedsPair Index
                pair_roiDistMatIdx  = [ii,jj];
                pair_roiOrisfId     = [R.roi_ID_Tun(ii), R.roi_ID_Tun(jj)];
                pair_roiHoughId     = [R.roi_ID_Ret(ii), R.roi_ID_Ret(jj)];        
                numPairs            = size(pair_roiDistMatIdx,1);

                % Create Table
                mouse       = repmat({mouse_id}, numPairs ,1);
                uniqImFieldNum         = repmat(uif_id, numPairs ,1);
                
                
                D = table(mouse, uniqImFieldNum, pair_roiDistMatIdx, pair_roiOrisfId, pair_roiHoughId,... dRetKern_c2cDist,...
                                dRetKern_corr, dRetKern_corrStrong, dRetKern_corrWeak,...
                                dRetKern_Overlap, dRetKern_OverlapFiftyPrc, dRetKern_OverlapYes, dRetKern_OverlapNo,...
                                dTunKern_cos, dTunKern_corr,...
                                dSfResp,...
                                dLogSfEst, dLogSfEst_lessThanHalfOct, dLogSfEst_lessThanOneOct , dLogSfEst_lessThanTwoOct ,dLogSfEst_lessThanThreeOct,dLogSfEst_AnyDiff,...
                                           dLogSfEst_moreThanHalfOct, dLogSfEst_moreThanOneOct, dLogSfEst_moreThanTwoOct, dLogSfEst_moreThanThreeOct,...
                                dOriResp_corr, dOriResp_corrStrong, dOriResp_corrWeak,...
                                dOriEst,    dOriEst_Same ,  dOriEst_Different,...
                                dSegRoi);          
                        
                
            end
        end
        
    end
    
    %% ANALYSIS    
    properties
        ranksum_stack
        lmfit_stack
        
    end
    
    methods
        %% - - - RUN A RANK SUM TEST
        function [rnkTest, binCuttoffs, binMedians] = test_and_plot_ranksum(obj, uif, dY_name, dX_name, options)
            
            arguments
                obj
                uif                     
                dY_name                 = 'dTunKern_corr'
                dX_name                 = 'dRetKern_corr'
                options.tail            = 'both'
                options.axis_handle     = nexttile
            end            
            cla(options.axis_handle)                
   
            % -----------
            % INITIALIZE
            % Get the distance data table of all imaging fields
            imfield_idx =  ismember(obj.dist_stack.uif, uif);            
            d           = obj.dist_stack(imfield_idx,:);
              
            % Stop exectution if imaging fields ids do not exist or if imaging fields 
            % do not have enough pairs of rois to analyze
            if ~all(ismember(uif, obj.dist_stack.uif))
                text(options.axis_handle, .5,.5, sprintf('Distance table for one of these imaging fields has not been created'),...
                'units', 'normalized', 'horizontalalignment', 'center', 'fontsize', 8) 
                axis off square
                set(gca,'color','none') 
                return
            end
            

            if height(d) < 2
                text(options.axis_handle, .5,.5, sprintf('One of these imaging fields does not have enough roi pairs to analyze'),...
                'units', 'normalized', 'horizontalalignment', 'center', 'fontsize', 8) 
                axis off square
                set(gca,'color','none')
                return
            end     
            % -----------
            
            % -----------
            % GET THE DATA BINS
            upperBinThresh = 75;
            lowerBinThresh = 25;
            isInUpperBin = @(x) x(:)>=prctile(x(:), upperBinThresh);
            isInLowerBin = @(x) x(:)<=prctile(x(:), lowerBinThresh);
            isInMidFifty  = @(x) x(:)>=prctile(x(:), 25) & x(:)<=prctile(x(:), 75);


            dY_lowerBin   = d{isInLowerBin(d{:,dX_name}), dY_name}; %dY, for observations that fall in the lower bin of dX
            dY_upperBin   = d{isInUpperBin(d{:,dX_name}), dY_name}; %dY, for observations that fall in the upper bin of dX
            num_LowerBin = numel(dY_lowerBin);
            num_UpperBin = numel(dY_upperBin);

            % -----------
             % RUN TEST AND PLOT
            if ~isempty(dY_lowerBin) && ~isempty(dY_upperBin)
                rnkTest             = ranksum(dY_lowerBin, dY_upperBin, 'tail', options.tail);  
                dY_lowerBin_median  = median(dY_lowerBin);
                dY_upperBin_median  = median(dY_upperBin);
                dMedian             = abs(dY_lowerBin_median - dY_upperBin_median);


                binCuttoffs = [prctile(d{:,dX_name}, lowerBinThresh), prctile(d{:,dX_name}, upperBinThresh)];
                binMedians = [dY_lowerBin_median dY_upperBin_median];

                H(1) = histogram(options.axis_handle , dY_lowerBin); axis square; hold on;
                H(2) = histogram(options.axis_handle , dY_upperBin); axis square;   

                box off; axis tight;
                set(gca,'color','none', 'Box', 'off')


                % Fit the histogram        
                nrm = @(x) (x - min(x(:)))/(max(x(:)-min(x(:))));
                x               = linspace(min(xlim), max(xlim), 100);
                pdf_dY_lower    = fitdist(dY_lowerBin,'Kernel', 'width', H(1).BinWidth);
                pdf_dY_lower    = nrm(pdf(pdf_dY_lower,x))*max(H(1).Values);
                pdf_dY_upper    = fitdist(dY_upperBin,'Kernel','width', H(1).BinWidth);
                pdf_dY_upper    = nrm(pdf(pdf_dY_upper,x))*max(H(2).Values);  
                plot(options.axis_handle, x,pdf_dY_lower,'LineWidth',.5, 'Color', 'b')
                plot(options.axis_handle, x,pdf_dY_upper,'LineWidth',.5, 'Color', 'r')

             % provide info   
                xlabel([sprintf('LOWER BIN BELOW %.02f | Med =%.02f| n=%d)\n',binCuttoffs(1), dY_lowerBin_median, num_LowerBin),...
                        sprintf('UPPER BIN ABOVE %.02f | Med =%.02f| n=%d)\n',binCuttoffs(2), dY_upperBin_median, num_UpperBin),...
                        sprintf('dMed=%.02f | P(LOW is *%s* of UP) = %.01d)', dMedian, options.tail, rnkTest)],...
                        'FontSize', 8)
            else
                rnkTest = [];    
                binCuttoffs = [nan nan];
                binMedians = [nan nan];
                text(AXS, .5,.5, sprintf('Not Enough\nData Points'),...
                    'units', 'normalized', 'horizontalalignment', 'center', 'fontsize', 6)        
            end
 
            text(.1, .8, ['fields: ', sprintf('%d ', uif)],...
                'Units', 'normalized', 'BackgroundColor', [0,0,0, .1])  

            figure(gcf)
        end
            
        %% - - - Regress receptive field overlap on tuning similarity
        function P = fit_plot_lm(obj, uif_list, dY_name, dX_name, options)
             arguments
                obj
                uif_list
                dY_name                     = 'dTunKern_corr'
                dX_name                     = 'dRetKern_corr'
                options.binary_var_names    = ''
                options.axis_handle         = nexttile
             end                
             
             
            % Get the distance data table of all imaging fields specified
            uifs_available = obj.dist_stack.uniqImFieldNum;
            imfield_idx =  ismember(uifs_available, uif_list);            
            d           = obj.dist_stack(imfield_idx,:);
            dY          = d{:, dY_name};
            dX          = d{:, dX_name};            

            % Stop exectution if imaging field ids do not exist 
            cla(options.axis_handle)
            if ~all(ismember(uif_list, uifs_available))
                text(options.axis_handle, .5,.5, sprintf('Distance table for one of these imaging fields has not been created'),...
                'units', 'normalized', 'horizontalalignment', 'center', 'fontsize', 8) 
                axis off square
                set(gca,'color','none') 
                return
            end

            % Stop exectution  if imaging fields do not have enough roi pairs to analyze
            if height(d) < 2
                text(options.axis_handle, .5,.5, sprintf('One of these imaging fields does not have enough roi pairs to analyze'),...
                'units', 'normalized', 'horizontalalignment', 'center', 'fontsize', 8) 
                axis off square
                set(gca,'color','none')
                return
            end     

            % Fit Model 
            lm   = fitlm(dX, dY, 'ResponseVar', dY_name, 'PredictorVars', dX_name );

            beta       = lm.Coefficients.Estimate(2);
            p_val       = lm.Coefficients.pValue(2);
            r           = sign(beta)*sqrt(lm.Rsquared.Ordinary);
            n_pairs     = lm.NumObservations;

           % Plot Model Fit
            axis(options.axis_handle);
            P               = lm.plotAdded; 
            P(1).Marker     = '.';
            P(1).MarkerSize = 8;        
            
            title(sprintf('r=%.02f | b=%.02f | p=%.01d\nn=%d', r, beta, p_val, n_pairs),...
                'fontweight', 'normal', 'fontsize', 8)%, 'BackgroundColor', [0,0,0, .1])
            
            text(.1, .8, ['fields: ', sprintf('%d ', uif_list)],...
                'Units', 'normalized')%, 'BackgroundColor', [0,0,0, .1])

            hold on;
            legend off;  
            axis square;   
            set(gca,'color','none', 'Box', 'off')

            
            % Show the condition the pair of boutons
            if ~isempty(options.binary_var_names)
                dC0 =  logical(d{:,options.binary_var_names{1}});
                dC1 =  logical(d{:,options.binary_var_names{2}});

                dX_given_dC0 = dX(dC0,:);
                dY_given_dC0 = dY(dC0,:);

                dX_given_dC1 = dX(dC1,:);
                dY_given_dC1 = dY(dC1,:);                


                plot(options.axis_handle, dX_given_dC0, dY_given_dC0, 'bo')
                plot(options.axis_handle, dX_given_dC1, dY_given_dC1, 'go')
                legend(options.binary_var_names)
            end

            drawnow;
            figure(gcf)
        end
        
    end
    
    methods (Static)
        %% - - - Plot Histogram 
        function plot_univar_hist(data_table, variable_name, options)
            arguments
                data_table           
                variable_name           = 'dTunKern_corr'
                options.axis_handle     = nexttile
            end
            
            X   = data_table{:, variable_name};
            uifs_in_data_table = unique(data_table.uniqImFieldNum);
                    
            cla(options.axis_handle) 

            % Stop exectution if imaging fields do not have enough rois (or pairs of rois) to analyze            
            if height(X) < 2
                text(options.axis_handle, .5,.5, sprintf('Imaging field does not have rois'),...
                'units', 'normalized', 'horizontalalignment', 'center', 'fontsize', 8) 
                axis off square
                set(gca,'color','none')
                return
            end        


            
            number_of_units     = numel(X);
            mean_of_variable    = mean(X);

            histogram(options.axis_handle, X); 
            xlabel(variable_name, 'Interpreter', 'none')
            legend off;  
            axis square;   
            line([mean_of_variable mean_of_variable], ylim)                
            set(gca,'color','none', 'Box', 'off')
            title(sprintf('N: %d', number_of_units), 'fontsize', 10);
            text(.1,.9, ['imaging fields: ' sprintf('%d ', uifs_in_data_table)], 'Units', 'normalized')
            

        end        
        
    end



end

