classdef v1anz < handle
    %V1ANZ Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Constant)
        path_raw_data = '/Users/luis/Box/prjLGNTB/v1DATA/'
        fname_distance_stack = [v1anz.path_raw_data, 'database_tuning_on_off.mat'];
        
    end
    
    
    %% GET TABLES
    properties
        roi_stack 
        roi_stack_uif_list 
        dist_stack 
        dist_stack_uif_list
        dist_stack_w_both_subreg
        dist_stack_w_both_subreg_uif_list
        
        
    end 
    
    methods 
        %% - - - Load Roi/single neuron table 
        function obj = get_roi_table(obj)
            roi_stack               = load(obj.fname_distance_stack);
            obj.roi_stack           = roi_stack.X;
            obj.roi_stack  = add_unique_imaging_fields(obj.roi_stack);
            obj.roi_stack_uif_list = unique(obj.roi_stack.uniqueImField);
            head(obj.roi_stack)
            tail(obj.roi_stack)
            
            % Reset Distance Table
            obj.dist_stack = [];
            obj.dist_stack_uif_list = [];
            
            function roi_stack = add_unique_imaging_fields(roi_stack)
                A = categorical(categorical(roi_stack.mouseID));
                B = categorical(categorical(roi_stack.imageDepth));
                C = A.*B;
                uniqueImField = grp2idx(string(C));
                roi_stack = [table(uniqueImField), roi_stack];
            end
        
        end
          
  
        %% - - - Generate a distance table from each imaging field, and stack. 
        %% - - - NOTE: double check if pdist is actially performing proper correlation
        function dist_stack = get_distance_table(obj,roi_stack_uif_list)
            arguments
                obj
                roi_stack_uif_list = obj.roi_stack_uif_list % default run all imaging fields with an roi table 
            end
            dist_stack  = [];
            
            for uif_val = roi_stack_uif_list(:)'
                r           = obj.roi_stack(obj.roi_stack.uniqueImField == uif_val,:);
                mouse_id    = r.mouseID{1}; % lgn names different
                [distances] = get_one_dist_table(r,mouse_id,uif_val);  
                dist_stack = [dist_stack; distances];        
            end

            obj.dist_stack          = dist_stack;
            
            % bad rois from imaging fields may* be tossed when creating distance table,
            % so tables cannot be created for some imaging fields. therefore keep track of those that
            % are left
            
            obj.dist_stack_uif_list = unique(dist_stack.uniqImFieldNum);
            
            
            function D = get_one_dist_table(roi_table_from_one_uif, mouse_id, uif_id)
                R = roi_table_from_one_uif;
                % Compute retinotopic Overlap/similarity
                dON_corr    = droi(R.onSubfield, 'kerncorr');   
                dOFF_corr   = droi(R.offSubfield, 'kerncorr'); 

                % Compute Tuning Similarity Measures
                dTunKern_corr   = droi(R.tuningKernel, 'kerncorr'); % Tuning Kernel  Similirity
                [~, ii, jj] = droi(R.tuningKernel, 'cosangle' );% Tuning Kernel  cosign angle

                % Generate table...
                % Get Roi Pairing Info
                pair_roiDistMatIdx  = [ii,jj];   
                numPairs            = size(pair_roiDistMatIdx,1);
                
                % Create Table
                mouse       = repmat({mouse_id}, numPairs ,1);
                uniqImFieldNum         = repmat(uif_id, numPairs ,1);
                
                D = table(mouse, uniqImFieldNum, pair_roiDistMatIdx,... 
                    dON_corr, dOFF_corr,...
                    dTunKern_corr);  
                
                %% Create a distance table but only for pairs with both
                % On AND off subregions
                
                function d = compute_dist(R)

                    
                   
                    
                end

                
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
                dY_name                 = ''
                dX_name                 = ''
                options.tail            = 'both'
                options.axis_handle     = nexttile
            end            
            cla(options.axis_handle)                
   
            % -----------
            % INITIALIZE
            % Get the distance data table of all imaging fields
            imfield_idx =  ismember(obj.dist_stack.uniqImFieldNum, uif);            
            d           = obj.dist_stack(imfield_idx,:);
              
            % Stop exectution if imaging fields ids do not exist or if imaging fields 
            % do not have enough pairs of rois to analyze
            if ~all(ismember(uif, obj.dist_stack.uniqImFieldNum))
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
                dY_name                    
                dX_name                     = 'dTunKern_corr'
                options.binary_var_names    = ''
                options.axis_handle         = nexttile
             end                
             
             
            % Get the distance data table of all imaging fields specified
            uifs_available = unique(obj.dist_stack.uniqImFieldNum);
            imfield_idx =  ismember(uif_list, uifs_available);            
            d           = obj.dist_stack(imfield_idx,:);
            dY          = d{:, dY_name};
            dX          = d{:, dX_name};            

            % Stop exectution if imaging field ids do not exist 
            cla(options.axis_handle)
            if ~all(ismember(uifs_available, uif_list))
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
                'fontweight', 'normal', 'fontsize', 8, 'BackgroundColor', [0,0,0, .1])
            
            text(.1, .8, ['fields: ', sprintf('%d ', uif_list)],...
                'Units', 'normalized',  'BackgroundColor', [0,0,0, .1])

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
    
    
end
        
    


