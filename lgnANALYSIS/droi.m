function [D, ii, jj] = droi(X, distanceMeasure)
%UNTITLED Compute distances between properties of rois
%   Detailed explanation goes here
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
                D = D'; % transpose into to column
            case 'absdiff'
                D  = abs(takeDiff(X));
                [~, ii, jj] = takeDiff(X); 
                
            case 'distance2dNorm' 
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
%                     sprintf('%d/%d', p, numPairs)
                end   

        end

        % Kernel Pixel Overlap
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

%                 sprintf('%d/%d', p, numPairs)
            end

        end

        % Differences
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

        % Normalized Distance Between RF Centers
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

%                     sprintf('%d/%d', p, numPairs)
                end 

        end
end

