function roiPropHistogram(ROISCAT, property, varargin)

ImFldID_list = unique(ROISCAT.uif);
numImFlds = numel(ImFldID_list);
plotPropValPairs = varargin;
figure(1);  
AxProps = properties(axes);
HistProps = properties(histogram);
figure(1);
   

for j = 1:numImFlds+1

    if j <= numImFlds
        imFldID   = ImFldID_list(j);        
        X           = ROISCAT(ROISCAT.uif == imFldID,:) ;        
        mouseName   = X.mouse{1};
        titleString = sprintf('%d (%s)',imFldID, mouseName )   ;     
        A           = subplot(2, numImFlds, j); 
        
    elseif j >numImFlds
        X           = ROISCAT;        
        A           = subplot(2, 1, 2); 
        if iscell(property)
            p = property{:};
        else
            p = property;
        end
          
        titleString = sprintf(sprintf('%s \nALL FIELDS', p))       ; 
    end
    
    if ~iscell(property)
        roiProp = X(:, property);   
        H       = histogram(A, roiProp.Variables); 
    elseif iscell(property)
        roiPropX = X(:, property{1});  
        roiPropY = X(:, property{2});           
        H       = histogram2(A, roiPropX.Variables,roiPropY.Variables);   
        set(gca,'YDir', 'reverse')
    end
        set(gca,'color','none')    

    
    for k = 1:2:numel(plotPropValPairs)
        plotProp = plotPropValPairs{k};
        plotPropVal =  plotPropValPairs{k+1};
       
        if contains(plotProp, AxProps)
            set(A, plotProp, plotPropVal)  
        elseif contains(plotProp, HistProps)
            set(H, plotProp, plotPropVal)
        else
            error('Unknown hist or axes property')
        end
        
    end    
    
    title(titleString, 'BackgroundColor', 'w');
    axis square;
    box off

end

figure(gcf)

end