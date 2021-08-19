function M = loj_dX_vs_dY_by_dC(DISTCAT, biVarFormula, groupVarList, tailTests, varargin)
% Models of dY ~ dX, for each group in groupVarList
% One For Each Imaging Field, ignoring Spat Freq Sim
% One For Each Imaging Field And For Each Spatial Freq Similarity Level

%% INITIALIZE
% DESCRIPTIVE TEXTINITIALZE 

plotPropValPairs = varargin{1};


close all; 
FIG1 = figure(1);
FIG2 = figure(2);

numelUif    = max(DISTCAT.uif);% number of unique imaging fields,  last element is the total number 
numGroups   = numel(groupVarList);
M       = []; % Model Table of imaging fields


%% RUN

for i = 1:numGroups + 1

    for uif = 1:numelUif
        d = DISTCAT(DISTCAT.uif == uif, :);
        if i <= numGroups
            dC  = groupVarList{i};
            d =  d(d{:, dC} == 1,:);
        elseif i == numGroups + 1
            dC = 'ALLPAIRS';
        end
        
        uif_strName = num2str(uif)';
        uif_strName = uif_strName(~isspace(uif_strName)); 
        mouse       = unique(d.mouse); 
        if numel(mouse) >1
            error('you got more than one mouse goin on here som up')
        elseif isempty(mouse) 
            mouse = 'no data';
        else
            mouse = mouse{1};
        end
        
        k = numelUif*(i-1)+uif;
        figure(FIG1);
        AXS = subplot(numGroups+1,  numelUif, k);     
        model = fitAndPlotModel(d, biVarFormula, groupVarList, AXS, plotPropValPairs); 

        
        if i == 1
            text(.5, 1.75, sprintf('Field%s(%s)', uif_strName, mouse),...
                'units', 'normalized', 'HorizontalAlignment', 'center',...
                'fontweight', 'bold', 'fontsize', 8)
        end

        if uif == 1
            text(-1.25, .5,sprintf('%s', dC), 'units', 'normalized',...
                'fontweight', 'bold', 'fontsize', 8)
        else      
            ylabel('');
        end           

        M = [M; table({dC}, {uif_strName}, {uif},  {d}, height(d), {model},...
            'variablenames', {'group', 'imFieldName', 'imageFieldId', 'data', 'numPairs',   'model'})];
        drawnow;        
    end



        

    
end

for i = 1:numGroups + 1
        
    uif      = unique(DISTCAT.uif);
    uif_strName = num2str(uif)';
    uif_strName = uif_strName(~isspace(uif_strName)); 
    mouse           = 'ALLMICE';    

    d = DISTCAT;
    if i <= numGroups
        dC  = groupVarList{i};
        d =  d(d{:, dC} == 1,:);
    elseif i == numGroups + 1
        dC = 'ALLPAIRS';
    end
    
    figure(FIG2);
    AXS = subplot(2,numGroups+1,  i);      hold on;
    model = fitAndPlotModel(d, biVarFormula, groupVarList, AXS, plotPropValPairs);
    
    AXS_rnk = subplot(2,numGroups+1,  numGroups+1+i);    
    [rnkTest, binCuttoffs, binMedians] = runAndPlotRankTest(d,biVarFormula, AXS_rnk,  tailTests{i}, plotPropValPairs);   

    % DRAW MEDIAN AND BINZZZ
        lwdth = 2;
        lstyl = ':';
        line(AXS, binCuttoffs(1)*[1 1], get(AXS, 'ylim'), 'color', 'b', 'linewidth', lwdth, 'linestyle', lstyl);
        line(AXS, binCuttoffs(2)*[1 1], get(AXS, 'ylim'), 'color', 'r','linewidth', lwdth, 'linestyle', lstyl);         
        line(AXS, [AXS.XLim(1) binCuttoffs(1)], binMedians(1)*[1 1], 'color', 'b','linewidth', lwdth, 'linestyle', lstyl); 
        line(AXS, [binCuttoffs(2) AXS.XLim(2)], binMedians(2)*[1 1], 'color', 'r','linewidth', lwdth, 'linestyle', lstyl);   
    
    % 
    text(AXS, .5, 1.35, sprintf('%s\nField%s(%s)', dC, uif_strName, mouse),...
        'units', 'normalized', 'HorizontalAlignment', 'center', 'fontweight', 'bold', 'fontsize', 8)    
    
    
    M = [M; table({dC}, {uif_strName}, {uif},  {d}, height(d), {model},...
        'variablenames', {'group', 'imFieldName', 'imageFieldId', 'data', 'numPairs',   'model'})];
    drawnow;        

    
end   
close(F3)
 

end

%% FIT AND PLOT LINEAR MODEL
function model = fitAndPlotModel(d, biVarFormula, groupVarList, AXS, plotPropValPairs)

if height(d) > 2
    
        if contains(biVarFormula, '~') % run biavariate model
                model = fitlm(d, biVarFormula) ;
                P = plotLM(model, AXS, groupVarList, plotPropValPairs);                
        else % Plot Univariate Histogram 
            plotUnivarHist(d, biVarFormula, AXS, plotPropValPairs);
        end
        set(gca,'color','none')
        
else    
    model = [];
    text(AXS, .5,.5, sprintf('Error Not Enough\nData Points'),...
    'units', 'normalized', 'horizontalalignment', 'center', 'fontsize', 8) 
    axis square
    set(gca,'color','none')
end        
      
        

 %% PLOT BIVARIATE MODEL
function P = plotLM(model, AXS, groupVarList, plotPropValPairs)

        pval    = model.coefTest; 
        bCoeff  = model.Coefficients.Estimate(2);
        r       = sign(bCoeff)*sqrt(model.Rsquared.Ordinary);
        n       = model.NumObservations;
       % Plot
        axis(AXS);
        P       = model.plotAdded; hold on;
        title(sprintf('r=%.02f|b=%.02f\np=%.01d\n(n=%d)', r, bCoeff, pval, n), 'fontweight', 'normal', 'fontsize', 8);
        
        P(1).Marker    = '.';
        P(1).MarkerSize= 8;
        
        
        AxProps = properties(AXS);       
        LineProps = properties(P(1));       
    for k = 1:2:numel(plotPropValPairs)
        plotProp    = plotPropValPairs{k};
        plotPropVal =  plotPropValPairs{k+1};
       
        if contains(plotProp, AxProps)
            set(AXS, plotProp, plotPropVal)  
        elseif contains(plotProp, LineProps)
            set(P(1), plotProp, plotPropVal)
        end
        
    end 
    
        legend off;  
        axis square;   
        box off;          
        set(gca,'color','none')
        

        dC1 = model.Variables{:, groupVarList{1}};
        dC2 = model.Variables{:, groupVarList{2}};
       
        if all(dC1)
                P(1).MarkerEdgeColor = 'b';                                
        elseif all(dC2)
                P(1).MarkerEdgeColor = 'g';      
        else
            dY = model.ResponseName;
            dX = model.PredictorNames{1};
            
            delete(AXS.Children(end))
            plot(AXS, model.Variables{logical(dC1),dX} , model.Variables{logical(dC1),dY}, 'b.')
            plot(AXS, model.Variables{logical(dC2), dX} , model.Variables{logical(dC2),dY}, 'g.')
                          
        drawnow;
        end

end

%% PLOT UNIVARIATE MODEL
function plotUnivarHist(d, formula, AXS, plotPropValPairs)
    x = d{:, formula};
    
    n = numel(x);
    mu = mean(x);

    H = histogram(AXS, x); 
    xlabel(formula)
 
    legend off;  
    axis square;   
    box off;
    set(gca,'color','none')
        
    title(sprintf('n=%d', n), 'fontweight', 'normal', 'fontsize', 8);

        AxProps     = properties(AXS);       
        histProps   = properties(H);       
    for k = 1:2:numel(plotPropValPairs)
        pltProp     = plotPropValPairs{k};
        pltPropVal  =  plotPropValPairs{k+1};

        if contains(pltProp, AxProps)
            set(AXS, pltProp, pltPropVal)  
        elseif contains(pltProp, histProps)
            set(H, pltProp, pltPropVal)
        else
            error('no plot prop named %s',pltProp);
        end

    end 
    
        
    
end


end

%% RANK TEST FOR dY FOR UPPER AND LOWER BINS of dX

function [rnkTest, binCuttoffs, binMedians] = runAndPlotRankTest(d, formula, AXS, tail, plotPropValPairs)
upperBinThresh = 75;
lowerBinThresh = 25;
isInUpperBin = @(x) x(:)>=prctile(x(:), upperBinThresh);
isInLowerBin = @(x) x(:)<=prctile(x(:), lowerBinThresh);
isInMidFifty  = @(x) x(:)>=prctile(x(:), 25) & x(:)<=prctile(x(:), 75);


    f           = formula(~isspace(formula));
    dYname  = f(1:strfind(f, '~')-1);
    dXname = f(strfind(f, '~')+1:end);

    dY_lowerBin   = d{isInLowerBin(d{:,dXname}), dYname}; %dY, for observations that fall in the lower bin of dX
    dY_upperBin   = d{isInUpperBin(d{:,dXname}), dYname}; %dY, for observations that fall in the upper bin of dX
    num_LowerBin = numel(dY_lowerBin);
    num_UpperBin = numel(dY_upperBin);
    
    if ~isempty(dY_lowerBin) && ~isempty(dY_upperBin)
        rnkTest = ranksum(dY_lowerBin, dY_upperBin, 'tail',tail);  
        dY_lowerBin_median = median(dY_lowerBin);
        dY_upperBin_median = median(dY_upperBin);
        dMedian = abs(dY_lowerBin_median - dY_upperBin_median);
        

        binCuttoffs = [prctile(d{:,dXname}, lowerBinThresh), prctile(d{:,dXname}, upperBinThresh)];
        binMedians = [dY_lowerBin_median dY_upperBin_median];
        
        H(1) = histogram(AXS, dY_lowerBin); axis square; hold on;
        H(2) = histogram(AXS, dY_upperBin); axis square;   
           
        box off; axis tight;
        set(gca,'color','none')
           
        axsProps   = properties(AXS);  
        histProps  = properties(H(1));  
        
    for k = 1:2:numel(plotPropValPairs)
        pltProp     = plotPropValPairs{k};
        pltPropVal  =  plotPropValPairs{k+1};
     
        if contains(pltProp, axsProps)   
            if strcmp(pltProp, 'YLim')
                pltProp = 'XLim';
                set(AXS, pltProp, pltPropVal)
            elseif strcmp(pltProp, 'XLim')
                continue;  
            end
        end
            
        if contains(pltProp, histProps)
            set(H(1), pltProp, pltPropVal)
            set(H(2), pltProp, pltPropVal)
        end
        

    end 
        

        
% FIT KERN HISTORGRAM        
%     sp = [-10 100];
        nrm = @(x) (x - min(x(:)))/(max(x(:)-min(x(:))));
        x = linspace(min(xlim), max(xlim), 100);
        pdf_dY_lower = fitdist(dY_lowerBin,'Kernel', 'width', H(1).BinWidth);
        pdf_dY_lower = nrm(pdf(pdf_dY_lower,x))*max(H(1).Values);
        pdf_dY_upper = fitdist(dY_upperBin,'Kernel','width', H(1).BinWidth);
        pdf_dY_upper = nrm(pdf(pdf_dY_upper,x))*max(H(2).Values);  
        plot(AXS, x,pdf_dY_lower,'LineWidth',.5, 'Color', 'b')
        plot(AXS, x,pdf_dY_upper,'LineWidth',.5, 'Color', 'r')
        
     % PROVIDE INFORMATION   
        xlabel([sprintf('LOWER BIN BELOW %.02f |Med =%.02f|n=%d)\n',binCuttoffs(1), dY_lowerBin_median, num_LowerBin),...
                sprintf('UPPER BIN ABOVE %.02f |Med =%.02f|n=%d)\n',binCuttoffs(2), dY_upperBin_median, num_UpperBin),...
                sprintf('dMed=%.02f | P(LOW is *%s* of UP) =%.01d)', dMedian, tail, rnkTest)],...
                'FontSize', 6)
                    

    
    else
        rnkTest = [];    
        binCuttoffs = [nan nan];
        binMedians = [nan nan];
        text(AXS, .5,.5, sprintf('Not Enough\nData Points'),...
            'units', 'normalized', 'horizontalalignment', 'center', 'fontsize', 6)        
    end
    
 

    figure(gcf)
end


