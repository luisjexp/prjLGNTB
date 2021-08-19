%% ALIGNMENT EXAMPLE
uifID = 1;
ImAlign = UIFCAT.alignImage{uifID};
ImSegment = UIFCAT.segMask{uifID}>0;


clf
imagesc(ImAlign); axis equal; axis off 
imcontrast

colormap('winter')
figure(gcf)

%% ORI PREF DISTRIBUTIONS
cst = ROISCAT.signif_tun == 1;
R = ROISCAT(cst,:);

roiPropHistogram(R, 'oriest',...
    'BinWidth', 20, 'XLim', [0 180], 'XTick', [0 45 90 135 180])

%% SPFRQ PREF DISTRIBUTIONS
cst = ROISCAT.signif_tun == 1;
sf = log2([0.0079 0.0104 0.0136 0.0179 0.0234 0.0306 0.0401 0.0526 0.0689 0.0902 0.1182 0.1549]);


R = ROISCAT(cst,:);
% R.sfPeak = log2(R.sfPeak);
roiPropHistogram(R, 'sfPeak')
% 'BinWidth', 0.025, 'BinLimits', [min(R.sfPeak) max(R.sfPeak)] )
    
%     'BinWidth', 0.1, 'BinLimits', [sf(1) sf(end)] )


%% SPFRQ PREFxORIPREF DISTRIBUTIONS
clf
cst = ROISCAT.signif_tun == 1;
R = ROISCAT(cst,:);
roiPropHistogram(R, {'sfPeak','oriest'},...
    'YBinEdges', [0:30:180],...    
    'EdgeAlpha',.5, 'FaceAlpha', .25)


%% SPFRQ PREFx ORISelect DISTRIBUTIONS
clf
cst = ROISCAT.signif_tun == 1;
R = ROISCAT(cst,:);
roiPropHistogram(R, {'sfPeak','osi'},...
    'YBinEdges', [0:.1:1],...    
    'EdgeAlpha',.5, 'FaceAlpha', .25)




%%
nrm = @(x) (x - min(x(:)))/(max(x(:)-min(x(:))));
clf;
cst = ROISCAT.signif_tun == 1 & ROISCAT.kurt_ret>13;
R = ROISCAT(cst,:);
kTunAggregate = zeros(size(R.kern_tun{1})) ;
for i = 1:height(R)
   kTunAggregate = kTunAggregate + double(nrm(R.kern_tun{i} )>.95);  
   imagesc(kTunAggregate)
pause(.1)
figure(gcf)

end
figure(gcf)

[sel, ptheta] = lojResultant(resp, f, thvals_deg);


%%  Kernel(unfilt)CORR VS RF CORR, r = .19
plotProps = {'YLim', [min(DISTCAT.dTunKern_corr) max(DISTCAT.dTunKern_corr)],...
    'XLim', [min(DISTCAT.dRetKern_corr) 1],...
    'BinWidth', .25,...
    'Normalization', 'count',...
    'ShowOthers', 'On'}; 

M = loj_dX_vs_dY_by_dC(  DISTCAT,...
    'dTunKern_corr ~ dRetKern_corr',...
    {'dLogSfEst_lessThanOneOct', 'dLogSfEst_moreThanOneOct'},...
    {'both', 'both', 'both'},...    
    plotProps);

%%

plotProps = {'YLim', [0 max(DISTCAT.dOriEst)],...
    'XLim', [min(DISTCAT.dRetKern_corr) 1],...
    'BinWidth', 10,...
    'Normalization', 'count',...
    'ShowOthers', 'On'}; 

M = loj_dX_vs_dY_by_dC(  DISTCAT,...
    'dOriEst ~ dRetKern_corr',...
    {'dLogSfEst_moreThanOneOct','dLogSfEst_lessThanOneOct' },...
    {'both', 'both', 'both'},...    
    plotProps);


%%  Kernel(unfilt)CORR VS RF CORR, r = .19
plotProps = {'YLim', [min(DISTCAT.dTunKern_corr) max(DISTCAT.dTunKern_corr)],...
    'XLim', [min(DISTCAT.dRetKern_c2cDist) max(DISTCAT.dRetKern_c2cDist)],...
    'Normalization', 'count',...
    'ShowOthers', 'On'}; 

M = loj_dX_vs_dY_by_dC(  DISTCAT,...
    'dTunKern_corr ~ dRetKern_c2cDist',...
    {'dLogSfEst_lessThanOneOct', 'dLogSfEst_moreThanOneOct'},...
    {'both', 'both', 'both'},...    
    plotProps);




