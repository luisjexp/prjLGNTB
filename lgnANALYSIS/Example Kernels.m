%% Example Figs
vars = {'pair_roiMasterEntry', 'dRetKern_corr','dTunKern_corr', 'dLogSfEst_lessThanOneOct', 'dOriEst'}

dIJ =['ij| ', num2str(table2array(DISTCAT(232,vars))) ]; %
dJK =['jk| ', num2str(table2array(DISTCAT(230,vars))) ]; %k157 vs j127
dIK =['ik| ', num2str(table2array(DISTCAT(252,vars))) ]; %k157 vs j127


i = 194;
j = 127;
k = 157;
subplot(3,2,1); imagesc(ROISCAT(i,:).kern_tun{:} ); axis off square;
subplot(3,2,2); imagesc(ROISCAT(i,:).kern_ret{:});axis off equal;
subplot(3,2,3); imagesc(ROISCAT(j,:).kern_tun{:} ); axis off square;
subplot(3,2,4); imagesc(ROISCAT(j,:).kern_ret{:}); axis off equal;
subplot(3,2,5); imagesc(ROISCAT(k,:).kern_tun{:} ); axis off square;
subplot(3,2,6); imagesc(ROISCAT(k,:).kern_ret{:}); axis off equal

text(-.9,-.04, sprintf('%s  |  ', vars{:}), 'Units', 'normalized', 'FontSize', 8, 'BackgroundColor', 'w')
text(-.8,-.3, [dIJ; dJK; dIK], 'Units', 'normalized', 'FontSize', 8, 'BackgroundColor', 'w')


set(gcf, 'color', 'none')


figure(gcf)