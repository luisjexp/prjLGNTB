for k = 1:height(DISTCAT)
    
entry = k;
i = DISTCAT(k,:).pair_roiMasterEntry(1);
j = DISTCAT(k,:).pair_roiMasterEntry(2);

dIJ =['ij| ', num2str(table2array(DISTCAT(entry,vars))) ]; %
subplot(2,2,1); imagesc(ROISCAT(i,:).kern_tun{:} ); axis off square;
subplot(2,2,2); imagesc(ROISCAT(i,:).kern_ret{:});axis off equal;
subplot(2,2,3); imagesc(ROISCAT(j,:).kern_tun{:} ); axis off square;
subplot(2,2,4); imagesc(ROISCAT(j,:).kern_ret{:});axis off equal;

text(-.9,-.04, sprintf('%s  |  ', vars{:}), 'Units', 'normalized', 'FontSize', 8, 'BackgroundColor', 'w')
text(-.8,-.3, [dIJ], 'Units', 'normalized', 'FontSize', 8, 'BackgroundColor', 'w')
set(gcf, 'color', 'none')

pause;
clf;

end


