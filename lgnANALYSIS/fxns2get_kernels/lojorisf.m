function [r,stat] = lojorisf(fname)
% Orientation/Sf Kernel stats from randorisf experiment signals
% Modified version, by Luis
    

    % modify file name for reading logfile
    if contains(fname,'rigid') % search for rigid in filename
        si = strfind(fname,'_'); 
        fnameLog = fname( 1:si(end)-1); % remove it
    else 
        fnameLog = fname;
    end
    

%%
log = sbxreadorisflog(fnameLog); % read log
log = log{1};       % assumes only 1 trial
max_ori = max(abs(log.ori))+1;
max_sphase = max(abs(log.sphase))+1;
max_sper = max(abs(log.sper))+1;

load([fname, '.signals'],'-mat');    % load signals
if numel(size(spks)) == 3
    spks = squeeze(spks(1,:,:));
    sig = squeeze(sig(1,:,:));
end



dsig = spks ;

ncell = size(dsig,2);
nstim = size(log,1);


%%
tauwindow = [-2 17];
ntau = diff(tauwindow)+1;
tauvals = tauwindow(1):tauwindow(2);

r = zeros(max_ori,max_sphase,max_sper,ntau,ncell);
N = zeros(max_ori,max_sphase,max_sper);

disp('Processing...');
for(i=1:nstim)
        r(log.ori(i)+1,log.sphase(i)+1,log.sper(i)+1,:,:) =  ... 
            squeeze(r(log.ori(i)+1,log.sphase(i)+1,log.sper(i)+1,:,:)) + ...
            dsig(log.sbxframe(i)+tauwindow(1):log.sbxframe(i) + tauwindow(2),:);
        
         N(log.ori(i)+1,log.sphase(i)+1,log.sper(i)+1) = ...
            N(log.ori(i)+1,log.sphase(i)+1,log.sper(i)+1) + 1; 
end

% Normalize by stim count
for(idx_stdmax_tune=1:ntau)
    for(n = 1:ncell)
        r(:,:,:,idx_stdmax_tune,n) = r(:,:,:,idx_stdmax_tune,n)./N;
    end
end


% average across spatial phase

R = r;
r = squeeze(nanmean(r,2));

% now r = r(ori,sper,time,cell)

h = fspecial('gauss',5,1);
disp('Filtering');
k = 0;
for(idx_stdmax_tune=1:ntau)
    for(n = 1:ncell)
        rf = squeeze(r(:,:,idx_stdmax_tune,n));
        rf2 = [rf(end-1,:); rf(end,:); rf ;rf(1,:) ; rf(2,:)];
        rf3 = filter2(h,rf2,'same');
        r(:,:,idx_stdmax_tune,n) = rf3(3:end-2,:);
        k = k+1;
    end
end


% param values
ori = 0:10:170;
sper = 1920./(1.31.^(0:11));   % spatial period in pixels
sper_pix = sper;
sper = sper / 15.25;           % 15.25 pixels per degree 
sf = 1./sper;                  % cycles/deg
sphase = linspace(0,360,9);
sphase = sphase(1:end-1);


disp('Collecting stats...')

% [xx,yy] = meshgrid(-12:12,-12:12);
% zz = xx+1i*yy;
% zz = abs(zz).*exp(1i*angle(zz)*2);

% calculate for each case...


for i = 1:size(r,4)
    
    z = squeeze(r(:,:,:,i));
    
    if any(isnan(z(:)))
        fprintf('\n %dth kernel not processed\n', i)        
        continue;
    end
    
    q = reshape(z,max_ori*max_sper,[]);
    stdlist_tune = std(q);
    [stdmax_tune,idx_stdmax_tune] = max(stdlist_tune);
    
    tmax = idx_stdmax_tune;    
    stat(i).stdlist_tune = stdlist_tune;
    stat(i).idx_stdmax_tune = idx_stdmax_tune;
    stat(i).tau_maxstd_tune = tauvals(idx_stdmax_tune);
    
    stat(i).stdmax_tune = stdmax_tune;
    stat(i).kern_tune = squeeze(z(:,:,idx_stdmax_tune));

    imagesc(log10(sf),ori,stat(i).kern_tune);
    xlabel('Spatial Frequency (cycles/deg)');
    ylabel('Orientation (deg)');
    xval = get(gca,'xtick');
    l = cell(length(xval),1);
    for k = 1:length(l)
        l{k} = sprintf('%.2f',10^xval(k));
    end
    set(gca,'xticklabel',l);
    
% separability measure
    [u,s,v] = svd(stat(i).kern_tune);    
    stat(i).lambda12 = s(1,1)/s(2,2);
    
    % energy measure
    
    stat(i).delta = max(stat(i).stdlist_tune)/stat(i).stdlist_tune(1); 
    stat(i).sig_tune = stat(i).delta>1.75;
    
    % estimate preferred ori and sf    
    q                       = s(1,1);
    s                       = zeros(size(s));
    s(1,1)                  = q;
    kern_tune_smth          = u*s*v';
    stat(i).kern_tune_smth  = kern_tune_smth;
    
%     resp_sf = mean(kern_tune_smth);
%     resp_ori = mean(kern_tune_smth');
    
    [ii,jj] = find(stat(i).kern_tune == max(stat(i).kern_tune(:)));   % take slices through max
    
    sfresp  = stat(i).kern_tune(ii(1),:);
    oriresp = stat(i).kern_tune(:,jj(1))';

    stat(i).sfresp_tune  = sfresp;
    stat(i).oriresp_tune = oriresp;
    
    sfresp          = sfresp - max(sfresp)*.75; % clip tails...
    sfresp(sfresp<0)= 0;
    
    stat(i).oriest_tune = rad2deg( angle(sum(oriresp.*exp(1i*ori*2*pi/180)))/2 );
    stat(i).osi_tune    = resultant(oriresp, 2, ori);     % selectivity
    
    if(stat(i).oriest_tune<0)
        stat(i).oriest_tune = stat(i).oriest_tune+180;
    end
    stat(i).sfest_tune  = 10^(sum(sfresp.*log10(sf))/sum(sfresp));
    
    hold on
    plot(log10(stat(i).sfest_tune),stat(i).oriest_tune,'k.','markersize',15);
    hold off;
    
    title(sprintf('Cell #%d',i));
        
    drawnow;
%     pause(.05)
    % try to reconstruct linear rf - decimate by 4
%     [xx,yy] = meshgrid((1:1920/4)-1920/8,(1:1080/4)-1080/8);
%     rf = zeros(size(xx));
%     sp = sper_pix/4;
%     for m=1:max_ori
%         for j=1:max_sper
%             for k=1:max_sphase
%                 if(N(m,k,j)>0)
%                     sth = sind(ori(m));
%                     cth = cosd(ori(m));
%                     stim = cos (2*pi * (cth*xx + sth*yy) / sp(j) + sphase(k));
%                     rf = rf + R(m,k,j,stat(i).idx_stdmax_tune,i)*stim;
%                 end
%             end
%         end
%     end
% 
%     subplot(1,2,2)
%     imagesc(rf)
%     axis off;
    
%    pause;
end


%% Savee Configuration
config.ori      = ori;
config.sf       = sf;         % cycles/deg
config.ph       = sphase(1:end-1);
config.tauvals  = tauvals;


%% Stim Response for Noise Correlations
nc_r = cell(max_ori,max_sphase,max_sper,ncell);

disp('Computing Noise Correlations...');
for(i=1:nstim)
    for w = 1:ncell
        a = nc_r{log.ori(i)+1,log.sphase(i)+1,log.sper(i)+1,w};
        nc_r{log.ori(i)+1,log.sphase(i)+1,log.sper(i)+1,w} = [a dsig(log.sbxframe(i)+5,w)];
    end
end


for i = 1:ncell
    stat(i).stimresp_tune = cell(max_ori,max_sper);
    for o = 1:max_ori
        for s = 1:max_sper
            li = [];
            for p = 1:max_sphase
                li = [li(:) ; nc_r{o,p,s,i}(:)];
            end
            
            stat(i).stimresp_tune{o,s} = li;
        end
    end
end





%% get cell body positions
if exist([fnameLog, '_rigid.segment'],'file')
    disp('Finding cell body properties...');
    m   = load([fnameLog, '_rigid.segment'], '-mat');
    p   = regionprops(m.mask, 'centroid', 'Area'); 

    for i = 1:size(p,1)
    stat(i).xy_soma     = p(i).Centroid;
    stat(i).area_soma 	= p(i).Area;   
    end
    
else
    
     for i = 1:size(p,1)
         stat(i).xy_soma    = [nan nan];
         stat(i).area_soma 	= nan;   
    end   
end

%% Save Spikes 
disp('Storing Spikes...');

for i = 1:size(spks,2)
stat(i).spks_orisf     = spks(:,i);
stat(i).signal_orisf   = sig(:,i);
stat(i).K = squeeze(r(:,:,:,i));
end

%% Load and save running data

if exist([fnameLog '_quadrature.mat'], 'file')
   quad_data = load([fnameLog '_quadrature.mat']);
   runsig = quad_data.quad_data;
else 
    runsig = nan;
end

if exist([fnameLog '_eye.mat'], 'file')
   eyesig = quad_data;
else 
    eyesig = nan;
end

%% Save
disp('Saving...');

save([fname '.orisf'],'stat','config','runsig', '-v7.3');    

disp('Done!');

%% 

k = 0;
for i = 1:size(stat,2)
    if mod(i,64) == 0 
        F = figure;
        k = 0;
    end
    k = k+1;
    axs = subplot(8,8,k);
    imagesc(stat(i).kern_tune_smth); 
    set(axs, 'Units', 'normalized');
    axs.Position(1:2) = axs.Position(1:2) + [-.01 .01];
    axs.Position(3:4) = axs.Position(3:4) + .01;
    
    axis square off;
    drawnow;
    
end

end


%% Get OSI
function [sel, ptheta] = resultant(resp, f, thvals_deg)
    % resp: response vector

    % resultant
    result  = sum(resp.*exp(1i*f*thvals_deg*pi/180)); 
    
    % selectivity index
    sel     = abs(result)/sum(abs(resp)); 
    
    % est preffered ori
    ptheta    = rad2deg(angle(result)/f); 
    ptheta(ptheta < 0) = ptheta(ptheta < 0) +  360/f;

end



