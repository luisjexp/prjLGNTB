function S = lojhough_getKernStats(fname)
%% returns stats of hough kernels (need to process using hough_getKerns first)


%% I. INITIALIZE

% Load ON and OFF Maps
% If kernels exist, aka if have been extracted previously
if exist([fname '.houghkernels'], 'file') 
    fprintf('\n\nLoading kernels (if there are a large number of rois, this may take a while (>30sec) )...\n')
    K   = load([fname '.houghkernels'], '-mat'); % load kernels
    K   = K.K;
    SIG = load([fname '.signals'], '-mat'); % load spike rate time series
    if ndims(SIG.sig) == 3 
        sig = squeeze(SIG.sig(1,:,:)); % load green channel
        spks = squeeze(SIG.spks(1,:,:)); 
    else
        sig = SIG.sig; % load green channel
        spks = SIG.spks;         
    end
    nroi = size(K,3);
    ntau = size(K,2);
    h = fspecial('gauss',250,15);
    
% Otherwise if not extracted call error
elseif ~exist([fname '.houghkernels'], 'file')
    error('No hough kernels found for this filename. Use the lojhough_getKerns fxn to extract kernels first')
end    
    
Fig = figure;


for i = 1:nroi  
    k = cell(ntau,1);

    for j = 1:ntau
        k{j} = reshape(K(:,j,i),  1081, 3841 ); %3841 x 1081
        k{j} = imresize( filter2(h,k{j}, 'same'),0.5 );
    end
       
    % Variance, Kurtosis, and Magnitute of Kernels at each tau
    S(i).varlist      = cellfun(@(k) var(k(:)), k) ;    
    S(i).kurtlist     = cellfun(@(k) kurtosis(k(:)), k) ;
    S(i).amplist      = cellfun(@(k) max(k(:)), k);
    
    % Max Variance, Max Kurtosis, and Max Magnitute of Kernels
    [S(i).varmax, S(i).maxvar_tau]       = max(S(i).varlist);    
    S(i).kurtmax   = max(S(i).kurtlist);    
    S(i).ampmax    = max(S(i).amplist);     

    % Compute On and Off SNR
	S(i).snr      = S(i).varlist(S(i).maxvar_tau)/S(i).varlist(1);
         
    % Kernels at their  optimal taus 
    S(i).kern        = k{S(i).maxvar_tau};  
   
    % Find kernel Peaks Location       
    [y0, x0]  = find(S(i).kern == max(S(i).kern(:)) );
    S(i).xy      = [x0(1)    y0(1)];

    % Spikes
    S(i).spks     = spks(:,i)';
    S(i).sig  = sig(:,i)';

    
    subplot(3,1,2); cla;
    subplot(3,1,3); cla;
    
    S(i).gausFit = [];
    
%     if S(i).kurtmax > 5
%         try
% 
%             S(i).gausFit = fitGaussian(S(i).kern, S(i).xy);
% 
%             subplot(3,1,2); cla;
%                 plot(S(i).gausFit); view(0,90); grid off; axis equal off; 
%                 hold on
%                 plot(S(i).xy(1), S(i).xy(2) , 'r*');    
% 
%              subplot(3,1,3); cla;
%                     plot(S(i).gausFit);            
%         catch
% 
%         end
%     end
    
    % Plot
    subplot(3,1,1); cla;
        imagesc(S(i).kern); 
        axis equal off
        title(sprintf('%d/%s',i,nroi ));
        
    drawnow;
    
end
close(Fig)

S = orderfields(S);

save([fname '.houghkernstats'],'S', '-v7.3');
fprintf('\nSaved...\n %s\n', [fname '.houghkernstats'])        
disp('Done!');
    

end



function G = fitGaussian(K, cXcY)

    %% Fit Gaussian to Retino kernels
    [x, y] =meshgrid(1:size(K,2),1:size(K,1));
    x = x(:);
    y = y(:);
    G  = @(cx,cy, sigma, A, b0, x, y)  b0 + A*exp(-(   (((cx-x).^2)+((cy-y).^2)) /(2 * sigma^2) ) );
    nrm = @(s) (s - min(s(:)))/ (max(s(:)) - min(s(:)));
    
    kProps = regionprops((reshape(nrm(K), size(K)) > .9), 'EquivDiameter', 'boundingbox');
    z = K(:);

    startSigma  = round(kProps(1).EquivDiameter/4);
    startx      = cXcY(1);
    starty      = cXcY(2);
    startA      = max(K(:));
    startb0     = max(K(:))*.1;

    try
        G = fit([x,y],z,G, 'StartPoint',[startx starty startSigma startA startb0]);
    catch

    end
end



