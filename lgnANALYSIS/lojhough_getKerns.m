function lojhough_getKerns(fname)
%% returns hough kernels at multuple time delays 
clearvars -except fname

%%  Allow for rigid and nonrigid .signals files as input..
    if contains(fname,'rigid') 
        si = strfind(fname,'_'); 
        fnamelog = fname( 1:si(end)-1); 
    else 
        fnamelog = fname;
    end
    
    
%%
log = sbxreadhoughlog(fnamelog); % read log
log = log{1};       % assumes only 1 trial

load([fname, '.signals'],'-mat');    % load signals
if(ndims(spks)>2)
    spks = squeeze(spks(1,:,:));
end


bad = find(isnan(sum(spks))); % some bad ROIs?
spks(:,bad) = abs(randn(size(spks,1),length(bad))); % replace with white noise

ncell = size(spks,2);
nstim = size(log,1);
tauWindow = [2 13]; 
ntau = numel(tauWindow(1):tauWindow(end));

clear np sig si bad
fprintf('\nProcessing...\n%s\n', fname)        

%%
[xx,yy] = meshgrid(-1920:1920,-540:540); %3841 x 1081
P = single([xx(:) yy(:)]);
K = zeros(numel(xx),ntau,ncell,'single');
spks = single(spks);
clear xx yy

for m=1:nstim
    tic
    fprintf('\nStim #%d/%d\n \t', m,nstim)   
    
    z = abs(P * [cos(log.angle(m)) sin(log.angle(m))]' - log.dist(m)) < 90; % image
    B = reshape(spks(log.sbxframe(m)+tauWindow(1):log.sbxframe(m)+tauWindow(end),:),1,ntau*ncell); % weight
    K = K + reshape(z * B,size(K));
    
    fprintf('DONE!(dur:%.02f)', toc);

end
    
  
% for m=1:nstim
%     z = abs(P * [cos(log.angle(m)) sin(log.angle(m))]' - log.dist(m)) < 90;
%     km = zeros(size(xx));
%     km = reshape(z,size(km));
% 
%     fprintf('\nStim #%d/%d\n \t', m,nstim)       
%     for i = 1:ncell
%         fprintf(' .');
%         for tau = 1:ntau
%             y = dsig(log.sbxframe+tau-2,i);
%             K{tau, i} = K{tau, i}+y(m)*km;
%         end 
%     end
%     fprintf(' DONE!');
%     
% end

disp('Saving...');
extention = '.houghkernels';
save([fname extention],'K', '-v7.3');
fprintf('\nSaved as...\n%s\n', [fname extention])        
disp('Done!');

end