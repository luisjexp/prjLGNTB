function L = EXPERIMENTLOG

luisDir = cd;
main2pdataDir = lojGetDataPath; 

% %% TO DO...

% sbxpullsignals('c:\2pdata\lgn05\lgn03_MHO_003_001_WITH_003_000_rigid')
% lojhough_getKerns('c:\2pdata\lgn03\lgn03_003_001_rigid')

% lojhough_getKerns('c:\2pdata\lgn05\lgn05_006_002_rigid')
% lojhough_getKerns('c:\2pdata\lgn05\lgn05_007_000_rigid')


% lojhough_getKerns
%% Hard Coded
                                %[HOUGH EXP, ORISF EXP]
% Lgn02 mouse | all sessions processed
mouseNameList{1}            = 'lgn02';
imfield_expNames_All{1}     = {'_000_004', '_000_003';...%merged|aligned|segmented|pulled|split|estimated+ready
                              '_001_000', '_001_001';... %merged|aligned |segmented|pulled|split|estimated(lots of good ret maps)|ready 
                              '_002_000', '_002_002';... %merged|BADALIGNMENT
                              '_003_000', '_003_001'};   %merged|aligned|segmented|pulled|split|estimated+ready
                          
                          
% Lgn03 mouse | all potential sessions included 
mouseNameList{2}            = 'lgn03';
imfield_expNames_All{2}     = {'_001_001', '_001_000';...%merged|aligned|segmented|pulled|split|estimated|ready 
                               '_001_003', '_001_002';...%merged|aligned|segmented|pulled|split|estimated|ready 
                               '_002_001', '_002_000';...%merged|aligned|segmented|pulled|split|estimated(no good ret maps)|ready (Consider removing)
                               '_003_001', '_003_000';...%merged|aligned|segmented|needPulling
                               '_004_000', '_004_002';...%(damageduring imaging but try merging)|NEEDSALIGN                                
                               '_005_000', '_005_004'};  %merged|aligned|segmented|pulled|split|estimated(no good ret maps)|ready (Consider removing)
                           
% Lgn04 mouse | all potential sessions included       
mouseNameList{3}            = 'lgn04'; 
imfield_expNames_All{3}     = {'_000_000', '_000_005';... %needsMerging
                                '_001_000', '_001_001';...%merged|aligned|segmented|pulled|isEstimating
                               '_002_000', '_002_003'};   %merged|BADALIGNMENT

% Lgn05 mouse | all potential sessions included   
mouseNameList{4}            = 'lgn05'; 
imfield_expNames_All{4}     = { '_000_002', '_000_000';...% merged|ALIGNERROR
                                '_001_001', '_001_000';...% merged|aligned|segmented|pulled|split|estimated|ready 
                                '_002_001', '_002_000';...% merged|aligned|SEGMENTINGNOTDONE
                                '_003_004', '_003_005';...% merged|aligned|segmented|pulled|split|estimated(poor ret maps)|ready
                                '_004_001', '_004_000';...% merged|aligned|segmented|pulled|split|HOUGH_READLOG_ERROR
                                '_006_002', '_006_001';...% merged|aligned|segmented|pulled|split|ERROREEST
                                '_007_000', '_007_001';...% merged|aligned|segmented|pulled|split|ERROREEST                                
                                '_008_001', '_008_000';...% merged|aligned|segmented|pulled|split|estimated|ready                                
                                '_010_001', '_010_000';...% merged|aligned|segmented|pulled|split|isEstimating
                                '_011_000', '_011_001'};  % merged|aligned|segmented|pulled|split|estimated|ready
       
mouseNameList{4}            = 'lgn06'; 
imfield_expNames_All{4}     = { '_000_000', '_000_001';...% merged|alignmentHorrible|
                                '_001_000', '_001_001';...% merged|alignmentHorrible|
                                '_002_000', '_002_001';...% merged|alignmentHorrible|
                                '_003_000', '_003_001';...% merged|alignedBad|Segmented12|pulled|split|estimated|ready                         
                                '_004_000', '_004_002'};  % merged|alignedBad||Segmented26|pulled|split|estimatedOrisfOnly
                            
                            
%%
clc;

varNames = {'uif', 'mouse', 'houghExpId', 'orisfExpId', 'mergeSessId', 'houghExp_fName', 'orisfExp_fName', 'mergeExp_fName',...
            'hasMerged','hasAligned','hasSegment','hasSignal','hasSplit',...
            'hasHoughKern','hasHoughKernStats', 'hasOrisfKern', 'isReady'};

L = cell2table(cell(0,numel(varNames)));

uif = 0;
for i = 1:numel(mouseNameList)
    for j = 1:size(imfield_expNames_All{i},1)
        uif = uif+1;
        mouse = mouseNameList{i};
        dataPath = sprintf('%s%s/',main2pdataDir, mouse);    
        cd(dataPath)
        
        % File IDs
        houghExpId = imfield_expNames_All{i}{j,1};
        orisfExpId = imfield_expNames_All{i}{j,2};
        if strcmp(mouse,'lgn02') && strcmp('_000_004', houghExpId) && strcmp('_000_003', orisfExpId)
            mergeSessId = '_m';
        elseif  strcmp(mouse,'lgn05') && strcmp('_008_001', houghExpId) && strcmp('_008_000', orisfExpId)
            mergeSessId = '_m_008';
        elseif strcmp(mouse,'lgn05') && strcmp('_003_004', houghExpId) && strcmp('_003_005', orisfExpId)% flipped names on accident
            mergeSessId = sprintf('_MHO%s_WITH%s',houghExpId,orisfExpId); % keep incorrect merge file name            
            houghExpId = imfield_expNames_All{i}{j,2}; % just switch individual file names
            orisfExpId = imfield_expNames_All{i}{j,1};            
        else
            mergeSessId = sprintf('_MHO%s_WITH%s',houghExpId,orisfExpId);
        end
        
        % File Names
        mergeExp_fName  = [mouse, mergeSessId];        
        houghExp_fName  = sprintf('%s%s%s',mouse,houghExpId);
        orisfExp_fName  = sprintf('%s%s%s',mouse,orisfExpId);       
        
        hasMerged   = exist([mergeExp_fName, '.sbx'],'file');
        hasAligned  = exist([mergeExp_fName, '_rigid.sbx'],'file');        
        hasSegment  = exist([mergeExp_fName, '_rigid.segment'],'file');
        hasSignal   = exist([mergeExp_fName, '_rigid.signals'],'file');
        hasSplit    = hasSignal & exist([houghExp_fName, '_rigid.signals'],'file') & exist([orisfExp_fName, '_rigid.signals'],'file');
        hasHoughKern= hasSplit & exist([houghExp_fName, '_rigid.houghkernels'],'file');
        hasHoughKernStats= hasHoughKern & exist([houghExp_fName, '_rigid.houghkernstats'],'file');
        hasOrisfKern= hasSplit & exist([orisfExp_fName, '_rigid.orisf'],'file');
        
        isReady = hasOrisfKern&hasHoughKernStats;
        
        L = [L; cell2table({uif, mouse, houghExpId, orisfExpId, mergeSessId, houghExp_fName, orisfExp_fName, mergeExp_fName ,...
            hasMerged,hasAligned,hasSegment,hasSignal,hasSplit,...
            hasHoughKern,hasHoughKernStats,hasOrisfKern, isReady})];
        
%         if ~hasMerged
%             sbxmergesessions(mergeExp_fName,{houghExp_fName,orisfExp_fName})
%         else
% %             fprintf('\nDid Not Merge %s (already exists)\n', mergeExp_fName)
%         end
        
        
    end
end

L.Properties.VariableNames = varNames;
L(:,[1 6 7 8 9:end])
cd(luisDir);

end


