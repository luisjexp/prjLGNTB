%% scratch pad for project: local tuning biases in mouse LGN
clear;
cd /Users/luis/Box/prjLGNTB/
addpath /Users/luis/Box/prjLGNTB/lgnDATA/   
addpath /Users/luis/Box/prjLGNTB/lgnANALYSIS/
addpath /Users/luis/Box/boxPHD/Toolbox/sbx/
%%

clc;
lgn = lgnanz
lgn.uif_tbl_base

[roi_table, uif_info] = lgn.get_stacked_roi_tables();


dist_stack = lgn.get_stacked_distance_table()
%% CREATE TABLE  
ROISCAT     = [];
DISTANCECAT = [];
numberof_uif = height(tbl_basic_project_info);

for uif = 1:numberof_uif
    mouse_id = tbl_basic_project_info.mouse_id{uif};
    houghExpId = tbl_basic_project_info.fname_base_hough{uif};
    orisfExpId = tbl_basic_project_info.fname_base_hough{uif};
    

    [rois, uifInfo] = getRoiTable(mouse, houghExpId, orisfExpId, unique_Im_field_number);
    ROISCAT     = [ROISCAT;rois];   

end



ROISCAT = [table((1:height(ROISCAT))', 'VariableNames', {'masterEntry'} ), ROISCAT];



for uif = 1:numel(unique(ROISCAT.uniqImFieldNum))
    ROIS = ROISCAT(ROISCAT.uniqImFieldNum == uif,:);
    [distances, distSelection] = getDistanceTable(ROIS);  
    DISTANCECAT = [DISTANCECAT; distances];        
end

close all
DISTANCECAT = [table((1:height(DISTANCECAT))', 'VariableNames', {'masterEntry'} ), DISTANCECAT];   

    
    
    
    
