%% scratch pad for project: local tuning biases in mouse LGN
clear;
cd /Users/luis/Box/prjLGNTB/
addpath /Users/luis/Box/prjLGNTB/lgnDATA/   
addpath /Users/luis/Box/prjLGNTB/lgnANALYSIS/
addpath /Users/luis/Box/boxPHD/Toolbox/sbx/

%%
clc;
lgn = lgnanz;
lgn.uif_tbl_base


%% GENERATE ROI DATA TABLE FOR IMAGING FIELDS, and stack them 






uifs_to_process = [3 4 5];
lgn.get_stacked_roi_tables(uifs_to_process);

head(lgn.roi_stack) % now you have a table where each row corresponds to an roi in the imaging field,

%% generate a distance data table from each imaging field's ROI table, and stack them
lgn.get_stacked_distance_table(); % 
head(lgn.dist_stack)    % now you have a table where each row corresponds to a pair of rois, 
summary(lgn.dist_stack)


%% Linear Models
figure(1)
clf;
%   Run the model for data in one imaging field
uif  = lgn.roi_stack_uif_list(1);
lgn.fit_plot_lm(uif);

%   Run the test for combined data accross imaging fields
uif_stack       = lgn.roi_stack_uif_list;
lgn.fit_plot_lm(uif_stack);



%% Rank Sum Tests
figure(2)
clf;
%   Run the test for data in one imaging field
uif     = lgn.roi_stack_uif_list(1);
lgn.test_and_plot_ranksum(4);

%   Run the test for combined data accross imaging fields
uif_stack   = lgn.roi_stack_uif_list;
lgn.test_and_plot_ranksum(uif_stack);


%%
dist_var_name = 'dTunKern_corr';
plotProps = {}; 

axis_handle = subplot(1,1,1);
uif = 4;


lgn.plot_univar_hist(uif, dist_var_name)
lgn.fit_plot_lm()

clf
lgn.test_and_plot_ranksum()

