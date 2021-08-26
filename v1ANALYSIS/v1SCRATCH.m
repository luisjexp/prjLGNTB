clear;
clc;
close all;
addpath /Users/luis/Box/prjLGNTB/v1DATA/   
addpath /Users/luis/Box/prjLGNTB/v1ANALYSIS/
cd /Users/luis/Box/prjLGNTB/v1ANALYSIS

v1  = v1anz();
v1.get_roi_table

%%
clc
v1.get_distance_table()

%% - - - 
clc
uif_list    = v1.dist_stack.uniqImFieldNum;
dY_name     = 'dTunKern_corr';

dON_name    = 'dON_corr';
dOFF_name   = 'dOFF_corr';

close all
axis_handle_dON     = axis;

%% - - -  LINEAR MODEL
figure(1)
v1.fit_plot_lm(uif_list, dY_name, dON_name, axis_handle = axis_handle_dON);
% axis_handle_dOFF     = subplot(2,1,2);
figure(2)
v1.fit_plot_lm(uif_list, dY_name, dOFF_name, axis_handle = axis_handle_dON);
%% --- RANK SUM

v1.test_and_plot_ranksum(uif_list, dY_name, dOFF_name)

