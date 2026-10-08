%% Example 2: six-intermediate-node line graph
% Run this script to compare uniform and time-varying capacities.
% Numerical results and figures are written to this example's results/ folder.
clear; clc; close all;
here=fileparts(mfilename('fullpath'));
addpath(fullfile(here,'helpers'));
example2_compare('line',fullfile(here,'results'));
