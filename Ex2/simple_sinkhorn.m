%% Example 2: three paths with shared nodal capacities
% Run this script to compare uniform and time-varying capacities.
% Numerical results and figures are written to this example's results/ folder.
clear; clc; close all;
here=fileparts(mfilename('fullpath'));
addpath(fullfile(here,'helpers'));
example2_compare('shared',fullfile(here,'results'));
