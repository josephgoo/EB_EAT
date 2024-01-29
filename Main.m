%% FUNCTION NAME: Main
% Main entry point function for key rate calculation.

format long 
clear all;
close all;
clc;

%%%%%%%%%%%%%%%%%%%%% Setting MATLAB Library Path %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
restoredefaultpath
cd /Users/joseph/Documents
addpath(genpath('/Users/joseph/Documents/GitHub/EB_EAT'))
addpath(genpath('/Users/joseph/Documents/MATLAB/cvx'))
addpath(genpath('/Users/joseph/Documents/MATLAB/YALMIP-master'))

cvx_solver mosek

%%%%%%%%%%%%%%%%%%%%% Setting User Input %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
preset='BB84';
[protocolDescription,channelModel,leakageEC,parameters,solverOptions]=feval(preset);

%%%%%%%%%%%%%%%%%%%%% Run Main Iteration %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%call main iteration function
results=mainIteration(protocolDescription,channelModel,leakageEC,parameters,solverOptions);
results.upperBound

%automatically parse and plot the results (optional)
%the third optional argument is the plotting style
%available options for 1D data:
%1.'linear': plot x and y linearly
%2.'linear-log': plot x versus log10(y)
%3.'km-log': plot -log10(x)*10/0.2 (x is assumed to be transmittance and converted to km) versus log(y)
%4.'none': do not plot
plotResults(results,parameters,'linear')








