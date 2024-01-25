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

%%%%%%%%%%%%%%%%%%%%% Step 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[gval]=Step1(protocolDescription,channelModel,leakageEC,parameters,solverOptions);

%%%%%%%%%%%%%%%%%%%%% Step 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%user should supply the list of parameters used in this description/channel file
%this list varNames should be a subset of the full parameter list declared in the preset file
%parameters specified in varNames can be used below like any other MATLAB variables
varNames=["pz"];

%%%%%%%%%%%%%%%%%%%%% interfacing (please do not modify) %%%%%%%%%%%%%%%%%%%%%%%%%

%the functions findVariables and addVariables automatically search the input (names,p) for
%the parameter values based on varNames, and convert them to MATLAB variables.
varValues = findVariables(varNames,names,p);
addVariables(varNames,varValues);

%%%%%%%%%%%%%%%%%%%%% user-supplied description begin %%%%%%%%%%%%%%%%%%%%%%%%%
%maximum and minimum of g
gMax=max(gval);
gMin=min(gval);

%coefficient of q(0) and q(1) in g
gzero=gval(1);
gone=gval(2);

%coefficient of q(0) in f
fzero=gMax+(1/(1-pz))*(gzero-gMax)











