%% BDG 2020
        
%% 

%% 0.Preparation
% 0.1.clear space
clear all;
close all;
clc;


% 0.2.set work directory
cd F:\GitHub\SVARToolBox;


% 0.3.load data
Prep_ImportData('Data_BDG2020.xlsx')    
load ..\SVARToolBox\DATASET DATASET_VAR;