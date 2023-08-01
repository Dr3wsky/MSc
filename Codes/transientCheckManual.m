clc
close all
clear all

%% Load Data
start = 8;
nRuns = 9;

loc = 'C:\Users\drewm\Documents\2.0 MSc\2.0 Simulations\';
sim = 'RSM_LRR';
allFolder = dir(append(loc,sim,'\postProcessing-',string(start)));
allTimes = dir(append(loc, sim,'\postProcessing-',string(start),'\',allFolder(3).name));

rows = length(allTimes);
for i = start+1:nRuns
    rows = rows+1;
    tempDir = dir(append(loc, sim,'/postProcessing-',num2str(i),'/',allFolder(3).name));
    allTimes(rows).name = tempDir(3).name;
    allTimes(rows).folder = tempDir(3).folder;
    allTimes(rows).date = tempDir(3).date;
    allTimes(rows).bytes = tempDir(3).bytes;
    allTimes(rows).isdir = tempDir(3).isdir;
    allTimes(rows).datenum = tempDir(3).datenum;
end

%% Setup data for analysis