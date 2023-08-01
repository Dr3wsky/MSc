clc
close all
clear all

%% Load Folder Structure
start = 7;
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

latest_time = allTimes(end).name;

vars = [9, 13, 19];
data_headers = ["time", "jetOutletFlow", "tubeInletFlow","probe_k", "probe_Ma", ...
    "probe_p", "probe_rho", "probe_T", "probe_U", "probe_UMean", "probe_UPrime2Mean" ];

%% Load Data
for i = vars
    fold_name = allFolder(i).name;
    files = dir(append(append(loc, sim,'\postProcessing-',num2str(nRuns),'\',fold_name,'\',num2str(latest_time))));

    for j = 3:height(files)
        data_name_holder = append(append(loc, sim,'\postProcessing-',num2str(nRuns),'\',fold_name,'\',num2str(latest_time),'/',files(j).name));
        temp_data = readtable(data_name_holder);
        % Assign time var
        if i==9                                 
            postProcessing_data(:,1) = temp_data(:,1);
        end
        disp(files(j).name)
        cols = width(temp_data);
        % Fix for same header names
        if string(temp_data.Properties.VariableNames{2}) == "sum_phi_"
            temp_data.Properties.VariableNames{2} = fold_name;
        elseif (string(temp_data.Properties.VariableNames{2}) == "Var2") && (files(j).name(1)~='U')
            temp_data.Properties.VariableNames{2} = append('probe1250_',files(j).name);
        elseif (files(j).name(1)~='U')
            keyboard 
        end
        % Fix to concat tables of diff height
        if height(temp_data)>height(postProcessing_data)
            postProcessing_data = cat(2, postProcessing_data, temp_data(1:height(postProcessing_data), 2));
        else
            postProcessing_data = cat(2, postProcessing_data, temp_data(:, 2:end));
        end

    end


end

