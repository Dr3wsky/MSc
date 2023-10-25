clc
close all
clear all
format longg

%% Load Folder Structure
start = 9;
nRuns = 10;

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
        %Skip Prime2MeanData, Not important right now
        if string(files(j).name) == "UPrime2Mean"
            continue
        end

        % Get variable name and import data
        data_name_holder = append(append(loc, sim,'\postProcessing-',num2str(nRuns),'\',fold_name,'\',num2str(latest_time),'/',files(j).name));
        disp(files(j).name)
        temp_data = readtable(data_name_holder);

        % Assign time var
        if i==9
            pp_data(:,1) = temp_data(:,1);     
            pp_data.Properties.VariableNames{1} = 'time';
        end
        cols = width(temp_data);

        % Assign header names and handle datatype conversions
        if string(temp_data.Properties.VariableNames{2}) == "sum_phi_"
            temp_data.Properties.VariableNames{2} = fold_name;

        elseif (string(temp_data.Properties.VariableNames{2}) == "Var2") && (files(j).name(1)~='U' && files(j).name(1)~='R')
            temp_data.Properties.VariableNames{2} = append(files(j).name, '_probe1250');

        elseif (files(j).name(1)=='U') && (string(files(j).name)~="UPrime2Mean")
            if string(files(j).name)=="UMean"
                temp_data_convU(1:height(temp_data),2) = str2double(erase(string(temp_data{:,2}), '('));
                temp_data_convU(1:height(temp_data),3) = temp_data{:,3};
                temp_data_convU(1:height(temp_data),4) = str2double(erase(string(temp_data{:,4}), ')'));
            else
                temp_data_convU(:,2) = str2double(erase(string(temp_data{:,2}), '('));
                temp_data_convU(:,3) = temp_data{:,3};
                temp_data_convU(:,4) = str2double(erase(string(temp_data{:,4}), ')'));
            end

        elseif (files(j).name(1)=='R')
            if string(files(j).name)=="RMean"
                temp_data_convR(1:height(temp_data),2) = str2double(erase(string(temp_data{:,2}), '('));
                temp_data_convR(1:height(temp_data),3:6) = temp_data{:,3:6};
                temp_data_convR(1:height(temp_data),7) = str2double(erase(string(temp_data{:,7}), ')'));
            else
                temp_data_convR(:,2) = str2double(erase(string(temp_data{:,2}), '('));
                temp_data_convR(:,3:6) = temp_data{:,3:6};
                temp_data_convR(:,7) = str2double(erase(string(temp_data{:,7}), ')'));
            end
        end

        % Fix to concat tables and handle array length di8screpancies
        if (height(temp_data)>height(pp_data)) || exist('temp_data_convU') || exist('temp_data_convR')
            if (files(j).name(1)=='U')
                pp_data = cat(2, pp_data, table(temp_data_convU(1:height(pp_data), 2)));
                pp_data.Properties.VariableNames{end} = append(files(j).name(:)', '_1_','probe1250');
                pp_data = cat(2, pp_data, table(temp_data_convU(1:height(pp_data), 3)));
                pp_data.Properties.VariableNames{end} = append(files(j).name(:)', '_2_', 'probe1250');
                pp_data = cat(2, pp_data, table(temp_data_convU(1:height(pp_data), 4)));
                pp_data.Properties.VariableNames{end} = append(files(j).name(:)', '_3_','probe1250');
            elseif (files(j).name(1)=='R')
                pp_data = cat(2, pp_data, table(temp_data_convR(1:height(pp_data), 2)));
                pp_data.Properties.VariableNames{end} = append(files(j).name(:)', '_xx_', 'probe1250');
                pp_data = cat(2, pp_data, table(temp_data_convR(1:height(pp_data), 3)));
                pp_data.Properties.VariableNames{end} = append(files(j).name(:)', '_yy_','probe1250');
                pp_data = cat(2, pp_data, table(temp_data_convR(1:height(pp_data), 4)));
                pp_data.Properties.VariableNames{end} = append(files(j).name(:)', '_zz_','probe1250');
                pp_data = cat(2, pp_data, table(temp_data_convR(1:height(pp_data), 5)));
                pp_data.Properties.VariableNames{end} = append(files(j).name(:)', '_xy_','probe1250');
                pp_data = cat(2, pp_data, table(temp_data_convR(1:height(pp_data), 6)));
                pp_data.Properties.VariableNames{end} = append(files(j).name(:)', '_yz_','probe1250');
                pp_data = cat(2, pp_data, table(temp_data_convR(1:height(pp_data), 7)));
                pp_data.Properties.VariableNames{end} = append(files(j).name(:)', '_xz_','probe1250');

            else
                pp_data = cat(2, pp_data, temp_data(1:height(pp_data), 2));
            end
        else
            pp_data = cat(2, pp_data, temp_data(:, 2:end));
        end
    end
end

%% Transient Data Analysis
format longg

% % Variable of interest
% time_var = pp_data.time;
% y_vars = [pp_data.probe1250_U_1, pp_data.probe1250_U_2, pp_data.probe1250_UMean_2, pp_data.probe1250_UMean_1, pp_data.probe1250_k, pp_data.jetOutletFlow];
% 
% for i = 1:length(y_vars) 
%     % Find periodicity
%     [pks,locs] = findpeaks(y_vars(:,i), time_var);
%     period(i) = max(diff(locs));
% end









