%% ParaView Data Plotting

% The following script is meant to interpret data from .csv exported from
% Paraview and plot desired variable(s) for the fgiven domain. 

% The idea is to eentually add some user-inout functionality but for now
% will be a brute-force ploitting script. 

clc
close all
clear all

%% Define Working Space and Read Data

% load exp Data
load('exp_data.mat')

% load k-w data
% for i = 1:4
mesh = 4;
loc = strcat('C:\Users\drewm\Documents\2.0 MSc\2.0 Simulations\PIV_JetSim_M', string(mesh),'\');
type = sprintf("axial_data");
eval(sprintf("data_kw_M%s = readtable(strcat(loc, type, '.csv'), 'Format','auto');", string(mesh)))
eval(sprintf("kw_M%s_tke_dimless = data_kw_M%s.kMean./(data_kw_M%s.UMean_0.^2);", string(mesh), string(mesh), string(mesh)))
% end

% load rsm data
model = 'LRR';
loc2 = strcat('C:\Users\drewm\Documents\2.0 MSc\2.0 Simulations\RSM_', string(model),'\');
eval(sprintf("data_RSM_%s = readtable(strcat(loc2, type, '.csv'), 'Format','auto');", string(model)))
eval(sprintf("RSM_%s_tke_dimless = data_RSM_%s.kMean./(data_RSM_%s.UMean_0.^2);", string(model), string(model), string(model)))
eval(sprintf("RSM_%s_tke_dimless_unsteady = data_RSM_%s.k./(data_RSM_%s.U_0.^2);", string(model), string(model), string(model)))

% Other properties and data handling
R_stress_legend = ['Rxx', 'Rxy', 'Rxz', 'Ryy', 'Ryz', 'Rzz'];
eval(sprintf("names = data_RSM_%s.Properties.VariableNames;", string(model)));
eval(sprintf("xd = (data_RSM_%s.Points_0-0.33533)/0.0064;", string(model)));
start = find(xd>0,1);
fin = find(xd>40, 1);
var = 'UMean_1';
% assign start-point to exp data so do not plot zeros.
exp_xd = exp_x(1,:)/6.4;
exp_start = find(exp_u_centreline == 0);

% Non-dimensionalize variable of choice
exp_u_centreline_nonDim = exp_u_centreline/max(exp_u_centreline);
eval(sprintf("kw_M4_U_nonDim(1,:) = data_kw_M%s.%s/mean(data_kw_M%s.%s(1:37));", string(mesh), var, string(mesh), var)); 
eval(sprintf("RSM_LRR_U_nonDim(1,:) = data_RSM_%s.%s/mean(data_RSM_%s.%s(1:37));", string(model), var, string(model), var)); 
% eval(sprintf("RSM_LRR_U_nonDim(1,:) = data_RSM_%s.%s/mean(data_RSM_%s.%s(1:37));", string(model), "U_1", string(model), "U_1")); 

    
%% Plotting
% Display options
markers = ["o" "s" "^" ">" "<" "*"];
colormap = ["#e74c3c"; "#e67e22";"#2ecc71"; "#8e44ad";  "#3498db"; "#f1c40f"];
linetype = ["-" "--" ":" "-." "-" "--" ":" "-."];

% Plot data
figure(1)
hold on 
eval(sprintf("plot(xd, kw_M%s_tke_dimless, linetype(1), 'Color', '%s', 'LineWidth', 1.25)", string(mesh),  colormap(4)))
eval(sprintf("plot(xd, RSM_%s_tke_dimless, linetype(2), 'Color', '%s', 'LineWidth', 1.25)", string(model),  colormap(5)))
eval(sprintf("plot(xd, RSM_%s_tke_dimless_unsteady, linetype(2), 'Color', '%s', 'LineWidth', 1.25)", string(model),  colormap(2)))
% plot(exp_xd(11:end-11), exp_v_centreline(11:end-7), 'LineWidth', 1.15, Color=colormap(5)) 

% plot(exp_xd(11:455), exp_v_centreline(11:end), 'LineWidth', 1.15, Color=colormap(5))

grid on
grid minor
xlim([0 100])
ylim auto
ylabel('$k/{u^2}_{jet}$' ,'interpreter','latex','FontSize', 14)
xlabel('$x/D_{jet}$','interpreter','latex', 'FontSize', 14)
title('Turbulent Kinetic Energy Profile', 'FontSize', 14)
legend("k-w SST, M4", "7-equation RSM, LRR", "7-equation RSM-unsteady, LRR")
% legend("M1: k-w SST", "M2: k-w SST","M3: k-w SST", "M4: k-w SST","Experiment: PIV")
% legend("M4: k-w SST","Experiment: PIV")
lgd.FontSize = 14;

% "RANS CFD, M1",

% legend('CFD RANS M1', 'CFD RANS M2', 'CFD RANS M3')

% figure(2)
% eval(sprintf("plot(xd(start:fin), data.%s(start:fin));", x_var))
% ylim auto
% xlim auto
% ylabel('U [m/s]')
% xlabel('x/d')