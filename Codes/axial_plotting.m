%% ParaView Data Plotting

% The following script is meant to interpret data from .csv exported from
% Paraview and plot desired variable(s) for the fgiven domain. 

% The idea is to eentually add some user-inout functionality but for now
% will be a brute-force ploitting script. 

clc
close all
% clear all

%% Define Working Space and Read Data

load('exp_data.mat')

for i = 1:4
    mesh = i;
    loc = strcat('C:\Users\drewm\Documents\2.0 MSc\2.0 Simulations\PIV_JetSim_M', string(mesh),'\');
    type = sprintf("axial_data");
    eval(sprintf("data_M%s = readtable(strcat(loc, type, '.csv'), 'Format','auto');", string(mesh)))
    eval(sprintf("M%s_tke_dimless = data_M%s.kMean./(data_M%s.UMean_0.^2);", string(mesh), string(mesh), string(mesh)))
end

R_stress_legend = ['Rxx', 'Rxy', 'Rxz', 'Ryy', 'Ryz', 'Rzz'];
eval(sprintf("names = data_M%s.Properties.VariableNames;", string(mesh)));
eval(sprintf("xd = (data_M%s.Points_0-0.33533)/0.0064;", string(mesh)));
start = find(xd>0,1);
fin = find(xd>40, 1);
var = 'UMean_0';


exp_xd = exp_x(1,:)/6.4;
exp_start = find(exp_u_centreline == 0);
% Add func to normalize u_centrline by u_centreline_max

for j = 1:4
    eval(sprintf("U_nonDim(%s,:) = data_M%s.%s/mean(data_M%s.%s(1:37));", string(j), string(j), var, string(j), var)); 
end
exp_u_centreline_nonDim = exp_u_centreline/max(exp_u_centreline);
    
%% Plotting
markers = ["o" "s" "^" ">" "<" "*"];
colormap = ["#e74c3c"; "#e67e22";"#2ecc71"; "#8e44ad";  "#3498db"; "#f1c40f"];

figure(1)
hold on 
for j = 1:4
eval(sprintf("plot(xd, M%s_tke_dimless, 'Color', '%s', 'LineWidth', 1.25)", string(j), colormap(j)))
end
% plot(exp_xd(11:end-11), exp_u_centreline(11:end-11), 'LineWidth', 1.15, Color=colormap(5)) 
% plot(exp_xd(11:455), exp_v_centreline(11:end), 'LineWidth', 1.15, Color=colormap(5))

grid on
grid minor
xlim([0 100])
ylim([-0.01 0.06])
ylabel('$Turbulent\ Kinetic\ Energy:\ k/{u^2}_{jet}$' ,'interpreter','latex','FontSize', 14)
xlabel('$Axial\ Position:\ x/D_{jet}$','interpreter','latex', 'FontSize', 14)
% title('Centreline Velocity Profile', 'FontSize', 14)
legend("M1: k-w SST", "M2: k-w SST","M3: k-w SST", "M4: k-w SST")
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