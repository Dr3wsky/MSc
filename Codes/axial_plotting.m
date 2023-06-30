%% ParaView Data Plotting

% The following script is meant to interpret data from .csv exported from
% Paraview and plot desired variable(s) for the fgiven domain. 

% The idea is to eentually add some user-inout functionality but for now
% will be a brute-force ploitting script. 

clc
close all

%% Define Working Space and Read Data

load('exp_data.mat')

for i = 2:4
    mesh = i;
    loc = strcat('C:\Users\drewm\Documents\2.0 MSc\2.0 Simulations\PIV_JetSim_M', string(mesh),'\');
    type = sprintf("axial_data");
    eval(sprintf("data_M%s = readtable(strcat(loc, type, '.csv'), 'Format','auto');", string(mesh)))
    eval(sprintf("M%s_tke_dimless = data_M%s.k./(data_M%s.U_0.^2);", string(mesh), string(mesh), string(mesh)))
end

R_stress_legend = ['Rxx', 'Rxy', 'Rxz', 'Ryy', 'Ryz', 'Rzz'];
eval(sprintf("names = data_M%s.Properties.VariableNames;", string(mesh)));
eval(sprintf("xd = (data_M%s.Points_0-0.33533)/0.00635;", string(mesh)));
start = find(xd>0,1);
fin = find(xd>40, 1);
var = 'UMean_0';

exp_xd = exp_x(1,:)/6.35;
exp_start = find(exp_u_centreline == 0);

    
%% Plotting
markers = ["o" "s" "^" ">" "<" "*"];
colormap = ["#3498db"; "#2ecc71";  "#8e44ad"; "#e74c3c"; "#e67e22"; "#f1c40f"];

figure(1)
hold on 
for j = 2:4
eval(sprintf("plot(xd, data_M%s.%s)", string(j), var))
end
plot(exp_xd(11:end), exp_u_centreline(11:end)) % don't have data for initial jet exit due to PIV reflections

% plot(xd(start:fin), data_M2.U_0(start:fin))
% plot(xd(start:fin), data_M3.U_0(start:fin))
% plot(xd(start:fin), data_M4.U_0(start:fin))
% plot(data_Exp.Var1, data_Exp.Var3)
xlim([0 25])
ylim auto
ylabel('$Axial Centreline Velocity: U [m/s] $','interpreter','latex')
xlabel('$Non-Dimensional Axial Position: x/D_{jet}$','interpreter','latex')
title('Centreline Velocity Decay')
legend("RANS CFD, M2","RANS CFD, M3", "RANS CFD, M4","PIV Experimental")

% "RANS CFD, M1",

% legend('CFD RANS M1', 'CFD RANS M2', 'CFD RANS M3')

% figure(2)
% eval(sprintf("plot(xd(start:fin), data.%s(start:fin));", x_var))
% ylim auto
% xlim auto
% ylabel('U [m/s]')
% xlabel('x/d')