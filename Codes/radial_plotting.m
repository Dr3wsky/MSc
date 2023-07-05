%% ParaView Data Plotting

% The following script is meant to interpret data from .csv exported from
% Paraview and plot desired variable(s) for the fgiven domain. 

% The idea is to eentually add some user-inout functionality but for now
% will be a brute-force ploitting script. 

clc
close all
% clear all 

%% Define Working Space and Variables 
mesh = 4;
loc = strcat('C:\Users\drewm\Documents\2.0 MSc\2.0 Simulations\PIV_JetSim_M', string(mesh),'\');
type = "radial_data_xd_";
R_stress_legend = ['Rxx', 'Rxy', 'Rxz', 'Ryy', 'Ryz', 'Rzz'];
r_jet = 0.0032;
r_tube = .0254;

%% Read Data 
j=0;
i=0;
data_0 = readtable(strcat(loc, type, '0.0', '.csv'), 'Format','auto');
names = data_0.Properties.VariableNames;
r_dimless = data_0.Points_1/data_0.Points_1(end);
x_var = "UMean_0";

for i =  10:2.5:25
    j = j + 1;
    xd(j) = i;
    data = readtable(strcat(loc, type, sprintf('%0.1f',xd(j)), '.csv'), 'Format','auto');
%     u_mag = sqrt((data.U_0.^2) + (data.U_1.^2) + (data.U_2.^2));
    eval(sprintf("%s(:,%d) = data.%s;", x_var, j, x_var))


%% Radial Profile 

    figure(1)
    hold on
%     plot(u_mag, r_dimless)
    eval(sprintf("plot(%s(:,%d), r_dimless);",x_var, j)) 
%     plot(b_vel,b_dimless, 'Marker',  '.', 'Color', 'k', 'MarkerSize', 9)

end
hold off
    ylabel(sprintf('%s [m/s]', x_var), 'FontSize', 14)
    xlabel('Non-DImensional Radial Position: r/R', 'FontSize', 14)
    title('Radial Profile', 'FontSize', 14)
    xlim auto
    ylim auto
    legend ('x/d = 12.5', 'x/d = 15', 'x/d = 17.5', 'x/d = 20', 'x/d = 22.5')
% legend ('x/d = 15', 'x/d = 30', 'x/d = 45', 'x/d = 60', 'x/d = 75', 'x/d = 90','x/d = 105', 'x/d = 120', 'x/d = 135')