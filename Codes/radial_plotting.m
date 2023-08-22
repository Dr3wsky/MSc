 %% ParaView Data Plotting

% The following script is meant to interpret data from .csv exported from
% Paraview and plot desired variable(s) for the fgiven domain. 

% The idea is to eentually add some user-inout functionality but for now
% will be a brute-force ploitting script. 

clc
close all
clear all 

%% Define Working Space and Variables 
mesh = 4;
model = 'LRR';
loc = strcat('C:\Users\drewm\Documents\2.0 MSc\2.0 Simulations\RSM_', string(model),'\');
type = "radial_data_xd_";
Rstress_legend = ['Rxx', 'Rxy', 'Rxz', 'Ryy', 'Ryz', 'Rzz'];
UPrime2Mean_legend = ['UPrime2Mean_xx', 'UPrime2Mean_xy', 'UPrime2Mean_xz', 'UPrime2Mean_yy', 'UPrime2Mean_yz', 'UPrime2Mean_zz'];

r_jet = 0.0032;
r_tube = 0.025396;

%% Read Data 
j=0;
i=0;
data_0 = readtable(strcat(loc, type, '0.0', '.csv'), 'Format','auto');
names = data_0.Properties.VariableNames;
r_dimless = data_0.Points_1/data_0.Points_1(end);
x_var = "RMean_1";

markers = ["o" "s" "^" ">" "<" "*"];
colormap = ["#f1c40f"; "#e67e22"; "#f51818"; "#f03acb"; "#8e44ad";  "#3498db"; "#2ecc71";];

for i =  10:2.5:25
    j = j + 1;
    xd(j) = i;
    data = readtable(strcat(loc, type, sprintf('%0.1f',xd(j)), '.csv'), 'Format','auto');
    eval(sprintf("%s(:,%d) = data.%s;", x_var, j, x_var))
    %     u_mag = sqrt((data.U_0.^2) + (data.U_1.^2) + (data.U_2.^2));


%% Radial Profile 

    figure(1)
    hold on
    eval(sprintf("plot(r_dimless, %s(:,%d), 'Color', '%s', 'LineWidth', 1);",x_var, j, colormap(j))) 
%     plot(b_vel,b_dimless, 'Marker',  '.', 'Color', 'k', 'MarkerSize', 9)
%     plot(u_mag, r_dimless)
    eval(sprintf('legend_array(%d) = "x/D_{jet} = %s";', j, string(i))) 

end
hold off
    grid on
    xlim auto
    ylim auto
    legend(legend_array)
    ylabel('$U_x\ [m/s]$','interpreter','latex', 'FontSize', 14)
    xlabel('$Radial\ Position:\ r/R_{tube}$','interpreter','latex', 'FontSize', 14)
    title('Radial Velocity Profiles', 'FontSize', 14)
    xlim auto
    ylim auto
    legend(legend_array)