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
loc = strcat('C:\Users\drewm\Documents\2.0 MSc\2.0 Simulations\PIV_JetSim_M', string(mesh),'\');
type = "radial_data_xd_";
R_stress_legend = ["Rxx", "Ryy", "Rzz", "Rxy", "Ryz", "Rxz"];
UPrime_2_Mean_legend = ["UPrime2Mean_xx", "UPrime2Mean_xy", "UPrime2Mean_xz", "UPrime2Mean_yy", "UPrime2Mean_yz", "UPrime2Mean_zz"];
r_jet = 0.0032;
r_tube = 0.0254;

%% Read Data 
j=0;
i=0;
data_0 = readtable(strcat(loc, type, '0.0', '.csv'), 'Format','auto');
names = data_0.Properties.VariableNames;
r_dimless = data_0.Points_1/data_0.Points_1(end);
x_var = "UMean_0";

for i = 17.375:1.625:23.875                                                 
    j = j + 1;
    xd(j) = i;                                                                                  % Assign radial position array
    data = readtable(strcat(loc, type, sprintf('%0.3f',xd(j)), '.csv'), 'Format','auto');       % Non-dimensionalize radial position
%     u_mag = sqrt((data.U_0.^2) + (data.U_1.^2) + (data.U_2.^2));                              % Calculate overall U_magnitude

    eval(sprintf("%s(:,%d) = data.%s;", x_var, j, x_var))                                       % Plot line
    
    %% Find location and calculate secondary stream for radial pos
    u_m(j) = max(data.UMean_0(j));                                                                 % Absolute max jet velocity
    counter = 0;
    if (xd(j)<20.5)
        for k = find(r_dimless>=(r_jet/r_tube),1):find(r_dimless>=0.9,1)
            counter = counter+1;
            UMean_0_trimmed(k-126,j) = UMean_0(k-126,j);  
            UMean_0_trimmed(UMean_0_trimmed<= 0.5) = NaN;                   % Trim jet and wall effects
            r_trimmed(counter) = k/length(r_dimless);                       % Remove zeros for histogram
        end
    else 
         for k = find(r_dimless>=3*(r_jet/r_tube),1):find(r_dimless>=0.95,1)
            counter = counter+1;
            UMean_0_trimmed(k,j) = UMean_0(k,j);                            % Trim jet and wall effects
            UMean_0_trimmed(UMean_0_trimmed<= 0.5) = NaN;                   % Remove zeros for histogram
            r_trimmed(counter) = k/length(r_dimless);
         end
    end
    
    % Use histogram to find average secondary velocity
    h = histogram(UMean_0_trimmed(:,j), 'BinWidth', 0.1);
    holder = h.Values;
    idx = find(holder == max(h.Values));                                    % index tallest bin (y value)
    U2(j) = mean(h.BinEdges(idx));                                          % param value of tallex bin
    if (xd(j)==25)                                                          % Manually asssign U2 for xd = 25
        U2(j) = 32;
    end
    U_excess(:,j) = UMean_0(:,j) - U2(j);
    U_excess_centreline(j) = max(U_excess(:,j));
    close(figure)

    %% Find half-width and non-dimensionalize
    b_vel(j) = U_excess_centreline(j)/2;
    b_ind(j) = find(U_excess(:,j) < b_vel(j), 1, "first");
    b(j) = data.Points_1(b_ind(j));
    b_dimless(j) = b_ind(j)/length(r_dimless);
    % b(j) = b_dimless(j)*r_tube;
    zeta(:,j) = data.Points_1/b(j);
    f_zeta(:,j) = U_excess(:,j)/U_excess_centreline(j);
    

    %% Pull Reynolds Stresses & Normalize
    Rxx(:,j) = data.turbulenceProperties_R_0;
    Ryy(:,j) = data.turbulenceProperties_R_1;
    Rxy(:,j) = data.turbulenceProperties_R_3;
end




% % To normalize Reynolds Stresses
% for j = 2:length(xd)-1
%     b_prime(j) = ((b(j-1)+b(j+1))*.5)/(((xd(j-1)+xd(j+1))*.5)*2*r_jet);
%     Rxx_norm(:,j) = Rxx(:,j)/((U_excess_centreline(j)^2)*b_prime(j));
%     Rxy_norm(:,j) = Rxy(:,j)/((U_excess_centreline(j)^2)*b_prime(j));
%     Ryy_norm(:,j) = Ryy(:,j)/((U_excess_centreline(j)^2)*b_prime(j));
% end


%% Collect comparison data & reorg for plotting
load exp_data.mat;

% Assign xd array and find locations
scale = exp_x_raw(1,2,1)-exp_x_raw(1,1,1);
exp_xd = linspace(12.5, 19, 5);
exp_x_loc= round((exp_xd*(r_jet*2*1000)/scale));

% % Normalize y-component velocity
% for i = 1:length(exp_v_centreline)
% exp_Vnorm(:,i) = exp_v(:,i)/exp_v_centreline(i);
% end


% % To inverse sign direction for exp data
% exp_Rxy = -1*(exp_Rxy);

% % To calculate b_prime for exp data and normalize Reynolds Stresses
% exp_b_prime_fine(1) = 0;
% for i = 2:length(exp_x(1,:))-1
%     exp_b_prime_fine(i) = ((exp_b(i+1)-exp_b(i-1))*.5)/((exp_x(1,i+1)-exp_x(1,i-1))*.5);
% end
% exp_b_prime_fine(459) = 0;
% 
% for i = 2:length(xd)-1    
%     x_ind(i) = exp_x_loc(find(exp_xd==xd(i)));              % Parsing index for xd location 
%     x_ind(i-1) = exp_x_loc(find(exp_xd==xd(i-1)));
%     x_ind(i+1) = exp_x_loc(find(exp_xd==xd(i+1)));
%     exp_b_prime(i) = ((exp_b(x_ind(i-1))+exp_b(x_ind(i+1)))*.5)/((exp_x(1,x_ind(i-1))+exp_x(1,x_ind(i+1)))*.5);
%     exp_Rxx_norm(:,i) = exp_Rxx(:,x_ind(i))/((exp_Uexcess_centreline(x_ind(i))^2)*exp_b_prime(i));
%     exp_Rxy_norm(:,i) = exp_Rxy(:,x_ind(i))/((exp_Uexcess_centreline(x_ind(i))^2)*exp_b_prime(i));
%     exp_Ryy_norm(:,i) = exp_Ryy(:,x_ind(i))/((exp_Uexcess_centreline(x_ind(i))^2)*exp_b_prime(i));
% end





markers = ["o" "s" "^" ">" "<" "*"];
linetype = ["-" "--" ":" "-." "-" "--" ":" "-."];
colormap = ["#e60049", "#0bb4ff", "#50e991", "#9b19f5", "#ffa300", "#dc0ab4", "#b3d4ff", "#00bfa0"];

%% Plot Results and Figures
close all

% % % Plot Ken's data to recreate figs from manuscript
% % Initialize Legend array
% legend_array(1) = "Gaussian Shape Function";
% legend_array(2) = "Schlicting Jet Equation";
% 
% figure (1)
% hold on
% % Plot reference relations
% fplot(@(x) exp(-log(2)*x^2),'--k',[-2.5 2.5], 'LineWidth', 2.5);
% fplot(@(x) (1+x^2/2)^(-2),'b',[-1.5 1.5], 'LineWidth', 2);
% 
% % Loop through each xd position
% for j = 1:length(exp_xd)
%     legend_count = 2+j;
%     % plot(zeta(:,j), f_zeta(:,j), 'Color', colormap(1), 'LineWidth', 1.5);
%     plot(exp_zeta(:,exp_x_loc(j)), exp_f_zeta(:,exp_x_loc(j)),markers(j),'MarkerSize',8,'MarkerEdgeColor', colormap(j));
%     eval(sprintf('legend_array(%d) = "x/D_{jet} = %0.3f";', legend_count, exp_xd(j)));
% end
% legend(legend_array)
% xlabel('$\eta = r/b$','interpreter','latex', 'FontSize', 14)
% ylabel('$\frac{\overline{U}}{\overline{U}_{m}}\ \ \ \ \ $','interpreter','latex', 'FontSize', 18)
% hYLabel = get(gca,'YLabel');
% set(hYLabel,'rotation',0,'VerticalAlignment','middle')
% grid on
% ylim([0 1])
% xlim([-2.5 2.5])


% Plot non-dimaensionalized velocity comparison
legend_array(1) = "Gaussian Shape Function";
legend_array(2) = "Schlicting Jet Equation";

figure(1)
% subplot(2,3,j)          % plots on all one figure
hold on
fplot(@(x) exp(-log(2)*x^2),'--k',[-2.5 2.5], 'LineWidth', 2.5);
fplot(@(x) (1+x^2/2)^(-2),'b',[-1.5 1.5], 'LineWidth', 2);
legend_offset = 5;
for j = 1:length(xd)
    legend_count = 2+j;
    plot(zeta(:,j), f_zeta(:,j), markers(j),'MarkerSize',8,'MarkerEdgeColor', colormap(j), 'MarkerIndices',1:25:length(f_zeta));
    % plot(-zeta(:,j), f_zeta(:,j), linetype(j),'Color', colormap(j), 'LineWidth', 2);
    eval(sprintf('legend_array(%d) = "CFD k-w SST, x/D_{jet} = %0.3f";', legend_count, xd(j)));
    % plot(exp_zeta(:,exp_x_loc(j)), exp_f_zeta(:,exp_x_loc(j)), markers(j),'MarkerSize',8,'MarkerEdgeColor', colormap(j));
    % eval(sprintf('legend_array(%d) = "PIV Experiment, x/D_{jet} = %0.3f";', legend_offset, exp_xd(j)));
end
%Plot reverse side as well
for j = 1:length(xd)
    plot(-zeta(:,j), f_zeta(:,j), markers(j),'MarkerSize',8,'MarkerEdgeColor', colormap(j), 'MarkerIndices',1:25:length(f_zeta));
end
hold off
legend(legend_array(1:7))
xlabel('$\eta = r/b$','interpreter','latex', 'FontSize', 14)
ylabel('$\frac{\overline{U}}{\overline{U}_{m}}\ \ \ \ \ $','interpreter','latex', 'FontSize', 18)
hYLabel = get(gca,'YLabel');
set(hYLabel,'rotation',0,'VerticalAlignment','middle')
grid on
ylim([0 1.05])
xlim([-2.5 2.5])

%% Plot b_relation


% % Plot Reynolds Stresses
% figure(1)
% legend_array(1) = "k-w SST M4";
% legend_array(2) = "PIV Experimenal";
% 
% for j = 2:length(xd)-1
%     subplot(2,3,j-1)
%     hold on
%     plot(xd(j), b(j),'LineWidth', 1.5, 'Color', colormap(2))
%     plot(zeta(j), Ryy_norm(:,j), 'Color', colormap(2), 'LineWidth', 1.5)
%     plot(zeta(j), Rxy_norm(:,j), 'Color', colormap(3), 'LineWidth', 1.5)
%     plot(exp_xd(find(exp_xd==xd(j))), exp_b(exp_x_loc(find(exp_xd==xd(j)))), 'LineWidth', 1.5, 'Color', colormap(3))
% end
%     xlabel('$x/D_{jet}','interpreter','latex', 'FontSize', 14')
%     ylabel("$Half-Width:/ b\ \ \ \ $",'interpreter','latex', 'FontSize', 16)
%     ylabel({'Reynolds\ \ \ \ \ \ \ \ \ \ \ '; 'Stresses\ \ \ \ \ \ \ \ \ \ \ '; 'Normalized\ \ \ \ \ \ \ \ \ \ \ '}, 'interpreter','latex', 'FontSize', 12')
%     hYLabel = get(gca,'YLabel');
%     set(hYLabel,'rotation',0,'VerticalAlignment','middle')
%     eval(sprintf('title("x/D_{jet} = %s");', string(xd(j))));
%      legend(legend_array)
%     legend("$\overline{u'u'}$","$\overline{v'v'}$", "$\overline{u'v'}$",'interpreter','latex', 'FontSize', 14')
%     ylim auto
%     xlim([0 27.5])
%     grid on
%     grid minor


% For multi Reynolds Plotting: 
% plot(r_dimless, Rxx(:,j), '*k', 'MarkerIndices',1:15:length(r_dimless), 'MarkerSize',4,'MarkerEdgeColor',colormap(3))
% legend("$\overline{u'u'}$","$\overline{v'v'}$", "$\overline{u'v'}$",'interpreter','latex'






    
    %% Future Calcs
    % Convective Mach No.
    % Vorticity thickness growth rate
