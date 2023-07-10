%% ParaView Data Plotting

% The following script is meant to interpret data from .csv exported from
% Paraview and plot desired variable(s) for the fgiven domain. 

% The idea is to eentually add some user-inout functionality but for now
% will be a brute-force ploitting script. 

clc
% close all
clear all 

%% Define Working Space and Variables 
mesh = 4;
loc = strcat('C:\Users\drewm\Documents\2.0 MSc\2.0 Simulations\PIV_JetSim_M', string(mesh),'\');
type = "radial_data_xd_";
R_stress_legend = ['Rxx', 'Ryy', 'Rzz', 'Rxy', 'Ryz', 'Rxz'];
UPrime_2_Mean_legend = ['UPrime2Mean_xx', 'UPrime2Mean_xy', 'UPrime2Mean_xz', 'UPrime2Mean_yy', 'UPrime2Mean_yz', 'UPrime2Mean_zz'];
r_jet = 0.0032;
r_tube = 0.0254;

markers = ["o" "s" "^" ">" "<" "*"];
colormap = ["#f1c40f"; "#e67e22"; "#f51818"; "#f03acb"; "#8e44ad";  "#3498db"; "#2ecc71"];

%% Read Data 
j=0;
i=0;
data_0 = readtable(strcat(loc, type, '0.0', '.csv'), 'Format','auto');
names = data_0.Properties.VariableNames;
r_dimless = data_0.Points_1/data_0.Points_1(end);
x_var = "UMean_0";

for i =  10:2.5:25                                                   
    j = j + 1;
    xd(j) = i;                                                                                  % Assign radial position array
    data = readtable(strcat(loc, type, sprintf('%0.1f',xd(j)), '.csv'), 'Format','auto');       % Non-dimensionalize radial position
%     u_mag = sqrt((data.U_0.^2) + (data.U_1.^2) + (data.U_2.^2));                              % Calculate overall U_magnitude

    eval(sprintf("%s(:,%d) = data.%s;", x_var, j, x_var))                                       % Plot line
    
    %% Find location and calculate secondary stream for radial pos
    u_m(j) = max(UMean_0(:,j));                                                                 % Absolute max jet velocity
    counter = 0;
    for k = find(r_dimless>=(r_jet/r_tube),1):find(r_dimless>=0.9,1)
        counter = counter+1;
        UMean_0_trimmed(k-126,j) = UMean_0(k,j);                                                % Trim jet and wall effects
        r_trimmed(counter) = k/length(r_dimless);
    end
    
    % Use histogram to find average secondary velocity
    h = histogram(UMean_0_trimmed(:,j), 'BinWidth', 0.1);
    holder = h.Values;
    idx = find(holder == max(h.Values));                                    % index tallest bin (y value)
    U2(j) = mean(h.BinEdges(idx));                                          % param value of tallex bin
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
    

    %% Pull Reynolds Stresses & Normalize                                   % Need the b_prime, so cannot plot xd = 1, as is okay because it it not self-similat in this region anyways
    Rxx(:,j) = data.turbulenceProperties_R_0;
    Ryy(:,j) = data.turbulenceProperties_R_1;
    Rxy(:,j) = data.turbulenceProperties_R_3;

    % if j >= 2                                                             
    %     b_prime(j) = (b(j)-b(j-1))/(2.5*2*r_jet);
    %     Rxx_norm(:,j) = Rxx(:,j)/((U_excess_centreline(j)^2)*b_prime(j));
    %     Rxy_norm(:,j) = Rxy(:,j)/((U_excess_centreline(j)^2)*b_prime(j));
    %     Ryy_norm(:,j) = Ryy(:,j)/((U_excess_centreline(j)^2)*b_prime(j));
    % end

end


%% Collect comparison data & reorg for plotting
load exp_data.mat;


% exp_Rxy = -1*(exp_Rxy);

% count = 0;
% for i = 2:length(exp_x_loc)    
%     count = count +1;
%     % WRONG - ussue here with finding b_prime
%     exp_b_prime(count) = (exp_half_width(exp_x_loc(i))-exp_half_width(exp_x_loc(i-1)))/(2.5*2*r_jet);
%     exp_Rxx_norm(:,i) = exp_Rxx_g(:,i)/((exp_u_max(i)^2)*exp_b_prime(count));
%     exp_Rxy_norm(:,i) = exp_Rxy_g(:,i)/((exp_u_max(i)^2)*exp_b_prime(count));
%     exp_Ryy_norm(:,i) = exp_Ryy_g(:,i)/((exp_u_max(i)^2)*exp_b_prime(count));
% end




markers = ["o" "s" "^" ">" "<" "*"];

colormap = ["#e60049", "#0bb4ff", "#50e991", "#9b19f5", "#ffa300", "#dc0ab4", "#b3d4ff", "#00bfa0"];

%% Plot Results and Figures
close all

% % Velocity data and non-dimensionalized
% legend_array(1) = "Gaussian Shape Function";
% legend_array(2) = "Schlicting Jet Equation";
% legend_array(3) = "k-w SST M4";
% legend_array(4) = "PIV Experimenal";
% for j = 1:length(xd)-1
%     figure(j)
%     hold on
%     fplot(@(x) exp(-log(2)*x^2),'--k',[-2.5 2.5], 'LineWidth', 1);
%     fplot(@(x) (1+x^2/2)^(-2),'*k',[-1.5 1.5], 'LineWidth', 0.8);
%     plot(zeta(:,j), f_zeta(:,j), 'Color', colormap(j), 'LineWidth', 1.5);
%     plot(exp_zeta(:,exp_x_loc(find(exp_xd==xd(j)))), exp_Unorm(:,exp_x_loc(find(exp_xd==xd(j)))),markers(2),'MarkerSize',6,'MarkerEdgeColor', 'k', 'MarkerFaceColor',colormap(j));
%     eval(sprintf('title("x/D_{jet} = %s");', string(xd(j))));
%     legend(legend_array)
%     xlabel('$\eta = r/b$','interpreter','latex', 'FontSize', 14)
%     ylabel('$\frac{\overline{U}}{\overline{U}_{jet}}\ \ \ \ \ $','interpreter','latex', 'FontSize', 18)
%     hYLabel = get(gca,'YLabel');
%     set(hYLabel,'rotation',0,'VerticalAlignment','middle')
%     grid on
%     ylim([0 1])
%     xlim([0 2])
%     hold off
% end


% Plot Velocity data vs. exp results
% figure(5)
% for j = 1:length(b_vel)
%     subplot(2,3,j)
%     hold on
%     plot(r_dimless, U(:,j), '-b')
%     plot(r_dimless, Ryy(:,j), '*k', 'MarkerIndices',1:15:length(r_dimless), 'MarkerSize',4,'MarkerEdgeColor',colormap(3))
%      plot(r_dimless, Rxy(:,j), 'or', 'MarkerIndices',1:15:length(r_dimless), 'MarkerSize',4,'MarkerEdgeColor',colormap(5))
%     plot(exp_r_norm(:,exp_x_loc(j)), exp_Ryy_g(:,exp_x_loc(j)),markers(j),'MarkerSize',4,'MarkerEdgeColor',colormap(j))
%     hold off
%     ylim([0 2500])
%     xlim([0 1])
%     xlabel('$r/R$','interpreter','latex')
%     ylabel("$Reynolds Stress [m^2/s^2]$",'interpreter','latex')
%     hYLabel = get(gca,'YLabel');
%     set(hYLabel,'rotation',0,'VerticalAlignment','middle')
%     legend("$\overline {u'u'}$","$\overline {v'v'}$", "$\overline {u'v'}$",'interpreter','latex')
%     eval(sprintf("title('Reynold Stresses x/d = %0.1f')", xd(j)))
% end







    
    %% Future Calcs
    % Convective Mach No.
    % Vorticity thickness growth rate
