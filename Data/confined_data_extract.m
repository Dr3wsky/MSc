% close all;
clear;
clc;

%Characteristic Velocity Profile

rawI = 459; % Number of points in x-direction in raw data (found in header of .dat DaVis files)
rawJ = 130; % Number of points in y-direction in raw data (found in header of .dat DaVis files)

jet_diameter = 6.4;

% 
data = dlmread('B0001.dat',' ',4, 0);

x_g1=data(:,1);
y_g1=data(:,2);
u_g1=data(:,3);
v_g1=data(:,4);
Rxx_g1=data(:,6);
Rxy_g1=data(:,7);
Ryy_g1=data(:,8);

x_g=permute(reshape(x_g1,rawI,rawJ,[]),[2 1 3]);
y_g=permute(reshape(y_g1,rawI,rawJ,[]),[2 1 3]);
u_g=permute(reshape(u_g1,rawI,rawJ,[]),[2 1 3]);
v_g=permute(reshape(v_g1,rawI,rawJ,[]),[2 1 3]);
Rxx_g=permute(reshape(Rxx_g1,rawI,rawJ,[]),[2 1 3]);
Rxy_g=permute(reshape(Rxy_g1,rawI,rawJ,[]),[2 1 3]);
Ryy_g=permute(reshape(Ryy_g1,rawI,rawJ,[]),[2 1 3]);


scale = x_g(1,2,1)-x_g(1,1,1);
%% Locate Secondary Stream

u_diff = zeros(69,rawI);
r_peak_location = zeros(rawI,1);
u_peak_location = zeros(rawI,1);
y_centered = zeros(rawJ,rawI);
y_guess = zeros(rawJ,rawI);
sigma_udiff = zeros(rawI,1);
u_thickness = zeros(rawI,1);
u_secondary_index = zeros(rawI,1);
u_secondary = zeros(rawI,1);
u_excess = zeros(rawJ,rawI);
u_centerline = zeros(rawI,1);
v_centerline = zeros(1,rawI);
u_normalized = zeros(rawJ,rawI);
u_max = zeros(rawI,1);
exclude = zeros(rawI,1); 
exclude(1) = 20; %starting value of exclude
u_excess_peak_location = zeros(rawI);
vorticity_thickness = zeros(rawI,1);

%hold on
for i = 1:rawI
   options = fitoptions('gauss1');
   options.Lower = [10 23.5 0];
   options.Upper = [500 26.5 inf];
     
   %fit Gaussian to velocity profile to estimate center for jet
   f_1 = fit(-y_g(:,i),u_g(:,i),'gauss1',options);
   coeff_1 = coeffvalues(f_1);
   u_peak_location(i) = coeff_1(2);
   y_guess(:,i) = -y_g(:,i) - u_peak_location(i);
   
    %figure(1)
    % plot(f_1, -y_g(:,i),u_g(:,i))   
 
   options = fitoptions('gauss1');
   options.Lower = [.1 0 0];
   options.Upper = [300 15 20];
   options.Exclude = y_guess(floor(rawJ/2):end,i) > exclude(i);
   
   u_diff(:,i) = gradient(u_g(floor(rawJ/2-rawJ/64):end,i),-y_g(floor(rawJ/2-rawJ/64):end,i));
   f_2 = fit(y_guess(floor(rawJ/2-rawJ/64):end,i),-u_diff(:,i),'gauss1',options);
   coeff_2 = coeffvalues(f_2);
   sigma_udiff(i) = coeff_2(3)/sqrt(2); %3 S.D. 
   r_peak_location(i) = coeff_2(2); %location from center
% %      
   %find location of secondary stream
   exclude(i+1) =  r_peak_location(i) + 4*sigma_udiff(i); %excluse the boundary layer at the mixing tube wall
   u_thickness(i) = r_peak_location(i) + 3*sigma_udiff(i); %
% % 

%     figure(2)
%     plot(f_2, y_guess(floor(rawJ/2-rawJ/64):end,i),-u_diff(:,i))
%     xline(u_thickness(i))
   
   %Remove secondary stream
   u_secondary_index(i) = ceil(u_thickness(i)/scale) + floor(rawJ/2); %find Index of secondary stream
   if u_secondary_index(i) > rawJ-1
       u_secondary_index(i) = rawJ-1;
   else
   end
   
   u_secondary(i) = u_g(u_secondary_index(i),i);
   u_excess(:,i) = u_g(:,i)-u_secondary(i); %calculating the excess velocity
   u_normalized(:,i) = max(u_excess(:,i),0); %removing negatives
   u_max(i) = max(u_normalized(:,i)); %calculating maximum excess velocity
   u_normalized(:,i) = u_excess(:,i)/u_max(i); %normalizing by excess velocity
   u_normalized(isnan(u_normalized))=0; %replaceing nan values with zero
   
   options = fitoptions('gauss1');
   options.Lower = [.9 -25 -inf];
   options.Upper = [1.2 25 inf];
   
   f_3 = fit(y_guess(:,i),u_normalized(:,i),'gauss1');

%  figure(3)
   % plot(f_3, y_guess(:,i),u_normalized(:,i))  
   % xline(u_thickness(i))
   
   coeff_3 = coeffvalues(f_3);
   u_excess_peak_location(i) = coeff_3(2);
   y_centered(:,i) = y_guess(:,i) - u_excess_peak_location(i);
   %u_normalized = max(u_normalized,0);

   %find centerline peak velocity
   if length(find(u_normalized(:,i)==1)) == 1
       u_centerline(i) = u_g(find(u_normalized(:,i)==1),i);
       v_centreline(i) = v_g(find(u_normalized(:,i)==1),i);
   else
       u_centerline(i) = 0;
       v_centerline(i) = 0;
   end        
end

%hold off
%u_normalized = max(u_normalized,0);
%% Normalizing

%find half width
half_width = zeros(rawI,1);
peak_location = zeros(rawI,1);
y_normalized = zeros(rawJ, rawI);
y_normalized_1 = zeros(rawJ, rawI);

% 
% figure
% hold on
 for i = 1:rawI
   [pks,locs,widths,proms] =  findpeaks(u_normalized(:,i),y_centered(:,i),'MinPeakProminence',.99,'Annotate','extents','WidthReference','halfheight');
   if isempty(widths)
       half_width(i) = 0;
       peak_location(i) = 0;
   else
   half_width(i) = .5*widths;

   peak_location(i) = locs;
   end
   y_normalized(:,i) = y_centered(:,i)/(half_width(i));

   %Calculate 1 SD
   sd = half_width(i)*2/(2*sqrt(2*log(2)));
   jet_indicies = round((rawJ-u_secondary_index(i))/2) : u_secondary_index(i);

 end 

 %%

x_locations_dx = linspace(12.5, 19, 5);
x_locations = round((x_locations_dx*jet_diameter)/scale);
markers = ["o" "s" "^" ">" "<" "*"];
colormap = ["#003f5c"; "#665191";  "#d45087"; "#f95d6a"; "#ff7c43"; "#ffa600"];

figure(1)
hold on
fplot(@(x) exp(-log(2)*x^2),'--k',[-2.5 2.5]);
fplot(@(x) (1+x^2/2)^(-2),'-b',[-1.5 1.5]);
for i = 1:6
    plot(y_normalized(:,x_locations(i)), u_normalized(:,x_locations(i)),markers(i),'MarkerSize',4,'MarkerEdgeColor',colormap(i))
end
hold off
legend ('Gaussian Shape Function', 'Schlicting Jet eqn', 'x/d = 10', 'x/d = 12.5', 'x/d = 15', 'x/d = 17.5', 'x/d = 20', 'x/d = 22.5')
xlabel('$\eta = r/b$','interpreter','latex')
ylabel('$U/U_{m}$','interpreter','latex')
axis([-2.5 2.5 0 1])
yticks(0:.25:1);
box on

centerline_data= [transpose(x_g(1,:)) u_centerline];
writematrix(centerline_data,'centerline_data_2_4_confined.xlsx')