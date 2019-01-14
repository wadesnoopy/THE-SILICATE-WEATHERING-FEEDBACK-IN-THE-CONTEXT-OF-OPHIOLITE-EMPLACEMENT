%----------------------------------------------------------------------
% The code here is Model 1 for Os system in the appendix with temperature effect
% integrated

% With Monte Carlo error analysis


% !!! Be aware that some normalization will generate negative number, use
% A(A<0) = 0 to transform or just treat it as fail. I treat it as fail.

clear all
close all

fprintf('\n\n')

cpu_start = cputime;

%----------------------------------------------------------------------

% Load the equation solving

load('equi_simple.mat');

% % Normalize ratio
% ratio_ht = 0.26;
% ratio_lt = 0.11;
% 
% ratio_ht = ratio_ht / (7.4 + ratio_ht);
% ratio_lt = ratio_lt / (7.4 + ratio_lt);
% ratio_cosm = ratio_cosm / (7.4 + ratio_cosm);
% ratio_dust = ratio_dust / (7.4 + ratio_dust);
% ratio_river = ratio_river / (7.4 + ratio_river);

% ------------------------------------
% Print finishing loading the equations
% ------------------------------------

fprintf('Finish the equation loading\n')

e = cputime-cpu_start;
fprintf('%f mins passed...\n\n', e /60)


OS_MOL = 190.2;

% Control the seeding the same !!!!!!!!!!!!!!!!!!
% rng(1);

% This is the time range of each segment
t_interval=1*10^4;

% Set monte carlo number
num_monte = 1000;

%----------------------------------------------------------------------

% Load data_record from Li (already smooth and interpolate them)
load('data_record_smooth.mat')

% ------------------------------------
% Print finishing loading the data_record
% ------------------------------------

fprintf('Finish the data_record loading\n')

e = cputime-cpu_start;
fprintf('%f mins passed...\n\n', e /60)

% Total points
num_points = floor((AGE_OLD - AGE_YOUNG) * 10^6 / t_interval) + 1;

% Construct vectors to store parameters (Some are zeros after bounding)
K_total = ones(1,num_monte) * -1;


% Parameters from Jagoutz

Jagoutz_age = xlsread('Jagoutz_results', '250', 'A2:A1002');
Jagoutz_area = xlsread('Jagoutz_results', '250', 'B2:B1002');

% ------------------------- Interpolating the data record


Jagoutz_age = (AGE_OLD - Jagoutz_age) * 10^6;
% No need to smooth since already smoothed
Area_interp = @(t)interp1(Jagoutz_age, Jagoutz_area, t, 'PCHIP');

% ------------------------- All parameters for late pleistocene from Li

% From Li 2013
% For figure 1
e_dOs_ocean = 0.05;

% ------------------------- All parameters for ophiolite from Jagoutz (all
% mean value)

mafic = 0.49;
Os_conc_mafic = 0.021 * 10^-9;
Os_conc_ultra = 3.60 * 10^-9;
dOs_mafic = 0.143;
dOs_ultra = 0.129;
sw_slope = 0.0553;
sw_factor = 18.41;
co2_slope = 0.0642;
co2_factor = 323.44;
sw_runoff_SE = 1372; % mm/year
sw_T_SE = 25;

% add temperature
bottom_T_ini = T_interp(0);
sw_T = @(t)T_interp(t) + sw_T_SE - bottom_T_ini;

% add runoff
RT = 0.037; % from GEOCARB
sw_runoff = @(t)(1 + (T_interp(t) - T_interp(0)) * RT) * sw_runoff_SE;

% 1 sigma
e_mafic = 0;
e_Os_conc_mafic = Os_conc_mafic * 0.05;
e_Os_conc_ultra = Os_conc_ultra * 0.05;
e_dOs_mafic = 0.02;
e_dOs_ultra = 0.005;


e_sw_slope = 0;
e_sw_factor = 0;
e_co2_slope = 0;
e_co2_factor = 0;
e_sw_runoff_SE = 0;
e_sw_T_SE = 0;

% -------------------------From Jagoutz
mafic_monte = normrnd(mafic,e_mafic,1,num_monte);
ultra_monte = 1 - mafic_monte;
Os_conc_mafic_monte = normrnd(Os_conc_mafic,e_Os_conc_mafic,1,num_monte);
Os_conc_ultra_monte = normrnd(Os_conc_ultra,e_Os_conc_ultra,1,num_monte);
dOs_mafic_monte = normrnd(dOs_mafic/(7.4+dOs_mafic),7.4*e_dOs_mafic/(7.4+dOs_mafic)^2,1,num_monte);
dOs_ultra_monte = normrnd(dOs_ultra/(7.4+dOs_ultra),7.4*e_dOs_ultra/(7.4+dOs_ultra)^2,1,num_monte);


% Silicate weathering for SE Asia/Indonesia

% silicate dissolving flux
sw_slope_monte = normrnd(sw_slope,e_sw_slope,1,num_monte);
sw_factor_monte = normrnd(sw_factor,e_sw_factor,1,num_monte);

% co2 consumption flux (Note HCO3- flux is the same)
co2_slope_monte = normrnd(co2_slope,e_co2_slope,1,num_monte);
co2_factor_monte = normrnd(co2_factor,e_co2_factor,1,num_monte);
sw_runoff_SE_monte = normrnd(sw_runoff_SE,e_sw_runoff_SE,1,num_monte);
sw_T_SE_monte = normrnd(sw_T_SE,e_sw_T_SE,1,num_monte);





% Constitute flux

f_sw_monte = sw_factor_monte / 1000; % t/km^2/yr (need to divided by 1000) the original paper is wrong.

f_Os_mafic_monte = @(t)Area_interp(t) .* mafic_monte .* f_sw_monte * sw_runoff(t) .* exp(sw_slope_monte .* sw_T(t)) * 10^6 .* Os_conc_mafic_monte / OS_MOL;

f_Os_ultra_monte = @(t)Area_interp(t) .* ultra_monte .* f_sw_monte * sw_runoff(t) .* exp(sw_slope_monte .* sw_T(t)) * 10^6 .* Os_conc_ultra_monte / OS_MOL;

f_A_Jagoutz_monte = @(t)Area_interp(t) .* sw_runoff(t) .* co2_factor_monte .* exp(co2_slope_monte .* sw_T(t));

% Calculate Jagoutz Os flux
Os_jagoutz = zeros(num_points,num_monte);
A_jagoutz = zeros(num_points,num_monte);

for j=0:num_points-1
    t = t_interval * j;
    Os_jagoutz(j+1,:) = f_Os_mafic_monte(t) + f_Os_ultra_monte(t);
    A_jagoutz(j+1, :) = f_A_Jagoutz_monte(t);
end
   
Os_jagoutz_origin = Os_jagoutz;
Os_jagoutz(Os_jagoutz<0) = NaN;
Os_jagoutz_ave = nanmean(Os_jagoutz, 2);
Os_jagoutz_std = nanstd(Os_jagoutz, 0, 2);
Os_jagoutz_low = Os_jagoutz_ave - Os_jagoutz_std;
Os_jagoutz_high = Os_jagoutz_ave + Os_jagoutz_std;

A_jagoutz_origin = A_jagoutz;
A_jagoutz(A_jagoutz<0) = NaN;
A_jagoutz_ave = nanmean(A_jagoutz, 2);
A_jagoutz_std = nanstd(A_jagoutz, 0, 2);


% ----------------------------------------------- 55 Ma initialization
% started

% Calculate initial value at 55 Ma (steady state for river flux change and
% no change

t = 0;

dOs_ocean_now = dOs_ocean_interp(t);

dOs_ocean_now_monte = normrnd(dOs_ocean_now/(7.4+dOs_ocean_now),7.4*e_dOs_ocean/(7.4+dOs_ocean_now)^2,1,num_monte);

dOs_ocean_old = dOs_ocean_interp(t - t_interval);

dOs_ocean_old_monte = normrnd(dOs_ocean_old/(7.4+dOs_ocean_old),7.4*e_dOs_ocean/(7.4+dOs_ocean_old)^2,1,num_monte);

der_dOs_ocean_monte = (dOs_ocean_now_monte - dOs_ocean_old_monte) / t_interval;

der_dOs_ocean_monte = der_dOs_ocean_monte * 0;

Os_mafic_monte_now = f_Os_mafic_monte(t);
Os_ultra_monte_now = f_Os_ultra_monte(t);

for i=1:num_monte
    if (Os_mafic_monte_now(i) >= 0) && (Os_ultra_monte_now(i) >= 0)
        
        K_total(i) = F_K_total_norm(Os_mafic_monte_now(i), Os_ultra_monte_now(i), dOs_mafic_monte(i), dOs_ocean_now_monte(i), dOs_ultra_monte(i), der_dOs_ocean_monte(i));
        
        if (K_total(i) < 0)
            K_total(i) = -2;
        end
        
    end
    
end

dOs_ocean_old_monte = dOs_ocean_now_monte;

% ------------------------------------
% Print finishing the initialization
% ------------------------------------

fprintf('Finish the initialization (55 Ma values)\n')

e = cputime-cpu_start;
fprintf('%f mins passed...\n', e /60)

% ----------------------------------------------- 55 Ma initialization
% finished

% ----------------------------------------------- river flux no change
% (calculate Os ratio of the ocean) And for this value, no monte
% K_steady = K_total(1,:);
% K_steady(K_steady<0) = NaN;
% K_steady_ave = nanmean(K_steady);
% K_steady_std = nanstd(K_steady);
% K_steady_monte = normrnd(K_steady_ave,K_steady_std,1,num_monte);

dOs_total = ones(num_points, num_monte) * -1;

Os_mafic_monte_part =  mafic_monte .* f_sw_monte .* 10^6 .* Os_conc_mafic_monte / OS_MOL;
Os_ultra_monte_part =  ultra_monte .* f_sw_monte .* 10^6 .* Os_conc_ultra_monte / OS_MOL;

t_change = 0:t_interval:(AGE_OLD - AGE_YOUNG) * 10^6;

for i=1:num_monte
    if  (K_total(i) >= 0) && (Os_mafic_monte_part(i) >= 0) && (Os_ultra_monte_part(i) >= 0)
        flux_times_ratio = @(t,y)flux_ht * (ratio_ht_norm - y) + flux_lt * (ratio_lt_norm - y) + ...
            flux_cosm * (ratio_cosm_norm - y) + flux_dust * (ratio_dust_norm - y) + ...
            K_total(i) * flux_river * (ratio_river_norm - y)+ Area_interp(t) * Os_mafic_monte_part(i) * sw_runoff(t) * exp(sw_slope_monte(i) * sw_T(t)) *...
            (dOs_mafic_monte(i) - y) + Area_interp(t) * Os_ultra_monte_part(i) * sw_runoff(t) * exp(sw_slope_monte(i) * sw_T(t)) * (dOs_ultra_monte(i) - y);
%         flux_times_ratio = @(t,y)flux_ht * (ratio_ht - y) + flux_lt * (ratio_lt - y) + ...
%             flux_cosm * (ratio_cosm - y) + flux_dust * (ratio_dust - y) + ...
%             K_steady_monte(i) * flux_river * (ratio_river - y)+ Area_interp(0) * Os_mafic_monte_part(i) *...
%             (dOs_mafic_monte(i) - y) + Area_interp(0) * Os_ultra_monte_part(i) * (dOs_ultra_monte(i) - y);
        
        slope_ratio = @(t,y)flux_times_ratio(t,y) / Os_ocean;
        [~,Os_ratio] = ode45(slope_ratio, t_change, dOs_ocean_old_monte(i));
        dOs_total(:,i) = Os_ratio;
        fprintf('%d steps finished\n', i)
    end
end

dOs_total_origin = dOs_total;
dOs_total = 7.4 * dOs_total ./ (1 - dOs_total);
dOs_total(dOs_total<0) = NaN;
dOs_ave = nanmean(dOs_total, 2);
dOs_std = nanstd(dOs_total, 0, 2);
dOs_total_low = dOs_ave - dOs_std;
dOs_total_high = dOs_ave + dOs_std;


% Create file directory

file_dir = strcat('norm_simple_', num2str(ratio_river), '_nomatch_high_rate_change_T_runoff', num2str(t_interval / 1000), 'ky_', num2str(num_monte), 'monte');

if not(exist(file_dir, 'dir'))
    mkdir(file_dir)
end

% When concatenate dir (difference in pc and mac/unix)
if ispc
    book_name = strcat(file_dir, '\', 'output', '.xlsx');
    csv_name = strcat(file_dir, '\', 'output', '.csv');
    dOs_name = strcat(file_dir, '\', file_dir, '_dOs.jpg');
    Os_jagoutz_name = strcat(file_dir, '\', file_dir, '_s_jagoutz.jpg');
else
    book_name = strcat(file_dir, '/', 'output', '.xlsx');
    csv_name = strcat(file_dir, '/', 'output', '.csv');
    Os_jagoutz_name = strcat(file_dir, '/', file_dir, '_s_jagoutz.jpg');
end


% Begin to write data into file

age_all_plot_book = (AGE_OLD:-t_interval/10^6:AGE_YOUNG)'; % Must transform from 1 row to 1 col

% Title
headers = {'Age (Ma)', 'Os_ratio', 'Os_ratio_std', 'Os_total', 'Os_total_std', 'A_total', 'A_total_std'};

% Data
whole_data = [age_all_plot_book, dOs_ave, dOs_std, Os_jagoutz_ave, Os_jagoutz_std, A_jagoutz_ave, A_jagoutz_std];


% Use the online verstion to write to csv
csvwrite_with_headers(csv_name, whole_data, headers)

e = cputime-cpu_start;
fprintf('The whole program used %f mins\n', e /60)

% Plot all

% Note, don't reverse age_all_plot here, since when you use fill, it will
% be wrong.

age_all_plot = AGE_OLD:-t_interval/10^6:AGE_YOUNG;

label_fontsize = 14;
axis_fontsize = 12;

% Ocean Os ratio
figure1 = figure;
plot(age_all_plot, dOs_ave, 'r');
hold on

X=[age_all_plot,fliplr(age_all_plot)];
Y_dOs=[dOs_total_low',fliplr(dOs_total_high')];
fill(X,Y_dOs,'r', 'faceAlpha', 0.3, 'linestyle', 'none');

xlabel('Age (Ma)', 'FontSize', label_fontsize);
ylabel('dOs\_ocean', 'FontSize', label_fontsize);
set(gca,'fontsize',axis_fontsize);

% Save figure
saveas(figure1, dOs_name)

% Jagoutz Os flux

figure2 = figure;
plot(age_all_plot, Os_jagoutz_ave, 'r');
hold on

X=[age_all_plot,fliplr(age_all_plot)];
Y_Os_jagoutz=[Os_jagoutz_low',fliplr(Os_jagoutz_high')];
fill(X,Y_Os_jagoutz,'r', 'faceAlpha', 0.3, 'linestyle', 'none');

xlabel('Age (Ma)', 'FontSize', label_fontsize);
ylabel('Os\_jagoutz', 'FontSize', label_fontsize);
set(gca,'fontsize',axis_fontsize);

% Save figure
saveas(figure2, Os_jagoutz_name)

% saveas(gcf, dOs_name)
