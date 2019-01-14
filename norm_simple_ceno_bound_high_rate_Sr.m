%----------------------------------------------------------------------
% The code here is Model 2 for Sr system in the appendix


% With Monte Carlo error analysis

% !!! Be aware that some normalization will generate negative number, use
% A(A<0) = 0 to transform or just treat it as fail. I treat it as fail

clear all
close all

fprintf('\n\n')

cpu_start = cputime;

%----------------------------------------------------------------------

% Load the equation solving

load('equi_simple.mat');


% ------------------------------------
% Print finishing loading the equations
% ------------------------------------

fprintf('Finish the equation loading\n')

e = cputime-cpu_start;
fprintf('%f mins passed...\n\n', e /60)


SR_MOL = 87.62;

% Control the seeding the same !!!!!!!!!!!!!!!!!!
% rng(1);

% This is the time range of each segment
t_interval=5*10^5;

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
K_total = ones(num_points,num_monte) * -1;


% Parameters from Jagoutz

Jagoutz_age = xlsread('Jagoutz_results', '250', 'A2:A1002');
Jagoutz_area = xlsread('Jagoutz_results', '250', 'B2:B1002');

% ------------------------- Interpolating the data record


Jagoutz_age = (AGE_OLD - Jagoutz_age) * 10^6;
% No need to smooth since already smoothed
Area_interp = @(t)interp1(Jagoutz_age, Jagoutz_area, t, 'PCHIP');

% Area_interp = @(t)0;

% ------------------------- All parameters for late pleistocene from Li



% From Li 2013
% For figure 1
e_dSr_ocean = 0.00002;

% ------------------------- All parameters for ophiolite from Jagoutz (all
% mean value)

mafic = 0.49;
Sr_conc_mafic = 113 * 10^-6;
Sr_conc_ultra = 38 * 10^-6;
dSr_mafic = 0.7029;
dSr_ultra = 0.7056;
sw_slope = 0.0553;
sw_factor = 18.41;
co2_slope = 0.0642;
co2_factor = 323.44;
sw_runoff_SE = 1372; % mm/year
sw_T_SE = 25;

% sw_runoff_SE = 129; % mm/year This is for Ethiopia
% sw_T_SE = 21.3;

% 1 sigma
e_mafic = 0;
e_Sr_conc_mafic = Sr_conc_mafic * 0.24;
e_Sr_conc_ultra = 6.56 * 10^-6;
e_dSr_mafic = 0.0004;
e_dSr_ultra = 0.0008;
e_sw_slope = 0;
e_sw_factor = 0;
e_co2_slope = 0;
e_co2_factor = 0;
e_sw_runoff_SE = 0;
e_sw_T_SE = 0;




% -------------------------From Jagoutz
mafic_monte = normrnd(mafic,e_mafic,1,num_monte);
ultra_monte = 1 - mafic_monte;
Sr_conc_mafic_monte = normrnd(Sr_conc_mafic,e_Sr_conc_mafic,1,num_monte);
Sr_conc_ultra_monte = normrnd(Sr_conc_ultra,e_Sr_conc_ultra,1,num_monte);
dSr_mafic_monte = normrnd(dSr_mafic/(9.43+dSr_mafic),9.43*e_dSr_mafic/(9.43+dSr_mafic)^2,1,num_monte);
dSr_ultra_monte = normrnd(dSr_ultra/(9.43+dSr_ultra),9.43*e_dSr_ultra/(9.43+dSr_ultra)^2,1,num_monte);


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

f_sw_monte = sw_runoff_SE_monte .* sw_factor_monte .* exp(sw_slope_monte .* sw_T_SE_monte) / 1000; % t/km^2/yr

f_Sr_mafic_monte = @(t)Area_interp(t) .* mafic_monte .* f_sw_monte * 10^6 .* Sr_conc_mafic_monte / SR_MOL;

f_Sr_ultra_monte = @(t)Area_interp(t) .* ultra_monte .* f_sw_monte  * 10^6 .* Sr_conc_ultra_monte / SR_MOL;

% ----------------------------------------------- 55 Ma initialization
% started

% Calculate initial value at 55 Ma (steady state for river flux change and
% no change

t = 0;
% Solve Sr system

dSr_ocean_now = dSr_ocean_interp(t);
dSr_ocean_now_monte = normrnd(dSr_ocean_now/(9.43+dSr_ocean_now),9.43*e_dSr_ocean/(9.43+dSr_ocean_now)^2,1,num_monte);

dSr_ocean_old = dSr_ocean_interp(t - t_interval);
dSr_ocean_old_monte = normrnd(dSr_ocean_old/(9.43+dSr_ocean_old),9.43*e_dSr_ocean/(9.43+dSr_ocean_old)^2,1,num_monte);

der_dSr_ocean_monte = (dSr_ocean_now_monte - dSr_ocean_old_monte) / t_interval;

% der_dSr_ocean_monte = 0 * der_dSr_ocean_monte;


Sr_mafic_monte_now = f_Sr_mafic_monte(t);
Sr_ultra_monte_now = f_Sr_ultra_monte(t);


success_array = zeros(num_points, 1);
flag = 0;

for i=1:num_monte
    if (Sr_mafic_monte_now(i) >= 0) && (Sr_ultra_monte_now(i) >= 0)
        
        K_total(1,i) = F_K_total_Sr_norm(Sr_mafic_monte_now(i), Sr_ultra_monte_now(i), dSr_mafic_monte(i), dSr_ocean_now_monte(i), dSr_ultra_monte(i), der_dSr_ocean_monte(i));
        
        if (K_total(1,i) < 0)
            K_total(1,i) = -2;
        end
        
    end
    
end

success_array(1) = flag;

dSr_ocean_old_monte = dSr_ocean_now_monte;



% ------------------------------------
% Print finishing the initialization
% ------------------------------------

fprintf('Finish the initialization (55 Ma values)\n')

e = cputime-cpu_start;
fprintf('%f mins passed...\n', e /60)

% ----------------------------------------------- 55 Ma initialization
% finished


% ----------------------------------------------- river flux changes to
% match the record

for j=1:num_points-1
    t = t_interval * j;
    
    dSr_ocean_now = dSr_ocean_interp(t);
    dSr_ocean_now_monte = normrnd(dSr_ocean_now/(9.43+dSr_ocean_now),9.43*e_dSr_ocean/(9.43+dSr_ocean_now)^2,1,num_monte);
    der_dSr_ocean_monte = (dSr_ocean_now_monte - dSr_ocean_old_monte) / t_interval;
    
%     der_dOs_ocean_monte = der_dOs_ocean_monte * 0;
    
    Sr_mafic_monte_now = f_Sr_mafic_monte(t);
    Sr_ultra_monte_now = f_Sr_ultra_monte(t);
    
    
    flag = 0;
    
    for i=1:num_monte
        if (Sr_mafic_monte_now(i) >= 0) && (Sr_ultra_monte_now(i) >= 0)
        
            K_total(j+1, i) = F_K_total_Sr_norm(Sr_mafic_monte_now(i), Sr_ultra_monte_now(i), dSr_mafic_monte(i), dSr_ocean_now_monte(i), dSr_ultra_monte(i), der_dSr_ocean_monte(i));
        
            if (K_total(j+1, i) < 0)
                K_total(j+1, i) = -2;
            end
        
        end
    
    end
    
    success_array(j+1) = flag;

    dSr_ocean_old_monte = dSr_ocean_now_monte;
    % ------------------------------------
    % Print finishing t one time step
    % ------------------------------------
    
    fprintf('Finish by %f Ma\n', (AGE_OLD * 10^6 - t) / 10^6)
    
    e = cputime-cpu_start;
    fprintf('%f mins passed...\n', e /60)
end


% Transform to average with 95% confidence interval and row wise

% Transform all -1 to NaN. Then use nanmean and nanstd to calculate the
% rest.


K_total_origin = K_total;
K_total(K_total<0) = NaN;
K_ave = nanmean(K_total, 2);
K_std = nanstd(K_total, 0, 2);
K_low = K_ave - K_std;
K_high = K_ave + K_std;

Sr_total_ave = K_ave * flux_river_Sr;
Sr_total_std = K_std * flux_river_Sr;
Sr_total_low = Sr_total_ave - Sr_total_std;
Sr_total_high = Sr_total_ave + Sr_total_std;

% Create file directory

file_dir = strcat('norm_simple_', num2str(ratio_river), '_high_rate_Sr', num2str(t_interval / 1000), 'ky_', num2str(num_monte), 'monte');

if not(exist(file_dir, 'dir'))
    mkdir(file_dir)
end

% When concatenate dir (difference in pc and mac/unix)
if ispc
    book_name = strcat(file_dir, '\', 'output', '.xlsx');
    csv_name = strcat(file_dir, '\', 'output', '.csv');
    K_name = strcat(file_dir, '\', 'output', '_K.jpg');
    Os_name = strcat(file_dir, '\', 'output', '_Os.jpg');
    suc_name = strcat(file_dir, '\', 'output', '_success.jpg');
else
    book_name = strcat(file_dir, '/', 'output', '.xlsx');
    K_name = strcat(file_dir, '/', 'output', '_K.jpg');
    Os_name = strcat(file_dir, '/', 'output', '_Os.jpg');
    suc_name = strcat(file_dir, '/', 'output', '_success.jpg');
end


% Begin to write data into file

success_array = success_array / num_monte; % Transform it into percentage
age_all_plot_book = (AGE_OLD:-t_interval/10^6:AGE_YOUNG)'; % Must transform from 1 row to 1 col

% Title
headers = {'Age (Ma)', 'success percent', 'total K', 'total_K_std', 'total flux', 'total_flux_std'};

% Data
whole_data = [age_all_plot_book, success_array, K_ave, K_std, Sr_total_ave, Sr_total_std];


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

% Continental silicate weathering
figure1 = figure;
plot(age_all_plot, K_ave, 'r');
hold on

X=[age_all_plot,fliplr(age_all_plot)];
Y_K_total=[K_low',fliplr(K_high')];
fill(X,Y_K_total,'r', 'faceAlpha', 0.3, 'linestyle', 'none');

xlabel('Age (Ma)', 'FontSize', label_fontsize);
ylabel('Total flux', 'FontSize', label_fontsize);
set(gca,'fontsize',axis_fontsize);

% Save figure
saveas(figure1, K_name)

% Os total flux (not including Jagoutz)
figure2 = figure;
plot(age_all_plot, Sr_total_ave, 'r');
hold on


X=[age_all_plot,fliplr(age_all_plot)];
Y_Sr_total=[Sr_total_low',fliplr(Sr_total_high')];
fill(X,Y_Sr_total,'r', 'faceAlpha', 0.3, 'linestyle', 'none');

xlabel('Age (Ma)', 'FontSize', label_fontsize);
ylabel('Total flux', 'FontSize', label_fontsize);
set(gca,'fontsize',axis_fontsize);

% Save figure
saveas(figure2, Os_name)


% Success rate
figure3 = figure;
plot(age_all_plot, success_array, 'r');

% Save figure
saveas(figure3, suc_name)


ttt = 0:5e5:2e7;

xxx = dSr_ocean_interp(ttt);

plot(ttt, xxx)