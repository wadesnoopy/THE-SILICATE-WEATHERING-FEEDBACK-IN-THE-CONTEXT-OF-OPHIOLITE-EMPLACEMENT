% this is to smooth the data records

clear all
close all

cpu_start_data_record = cputime;

% Age interval (unit: Ma)

AGE_OLD = 55;
AGE_YOUNG = 30;

% -------------------------Data record


age_Sp = xlsread('shuang2016','Sp','A2:A13');
Sp_record = xlsread('shuang2016','Sp','B2:B13');

age_Sp_berner = xlsread('shuang2016','Sp_berner','A2:A59');
Sp_record_berner = xlsread('shuang2016','Sp_berner','B2:B59');


age_dC_carbb = xlsread('shuang2016','dC_carbb','A2:A241');
dC_carbb_record = xlsread('shuang2016','dC_carbb','B2:B241');

age_dC_orgb = xlsread('shuang2016','dC_orgb','A2:A115');
dC_orgb_record = xlsread('shuang2016','dC_orgb','B2:B115');

age_dSr_ocean = xlsread('shuang2016','dSr_ocean','A2:A305');
dSr_ocean_record = xlsread('shuang2016','dSr_ocean','B2:B305');

% % Old without burton but with ravizza
% age_dOs_ocean = xlsread('shuang2016','dOs_ocean','A2:A119');
% dOs_ocean_record = xlsread('shuang2016','dOs_ocean','B2:B119');

% New with burton but without ravizza
age_dOs_ocean = xlsread('shuang2016','dOs_ocean','A2:A146');
dOs_ocean_record = xlsread('shuang2016','dOs_ocean','B2:B146');

age_co2 = xlsread('shuang2016','CO2','A2:A332');
co2_record = xlsread('shuang2016','CO2','B2:B332');
co2_record_yerr_low = xlsread('shuang2016','CO2','C2:C332');
co2_record_yerr_high = xlsread('shuang2016','CO2','D2:D332');

age_T = xlsread('shuang2016','T','A2:A990');
T_record = xlsread('shuang2016','T','D2:D990');
T_record_yerr_low = T_record - xlsread('shuang2016','T','C2:C990');
T_record_yerr_high = xlsread('shuang2016','T','B2:B990') - T_record;

% Read caves et al. 2016

age_Mc = xlsread('shuang2016','scenario1','A2:A67');
Mc_record = xlsread('shuang2016','scenario1','B2:B118');
Mc_sd_record = xlsread('shuang2016','scenario1','C2:C118');

age_Ma = xlsread('shuang2016','scenario1','A2:A67');
Ma_record = xlsread('shuang2016','scenario1','D2:D118');
Ma_sd_record = xlsread('shuang2016','scenario1','E2:E118');

% Read Li et al. (2013) model output

age_Li_silw = xlsread('Li_model_output_origin', 'K_sili', 'A2:A73');
Li_silw = xlsread('Li_model_output_origin', 'K_sili', 'B2:B73');

age_Li_basw = xlsread('Li_model_output_origin', 'K_bas', 'A2:A77');
Li_basw = xlsread('Li_model_output_origin', 'K_bas', 'B2:B77');

age_Li_Mg = xlsread('Li_model_output_origin', 'Mg', 'A2:A64');
Li_Mg = xlsread('Li_model_output_origin', 'Mg', 'B2:B64');

age_Li_sedi = xlsread('Li_model_output_origin', 'K_sedi', 'A2:A62');
Li_sedi = xlsread('Li_model_output_origin', 'K_sedi', 'B2:B62');

% ------------------------- Interpolating the data record

% For spreading rate of seafloor
age_Sp = (AGE_OLD - age_Sp) * 10^6;
% Sp_record_smooth =  smooth(age_Sp, Sp_record,0.1,'loess');
Sp_interp = @(t)interp1(age_Sp, Sp_record, t, 'PCHIP');

age_Sp_berner = (AGE_OLD - age_Sp_berner) * 10^6;

% Following Cave's suggestion
% Sp_record_berner_smooth =  smooth(age_Sp_berner, Sp_record_berner,0.1,'loess');
Sp_berner_interp = @(t)interp1(age_Sp_berner, Sp_record_berner, t, 'PCHIP');

% For carbon isotope in carbonates
age_dC_carbb = (AGE_OLD - age_dC_carbb) * 10^6;
dC_carbb_record_smooth =  smooth(age_dC_carbb, dC_carbb_record,0.1,'loess');
dC_carbb_interp = @(t)interp1(age_dC_carbb, dC_carbb_record_smooth, t, 'PCHIP');

% For carbon isotope in organic carbon
age_dC_orgb = (AGE_OLD - age_dC_orgb) * 10^6;
dC_orgb_record_smooth =  smooth(age_dC_orgb, dC_orgb_record,0.1,'loess');
dC_orgb_interp = @(t)interp1(age_dC_orgb, dC_orgb_record_smooth, t, 'PCHIP');

% For Sr isotope in seawater
age_dSr_ocean = (AGE_OLD - age_dSr_ocean) * 10^6;
dSr_ocean_record_smooth =  smooth(age_dSr_ocean, dSr_ocean_record,0.05,'loess');
dSr_ocean_interp = @(t)interp1(age_dSr_ocean, dSr_ocean_record_smooth, t, 'PCHIP');

% For Os isotope in seawater
age_dOs_ocean = (AGE_OLD - age_dOs_ocean) * 10^6;
dOs_ocean_record_smooth =  smooth(age_dOs_ocean, dOs_ocean_record,0.1,'loess');
dOs_ocean_interp = @(t)interp1(age_dOs_ocean, dOs_ocean_record_smooth, t, 'PCHIP');

% For CO2 record
age_co2 = (AGE_OLD - age_co2) * 10^6;
co2_record_smooth =  smooth(age_co2, co2_record,0.1,'loess');
co2_interp = @(t)interp1(age_co2, co2_record_smooth, t, 'PCHIP');

co2_record_yerr_low_smooth =  smooth(age_co2, co2_record_yerr_low,0.1,'loess');
co2_yerr_low_interp = @(t)interp1(age_co2, co2_record_yerr_low_smooth, t, 'PCHIP');

co2_record_yerr_high_smooth =  smooth(age_co2, co2_record_yerr_high,0.1,'loess');
co2_yerr_high_interp = @(t)interp1(age_co2, co2_record_yerr_high_smooth, t, 'PCHIP');

% For T record

age_T = (AGE_OLD - age_T) * 10^6;
T_record_smooth =  smooth(age_T, T_record,0.1,'loess');
T_interp = @(t)interp1(age_T, T_record_smooth, t, 'PCHIP');

T_record_yerr_low_smooth =  smooth(age_T, T_record_yerr_low,0.1,'loess');
T_yerr_low_interp = @(t)interp1(age_T, T_record_yerr_low_smooth, t, 'PCHIP');

T_record_yerr_high_smooth =  smooth(age_T, T_record_yerr_high,0.1,'loess');
T_yerr_high_interp = @(t)interp1(age_T, T_record_yerr_high_smooth, t, 'PCHIP');

% % For Caves

% For Mc isotope in seawater
age_Mc = (AGE_OLD - age_Mc) * 10^6;
% No smooth
% Mc_record_smooth =  smooth(age_Mc, Mc_record,0.1,'loess');
Mc_record_smooth = Mc_record;
Mc_interp = @(t)interp1(age_Mc, Mc_record_smooth, t, 'PCHIP');


% No smooth

% Mc_sd_record_smooth =  smooth(age_Mc, Mc_sd_record,0.1,'loess'); 
Mc_sd_record_smooth =  Mc_sd_record;
Mc_sd_interp = @(t)interp1(age_Mc, Mc_sd_record_smooth, t, 'PCHIP');

age_Ma = (AGE_OLD - age_Ma) * 10^6;

% No smooth
% Ma_record_smooth =  smooth(age_Ma, Ma_record,0.1,'loess');
Ma_record_smooth =  Ma_record;
Ma_interp = @(t)interp1(age_Ma, Ma_record_smooth, t, 'PCHIP');

% No smooth
% Ma_sd_record_smooth =  smooth(age_Ma, Ma_sd_record,0.1,'loess');
Ma_sd_record_smooth = Ma_sd_record;
Ma_sd_interp = @(t)interp1(age_Ma, Ma_sd_record_smooth, t, 'PCHIP');



% For silicate weathering
age_Li_silw = (AGE_OLD - age_Li_silw) * 10^6;
Li_silw_interp = @(t)interp1(age_Li_silw, Li_silw, t, 'PCHIP');

% For basalt weathering
age_Li_basw = (AGE_OLD - age_Li_basw) * 10^6;
Li_basw_interp = @(t)interp1(age_Li_basw, Li_basw, t, 'PCHIP');

% For Mg concentration
age_Li_Mg = (AGE_OLD - age_Li_Mg) * 10^6;
Li_Mg_interp = @(t)interp1(age_Li_Mg, Li_Mg, t, 'PCHIP');

% For sediment weathering
age_Li_sedi = (AGE_OLD - age_Li_sedi) * 10^6;
Li_sedi_interp = @(t)interp1(age_Li_sedi, Li_sedi, t, 'PCHIP');


age_all = linspace(-5*10^6, 25*10^6, 301)'; % This is 60 to 30

% age_all = linspace(-5*10^6, 25*10^6, 31)'; % This is to create data for
% caves Rk

Sp_interp_all = Sp_interp(age_all);
Sp_berner_interp_all = Sp_berner_interp(age_all);
dC_carbb_interp_all = dC_carbb_interp(age_all);
dC_orgb_interp_all = dC_orgb_interp(age_all);
dSr_ocean_interp_all = dSr_ocean_interp(age_all);
dOs_ocean_interp_all = dOs_ocean_interp(age_all);
co2_interp_all = co2_interp(age_all);
co2_yerr_low_interp_all = co2_yerr_low_interp(age_all);
co2_yerr_high_interp_all = co2_yerr_high_interp(age_all);
T_interp_all = T_interp(age_all);
T_yerr_low_interp_all = T_yerr_low_interp(age_all);
T_yerr_high_interp_all = T_yerr_high_interp(age_all);

Mc_interp_all = Mc_interp(age_all);
Mc_sd_interp_all = Mc_sd_interp(age_all);

Ma_interp_all = Ma_interp(age_all);
Ma_sd_interp_all = Ma_sd_interp(age_all);

Li_silw_interp_all = Li_silw_interp(age_all);
Li_basw_interp_all = Li_basw_interp(age_all);
Li_Mg_interp_all = Li_Mg_interp(age_all);
Li_sedi_interp_all = Li_sedi_interp(age_all);


% Begin to write data into file

age_all_book = AGE_OLD - age_all / 10^6; % Must transform from 1 row to 1 col

% Title
headers = {'Age (Ma)', 'Sp', 'Sp_berner', 'dC_carbb', 'dC_orgb', 'dSr_ocean', 'dOs_ocean', 'co2', 'co2_yerr_low', 'co2_yerr_high', 'T', 'T_yerr_low', 'T_yerr_high'};

% Data
whole_data = [age_all_book, Sp_interp_all, Sp_berner_interp_all, dC_carbb_interp_all, dC_orgb_interp_all, dSr_ocean_interp_all, dOs_ocean_interp_all,...
    co2_interp_all, co2_yerr_low_interp_all, co2_yerr_high_interp_all, T_interp_all, T_yerr_low_interp_all, T_yerr_high_interp_all];

csv_name = 'data_record_smooth.csv';
% Use the online verstion to write to csv
csvwrite_with_headers(csv_name, whole_data, headers)


filename = 'data_record_smooth_60.mat';
save(filename)

% ------------------------------------
% Print finishing the equation solving.
% ------------------------------------

fprintf('Finish loading and interpolating the data\n')

label_fontsize = 14;
axis_fontsize = 12;

figure1 = figure;
scatter(AGE_OLD - age_Sp / 10^6, Sp_record, 'r')
hold on
scatter(AGE_OLD - age_Sp_berner / 10^6, Sp_record_berner, 'g')
hold on
plot(age_all_book, Sp_interp_all, 'b')
xlabel('Age (Ma)', 'FontSize', label_fontsize);
ylabel('Spreading rate', 'FontSize', label_fontsize);

xlim([0 55])

figure2 = figure;
scatter(AGE_OLD - age_dC_carbb / 10^6, dC_carbb_record, 'r')
hold on
plot(age_all_book, dC_carbb_interp_all, 'g')
xlabel('Age (Ma)', 'FontSize', label_fontsize);
ylabel('dC_carbb', 'FontSize', label_fontsize);

figure3 = figure;
scatter(AGE_OLD - age_dC_orgb / 10^6, dC_orgb_record, 'r')
hold on
plot(age_all_book, dC_orgb_interp_all, 'g')
xlabel('Age (Ma)', 'FontSize', label_fontsize);
ylabel('dC_orgb', 'FontSize', label_fontsize);

figure4 = figure;
scatter(AGE_OLD - age_dSr_ocean / 10^6, dSr_ocean_record, 'r')
hold on
plot(age_all_book, dSr_ocean_interp_all, 'g')
xlabel('Age (Ma)', 'FontSize', label_fontsize);
ylabel('dSr_ocean', 'FontSize', label_fontsize);

figure5 = figure;
scatter(AGE_OLD - age_dOs_ocean / 10^6, dOs_ocean_record, 'r')
hold on
plot(age_all_book, dOs_ocean_interp_all, 'g')
xlabel('Age (Ma)', 'FontSize', label_fontsize);
ylabel('dOs_ocean', 'FontSize', label_fontsize);

figure6 = figure;
scatter(AGE_OLD - age_co2 / 10^6, co2_record, 'r')
hold on
plot(age_all_book, co2_interp_all, 'g')
xlabel('Age (Ma)', 'FontSize', label_fontsize);
ylabel('co2', 'FontSize', label_fontsize);

figure7 = figure;
scatter(AGE_OLD - age_co2 / 10^6, co2_record_yerr_low, 'r')
hold on
plot(age_all_book, co2_yerr_low_interp_all, 'g')
xlabel('Age (Ma)', 'FontSize', label_fontsize);
ylabel('co2_record_yerr_low', 'FontSize', label_fontsize);

figure8 = figure;
scatter(AGE_OLD - age_co2 / 10^6, co2_record_yerr_high, 'r')
hold on
plot(age_all_book, co2_yerr_high_interp_all, 'g')
xlabel('Age (Ma)', 'FontSize', label_fontsize);
ylabel('co2_record_yerr_high', 'FontSize', label_fontsize);

figure9 = figure;
scatter(AGE_OLD - age_T / 10^6, T_record, 'r')
hold on
plot(age_all_book, T_interp_all, 'g')
xlabel('Age (Ma)', 'FontSize', label_fontsize);
ylabel('T', 'FontSize', label_fontsize);

figure10 = figure;
scatter(AGE_OLD - age_T / 10^6, T_record_yerr_low, 'r')
hold on
plot(age_all_book, T_yerr_low_interp_all, 'g')
xlabel('Age (Ma)', 'FontSize', label_fontsize);
ylabel('T_record_yerr_low', 'FontSize', label_fontsize);

figure11 = figure;
scatter(AGE_OLD - age_T / 10^6, T_record_yerr_high, 'r')
hold on
plot(age_all_book, T_yerr_high_interp_all, 'g')
xlabel('Age (Ma)', 'FontSize', label_fontsize);
ylabel('T_record_yerr_high', 'FontSize', label_fontsize);

figure12 = figure;
scatter(AGE_OLD - age_Mc / 10^6, Mc_record, 'r')
hold on
plot(age_all_book, Mc_interp_all, 'g')
xlabel('Age (Ma)', 'FontSize', label_fontsize);
ylabel('Mc', 'FontSize', label_fontsize);


figure13 = figure;
scatter(AGE_OLD - age_Mc / 10^6, Mc_sd_record, 'r')
hold on
plot(age_all_book, Mc_sd_interp_all, 'g')
xlabel('Age (Ma)', 'FontSize', label_fontsize);
ylabel('Mc_sd', 'FontSize', label_fontsize);

figure14 = figure;
scatter(AGE_OLD - age_Ma / 10^6, Ma_record, 'r')
hold on
plot(age_all_book, Ma_interp_all, 'g')
xlabel('Age (Ma)', 'FontSize', label_fontsize);
ylabel('Ma', 'FontSize', label_fontsize);


figure15 = figure;
scatter(AGE_OLD - age_Ma / 10^6, Ma_sd_record, 'r')
hold on
plot(age_all_book, Ma_sd_interp_all, 'g')
xlabel('Age (Ma)', 'FontSize', label_fontsize);
ylabel('Ma_sd', 'FontSize', label_fontsize);


figure16 = figure;
scatter(AGE_OLD - age_Li_silw / 10^6, Li_silw, 'r')
hold on
plot(age_all_book, Li_silw_interp_all, 'g')
xlabel('Age (Ma)', 'FontSize', label_fontsize);
ylabel('Ma_sd', 'FontSize', label_fontsize);

figure17 = figure;
scatter(AGE_OLD - age_Li_basw / 10^6, Li_basw, 'r')
hold on
plot(age_all_book, Li_basw_interp_all, 'g')
xlabel('Age (Ma)', 'FontSize', label_fontsize);
ylabel('Ma_sd', 'FontSize', label_fontsize);


figure18 = figure;
scatter(AGE_OLD - age_Li_Mg / 10^6, Li_Mg, 'r')
hold on
plot(age_all_book, Li_Mg_interp_all, 'g')
xlabel('Age (Ma)', 'FontSize', label_fontsize);
ylabel('Ma_sd', 'FontSize', label_fontsize);

figure19 = figure;
scatter(AGE_OLD - age_Li_sedi / 10^6, Li_sedi, 'r')
hold on
plot(age_all_book, Li_sedi_interp_all, 'g')
xlabel('Age (Ma)', 'FontSize', label_fontsize);
ylabel('Ma_sd', 'FontSize', label_fontsize);


e = cputime-cpu_start_data_record;
fprintf('Total data_record time is %f...\n', e /60)
