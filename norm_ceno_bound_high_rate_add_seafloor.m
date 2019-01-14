%----------------------------------------------------------------------
% The code here is Model 3 in the appendix with seafloor integrated

% With Monte Carlo error analysis


% !!! Be aware that some normalization will generate negative number, use
% A(A<0) = 0 to transform or just treat it as fail. I treat it as fail

clear all
close all

fprintf('\n\n')

cpu_start = cputime;

%----------------------------------------------------------------------

% Load the equation solving
load('equi_seafloor.mat')

% ------------------------------------
% Print finishing loading the equations
% ------------------------------------

fprintf('Finish the equation loading\n')

e = cputime-cpu_start;
fprintf('%f mins passed...\n\n', e /60)

OS_MOL = 190.2;
SR_MOL = 87.62;
MG_MOL = 24.31;

% Control the seeding the same !!!!!!!!!!!!!!!!!!
% rng(1);

% This is the time range of each segment
t_interval=5*10^5;

% Set monte carlo number
num_monte = 10;

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
K_csiliw_total = ones(num_points,num_monte) * -1;
K_basw_total = ones(num_points,num_monte) * -1;
K_sediw_total = ones(num_points,num_monte) * -1;
K_carbb_total = ones(num_points,num_monte) * -1;
K_orgb_total = ones(num_points, num_monte) * -1;
CMg_total = ones(num_points, num_monte) * -1;
co2_total = ones(num_points, 1) * -1;
co2_total_err = ones(num_points, 1) * -1;

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

Delta = 30;
dC_carbw = 1.8;
dC_orgw = -28.2;
dC_volc = -5.0;
C_carbw = 12.3 * 10^12;
C_ma = 4.7 * 10^12;
C_volc = 6.0 * 10^12;
A_csiliw = 8.7 * 10^12;
A_basw = 5.0 * 10^12;
Sr_ocean = 1.25 * 10^17;
Sr_cw = 3.4 * 10^10;
f_Sr_csiliw = 0.31;
f_Sr_basw = 0.73;
Sr_dia = 3.4 * 10^9;
dSr_cw = 0.7116;
dSr_carbw = 0.7077;
dSr_mantle = 0.7037;
dSr_dia = 0.7084;
Os_ocean = 7.20 * 10^7;
Os_S_cw = 0.34 * 10^-9;
S_cw = 3.3 * 10^12;
Os_Si_cw = 0.12 * 10^-9;
Si_cw = 4.7 * 10^12;
f_Os_cos = 0.14;
Os_Sr_basw = 0.61 * 10^-7;
dOs_cw = 1.54;
dOs_csiliw = 1.05;
dOs_mantle = 0.126;
CMg_plei = 54.4 * 10^-3; % mol/L
CMg_ini = 37 * 10^-3; % mol/L This needs to be changed for 55 Ma!!!
V_ocean = 1.37 * 10^21; % L
Mg_cw = 5.5 * 10^12;
Mg_C_csiliw = 0.13;
Mg_C_basw = 0.19;
Mg_hy = 4.2 * 10^12; % This also needs to be changed!!!

% 1 sigma

% For figure 1
e_dC_carbb = 0.2;
e_dC_orgb = 0.5;
e_dSr_ocean = 0.00002;
e_dOs_ocean = 0.05;
e_Sp = 0.04;

% Others
e_Delta = 0.5;
e_dC_carbw = 0.1;
e_dC_volc = 0.1;
e_C_carbw = C_carbw * 0.05;
e_C_ma = C_ma * 0.05;
e_C_volc = C_volc * 0.05;
e_A_csiliw = A_csiliw * 0.05;
e_A_basw = A_basw * 0.05;
e_Sr_ocean = Sr_ocean * 0.05;
e_Sr_cw = Sr_cw * 0.05;
e_f_Sr_csiliw = f_Sr_csiliw * 0.05;
e_f_Sr_basw = 0.1;
e_Sr_dia = Sr_dia * 0.2;
e_dSr_cw = 0.0002;
e_dSr_carbw = 0.0002;
e_dSr_mantle = 0.0002;
e_dSr_dia = 0.001;
e_Os_ocean = 0;
e_Os_S_cw = Os_S_cw * 0.05;
e_S_cw = S_cw * 0.05;
e_Os_Si_cw = 0;
e_Si_cw = Si_cw * 0.05;
e_f_Os_cos = f_Os_cos * 0.05; % Note change from 50% to 5%
e_Os_Sr_basw = Os_Sr_basw * 0.05;
e_dOs_cw = 0.05;
e_dOs_csiliw = 0.2;
e_dOs_mantle = 0.05;
e_CMg_plei = CMg_plei * 0.05;
e_CMg_ini = CMg_ini * 0.2;
e_V_ocean = V_ocean * 0.2;
e_Mg_cw = Mg_cw * 0.05;
e_Mg_C_csiliw = Mg_C_csiliw * 0.2;
e_Mg_C_basw = Mg_C_basw * 0.05;
e_Mg_hy = 0;

% Error big in Li 2013
e_Os_S_cw = Os_S_cw * 0.2;
e_S_cw = S_cw * 0.2;
e_Os_Si_cw = Os_Si_cw * 0.2;
e_Si_cw = Si_cw * 0.2;
e_Sr_ocean = Sr_ocean * 0.5;
e_f_Sr_csiliw = f_Sr_csiliw * 0.5;
e_f_Os_cos = f_Os_cos * 0.5;



% ------------------------- All parameters for ophiolite from Jagoutz (all
% mean value)


mafic = 0.49;
Os_conc_mafic = 0.021 * 10^-9;
Os_conc_ultra = 3.60 * 10^-9;
dOs_mafic = 0.143;
dOs_ultra = 0.129;
Sr_conc_mafic = 113 * 10^-6;
Sr_conc_ultra = 38 * 10^-6;
dSr_mafic = 0.7029;
dSr_ultra = 0.7056;
Mg_conc_mafic = 0.0436;
Mg_conc_ultra = 0.2308;
sw_slope = 0.0553;
sw_factor = 18.41;
co2_slope = 0.0642;
co2_factor = 323.44;
sw_runoff_SE = 1372; % mm/year
sw_T_SE = 25;


% sw_runoff_SE = 129; % mm/year
% sw_T_SE = 21.3;


% 1 sigma
e_mafic = 0;
e_Os_conc_mafic = Os_conc_mafic * 0.05;
e_Os_conc_ultra = Os_conc_ultra * 0.05;
e_dOs_mafic = 0.02;
e_dOs_ultra = 0.005;
e_Sr_conc_mafic = Sr_conc_mafic * 0.24;
e_Sr_conc_ultra = 6.56 * 10^-6;
e_dSr_mafic = 0.0004;
e_dSr_ultra = 0.0008;
e_Mg_conc_mafic = 0.0007;
e_Mg_conc_ultra = 0.014;

e_sw_slope = 0;
e_sw_factor = 0;
e_co2_slope = 0;
e_co2_factor = 0;
e_sw_runoff_SE = 0;
e_sw_T_SE = 0;


% Instead of check if the parameter is less than zero and then give it
% zero, just ignore it before computing the value

% -------------------------Carbon isotope value (Unit: per mil)
Delta_monte = normrnd(Delta,e_Delta,1,num_monte);
dC_carbw_monte = normrnd(dC_carbw,e_dC_carbw,1,num_monte);
dC_orgw_monte = dC_carbw_monte - Delta_monte;

dC_volc_monte = normrnd(dC_volc,e_dC_volc,1,num_monte);
dC_ma_monte = dC_volc_monte;
dC_plume_monte = dC_volc_monte;

% -------------------------Carbon flux (Unit: mol/year)
C_carbw_monte = normrnd(C_carbw,e_C_carbw,1,num_monte);
X_org_monte = (dC_carbw_monte - dC_volc_monte) ./ Delta_monte;
C_orgw_monte = X_org_monte ./ (1 - X_org_monte) .* C_carbw_monte;
C_carbw_monte(C_carbw_monte<0) = 0;


C_ma_monte = normrnd(C_ma,e_C_ma,1,num_monte);
C_volc_monte = normrnd(C_volc,e_C_volc,1,num_monte);
C_plume_monte = C_volc_monte - C_ma_monte;

% -------------------------Bicarbonate flux (Unit: mol/year)
A_csiliw_monte = normrnd(A_csiliw,e_A_csiliw,1,num_monte);

A_basw_monte = normrnd(A_basw,e_A_basw,1,num_monte);
A_carbw_monte = 2 * C_carbw_monte;

% -------------------------Sr flux (Unit: mol/year)
Sr_ocean_monte = normrnd(Sr_ocean,e_Sr_ocean,1,num_monte);
Sr_cw_monte = normrnd(Sr_cw,e_Sr_cw,1,num_monte);

f_Sr_csiliw_monte = normrnd(f_Sr_csiliw,e_f_Sr_csiliw,1,num_monte);
Sr_csiliw_monte = Sr_cw_monte .* f_Sr_csiliw_monte;
Sr_carbw_monte = Sr_cw_monte - Sr_csiliw_monte;

Sr_dia_monte = normrnd(Sr_dia,e_Sr_dia,1,num_monte);

% -------------------------Sr isotope value (Unit: per mil)

% Normalize
dSr_cw_monte = normrnd(dSr_cw/(9.43+dSr_cw),9.43*e_dSr_cw/(9.43+dSr_cw)^2,1,num_monte);
dSr_carbw_monte = normrnd(dSr_carbw/(9.43+dSr_carbw),9.43*e_dSr_carbw/(9.43+dSr_carbw)^2,1,num_monte);
dSr_csiliw_monte = (dSr_cw_monte - (1 - f_Sr_csiliw_monte) .* dSr_carbw_monte) ./ f_Sr_csiliw_monte;


% Normalize
dSr_mantle_monte = normrnd(dSr_mantle/(9.43+dSr_mantle),9.43*e_dSr_mantle/(9.43+dSr_mantle)^2,1,num_monte);
dSr_basw_monte = dSr_mantle_monte;
dSr_hy_monte = dSr_basw_monte;
dSr_dia_monte = normrnd(dSr_dia/(9.43+dSr_dia),9.43*e_dSr_dia/(9.43+dSr_dia)^2,1,num_monte);


% -------------------------Os flux (Unit: mol/year)
Os_ocean_monte = normrnd(Os_ocean,e_Os_ocean,1,num_monte);
Os_S_cw_monte = normrnd(Os_S_cw,e_Os_S_cw,1,num_monte);
S_cw_monte = normrnd(S_cw,e_S_cw,1,num_monte);
Os_sediw_monte = S_cw_monte .* Os_S_cw_monte;
Os_Si_cw_monte = normrnd(Os_Si_cw,e_Os_Si_cw,1,num_monte);
Si_cw_monte = normrnd(Si_cw,e_Si_cw,1,num_monte);
Os_csiliw_monte = Si_cw_monte .* Os_Si_cw_monte;

% -------------------------Os isotope value (Unit: per mil)
% Normalize
dOs_cw_monte = normrnd(dOs_cw/(7.4+dOs_cw),7.4*e_dOs_cw/(7.4+dOs_cw)^2,1,num_monte);
dOs_csiliw_monte = normrnd(dOs_csiliw/(7.4+dOs_csiliw),7.4*e_dOs_csiliw/(7.4+dOs_csiliw)^2,1,num_monte);
dOs_mantle_monte = normrnd(dOs_mantle/(7.4+dOs_mantle),7.4*e_dOs_mantle/(7.4+dOs_mantle)^2,1,num_monte);

dOs_basw_monte = dOs_mantle_monte;
dOs_hy_monte = dOs_mantle_monte;
dOs_cos_monte = dOs_mantle_monte;
dOs_sediw_monte = ((Os_sediw_monte + Os_csiliw_monte) .* dOs_cw_monte - Os_csiliw_monte .* dOs_csiliw_monte) ./ Os_sediw_monte;

% -------------------------Mg in the ocean
CMg_plei_monte = normrnd(CMg_plei,e_CMg_plei,1,num_monte);
CMg_ini_monte = normrnd(CMg_ini,e_CMg_ini,1,num_monte);
V_ocean_monte = normrnd(V_ocean,e_V_ocean,1,num_monte);

% -------------------------Mg flux (Unit: mol/year)
Mg_cw_monte = normrnd(Mg_cw,e_Mg_cw,1,num_monte);
Mg_C_csiliw_monte = normrnd(Mg_C_csiliw,e_Mg_C_csiliw,1,num_monte);
Mg_csiliw_monte = Mg_C_csiliw_monte .* A_csiliw_monte;
Mg_carbw_monte = Mg_cw_monte - Mg_csiliw_monte;
Mg_C_basw_monte = normrnd(Mg_C_basw,e_Mg_C_basw,1,num_monte);
Mg_basw_monte = Mg_C_basw_monte .* A_basw_monte;
% Value is tuned so that the resulting Mg concentration of modern ocean
% could match the observations, no error
Mg_hy_monte = normrnd(Mg_hy,e_Mg_hy,1,num_monte);

% -------------------------From Jagoutz
% Normalize
mafic_monte = normrnd(mafic,e_mafic,1,num_monte);
ultra_monte = 1 - mafic_monte;
Os_conc_mafic_monte = normrnd(Os_conc_mafic,e_Os_conc_mafic,1,num_monte);
Os_conc_ultra_monte = normrnd(Os_conc_ultra,e_Os_conc_ultra,1,num_monte);
dOs_mafic_monte = normrnd(dOs_mafic/(7.4+dOs_mafic),7.4*e_dOs_mafic/(7.4+dOs_mafic)^2,1,num_monte);
dOs_ultra_monte = normrnd(dOs_ultra/(7.4+dOs_ultra),7.4*e_dOs_ultra/(7.4+dOs_ultra)^2,1,num_monte);
Sr_conc_mafic_monte = normrnd(Sr_conc_mafic,e_Sr_conc_mafic,1,num_monte);
Sr_conc_ultra_monte = normrnd(Sr_conc_ultra,e_Sr_conc_ultra,1,num_monte);
dSr_mafic_monte = normrnd(dSr_mafic/(9.43+dSr_mafic),9.43*e_dSr_mafic/(9.43+dSr_mafic)^2,1,num_monte);
dSr_ultra_monte = normrnd(dSr_ultra/(9.43+dSr_ultra),9.43*e_dSr_ultra/(9.43+dSr_ultra)^2,1,num_monte);
Mg_conc_mafic_monte = normrnd(Mg_conc_mafic,e_Mg_conc_mafic,1,num_monte);
Mg_conc_ultra_monte = normrnd(Mg_conc_ultra,e_Mg_conc_ultra,1,num_monte);


% Silicate weathering for SE Asia/Indonesia

% silicate dissolving flux
sw_slope_monte = normrnd(sw_slope,e_sw_slope,1,num_monte);
sw_factor_monte = normrnd(sw_factor,e_sw_factor,1,num_monte);
% co2 consumption flux (Note HCO3- flux is the same)
co2_slope_monte = normrnd(co2_slope,e_co2_slope,1,num_monte);
co2_factor_monte = normrnd(co2_factor,e_co2_factor,1,num_monte);
sw_runoff_SE_monte = normrnd(sw_runoff_SE,e_sw_runoff_SE,1,num_monte);
sw_T_SE_monte = normrnd(sw_T_SE,e_sw_T_SE,1,num_monte);


% Add seafloor basalt weathering323.4



% Transform from Caves et al. (2016) supplement CO2 consumption to total
% basalt weathering amount

T_plei = 273.15 + 2; % From Jagoutz et al. (2016)

C_seafloor_plei = 1.75e12; % mol/yr

% bas_seafloor_plei = sw_factor * exp(sw_slope * T_plei) / (co2_factor * exp(co2_slope * T_plei)) * C_seafloor_plei;
% 
% Sr_seafloor_plei_monte = bas_seafloor_plei * 1e6 * Sr_conc_mafic_monte  / SR_MOL;
% 
% Os_seafloor_plei_monte = bas_seafloor_plei * 1e6 * Os_conc_mafic_monte  / OS_MOL;
% 
% Mg_seafloor_plei_monte = bas_seafloor_plei * 1e6 * Mg_conc_mafic_monte / MG_MOL;



% Constitute flux

f_sw_monte = sw_runoff_SE_monte .* sw_factor_monte .* exp(sw_slope_monte .* sw_T_SE_monte) / 1000; % t/km^2/yr

f_Os_mafic_monte = @(t)Area_interp(t) .* mafic_monte .* f_sw_monte * 10^6 .* Os_conc_mafic_monte / OS_MOL;

f_Os_ultra_monte = @(t)Area_interp(t) .* ultra_monte .* f_sw_monte  * 10^6 .* Os_conc_ultra_monte / OS_MOL;

f_A_Jagoutz_monte = @(t)Area_interp(t) .* sw_runoff_SE_monte .* co2_factor_monte .* exp(co2_slope_monte .* sw_T_SE_monte);

f_Sr_mafic_monte = @(t)Area_interp(t) .* mafic_monte .* f_sw_monte * 10^6 .* Sr_conc_mafic_monte / SR_MOL;

f_Sr_ultra_monte = @(t)Area_interp(t) .* ultra_monte .* f_sw_monte  * 10^6 .* Sr_conc_ultra_monte / SR_MOL;

f_Mg_mafic_monte = @(t)Area_interp(t) .* mafic_monte .* f_sw_monte * 10^6 .* Mg_conc_mafic_monte / MG_MOL;

f_Mg_ultra_monte = @(t)Area_interp(t) .* ultra_monte .* f_sw_monte  * 10^6 .* Mg_conc_ultra_monte / MG_MOL;


% ----------------------------------------------- Pleistocene setup started! (We don't include the ophiolite).



% Construct the Monte Carlo data to solve initial values

% Solve carbon system
dC_carbb = dC_carbb_interp(AGE_OLD * 10^6);
dC_carbb_monte = normrnd(dC_carbb,e_dC_carbb,1,num_monte);

dC_carbb_seafloor_monte = dC_carbb_monte - 1; % Here use zeebe and katz

dC_orgb = dC_orgb_interp(AGE_OLD * 10^6);
dC_orgb_monte = normrnd(dC_orgb,e_dC_orgb,1,num_monte);

C_carbb_monte = F_C_carbb(C_carbw_monte, C_ma_monte, C_orgw_monte, C_plume_monte, C_seafloor_plei, dC_ma_monte, ...
    dC_orgb_monte, dC_orgw_monte, dC_carbb_monte, dC_carbw_monte, dC_plume_monte, dC_carbb_seafloor_monte);

C_orgb_monte = F_C_orgb(C_carbw_monte, C_ma_monte, C_orgw_monte, C_plume_monte, C_seafloor_plei, dC_ma_monte, ...
    dC_orgb_monte, dC_orgw_monte, dC_carbb_monte, dC_carbw_monte, dC_plume_monte, dC_carbb_seafloor_monte);

A_rev_monte = F_A_rev(A_basw_monte, A_carbw_monte, A_csiliw_monte, C_carbw_monte, ...
    C_ma_monte, C_orgw_monte, C_plume_monte, C_seafloor_plei, dC_ma_monte, dC_orgb_monte, dC_orgw_monte, ...
    dC_carbb_monte, dC_carbw_monte, dC_plume_monte, dC_carbb_seafloor_monte);


% Solve Sr system

dSr_ocean = dSr_ocean_interp(AGE_OLD* 10^6);
dSr_ocean_monte = normrnd(dSr_ocean/(9.43+dSr_ocean),9.43*e_dSr_ocean/(9.43+dSr_ocean)^2,1,num_monte);

dSr_ocean_old = dSr_ocean_interp(AGE_OLD* 10^6 - t_interval);
dSr_ocean_old_monte = normrnd(dSr_ocean_old/(9.43+dSr_ocean_old),9.43*e_dSr_ocean/(9.43+dSr_ocean_old)^2,1,num_monte);

der_dSr_ocean_monte = (dSr_ocean_monte - dSr_ocean_old_monte) / t_interval;

% This can be changed by der_dSr_ocean_monte, really sensitive!!! Can be
% different from Li et al. (2013)
Sr_basw_hy_monte = F_Sr_basw_hy(Sr_dia_monte, Sr_ocean_monte, Sr_carbw_monte, Sr_csiliw_monte, dSr_dia_monte, ...
    dSr_ocean_monte, dSr_carbw_monte, dSr_mantle_monte, dSr_csiliw_monte, der_dSr_ocean_monte);

% Solve Os system

dOs_ocean = dOs_ocean_interp(AGE_OLD * 10^6);
dOs_ocean_monte = normrnd(dOs_ocean/(7.4+dOs_ocean),7.4*e_dOs_ocean/(7.4+dOs_ocean)^2,1,num_monte);

dOs_ocean_old = dOs_ocean_interp(AGE_OLD * 10^6 - t_interval);
dOs_ocean_old_monte = normrnd(dOs_ocean_old/(7.4+dOs_ocean_old),7.4*e_dOs_ocean/(7.4+dOs_ocean_old)^2,1,num_monte);

der_dOs_ocean_monte = (dOs_ocean_monte - dOs_ocean_old_monte) / t_interval;

% This can be changed by der_dOs_ocean_monte, really sensitive!!! Can be
% different from Li et al. (2013). In deed, Li et al. (2013) uses
% der_dOs_ocean_monte = 0

Os_basw_hy_cos_monte = F_Os_basw_hy_cos(Os_ocean_monte, Os_sediw_monte, Os_csiliw_monte, dOs_ocean_monte, dOs_sediw_monte, ...
    dOs_mantle_monte, dOs_csiliw_monte, der_dOs_ocean_monte);

% OK. Back to construct the remaining monte carlo parameters

% Bicarbonate
A_carbb_monte = 2 * C_carbb_monte;

% Sr island basw and hy

f_Sr_basw_monte = normrnd(f_Sr_basw,e_f_Sr_basw,1,num_monte);

Sr_basw_monte = Sr_basw_hy_monte .* f_Sr_basw_monte;

Sr_hy_monte = Sr_basw_hy_monte - Sr_basw_monte;


% Os cos
f_Os_cos_monte = normrnd(f_Os_cos,e_f_Os_cos,1,num_monte);

Os_cos_monte = Os_basw_hy_cos_monte .* f_Os_cos_monte;

% Os island basalt
Os_Sr_basw_monte = normrnd(Os_Sr_basw,e_Os_Sr_basw,1,num_monte);

Os_basw_monte = Sr_basw_monte .* Os_Sr_basw_monte;

%Os hy
Os_hy_monte = Os_basw_hy_cos_monte - Os_basw_monte - Os_cos_monte;


% ------------------------------------
% Print finishing the Pleistocene initialization
% ------------------------------------

fprintf('Finish the initialization (Pleistocene values)\n')

e = cputime-cpu_start;
fprintf('%f mins passed...\n', e /60)


% ----------------------------------------------- Pleistocene setup finished!

% ----------------------------------------------- 55 Ma initialization
% started

% Begin to solve all parameters (The reason is that CMg has its initial
% value, then others should also have their initial values. Afterwards,
% they can be solved following the time line.

% Calculate initial value at 55 Ma

CMg_total(1,:) = CMg_ini_monte;

t = 0;
Sp_now = Sp_interp(t);
e_Sp_now = Sp_now * e_Sp;
Sp_now_monte = normrnd(Sp_now,e_Sp_now,1,num_monte);

dC_carbb_now = dC_carbb_interp(t);
dC_carbb_now_monte = normrnd(dC_carbb_now,e_dC_carbb,1,num_monte);

dC_orgb_now = dC_orgb_interp(t);
dC_orgb_now_monte = normrnd(dC_orgb_now,e_dC_orgb,1,num_monte);

dOs_ocean_now = dOs_ocean_interp(t);
dOs_ocean_now_monte = normrnd(dOs_ocean_now/(7.4+dOs_ocean_now),7.4*e_dOs_ocean/(7.4+dOs_ocean_now)^2,1,num_monte);
dOs_ocean_old = dOs_ocean_interp(t - t_interval);
dOs_ocean_old_monte = normrnd(dOs_ocean_old/(7.4+dOs_ocean_old),7.4*e_dOs_ocean/(7.4+dOs_ocean_old)^2,1,num_monte);

der_dOs_ocean_monte = (dOs_ocean_now_monte - dOs_ocean_old_monte) / t_interval;

% der_dOs_ocean_monte = der_dOs_ocean_monte * 0;


dSr_ocean_now = dSr_ocean_interp(t);
dSr_ocean_now_monte = normrnd(dSr_ocean_now/(9.43+dSr_ocean_now),9.43*e_dSr_ocean/(9.43+dSr_ocean_now)^2,1,num_monte);

dSr_ocean_old = dSr_ocean_interp(t - t_interval);
dSr_ocean_old_monte = normrnd(dSr_ocean_old/(9.43+dSr_ocean_old),9.43*e_dSr_ocean/(9.43+dSr_ocean_old)^2,1,num_monte);

der_dSr_ocean_monte = (dSr_ocean_now_monte - dSr_ocean_old_monte) / t_interval;

% der_dSr_ocean_monte = der_dSr_ocean_monte * 0;

Os_mafic_monte_now = f_Os_mafic_monte(t);
Os_ultra_monte_now = f_Os_ultra_monte(t);
A_Jagoutz_monte_now = f_A_Jagoutz_monte(t);
Sr_mafic_monte_now = f_Sr_mafic_monte(t);
Sr_ultra_monte_now = f_Sr_ultra_monte(t);
Mg_mafic_monte_now = f_Mg_mafic_monte(t);
Mg_ultra_monte_now = f_Mg_ultra_monte(t);

% determine C_seafloor use berner's equation 2004
T_now = 273.15 + T_interp(t);

factor_T = exp(0.08 * T_now) / exp(0.08 * T_plei);

% C_seafloor_now = factor_T * C_seafloor_plei; % old

C_seafloor_now_monte = factor_T * C_seafloor_plei * Sp_now_monte;

dC_carbb_seafloor_now_monte = dC_carbb_now_monte - 1;

success_array = zeros(num_points, 1);
flag = 0;

Rsil_total = zeros(num_points, num_monte);
Rsil_jagoutz_total_land = zeros(num_points, num_monte);
Rsil_jagoutz_total = zeros(num_points, num_monte);
Imb_total = zeros(num_points, num_monte);
mafic_total = zeros(num_points, num_monte);
sil_total = zeros(num_points, num_monte);
seafloor_total = zeros(num_points,num_monte);

seafloor_total(1,:) = C_seafloor_now_monte;

mafic_seafloor_total = zeros(num_points, num_monte);

% Create fixed modern value monte check (reduce comparison)
monte_check = zeros(num_monte, 1);
for i=1:num_monte
    if (A_basw_monte(i) >= 0) && (A_carbb_monte(i) >= 0) && (A_carbw_monte(i) >= 0) && (A_csiliw_monte(i) >= 0) && ...
            (A_rev_monte(i) >= 0) && (CMg_plei_monte(i) >= 0) && (C_carbb_monte(i) >= 0) && (C_orgb_monte(i) >= 0) && (C_carbw_monte(i) >= 0) && ...
            (C_ma_monte(i) >= 0) && (C_orgw_monte(i) >= 0) && (C_plume_monte(i) >= 0) && (Mg_hy_monte(i) >= 0) && (Mg_basw_monte(i) >= 0) && ...
            (Mg_carbw_monte(i) >= 0)  && (Mg_csiliw_monte(i) >= 0) && (Os_hy_monte(i) >= 0) && (Os_cos_monte(i) >= 0) && (Os_basw_monte(i) >= 0) && ...
            (Os_ocean_monte(i) >= 0) && (Os_sediw_monte(i) >= 0)  && (Os_csiliw_monte(i) >= 0)  && (Sr_hy_monte(i) >= 0) && ...
            (Sr_dia_monte(i) >= 0) && (Sr_basw_monte(i) >= 0)  && (Sr_ocean_monte(i) >= 0) && (Sr_carbw_monte(i) >= 0) && ...
            (Sr_csiliw_monte(i) >= 0) && (V_ocean_monte(i) >= 0)
        monte_check(i) = 1;
    end
end


for i=1:num_monte
    if (monte_check(i) > 0) && (CMg_total(1,i) >= 0) && (Sp_now_monte(i) >= 0) && (A_Jagoutz_monte_now(i) >= 0) && (Mg_mafic_monte_now(i) >= 0) && (Mg_ultra_monte_now(i) >= 0) && (Os_mafic_monte_now(i) >= 0) && ...
            (Os_ultra_monte_now(i) >= 0) && (Sr_mafic_monte_now(i) >= 0) && (Sr_ultra_monte_now(i) >= 0)
        K_csiliw_total(1,i) = F_K_csiliw(A_Jagoutz_monte_now(i), A_basw_monte(i), A_carbb_monte(i), A_carbw_monte(i), A_csiliw_monte(i), A_rev_monte(i), ...
            CMg_total(1,i), CMg_plei_monte(i), C_carbb_monte(i), C_carbw_monte(i), C_ma_monte(i), C_orgw_monte(i), C_plume_monte(i), ...
            C_seafloor_now_monte(i), Mg_hy_monte(i), Mg_basw_monte(i), Mg_mafic_monte_now(i), Mg_carbw_monte(i), Mg_ultra_monte_now(i), Mg_csiliw_monte(i), Os_hy_monte(i), Os_cos_monte(i), Os_basw_monte(i), ...
            Os_mafic_monte_now(i), Os_ocean_monte(i), Os_sediw_monte(i), Os_ultra_monte_now(i), Os_csiliw_monte(i), Sp_now_monte(i), Sr_hy_monte(i), Sr_dia_monte(i), Sr_basw_monte(i), ...
            Sr_mafic_monte_now(i), Sr_ocean_monte(i), Sr_carbw_monte(i), Sr_ultra_monte_now(i), Sr_csiliw_monte(i), V_ocean_monte(i), dC_ma_monte(i), dC_orgb_now_monte(i),dC_orgw_monte(i), ...
            dC_carbb_now_monte(i), dC_carbw_monte(i), dC_plume_monte(i), dC_carbb_seafloor_now_monte(i), dOs_hy_monte(i), dOs_cos_monte(i), dOs_basw_monte(i), ...
            dOs_mafic_monte(i), dOs_ocean_now_monte(i), dOs_sediw_monte(i), dOs_ultra_monte(i), dOs_csiliw_monte(i), dSr_hy_monte(i), dSr_dia_monte(i), dSr_basw_monte(i), ...
            dSr_mafic_monte(i), dSr_ocean_now_monte(i), dSr_carbw_monte(i), dSr_ultra_monte(i), dSr_csiliw_monte(i), der_dOs_ocean_monte(i), der_dSr_ocean_monte(i), t_interval, factor_T);
        
        K_basw_total(1,i) = F_K_basw(A_Jagoutz_monte_now(i), A_basw_monte(i), A_carbb_monte(i), A_carbw_monte(i), A_csiliw_monte(i), A_rev_monte(i), ...
            CMg_total(1,i), CMg_plei_monte(i), C_carbb_monte(i), C_carbw_monte(i), C_ma_monte(i), C_orgw_monte(i), C_plume_monte(i), ...
            C_seafloor_now_monte(i), Mg_hy_monte(i), Mg_basw_monte(i), Mg_mafic_monte_now(i), Mg_carbw_monte(i), Mg_ultra_monte_now(i), Mg_csiliw_monte(i), Os_hy_monte(i), Os_cos_monte(i), Os_basw_monte(i), ...
            Os_mafic_monte_now(i), Os_ocean_monte(i), Os_sediw_monte(i), Os_ultra_monte_now(i), Os_csiliw_monte(i), Sp_now_monte(i), Sr_hy_monte(i), Sr_dia_monte(i), Sr_basw_monte(i), ...
            Sr_mafic_monte_now(i), Sr_ocean_monte(i), Sr_carbw_monte(i), Sr_ultra_monte_now(i), Sr_csiliw_monte(i), V_ocean_monte(i), dC_ma_monte(i), dC_orgb_now_monte(i),dC_orgw_monte(i), ...
            dC_carbb_now_monte(i), dC_carbw_monte(i), dC_plume_monte(i), dC_carbb_seafloor_now_monte(i), dOs_hy_monte(i), dOs_cos_monte(i), dOs_basw_monte(i), ...
            dOs_mafic_monte(i), dOs_ocean_now_monte(i), dOs_sediw_monte(i), dOs_ultra_monte(i), dOs_csiliw_monte(i), dSr_hy_monte(i), dSr_dia_monte(i), dSr_basw_monte(i), ...
            dSr_mafic_monte(i), dSr_ocean_now_monte(i), dSr_carbw_monte(i), dSr_ultra_monte(i), dSr_csiliw_monte(i), der_dOs_ocean_monte(i), der_dSr_ocean_monte(i), t_interval, factor_T);
        
        K_sediw_total(1,i) = F_K_sediw(A_Jagoutz_monte_now(i), A_basw_monte(i), A_carbb_monte(i), A_carbw_monte(i), A_csiliw_monte(i), A_rev_monte(i), ...
            CMg_total(1,i), CMg_plei_monte(i), C_carbb_monte(i), C_carbw_monte(i), C_ma_monte(i), C_orgw_monte(i), C_plume_monte(i), ...
            C_seafloor_now_monte(i), Mg_hy_monte(i), Mg_basw_monte(i), Mg_mafic_monte_now(i), Mg_carbw_monte(i), Mg_ultra_monte_now(i), Mg_csiliw_monte(i), Os_hy_monte(i), Os_cos_monte(i), Os_basw_monte(i), ...
            Os_mafic_monte_now(i), Os_ocean_monte(i), Os_sediw_monte(i), Os_ultra_monte_now(i), Os_csiliw_monte(i), Sp_now_monte(i), Sr_hy_monte(i), Sr_dia_monte(i), Sr_basw_monte(i), ...
            Sr_mafic_monte_now(i), Sr_ocean_monte(i), Sr_carbw_monte(i), Sr_ultra_monte_now(i), Sr_csiliw_monte(i), V_ocean_monte(i), dC_ma_monte(i), dC_orgb_now_monte(i),dC_orgw_monte(i), ...
            dC_carbb_now_monte(i), dC_carbw_monte(i), dC_plume_monte(i), dC_carbb_seafloor_now_monte(i), dOs_hy_monte(i), dOs_cos_monte(i), dOs_basw_monte(i), ...
            dOs_mafic_monte(i), dOs_ocean_now_monte(i), dOs_sediw_monte(i), dOs_ultra_monte(i), dOs_csiliw_monte(i), dSr_hy_monte(i), dSr_dia_monte(i), dSr_basw_monte(i), ...
            dSr_mafic_monte(i), dSr_ocean_now_monte(i), dSr_carbw_monte(i), dSr_ultra_monte(i), dSr_csiliw_monte(i), der_dOs_ocean_monte(i), der_dSr_ocean_monte(i), t_interval, factor_T);
        
        K_carbb_total(1,i) = F_K_carbb(A_Jagoutz_monte_now(i), A_basw_monte(i), A_carbb_monte(i), A_carbw_monte(i), A_csiliw_monte(i), A_rev_monte(i), ...
            CMg_total(1,i), CMg_plei_monte(i), C_carbb_monte(i), C_carbw_monte(i), C_ma_monte(i), C_orgw_monte(i), C_plume_monte(i), ...
            C_seafloor_now_monte(i), Mg_hy_monte(i), Mg_basw_monte(i), Mg_mafic_monte_now(i), Mg_carbw_monte(i), Mg_ultra_monte_now(i), Mg_csiliw_monte(i), Os_hy_monte(i), Os_cos_monte(i), Os_basw_monte(i), ...
            Os_mafic_monte_now(i), Os_ocean_monte(i), Os_sediw_monte(i), Os_ultra_monte_now(i), Os_csiliw_monte(i), Sp_now_monte(i), Sr_hy_monte(i), Sr_dia_monte(i), Sr_basw_monte(i), ...
            Sr_mafic_monte_now(i), Sr_ocean_monte(i), Sr_carbw_monte(i), Sr_ultra_monte_now(i), Sr_csiliw_monte(i), V_ocean_monte(i), dC_ma_monte(i), dC_orgb_now_monte(i),dC_orgw_monte(i), ...
            dC_carbb_now_monte(i), dC_carbw_monte(i), dC_plume_monte(i), dC_carbb_seafloor_now_monte(i), dOs_hy_monte(i), dOs_cos_monte(i), dOs_basw_monte(i), ...
            dOs_mafic_monte(i), dOs_ocean_now_monte(i), dOs_sediw_monte(i), dOs_ultra_monte(i), dOs_csiliw_monte(i), dSr_hy_monte(i), dSr_dia_monte(i), dSr_basw_monte(i), ...
            dSr_mafic_monte(i), dSr_ocean_now_monte(i), dSr_carbw_monte(i), dSr_ultra_monte(i), dSr_csiliw_monte(i), der_dOs_ocean_monte(i), der_dSr_ocean_monte(i), t_interval, factor_T);
        
        K_orgb_total(1,i) = F_K_orgb(A_Jagoutz_monte_now(i), A_basw_monte(i), A_carbb_monte(i), A_carbw_monte(i), A_csiliw_monte(i), A_rev_monte(i), ...
            CMg_total(1,i), CMg_plei_monte(i), C_carbb_monte(i), C_carbw_monte(i), C_ma_monte(i), C_orgb_monte(i), C_orgw_monte(i), C_plume_monte(i), ...
            C_seafloor_now_monte(i), Mg_hy_monte(i), Mg_basw_monte(i), Mg_mafic_monte_now(i), Mg_carbw_monte(i), Mg_ultra_monte_now(i), Mg_csiliw_monte(i), Os_hy_monte(i), Os_cos_monte(i), Os_basw_monte(i), ...
            Os_mafic_monte_now(i), Os_ocean_monte(i), Os_sediw_monte(i), Os_ultra_monte_now(i), Os_csiliw_monte(i), Sp_now_monte(i), Sr_hy_monte(i), Sr_dia_monte(i), Sr_basw_monte(i), ...
            Sr_mafic_monte_now(i), Sr_ocean_monte(i), Sr_carbw_monte(i), Sr_ultra_monte_now(i), Sr_csiliw_monte(i), V_ocean_monte(i), dC_ma_monte(i), dC_orgb_now_monte(i),dC_orgw_monte(i), ...
            dC_carbb_now_monte(i), dC_carbw_monte(i), dC_plume_monte(i), dC_carbb_seafloor_now_monte(i), dOs_hy_monte(i), dOs_cos_monte(i), dOs_basw_monte(i), ...
            dOs_mafic_monte(i), dOs_ocean_now_monte(i), dOs_sediw_monte(i), dOs_ultra_monte(i), dOs_csiliw_monte(i), dSr_hy_monte(i), dSr_dia_monte(i), dSr_basw_monte(i), ...
            dSr_mafic_monte(i), dSr_ocean_now_monte(i), dSr_carbw_monte(i), dSr_ultra_monte(i), dSr_csiliw_monte(i), der_dOs_ocean_monte(i), der_dSr_ocean_monte(i), t_interval, factor_T);
        
        if (K_csiliw_total(1,i) >= 0) && (K_basw_total(1,i) >= 0) && (K_sediw_total(1,i) >= 0) && (K_carbb_total(1,i) >= 0) && (K_orgb_total(1,i) >= 0)
            flag = flag + 1;
            Rsil_total(1,i) = (K_csiliw_total(1,i) * A_csiliw_monte(i) + K_basw_total(1,i) * A_basw_monte(i) + 2 * C_seafloor_now_monte(i)) / (A_csiliw_monte(i) + A_basw_monte(i) + 2 * C_seafloor_plei);
            Rsil_jagoutz_total(1,i) = (K_csiliw_total(1,i) * A_csiliw_monte(i) + K_basw_total(1,i) * A_basw_monte(i) + A_Jagoutz_monte_now(i) + 2 * C_seafloor_now_monte(i)) / (A_csiliw_monte(i) + A_basw_monte(i) + 2 * C_seafloor_plei);
            Rsil_jagoutz_total_land(1,i) = (K_csiliw_total(1,i) * A_csiliw_monte(i) + K_basw_total(1,i) * A_basw_monte(i) + A_Jagoutz_monte_now(i)) / (A_csiliw_monte(i) + A_basw_monte(i));
            Imb_total(1,i) = Sp_now_monte(i) * C_ma_monte(i) + C_plume_monte(i) + K_sediw_total(1,i) * C_orgw_monte(i) - K_orgb_total(1,i) * C_orgb_monte(i) - K_csiliw_total(1,i) * A_csiliw_monte(i) / 2 - K_basw_total(1,i) * A_basw_monte(i) / 2 - A_Jagoutz_monte_now(i) / 2 + Sp_now_monte(i) *  CMg_total(1,i) / CMg_plei * A_rev_monte(i) / 2 - C_seafloor_now_monte(i);
            mafic_total(1,i) =  K_basw_total(1,i) * A_basw_monte(i) + A_Jagoutz_monte_now(i);
            mafic_seafloor_total(1,i) = mafic_total(1,i) + seafloor_total(1,i)*2;
            sil_total(1,i) =  K_basw_total(1,i) * A_basw_monte(i) + K_csiliw_total(1,i) * A_csiliw_monte(i) + A_Jagoutz_monte_now(i) + 2 * C_seafloor_now_monte(i);
        else
            K_csiliw_total(1,i) = -2;
            K_basw_total(1,i) = -2;
            K_sediw_total(1,i) = -2;
            K_carbb_total(1,i) = -2;
            K_orgb_total(1,i) = -2;
        end
        
    end
    
end

success_array(1) = flag;


dOs_ocean_old_monte = dOs_ocean_now_monte;
dSr_ocean_old_monte = dSr_ocean_now_monte;

co2_total(1) = co2_interp(0);

if (co2_yerr_low_interp(0) > co2_yerr_high_interp(0))
    co2_total_err(1) = co2_yerr_low_interp(0);
else
    co2_total_err(1) = co2_yerr_high_interp(0);
end

% ------------------------------------
% Print finishing the initialization
% ------------------------------------

fprintf('Finish the initialization (55 Ma values)\n')

e = cputime-cpu_start;
fprintf('%f mins passed...\n', e /60)

% ----------------------------------------------- 55 Ma initialization
% finished


for j=1:num_points-1
    t = t_interval * j;
    
    Sp_now = Sp_interp(t);
    e_Sp_now = Sp_now * e_Sp;
    Sp_now_monte = normrnd(Sp_now,e_Sp_now,1,num_monte);
    
    dC_carbb_now = dC_carbb_interp(t);
    dC_carbb_now_monte = normrnd(dC_carbb_now,e_dC_carbb,1,num_monte);
    
    dC_orgb_now = dC_orgb_interp(t);
    dC_orgb_now_monte = normrnd(dC_orgb_now,e_dC_orgb,1,num_monte);
    
    dOs_ocean_now = dOs_ocean_interp(t);
    dOs_ocean_now_monte = normrnd(dOs_ocean_now/(7.4+dOs_ocean_now),7.4*e_dOs_ocean/(7.4+dOs_ocean_now)^2,1,num_monte);
    der_dOs_ocean_monte = (dOs_ocean_now_monte - dOs_ocean_old_monte) / t_interval;
    
    
    dSr_ocean_now = dSr_ocean_interp(t);
    dSr_ocean_now_monte = normrnd(dSr_ocean_now/(9.43+dSr_ocean_now),9.43*e_dSr_ocean/(9.43+dSr_ocean_now)^2,1,num_monte);
    der_dSr_ocean_monte = (dSr_ocean_now_monte - dSr_ocean_old_monte) / t_interval;
   
    CMg_total_old = CMg_total(j,:);
    
    CMg_total_old(CMg_total_old < 0) = NaN;

    CMg_total_old_ave = nanmean(CMg_total_old);
    CMg_total_old_std = nanstd(CMg_total_old);
    CMg_total_old_monte = normrnd(CMg_total_old_ave,CMg_total_old_std,1,num_monte);
    
    Os_mafic_monte_now = f_Os_mafic_monte(t);
    Os_ultra_monte_now = f_Os_ultra_monte(t);
    A_Jagoutz_monte_now = f_A_Jagoutz_monte(t);
    Sr_mafic_monte_now = f_Sr_mafic_monte(t);
    Sr_ultra_monte_now = f_Sr_ultra_monte(t);
    Mg_mafic_monte_now = f_Mg_mafic_monte(t);
    Mg_ultra_monte_now = f_Mg_ultra_monte(t);
    
    % determine C_seafloor
    T_now = 273.15 + T_interp(t);

    factor_T = exp(0.08 * T_now) / exp(0.08 * T_plei);
    
%     factor_T = 1;

%     C_seafloor_now = factor_T * C_seafloor_plei;
    
    C_seafloor_now_monte = factor_T * C_seafloor_plei * Sp_now_monte;

    dC_carbb_seafloor_now_monte = dC_carbb_now_monte - 1;

    seafloor_total(j+1,:) = C_seafloor_now_monte;
    
    flag = 0;
    
    for i=1:num_monte
        if (monte_check(i) > 0) && (CMg_total_old_monte(i) >= 0) && (Sp_now_monte(i) >= 0) && (A_Jagoutz_monte_now(i) >= 0) && (Mg_mafic_monte_now(i) >= 0) && (Mg_ultra_monte_now(i) >= 0) && (Os_mafic_monte_now(i) >= 0) && ...
                (Os_ultra_monte_now(i) >= 0) && (Sr_mafic_monte_now(i) >= 0) && (Sr_ultra_monte_now(i) >= 0)
            
            K_csiliw_total(j+1,i) = F_K_csiliw(A_Jagoutz_monte_now(i), A_basw_monte(i), A_carbb_monte(i), A_carbw_monte(i), A_csiliw_monte(i), A_rev_monte(i), ...
            CMg_total(1,i), CMg_plei_monte(i), C_carbb_monte(i), C_carbw_monte(i), C_ma_monte(i), C_orgw_monte(i), C_plume_monte(i), ...
            C_seafloor_now_monte(i), Mg_hy_monte(i), Mg_basw_monte(i), Mg_mafic_monte_now(i), Mg_carbw_monte(i), Mg_ultra_monte_now(i), Mg_csiliw_monte(i), Os_hy_monte(i), Os_cos_monte(i), Os_basw_monte(i), ...
            Os_mafic_monte_now(i), Os_ocean_monte(i), Os_sediw_monte(i), Os_ultra_monte_now(i), Os_csiliw_monte(i), Sp_now_monte(i), Sr_hy_monte(i), Sr_dia_monte(i), Sr_basw_monte(i), ...
            Sr_mafic_monte_now(i), Sr_ocean_monte(i), Sr_carbw_monte(i), Sr_ultra_monte_now(i), Sr_csiliw_monte(i), V_ocean_monte(i), dC_ma_monte(i), dC_orgb_now_monte(i),dC_orgw_monte(i), ...
            dC_carbb_now_monte(i), dC_carbw_monte(i), dC_plume_monte(i), dC_carbb_seafloor_now_monte(i), dOs_hy_monte(i), dOs_cos_monte(i), dOs_basw_monte(i), ...
            dOs_mafic_monte(i), dOs_ocean_now_monte(i), dOs_sediw_monte(i), dOs_ultra_monte(i), dOs_csiliw_monte(i), dSr_hy_monte(i), dSr_dia_monte(i), dSr_basw_monte(i), ...
            dSr_mafic_monte(i), dSr_ocean_now_monte(i), dSr_carbw_monte(i), dSr_ultra_monte(i), dSr_csiliw_monte(i), der_dOs_ocean_monte(i), der_dSr_ocean_monte(i), t_interval, factor_T);
            
            
            K_basw_total(j+1,i) = F_K_basw(A_Jagoutz_monte_now(i), A_basw_monte(i), A_carbb_monte(i), A_carbw_monte(i), A_csiliw_monte(i), A_rev_monte(i), ...
            CMg_total(1,i), CMg_plei_monte(i), C_carbb_monte(i), C_carbw_monte(i), C_ma_monte(i), C_orgw_monte(i), C_plume_monte(i), ...
            C_seafloor_now_monte(i), Mg_hy_monte(i), Mg_basw_monte(i), Mg_mafic_monte_now(i), Mg_carbw_monte(i), Mg_ultra_monte_now(i), Mg_csiliw_monte(i), Os_hy_monte(i), Os_cos_monte(i), Os_basw_monte(i), ...
            Os_mafic_monte_now(i), Os_ocean_monte(i), Os_sediw_monte(i), Os_ultra_monte_now(i), Os_csiliw_monte(i), Sp_now_monte(i), Sr_hy_monte(i), Sr_dia_monte(i), Sr_basw_monte(i), ...
            Sr_mafic_monte_now(i), Sr_ocean_monte(i), Sr_carbw_monte(i), Sr_ultra_monte_now(i), Sr_csiliw_monte(i), V_ocean_monte(i), dC_ma_monte(i), dC_orgb_now_monte(i),dC_orgw_monte(i), ...
            dC_carbb_now_monte(i), dC_carbw_monte(i), dC_plume_monte(i), dC_carbb_seafloor_now_monte(i), dOs_hy_monte(i), dOs_cos_monte(i), dOs_basw_monte(i), ...
            dOs_mafic_monte(i), dOs_ocean_now_monte(i), dOs_sediw_monte(i), dOs_ultra_monte(i), dOs_csiliw_monte(i), dSr_hy_monte(i), dSr_dia_monte(i), dSr_basw_monte(i), ...
            dSr_mafic_monte(i), dSr_ocean_now_monte(i), dSr_carbw_monte(i), dSr_ultra_monte(i), dSr_csiliw_monte(i), der_dOs_ocean_monte(i), der_dSr_ocean_monte(i), t_interval, factor_T);
            
            K_sediw_total(j+1,i) = F_K_sediw(A_Jagoutz_monte_now(i), A_basw_monte(i), A_carbb_monte(i), A_carbw_monte(i), A_csiliw_monte(i), A_rev_monte(i), ...
            CMg_total(1,i), CMg_plei_monte(i), C_carbb_monte(i), C_carbw_monte(i), C_ma_monte(i), C_orgw_monte(i), C_plume_monte(i), ...
            C_seafloor_now_monte(i), Mg_hy_monte(i), Mg_basw_monte(i), Mg_mafic_monte_now(i), Mg_carbw_monte(i), Mg_ultra_monte_now(i), Mg_csiliw_monte(i), Os_hy_monte(i), Os_cos_monte(i), Os_basw_monte(i), ...
            Os_mafic_monte_now(i), Os_ocean_monte(i), Os_sediw_monte(i), Os_ultra_monte_now(i), Os_csiliw_monte(i), Sp_now_monte(i), Sr_hy_monte(i), Sr_dia_monte(i), Sr_basw_monte(i), ...
            Sr_mafic_monte_now(i), Sr_ocean_monte(i), Sr_carbw_monte(i), Sr_ultra_monte_now(i), Sr_csiliw_monte(i), V_ocean_monte(i), dC_ma_monte(i), dC_orgb_now_monte(i),dC_orgw_monte(i), ...
            dC_carbb_now_monte(i), dC_carbw_monte(i), dC_plume_monte(i), dC_carbb_seafloor_now_monte(i), dOs_hy_monte(i), dOs_cos_monte(i), dOs_basw_monte(i), ...
            dOs_mafic_monte(i), dOs_ocean_now_monte(i), dOs_sediw_monte(i), dOs_ultra_monte(i), dOs_csiliw_monte(i), dSr_hy_monte(i), dSr_dia_monte(i), dSr_basw_monte(i), ...
            dSr_mafic_monte(i), dSr_ocean_now_monte(i), dSr_carbw_monte(i), dSr_ultra_monte(i), dSr_csiliw_monte(i), der_dOs_ocean_monte(i), der_dSr_ocean_monte(i), t_interval, factor_T);
            
            K_carbb_total(j+1,i) = F_K_carbb(A_Jagoutz_monte_now(i), A_basw_monte(i), A_carbb_monte(i), A_carbw_monte(i), A_csiliw_monte(i), A_rev_monte(i), ...
            CMg_total(1,i), CMg_plei_monte(i), C_carbb_monte(i), C_carbw_monte(i), C_ma_monte(i), C_orgw_monte(i), C_plume_monte(i), ...
            C_seafloor_now_monte(i), Mg_hy_monte(i), Mg_basw_monte(i), Mg_mafic_monte_now(i), Mg_carbw_monte(i), Mg_ultra_monte_now(i), Mg_csiliw_monte(i), Os_hy_monte(i), Os_cos_monte(i), Os_basw_monte(i), ...
            Os_mafic_monte_now(i), Os_ocean_monte(i), Os_sediw_monte(i), Os_ultra_monte_now(i), Os_csiliw_monte(i), Sp_now_monte(i), Sr_hy_monte(i), Sr_dia_monte(i), Sr_basw_monte(i), ...
            Sr_mafic_monte_now(i), Sr_ocean_monte(i), Sr_carbw_monte(i), Sr_ultra_monte_now(i), Sr_csiliw_monte(i), V_ocean_monte(i), dC_ma_monte(i), dC_orgb_now_monte(i),dC_orgw_monte(i), ...
            dC_carbb_now_monte(i), dC_carbw_monte(i), dC_plume_monte(i), dC_carbb_seafloor_now_monte(i), dOs_hy_monte(i), dOs_cos_monte(i), dOs_basw_monte(i), ...
            dOs_mafic_monte(i), dOs_ocean_now_monte(i), dOs_sediw_monte(i), dOs_ultra_monte(i), dOs_csiliw_monte(i), dSr_hy_monte(i), dSr_dia_monte(i), dSr_basw_monte(i), ...
            dSr_mafic_monte(i), dSr_ocean_now_monte(i), dSr_carbw_monte(i), dSr_ultra_monte(i), dSr_csiliw_monte(i), der_dOs_ocean_monte(i), der_dSr_ocean_monte(i), t_interval, factor_T);
            
            K_orgb_total(j+1,i) = F_K_orgb(A_Jagoutz_monte_now(i), A_basw_monte(i), A_carbb_monte(i), A_carbw_monte(i), A_csiliw_monte(i), A_rev_monte(i), ...
            CMg_total(1,i), CMg_plei_monte(i), C_carbb_monte(i), C_carbw_monte(i), C_ma_monte(i), C_orgb_monte(i), C_orgw_monte(i), C_plume_monte(i), ...
            C_seafloor_now_monte(i), Mg_hy_monte(i), Mg_basw_monte(i), Mg_mafic_monte_now(i), Mg_carbw_monte(i), Mg_ultra_monte_now(i), Mg_csiliw_monte(i), Os_hy_monte(i), Os_cos_monte(i), Os_basw_monte(i), ...
            Os_mafic_monte_now(i), Os_ocean_monte(i), Os_sediw_monte(i), Os_ultra_monte_now(i), Os_csiliw_monte(i), Sp_now_monte(i), Sr_hy_monte(i), Sr_dia_monte(i), Sr_basw_monte(i), ...
            Sr_mafic_monte_now(i), Sr_ocean_monte(i), Sr_carbw_monte(i), Sr_ultra_monte_now(i), Sr_csiliw_monte(i), V_ocean_monte(i), dC_ma_monte(i), dC_orgb_now_monte(i),dC_orgw_monte(i), ...
            dC_carbb_now_monte(i), dC_carbw_monte(i), dC_plume_monte(i), dC_carbb_seafloor_now_monte(i), dOs_hy_monte(i), dOs_cos_monte(i), dOs_basw_monte(i), ...
            dOs_mafic_monte(i), dOs_ocean_now_monte(i), dOs_sediw_monte(i), dOs_ultra_monte(i), dOs_csiliw_monte(i), dSr_hy_monte(i), dSr_dia_monte(i), dSr_basw_monte(i), ...
            dSr_mafic_monte(i), dSr_ocean_now_monte(i), dSr_carbw_monte(i), dSr_ultra_monte(i), dSr_csiliw_monte(i), der_dOs_ocean_monte(i), der_dSr_ocean_monte(i), t_interval, factor_T);
            
            
            CMg_total(j+1,i) = F_CMg(A_Jagoutz_monte_now(i), A_basw_monte(i), A_carbb_monte(i), A_carbw_monte(i), A_csiliw_monte(i), A_rev_monte(i), ...
            CMg_total(1,i), CMg_plei_monte(i), C_carbb_monte(i), C_carbw_monte(i), C_ma_monte(i), C_orgw_monte(i), C_plume_monte(i), ...
            C_seafloor_now_monte(i), Mg_hy_monte(i), Mg_basw_monte(i), Mg_mafic_monte_now(i), Mg_carbw_monte(i), Mg_ultra_monte_now(i), Mg_csiliw_monte(i), Os_hy_monte(i), Os_cos_monte(i), Os_basw_monte(i), ...
            Os_mafic_monte_now(i), Os_ocean_monte(i), Os_sediw_monte(i), Os_ultra_monte_now(i), Os_csiliw_monte(i), Sp_now_monte(i), Sr_hy_monte(i), Sr_dia_monte(i), Sr_basw_monte(i), ...
            Sr_mafic_monte_now(i), Sr_ocean_monte(i), Sr_carbw_monte(i), Sr_ultra_monte_now(i), Sr_csiliw_monte(i), V_ocean_monte(i), dC_ma_monte(i), dC_orgb_now_monte(i),dC_orgw_monte(i), ...
            dC_carbb_now_monte(i), dC_carbw_monte(i), dC_plume_monte(i), dC_carbb_seafloor_now_monte(i), dOs_hy_monte(i), dOs_cos_monte(i), dOs_basw_monte(i), ...
            dOs_mafic_monte(i), dOs_ocean_now_monte(i), dOs_sediw_monte(i), dOs_ultra_monte(i), dOs_csiliw_monte(i), dSr_hy_monte(i), dSr_dia_monte(i), dSr_basw_monte(i), ...
            dSr_mafic_monte(i), dSr_ocean_now_monte(i), dSr_carbw_monte(i), dSr_ultra_monte(i), dSr_csiliw_monte(i), der_dOs_ocean_monte(i), der_dSr_ocean_monte(i), t_interval, factor_T);
            
            if (K_csiliw_total(j+1,i) >= 0) && (K_basw_total(j+1,i) >= 0) && (K_sediw_total(j+1,i) >= 0) && (K_carbb_total(j+1,i) >= 0) && (K_orgb_total(j+1,i) >= 0) && (CMg_total(j+1,i) >= 0)
                flag = flag + 1;
                Rsil_total(j+1,i) = (K_csiliw_total(j+1,i) * A_csiliw_monte(i) + K_basw_total(j+1,i) * A_basw_monte(i) +  2 * C_seafloor_now_monte(i)) / (A_csiliw_monte(i) + A_basw_monte(i) +  2 * C_seafloor_plei);
                Rsil_jagoutz_total(j+1,i) = (K_csiliw_total(j+1,i) * A_csiliw_monte(i) + K_basw_total(j+1,i) * A_basw_monte(i) + A_Jagoutz_monte_now(i) +  2 * C_seafloor_now_monte(i)) / (A_csiliw_monte(i) + A_basw_monte(i) +  2 * C_seafloor_plei);           
                Rsil_jagoutz_total_land(j+1,i) = (K_csiliw_total(j+1,i) * A_csiliw_monte(i) + K_basw_total(j+1,i) * A_basw_monte(i) + A_Jagoutz_monte_now(i)) / (A_csiliw_monte(i) + A_basw_monte(i));
                
                Imb_total(j+1,i) = Sp_now_monte(i) * C_ma_monte(i) + C_plume_monte(i) + K_sediw_total(j+1,i) * C_orgw_monte(i) - K_orgb_total(j+1,i) * C_orgb_monte(i) - K_csiliw_total(j+1,i) * A_csiliw_monte(i) / 2 - K_basw_total(j+1,i) * A_basw_monte(i) / 2 - A_Jagoutz_monte_now(i) / 2 + Sp_now_monte(i) *  CMg_total(j+1,i) / CMg_plei * A_rev_monte(i) / 2 - C_seafloor_now_monte(i);
                mafic_total(j+1,i) =  K_basw_total(j+1,i) * A_basw_monte(i) + A_Jagoutz_monte_now(i);
                mafic_seafloor_total(j+1,i) = mafic_total(j+1,i) + seafloor_total(j+1,i) * 2;
                sil_total(j+1,i) =  K_basw_total(j+1,i) * A_basw_monte(i) + K_csiliw_total(j+1,i) * A_csiliw_monte(i) + A_Jagoutz_monte_now(i) +  2 * C_seafloor_now_monte(i); % wrong. need to add seafloor weathering
            else
                K_csiliw_total(j+1,i) = -2;
                K_basw_total(j+1,i) = -2;
                K_sediw_total(j+1,i) = -2;
                K_carbb_total(j+1,i) = -2;
                K_orgb_total(j+1,i) = -2;
                CMg_total(j+1,i) = -2;
            end
        end
        dOs_ocean_old_monte = dOs_ocean_now_monte;
        dSr_ocean_old_monte = dSr_ocean_now_monte;
        success_array(j+1) = flag;
        co2_total(j+1) = co2_interp(t);

        if (co2_yerr_low_interp(t) > co2_yerr_high_interp(t))
            co2_total_err(j+1) = co2_yerr_low_interp(t);
        else
            co2_total_err(j+1) = co2_yerr_high_interp(t);
        end
    end
    
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


K_csiliw_total_origin = K_csiliw_total;
K_csiliw_total(K_csiliw_total<0) = NaN;
K_csiliw_ave = nanmean(K_csiliw_total, 2);
K_csiliw_std = nanstd(K_csiliw_total, 0, 2);
K_csiliw_low = K_csiliw_ave - K_csiliw_std;
K_csiliw_high = K_csiliw_ave + K_csiliw_std;

F_csiliw_ave = K_csiliw_ave * A_csiliw / 2; % /2 means only silicate weathering
F_csiliw_std = sqrt(K_csiliw_ave.^2 * e_A_csiliw^2 + A_csiliw^2 * K_csiliw_std.^2) / 2;


K_basw_total_origin = K_basw_total;
K_basw_total(K_basw_total<0) = NaN;
K_basw_ave = nanmean(K_basw_total, 2);
K_basw_std = nanstd(K_basw_total, 0, 2);
K_basw_low = K_basw_ave - K_basw_std;
K_basw_high = K_basw_ave + K_basw_std;

F_basw_ave = K_basw_ave * A_basw / 2; % /2 means only silicate weathering
F_basw_std = sqrt(K_basw_ave.^2 * e_A_basw^2 + A_basw^2 * K_basw_std.^2) / 2;

Rsil_total_origin = Rsil_total;
Rsil_total(Rsil_total==0) = NaN;
Rsil_ave = nanmean(Rsil_total, 2);
Rsil_std = nanstd(Rsil_total, 0, 2);
Rsil_low = Rsil_ave - Rsil_std;
Rsil_high = Rsil_ave + Rsil_std;

Rsil_jagoutz_total_origin = Rsil_jagoutz_total;
Rsil_jagoutz_total(Rsil_jagoutz_total==0) = NaN;
Rsil_jagoutz_ave = nanmean(Rsil_jagoutz_total, 2);
Rsil_jagoutz_std = nanstd(Rsil_jagoutz_total, 0, 2);
Rsil_jagoutz_low = Rsil_jagoutz_ave - Rsil_jagoutz_std;
Rsil_jagoutz_high = Rsil_jagoutz_ave + Rsil_jagoutz_std;

Rsil_jagoutz_total_land_origin = Rsil_jagoutz_total_land;
Rsil_jagoutz_total_land(Rsil_jagoutz_total_land==0) = NaN;
Rsil_jagoutz_land_ave = nanmean(Rsil_jagoutz_total_land, 2);
Rsil_jagoutz_land_std = nanstd(Rsil_jagoutz_total_land, 0, 2);
Rsil_jagoutz_land_low = Rsil_jagoutz_land_ave - Rsil_jagoutz_land_std;
Rsil_jagoutz_land_high = Rsil_jagoutz_land_ave + Rsil_jagoutz_land_std;



mafic_total_origin = mafic_total;
mafic_total(mafic_total==0) = NaN;
mafic_total_ave = nanmean(mafic_total, 2);
mafic_total_std = nanstd(mafic_total, 0, 2);
mafic_total_low = mafic_total_ave - mafic_total_std;
mafic_total_high = mafic_total_ave + mafic_total_std;

mafic_seafloor_total_origin = mafic_seafloor_total;
mafic_seafloor_total(mafic_seafloor_total==0) = NaN;
mafic_seafloor_total_ave = nanmean(mafic_seafloor_total, 2);
mafic_seafloor_total_std = nanstd(mafic_seafloor_total, 0, 2);
mafic_seafloor_total_low = mafic_seafloor_total_ave - mafic_seafloor_total_std;
mafic_seafloor_total_high = mafic_seafloor_total_ave + mafic_seafloor_total_std;

sil_total_origin = sil_total;
sil_total(sil_total==0) = NaN;
sil_total_ave = nanmean(sil_total, 2);
sil_total_std = nanstd(sil_total, 0, 2);
sil_total_low = sil_total_ave - sil_total_std;
sil_total_high = sil_total_ave + sil_total_std;


K_sediw_total_origin = K_sediw_total;
K_sediw_total(K_sediw_total<0) = NaN;
K_sediw_ave = nanmean(K_sediw_total, 2);
K_sediw_std = nanstd(K_sediw_total, 0, 2);
K_sediw_low = K_sediw_ave - K_sediw_std;
K_sediw_high = K_sediw_ave + K_sediw_std;

F_carbw_ave = K_sediw_ave * C_carbw;
F_carbw_std = sqrt(K_sediw_ave.^2 * e_C_carbw^2 + C_carbw^2 * K_sediw_std.^2);

C_orgw = mean(C_orgw_monte);
e_C_orgw = std(C_orgw_monte);

F_orgw_ave = K_sediw_ave * C_orgw;
F_orgw_std = sqrt(K_sediw_ave.^2 * e_C_orgw^2 + C_orgw^2 * K_sediw_std.^2);

K_carbb_total_origin = K_carbb_total;
K_carbb_total(K_carbb_total<0) = NaN;
K_carbb_ave = nanmean(K_carbb_total, 2);
K_carbb_std = nanstd(K_carbb_total, 0, 2);
K_carbb_low = K_carbb_ave - K_carbb_std;
K_carbb_high = K_carbb_ave + K_carbb_std;

C_carbb = mean(C_carbb_monte);
e_C_carbb= std(C_carbb_monte);

F_carbb_ave = K_carbb_ave * C_carbb;
F_carbb_std = sqrt(K_carbb_ave.^2 * e_C_carbb^2 + C_carbb^2 * K_carbb_std.^2);

K_orgb_total_origin = K_orgb_total;
K_orgb_total(K_orgb_total<0) = NaN;
K_orgb_ave = nanmean(K_orgb_total, 2);
K_orgb_std = nanstd(K_orgb_total, 0, 2);
K_orgb_low = K_orgb_ave - K_orgb_std;
K_orgb_high = K_orgb_ave + K_orgb_std;

C_orgb = mean(C_orgb_monte);
e_C_orgb= std(C_orgb_monte);

F_orgb_ave = K_orgb_ave * C_orgb;
F_orgb_std = sqrt(K_orgb_ave.^2 * e_C_orgb^2 + C_orgb^2 * K_orgb_std.^2);

CMg_total_origin = CMg_total;
CMg_total(CMg_total<0) = NaN;
CMg_ave = nanmean(CMg_total, 2);
CMg_std = nanstd(CMg_total, 0, 2);
CMg_low = CMg_ave - CMg_std;
CMg_high = CMg_ave + CMg_std;

Imb_total_origin = Imb_total;
Imb_ave = nanmean(Imb_total, 2);
Imb_std = nanstd(Imb_total, 0, 2);
Imb_low = Imb_ave - Imb_std;
Imb_high = Imb_ave + Imb_std;


seafloor_total_origin = seafloor_total;
seafloor_total(seafloor_total<0) = NaN;
seafloor_ave = nanmean(seafloor_total, 2);
seafloor_std = nanstd(seafloor_total, 0, 2);
seafloor_low = seafloor_ave - seafloor_std;
seafloor_high = seafloor_ave + seafloor_std;

% Create file directory

file_dir = strcat('norm_high_rate_seafloor_', num2str(t_interval / 1000), 'ky_', num2str(num_monte), 'monte');

if not(exist(file_dir, 'dir'))
    mkdir(file_dir)
end

% When concatenate dir (difference in pc and mac/unix)
if ispc
    book_name = strcat(file_dir, '\', 'output', '.xlsx');
    csv_name = strcat(file_dir, '\', 'output', '.csv');
    sil_name = strcat(file_dir, '\', 'output', '_sil.jpg');
    bas_name = strcat(file_dir, '\', 'output', '_bas.jpg');
    sed_name = strcat(file_dir, '\', 'output', '_sed.jpg');
    car_name = strcat(file_dir, '\', 'output', '_car.jpg');
    oc_name = strcat(file_dir, '\', 'output', '_oc.jpg');
    cmg_name = strcat(file_dir, '\', 'output', '_cmg.jpg');
    suc_name = strcat(file_dir, '\', 'output', '_success.jpg');
    Rsil_name = strcat(file_dir, '\', 'output', '_Rsil.jpg');
    Rsil_jagoutz_name = strcat(file_dir, '\', 'output', '_Rsil_jagoutz.jpg');
else
    book_name = strcat(file_dir, '/', 'output', '.xlsx');
    csv_name = strcat(file_dir, '/', 'output', '.csv');
    sil_name = strcat(file_dir, '/', 'output', '_sil.jpg');
    bas_name = strcat(file_dir, '/', 'output', '_bas.jpg');
    sed_name = strcat(file_dir, '/', 'output', '_sed.jpg');
    car_name = strcat(file_dir, '/', 'output', '_car.jpg');
    oc_name = strcat(file_dir, '/', 'output', '_oc.jpg');
    Rsil_jagoutz_name = strcat(file_dir, '/', 'output', '_Rsil_jagoutz.jpg');
end


% Begin to write data into file

success_array = success_array / num_monte; % Transform it into percentage
age_all_plot_book = (AGE_OLD:-t_interval/10^6:AGE_YOUNG)'; % Must transform from 1 row to 1 col

% Title
headers = {'Age (Ma)', 'success percent', 'silicate_weather', 'silw_sd', 'basalt_weather', 'basw_sd', ...
    'sediment_weather', 'sedw_sd', 'carbonate_burial', 'carb_sd', 'oc_burial', 'ocb_sd', 'C_Mg', 'cmg_sd', ...
    'Rsil', 'Rsil_sd', 'Rsil_jagoutz', 'Rsil_jagoutz_sd',  'Rsil_jagoutz_land', 'Rsil_jagoutz_land_sd', 'CO2', 'CO2_err', 'Imb', 'Imb_err', 'F_csiliw_ave', ...
    'F_csiliw_sd', 'F_basw_ave', 'F_basw_sd', 'F_carbw_ave', 'F_carbw_sd', 'F_orgw_ave', 'F_orgw_sd',...
    'F_carbb_ave', 'F_carbb_sd', 'F_orgb_ave', 'F_orgb_sd', 'mafic_total_ave', 'mafic_total_sd', 'mafic_seafloor_total', 'mafic_seafloor_total_sd', 'sil_total_ave', 'sil_total_sd', 'seafloor_ave', 'seafloor_sd'};

% Data
whole_data = [age_all_plot_book, success_array, K_csiliw_ave, K_csiliw_std, K_basw_ave, K_basw_std, K_sediw_ave, ...
    K_sediw_std, K_carbb_ave, K_carbb_std, K_orgb_ave, K_orgb_std, CMg_ave, CMg_std, Rsil_ave, Rsil_std, ...
    Rsil_jagoutz_ave, Rsil_jagoutz_std, Rsil_jagoutz_land_ave, Rsil_jagoutz_land_std, co2_total, co2_total_err, Imb_ave, Imb_std, F_csiliw_ave, F_csiliw_std, ...
    F_basw_ave, F_basw_std, F_carbw_ave, F_carbw_std, F_orgw_ave, F_orgw_std, F_carbb_ave, F_carbb_std, F_orgb_ave, ...
    F_orgb_std, mafic_total_ave, mafic_total_std, mafic_seafloor_total_ave, mafic_seafloor_total_std, sil_total_ave, sil_total_std, seafloor_ave, seafloor_std];



% Use the online verstion to write to csv
csvwrite_with_headers(csv_name, whole_data, headers)

e = cputime-cpu_start;
fprintf('The whole program used %f mins\n', e /60)

% % Write data into excel
% 
% xlswrite(book_name, headers, 'sheet1', 'A1');
% 
% % No need to specify range
% % xlswrite(book_name, whole_data, 'sheet1', strcat('A2:N', num2str(line_end)));
% 
% xlswrite(book_name, whole_data, 'sheet1', 'A2');


% Plot all

% Note, don't reverse age_all_plot here, since when you use fill, it will
% be wrong.

age_all_plot = AGE_OLD:-t_interval/10^6:AGE_YOUNG;

label_fontsize = 14;
axis_fontsize = 12;

% Continental silicate weathering
figure1 = figure;
plot(age_all_plot, K_csiliw_ave, 'r');
hold on

X=[age_all_plot,fliplr(age_all_plot)];
Y_K_csiliw=[K_csiliw_low',fliplr(K_csiliw_high')];
fill(X,Y_K_csiliw,'r', 'faceAlpha', 0.3, 'linestyle', 'none');

xlabel('Age (Ma)', 'FontSize', label_fontsize);
ylabel('Continental silicate weathering', 'FontSize', label_fontsize);
set(gca,'fontsize',axis_fontsize);

% Save figure
saveas(figure1, sil_name)

% Island basalt weathering
figure2 = figure;
plot(age_all_plot, K_basw_ave, 'r');

hold on

X=[age_all_plot, fliplr(age_all_plot)];
Y_K_basw=[K_basw_low',fliplr(K_basw_high')];
fill(X,Y_K_basw,'r', 'faceAlpha', 0.3, 'linestyle', 'none');

xlabel('Age (Ma)', 'FontSize', label_fontsize);
ylabel('Island basalt weathering', 'FontSize', label_fontsize);
set(gca,'fontsize',axis_fontsize);

% Save figure
saveas(figure2, bas_name)

% Continental sediment weathering
figure3 = figure;
plot(age_all_plot, K_sediw_ave, 'r');

hold on

X=[age_all_plot,fliplr(age_all_plot)];
Y_K_sediw=[K_sediw_low',fliplr(K_sediw_high')];
fill(X,Y_K_sediw,'r', 'faceAlpha', 0.3, 'linestyle', 'none');

xlabel('Age (Ma)', 'FontSize', label_fontsize);
ylabel('Continental sediment weathering', 'FontSize', label_fontsize);
set(gca,'fontsize',axis_fontsize);

% Save figure
saveas(figure3, sed_name)


% Burial of carbonate
figure4 = figure;
plot(age_all_plot, K_carbb_ave, 'r');

hold on

X=[age_all_plot,fliplr(age_all_plot)];
Y_K_carbb=[K_carbb_low',fliplr(K_carbb_high')];
fill(X,Y_K_carbb,'r', 'faceAlpha', 0.3, 'linestyle', 'none');

xlabel('Age (Ma)', 'FontSize', label_fontsize);
ylabel('Burial of carbonate', 'FontSize', label_fontsize);
set(gca,'fontsize',axis_fontsize);

% Save figure
saveas(figure4, car_name)


% Burial of organic carbon
figure5 = figure;
plot(age_all_plot, K_orgb_ave, 'r');

hold on

X=[age_all_plot,fliplr(age_all_plot)];
Y_K_orgb=[K_orgb_low',fliplr(K_orgb_high')];
fill(X,Y_K_orgb,'r', 'faceAlpha', 0.3, 'linestyle', 'none');

xlabel('Age (Ma)', 'FontSize', label_fontsize);
ylabel('Burial of organic carbon', 'FontSize', label_fontsize);
set(gca,'fontsize',axis_fontsize);

% Save figure
saveas(figure5, oc_name)


% Seawater Mg concentration
figure6 = figure;
plot(age_all_plot, CMg_ave, 'r');

hold on

X=[age_all_plot,fliplr(age_all_plot)];
Y_CMg=[CMg_low',fliplr(CMg_high')];
fill(X,Y_CMg,'r', 'faceAlpha', 0.3, 'linestyle', 'none');

xlabel('Age (Ma)', 'FontSize', label_fontsize);
ylabel('Seawater Mg concentration (mmol/L)', 'FontSize', label_fontsize);
set(gca,'fontsize',axis_fontsize);
saveas(figure6, cmg_name)

% Success rate
figure7 = figure;
plot(age_all_plot, success_array, 'r');

% Save figure
saveas(figure7, suc_name)

% Rsil changes

figure8 = figure;
plot(age_all_plot, Rsil_ave, 'r');

hold on

X=[age_all_plot,fliplr(age_all_plot)];
Y_Rsil=[Rsil_low',fliplr(Rsil_high')];
fill(X,Y_Rsil,'r', 'faceAlpha', 0.3, 'linestyle', 'none');

xlabel('Age (Ma)', 'FontSize', label_fontsize);
ylabel('Rsil', 'FontSize', label_fontsize);
set(gca,'fontsize',axis_fontsize);

% Save figure
saveas(figure8, Rsil_name)

% Rsil_jagoutz changes

figure9 = figure;
plot(age_all_plot, Rsil_jagoutz_ave, 'r');

hold on

X=[age_all_plot,fliplr(age_all_plot)];
Y_Rsil_jagoutz=[Rsil_jagoutz_low',fliplr(Rsil_jagoutz_high')];
fill(X,Y_Rsil_jagoutz,'r', 'faceAlpha', 0.3, 'linestyle', 'none');

xlabel('Age (Ma)', 'FontSize', label_fontsize);
ylabel('Rsil_jagoutz', 'FontSize', label_fontsize);
set(gca,'fontsize',axis_fontsize);


% Imb changes

figure10 = figure;
plot(age_all_plot, Imb_ave, 'r');

hold on

X=[age_all_plot,fliplr(age_all_plot)];
Y_Imb=[Imb_low',fliplr(Imb_high')];
fill(X,Y_Imb,'r', 'faceAlpha', 0.3, 'linestyle', 'none');

xlabel('Age (Ma)', 'FontSize', label_fontsize);
ylabel('Imb', 'FontSize', label_fontsize);
set(gca,'fontsize',axis_fontsize);

% Save figure
saveas(figure8, Rsil_jagoutz_name)
