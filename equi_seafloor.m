% this is to create the symbolic equations for Model 3 with seafloor alteration
% integrated

%----------------------------------------------------------------------
clear all
close all

cpu_start_equi = cputime;

% ------------------------------Construct the parameters and equations


% This is the parameter for the whole run
syms Sp dC_carbb dC_orgb Sr_ocean dSr_ocean Os_ocean dOs_ocean K_csiliw K_basw K_sediw K_carbb ...
    K_orgb der_dSr_ocean der_dOs_ocean V_ocean CMg_plei CMg CMg_old dt...
    C_carbw C_orgw C_ma C_plume C_carbb C_orgb dC_carbw dC_orgw dC_ma dC_plume A_csiliw A_basw ...
    A_carbw A_carbb A_rev Sr_csiliw dSr_csiliw Sr_basw dSr_basw Sr_carbw dSr_carbw ...
    Sr_hy dSr_hy Sr_dia dSr_dia Os_csiliw dOs_csiliw Os_basw dOs_basw Os_sediw dOs_sediw ...
    Os_hy dOs_hy Os_cos dOs_cos Mg_csiliw Mg_basw Mg_carbw Mg_hy A_Jagoute Sr_mafic dSr_mafic ...
    Sr_ultra dSr_ultra Os_mafic dOs_mafic Os_ultra dOs_ultra Mg_mafic Mg_ultra

% This is the extra parameter for the modern solution (The ones need to be
% quantified)

syms Sr_basw_hy dSr_mantle Os_basw_hy_cos dOs_mantle

syms C_seafloor dC_seafloor factor_T

% -------------------------------------This is for whole run

f_C = K_sediw * C_carbw + K_sediw * C_orgw + Sp * C_ma + C_plume - K_carbb * C_carbb - K_orgb * C_orgb - C_seafloor == 0;

f_dC = K_sediw * C_carbw * dC_carbw + K_sediw * C_orgw * dC_orgw + Sp * C_ma * dC_ma + ...
    C_plume * dC_plume - K_carbb * C_carbb * dC_carbb - K_orgb * C_orgb * dC_orgb - C_seafloor * dC_seafloor == 0;

f_A = K_csiliw * A_csiliw + K_basw * A_basw + K_sediw * A_carbw + A_Jagoute - K_carbb * A_carbb - Sp * CMg / CMg_plei * A_rev == 0;

f_dSr = K_csiliw * Sr_csiliw * (dSr_csiliw - dSr_ocean) + K_basw * Sr_basw *(dSr_basw - dSr_ocean) + ...
    K_sediw * Sr_carbw * (dSr_carbw - dSr_ocean) + factor_T * Sp * CMg / CMg_plei * Sr_hy * (dSr_hy - dSr_ocean) + Sr_dia * ...
    (dSr_dia - dSr_ocean) + Sr_mafic * (dSr_mafic - dSr_ocean) + Sr_ultra * (dSr_ultra - dSr_ocean) - Sr_ocean * der_dSr_ocean == 0;

f_dOs = K_csiliw * Os_csiliw * (dOs_csiliw - dOs_ocean) + K_basw * Os_basw *(dOs_basw - dOs_ocean) + ...
    K_sediw * Os_sediw * (dOs_sediw - dOs_ocean) + factor_T * Sp * CMg / CMg_plei * Os_hy * (dOs_hy - dOs_ocean) + Os_cos * ...
    (dOs_cos - dOs_ocean) + Os_mafic * (dOs_mafic - dOs_ocean) + Os_ultra * (dOs_ultra - dOs_ocean) - Os_ocean * der_dOs_ocean == 0;

f_CMg = K_csiliw * Mg_csiliw + K_basw * Mg_basw + K_sediw * Mg_carbw + Mg_mafic + Mg_ultra - factor_T * Sp * CMg / CMg_plei * ...
    Mg_hy - V_ocean * (CMg - CMg_old) / dt == 0;

sol = solve([f_C, f_dC, f_A, f_dSr, f_dOs, f_CMg], K_csiliw, K_basw, K_sediw, K_carbb, K_orgb, CMg);


% Use simplify and symbol to function
F_K_csiliw = matlabFunction(simplify(sol.K_csiliw));
F_K_basw = matlabFunction(simplify(sol.K_basw));
F_K_sediw = matlabFunction(simplify(sol.K_sediw));
F_K_carbb = matlabFunction(simplify(sol.K_carbb));
F_K_orgb = matlabFunction(simplify(sol.K_orgb));
F_CMg = matlabFunction(simplify(sol.CMg));

% -------------------------------------This is for modern values

% Back calculated from equation 6-8


f_C_plei = C_carbw + C_orgw + C_ma + C_plume - C_carbb - C_orgb - C_seafloor == 0;

f_dC_plei = C_carbw * dC_carbw + C_orgw * dC_orgw + C_ma * dC_ma + ...
    C_plume * dC_plume - C_carbb * dC_carbb - C_orgb * dC_orgb - C_seafloor * dC_seafloor== 0;

f_A_plei = A_csiliw + A_basw + A_carbw - 2 * C_carbb - A_rev == 0;

f_dSr_plei = Sr_csiliw * (dSr_csiliw - dSr_ocean) + Sr_carbw * (dSr_carbw - dSr_ocean) + ...
    Sr_basw_hy *(dSr_mantle - dSr_ocean) + Sr_dia * (dSr_dia - dSr_ocean) - Sr_ocean * der_dSr_ocean == 0;

f_dOs_plei = Os_csiliw * (dOs_csiliw - dOs_ocean) + Os_sediw * (dOs_sediw - dOs_ocean) + ...
    Os_basw_hy_cos *(dOs_mantle - dOs_ocean) - Os_ocean * der_dOs_ocean == 0;


sol_plei = solve([f_C_plei, f_dC_plei, f_A_plei, f_dSr_plei, f_dOs_plei], C_carbb, C_orgb, A_rev, Sr_basw_hy, Os_basw_hy_cos);

% Use simplify and symbol to function
F_C_carbb = matlabFunction(simplify(sol_plei.C_carbb));
F_C_orgb = matlabFunction(simplify(sol_plei.C_orgb));
F_A_rev = matlabFunction(simplify(sol_plei.A_rev));
F_Sr_basw_hy = matlabFunction(simplify(sol_plei.Sr_basw_hy));
F_Os_basw_hy_cos = matlabFunction(simplify(sol_plei.Os_basw_hy_cos));

filename = 'equi_seafloor.mat';
save(filename)

% ------------------------------------
% Print finishing the equation solving.
% ------------------------------------

fprintf('Finish the equation solving\n')

e = cputime-cpu_start_equi;
fprintf('Total equation solving time is %f...\n', e /60)