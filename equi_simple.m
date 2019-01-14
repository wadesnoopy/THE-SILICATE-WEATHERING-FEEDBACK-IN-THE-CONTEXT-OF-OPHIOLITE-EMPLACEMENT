% This is to create the symbolic equations for Model 1 and 2

% flux mol/yr and ratio for Os

flux_ht = 10.5;
ratio_ht = 0.26;

flux_lt = 105.1;
ratio_lt = 0.11;


flux_cosm = 52.6;
ratio_cosm = 0.127;

flux_dust = 36.8;	
ratio_dust = 1.05;                   

flux_river = 1800;
ratio_river = 0.8;

ratio_ht_norm = ratio_ht / (7.4 + ratio_ht);
ratio_lt_norm = ratio_lt / (7.4 + ratio_lt);
ratio_cosm_norm = ratio_cosm / (7.4 + ratio_cosm);
ratio_dust_norm = ratio_dust / (7.4 + ratio_dust);
ratio_river_norm = ratio_river / (7.4 + ratio_river);

Os_ocean = 7.20 * 10^7;


% flux mol/yr and ratio for Sr

flux_dia_Sr = 3.4e9;
ratio_dia_Sr = 0.708;
ratio_dia_Sr_norm = ratio_dia_Sr / (ratio_dia_Sr + 9.43);

flux_hydro_Sr = 8.4e9;
ratio_hydro_Sr = 0.7025;
ratio_hydro_Sr_norm = ratio_hydro_Sr / (ratio_hydro_Sr + 9.43);

flux_river_Sr = 31e9;
ratio_river_Sr = 0.7103;
ratio_river_Sr_norm = ratio_river_Sr / (ratio_river_Sr + 9.43);

Sr_ocean = 1.25e17;

syms Sr_mafic dSr_mafic Sr_ultra dSr_ultra dSr_ocean der_dSr_ocean K_total_Sr


syms K_total Os_mafic dOs_mafic Os_ultra dOs_ultra dOs_ocean der_dOs_ocean ratio Os_ophio dOs_ophio K_ophio


% This is for norm Os
f_dOs_norm = flux_ht * (ratio_ht_norm - dOs_ocean) + flux_lt * (ratio_lt_norm - dOs_ocean) + ...
    flux_cosm * (ratio_cosm_norm - dOs_ocean) + flux_dust * (ratio_dust_norm - dOs_ocean) + ...
    K_total * flux_river * (ratio_river_norm - dOs_ocean) + Os_mafic * (dOs_mafic - dOs_ocean) + Os_ultra * (dOs_ultra - dOs_ocean) - Os_ocean * der_dOs_ocean == 0;

sol_norm = solve(f_dOs_norm, K_total);

% Use simplify and symbol to function
F_K_total_norm = matlabFunction(simplify(sol_norm));


% This is for norm Sr
f_dSr_norm = flux_hydro_Sr * (ratio_hydro_Sr_norm - dSr_ocean) + flux_dia_Sr * (ratio_dia_Sr_norm - dSr_ocean) + ...
    K_total_Sr * flux_river_Sr * (ratio_river_Sr_norm - dSr_ocean) + Sr_mafic * (dSr_mafic - dSr_ocean) + Sr_ultra * (dSr_ultra - dSr_ocean) - Sr_ocean * der_dSr_ocean == 0;

sol_Sr_norm = solve(f_dSr_norm, K_total_Sr);

% Use simplify and symbol to function
F_K_total_Sr_norm = matlabFunction(simplify(sol_Sr_norm));



filename = 'equi_simple.mat';
save(filename)