
%%%
%%% setExpname.m
%%%
%%% Convenience script to set the experiment name before loadexp.m is 
%%% called in other scripts.
%%%
gendir = '..';
% gendir = '/Volumes/Stewart-RAID1-B/UCLA/Julia/MITgcm_WS';
basedir = fullfile(gendir,'experiments');

expdir = basedir;
% expname = 'hires_seq_onethird';
% expname = 'hires_seq_onethird_strongwinds';
% expname = 'hires_seq_onethird_loGamT';
% expname = 'hires_seq_onethird_ECCObdry';
% expname = 'hires_seq_onethird_modbdry';
% expname = 'hires_seq_onethird_pyc1000';
% expname = 'hires_seq_onethird_SIbdryfix'a;
% expname = 'hires_seq_onethird_SIbdryfix_pyc1000_fixSIdrag';
% expname = 'hires_seq_onethird_SIbdryfix_pyc1000_fixSIdrag_fixSIsponge';
% expname = 'hires_seq_onethird_albFix_ERArad_pyc500';
% expname = 'hires_seq_onethird_noslip_loGamT';
% expname = 'hires_seq_onesixth';
% expname = 'hires_seq_onesixth_restart_hifreq';
% expname = 'hires_seq_onesixth_modbdries';
% expname = 'hires_seq_onesixth_albFix_ERArad_pyc500';
% expname = 'hires_seq_onesixth_loDeltaT';
% expname = 'hires_seq_onesixth_DeltaT60_2';
% expname = 'hires_seq_onesixth_fixNorthBdy_2';
% expname = 'hires_seq_onesixth_smoothtopog_tildes_6';
% expname = 'hires_seq_onesixth_iceflux2_iceAHconst_gammaFrict';
% expname = 'hires_seq_onesixth_iceflux2_iceAHconst_gammaFrict_ERAsw_tempoff';
% expname = 'hires_seq_onesixth_modbdrystrat_wind1.5';
% expname = 'hires_seq_onesixth_bdrysalt34';


% expname = 'hires_seq_onethird_RTOPO2_hifreq';
% expname = 'hires_seq_onethird_RTOPO2';
% expname = 'hires_seq_onethird';
% expname = 'hires_seq_onethird_notides_RTOPO2';
% expname = 'hires_seq_onesixth_modbdrystrat';
% expname = 'hires_seq_onesixth_notides';
% expname = 'hires_seq_onesixth_RTOPO2';
% expname = 'hires_seq_onetwelfth_notides_RTOPO2';
% % expname = 'hires_seq_onetwentyfourth_RTOPO2';
% expname = 'hires_seq_onetwentyfourth_notides_RTOPO2';
% expname = 'hires_seq_onetwentyfourth_notides_RTOPO2_SSH';

% expname = 'hires_nest_onethirtieth_notides_RTOPO2';

% expname = 'hires_nest_onethirtysecond_notides_RTOPO2';
% 
% expname = 'WC_seq_onethird_notides_RTOPO2';

% expname = 'WC_seq_onethird_RTOPO2_unmodEB_linHtFlx';
% expname = 'WC_seq_onethird_RTOPO2_unmodEB_alphaV0.5';

% expname = 'WC_seq_onethird_RTOPO2_unmodEB';
% expname = 'WC_seq_onethird_RTOPO2_unmodEB_alphaV0.5';
% expname = 'WC_seq_onethird_RTOPO2_unmodEB_alphaV0.5_restore34.7';

% expname = 'WC_seq_onethird_RTOPO2_dpyc200_alphaV0.5';
% expname = 'WC_seq_onethird_RTOPO2_dpyc200_SconstBCs2';

% expname = 'WC_seq_onethird_RTOPO2_dpyc100';
% expname = 'WC_seq_onethird_RTOPO2_dpyc100_alphaV0.5';
% expname = 'WC_seq_onethird_RTOPO2_dpyc100_alphaV0.25';
% % expname = 'WC_seq_onethird_RTOPO2_dpyc100_strat2e-5';
% expname = 'WC_seq_onethird_RTOPO2_dpyc100_Smin34';
% expname = 'WC_seq_onethird_RTOPO2_dpyc100_strat4e-5';
% expname = 'WC_seq_onethird_RTOPO2_dpyc100_gAlphaU0.5_gAlphaV0.5';
% expname = 'WC_seq_onethird_RTOPO2_dpyc100_offshorestrat2e-5';
% expname = 'WC_seq_onethird_RTOPO2_dpyc100_strat2e-5_alphaV0.5';
% expname = 'WC_seq_onethird_RTOPO2_dpyc0_alphaU0.25_alphaV0.25';
% expname = 'WC_seq_onethird_RTOPO2_dpyc100_alphaU0.25_alphaV0.25';
% expname = 'WC_seq_onethird_RTOPO2_dpyc100_Sf34.3_offshorestrat4e-6';
% expname = 'WC_seq_onethird_RTOPO2_dpyc100_Smin34_Sf34.3';
% expname = 'WC_seq_onethird_RTOPO2_dpyc0_Sf34.3_offshorestrat6e-6';
% expname = 'WC_seq_onethird_RTOPO2_dpyc100_Sf34.3_offshorestrat5e-6';
% expname = 'WC_seq_onethird_RTOPO2_dpyc200_Sf34.3_offshorestrat4e-6';

% expname = 'WC_onethird_ref';
% expname = 'WC_onethird_dpyc-150_strat1e-4';
% expname = 'WC_onethird_dpyc-150_strat2e-4';
% expname = 'WC_onethird_dpyc-150_strat2.25e-4';
% expname = 'WC_onethird_dpyc-150_strat2.5e-4';;
% expname = 'WC_onethird_dpyc-150_strat2.75e-4';
% expname = 'WC_onethird_dpyc-150_strat3e-4';
% expname = 'WC_onethird_dpyc-150_strat4e-4';
% expname = 'WC_onethird_dpyc-150_strat1e-4_ws';
% expname = 'WC_onethird_dpyc-150_strat1e-4_aU0.5_aV0.5_ws';
% expname = 'WC_onethird_dpyc-150_strat3e-4_aU0.5_aV0.5';
% expname = 'WC_onethird_dpyc-150_strat3e-4_bU0.5_bV0.5';

% expname = 'WC_onethird_dpyc-100_strat1e-4';
% expname = 'WC_onethird_dpyc-100_strat1.5e-4';
% expname = 'WC_onethird_dpyc-100_strat1.75e-4';
% expname = 'WC_onethird_dpyc-100_strat1.875e-4';
% expname = 'WC_onethird_dpyc-100_strat1.5e-4';
% expname = 'WC_onethird_dpyc-100_strat2e-4'

% expname = 'WC_onethird_dpyc-50_strat1e-4';
% expname = 'WC_onethird_dpyc-50_strat1.25e-4';
% expname = 'WC_onethird_dpyc-50_strat1.375e-4';
% expname = 'WC_onethird_dpyc-50_strat1.5e-4';
% expname = 'WC_onethird_dpyc-50_strat2e-4';

% expname = 'WC_onethird_dpyc0_strat5e-5';
% expname = 'WC_onethird_dpyc0_strat7.5e-5';
% expname = 'WC_onethird_dpyc0_strat8.75e-5';
% expname = 'WC_onethird_dpyc0_strat1e-4';
% expname = 'WC_onethird_dpyc0_strat1.125e-4';
% expname = 'WC_onethird_dpyc0_strat1.25e-4';
% expname = 'WC_onethird_dpyc0_strat1.5e-4';


% expname = 'WC_onethird_dpyc50_strat6.5e-5';
% expname = 'WC_onethird_dpyc50_strat7.5e-5';
% expname = 'WC_onethird_dpyc50_strat8e-5';
% expname = 'WC_onethird_dpyc50_strat1e-4';

% expname = 'WC_onethird_dpyc-100_strat1e-4';

% expname = 'WC_onethird_strat7e-5';
% expname = 'WC_onethird_strat8e-5';
% expname = 'WC_onethird_strat1e-4';
% expname = 'WC_onethird_strat4e-5';
% expname = 'WC_onethird_strat6e-5';
% expname = 'WC_onethird_strat6e-5';

% expname = 'WC_onethird_dpyc150_strat4e-5';
% expname = 'WC_onethird_dpyc150_strat5e-5';
% expname = 'WC_onethird_dpyc150_strat6e-5';

% expname = 'WC_onethird_dpyc200_strat4e-5';
% expname = 'WC_onethird_dpyc200_strat5e-5';
% expname = 'WC_onethird_dpyc200_strat6e-5';

% expname = 'WC_onethird_dpyc250_strat4e-5';
% expname = 'WC_onethird_dpyc250_strat5e-5';

% expname = 'WC_onethird_dpyc300_strat1e-5';
% expname = 'WC_onethird_dpyc300_strat3e-5';
expname = 'WC_onethird_dpyc300_strat5e-5';

% expname = 'WC_onethird_aU0.7_aV0.7_bU0.7_bV0.7_dT5_noA23';
% expname = 'WC_onethird_dpyc-150_strat1.5e-4_aU0.7_aV0.7_bU0.7_bV0.7_dT5_noA23';
% expname = 'WC_onethird_dpyc-150_strat1.75e-4_aU0.7_aV0.7_bU0.7_bV0.7_dT5_noA23';
% expname = 'WC_onethird_dpyc-150_strat2e-4_aU0.7_aV0.7_bU0.7_bV0.7_dT5_noA23';
% expname = 'WC_onethird_aU0.85_aV0.85_bU0.85_bV0.85_dT2.5_noA23';
% expname = 'WC_onethird_dpyc-150_strat2e-4_aU0.85_aV0.85_bU0.85_bV0.85_dT2.5_noA23';
% expname = 'WC_onethird_dpyc-150_strat2.25e-4_aU0.85_aV0.85_bU0.85_bV0.85_dT2.5_noA23';
% expname = 'WC_onethird_dpyc-150_strat2.5e-4_aU0.85_aV0.85_bU0.85_bV0.85_dT2.5_noA23';
% 

% expname = 'WC_onethird_aU0.55_aV0.55_bU0.55_bV0.55_dT7.5_noA23';
% expname = 'WC_onethird_aU0.4_aV0.4_bU0.4_bV0.4_dT10_noA23';
% expname = 'WC_onethird_aU0.5_aV0.5';
% expname = 'WC_onethird_dpyc-150_strat1e-4_aU0.5_aV0.5';

%%% Julia's control experiments
% expname = 'n_34452';

%%% wind experiments
% 
% expname = 'n_34452';
% expname = 'n_342';
% expname = 'a_34_20boundary';
% expname = 'a_3445_20boundary';

 
% expname = 's_349_20boundary';
%%%%% run above dies at 3750 days, going back to 300 seconds


% expname = 's_hist1_meridplus20';

% expname = 's_hist2_nowind';

% expname = 's_hist3_30colder';

% expname = 's_hist4_warm_meridplus20';
% 
% expname = 's_hist5_warm_windsplus20_30colder';

% expname = 's_hist6_warm_coldminus30';

% expname = 's_hist7_salty_nowind';

% expname = 's_hist8_nozonalwind';

% expname = 's_hist11_control_30warmer';

% expname = 's_hist12_warm_increasezonal20';

% expname = 's_hist13_control_meridminus20';
% 
% expname = 's_hist14_warm_meridplus10';
% 
% expname = 's_hist15_control_nomeridwind';
%  
% expname = 's_hist20_warm_windsplus5';
% 

% expname = 's_hist28_c_windsminus5';


% expname = 's_hist24_w_mplus15_2';  %(see above)
 
% expname = 's_hist31_c_meridplus10';

% expname = 's_hist32_c_meridplus5';

% expname = 's_hist32_w_meridminus10';

% expname = 's_hist32_w_meridminus20';

% expname = 's_hist32_w_meridminus5';

% expname = 's_hist19_cminus2_2';

% expname = 's_hist34_cminus5_2';
% 
% expname = 's_hist35_cminus10_2';

% expname = 's_hist36_c_minus20_2';

% expname = 's_hist37_c_minus15_2';

% expname = 's_hist38_c_minus50';

% expname = 's_hist39_c_0';

% expname = 's_hist41_c_minus30';

% expname = 's_hist40_w_plus10';

% expname = 'w_hist42_wminus50';

% expname = 's_hist43_w_meridplus5';

% expname = 's_hist44_w_meridminus30';
% 
% expname = 's_hist45_c_meridplus15';

% expname = 's_hist41_c_minus40';

% expname = 's_hist46_w_minus40';

% expname = 's_hist47_c_minus45';

% expname = 'z_warm_zonalplus10';

% expname = 'z_warm_zonalplus20';

% expname = 'z_control_zonalminus20';

% expname = 'z_control_zonalminus50';

% expname = 't_warm_tempminus5';

% expname = 't_warm_tempminus10';

% expname = 't_control_tempplus5';

% expname  = 't_control_tempplus10';

% expname = 'n_349';

% expname = 'minus30_conth';

%%%Triple check control runs 
% expname='n_342';

% expname = 'z_minus20_control';
% expname = 't_warm_tempminus5_h2';


% expname = 'control_plus10A';

