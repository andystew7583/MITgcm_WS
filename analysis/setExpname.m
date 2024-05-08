%%%
%%% setExpname.m
%%%
%%% Convenience script to set the experiment name before loadexp.m is 
%%% called in other scripts.
%%%
% gendir = '..';
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
expname = 'hires_seq_onethird_RTOPO2';
% expname = 'hires_seq_onethird';
% expname = 'hires_seq_onethird_notides_RTOPO2';
% expname = 'hires_seq_onesixth_modbdrystrat';
% expname = 'hires_seq_onesixth_notides';
% expname = 'hires_seq_onesixth_RTOPO2';
% expname = 'hires_seq_onetwelfth_notides_RTOPO2';
% % expname = 'hires_seq_onetwentyfourth_RTOPO2';
% expname = 'hires_seq_onetwentyfourth_notides_RTOPO2';

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

