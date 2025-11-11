%%%
%%% Checks closure of the surface heat budget in ice shelf cavities. 
%%%

%%% Options
expdir = '../experiments';
% expname = 'hires_seq_onetwentyfourth_notides_RTOPO2';
% expname = 'WC_seq_onethird_notides_RTOPO2_restore';
expname = 'WC_seq_onethird_notides_RTOPO2';
% expname = 'hires_seq_onethird_notides_RTOPO2';
loadexp;



tflux = rdmds(fullfile(resultspath,'TFLUX'),5478);
shifwflx = rdmds(fullfile(resultspath,'SHIfwFlx'),5478);
shihtflx = rdmds(fullfile(resultspath,'SHIhtFlx'),5478);
SIqsw = rdmds(fullfile(resultspath,'SIqsw'),5478);
latentflx = shifwflx*3.34e5; %%% Latent heat flux associated with meltwater input
diffhtflx = -shihtflx; %%% Downward diffusive heat flux from IOBL to resolved ocean
figure(1);pcolor(XC,YC,tflux+SIqsw-latentflx-diffhtflx);shading flat;colorbar;colormap redblue;caxis([-.01 .01]); %%% SHould be exactly zero in ice shelf caviti