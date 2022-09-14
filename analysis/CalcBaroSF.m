%%%
%%% Calculate and plot time-mean BSF
%%%

%%% Options
expdir = '../experiments';
% expname = 'hires_seq_onethird_RTOPO2';
% tmin = 18.05*86400*365;
% tmax = 27.05*86400*365;
% expname = 'hires_seq_onesixth_notides_RTOPO2';
% tmin = 9.05*86400*365;
% tmax = 18.05*86400*365;
% expname = 'hires_seq_onetwelfth_RTOPO2';
% tmin = 1.05*86400*365;
% tmax = 9.05*86400*365;
expname = 'hires_seq_onetwentyfourth_notides_RTOPO2';
tmin = 1.05*86400*365;
tmax = 4.05*86400*365;

%%% Load velocity
loadexp;

%%% Frequency of diagnostic output - should match that specified in
%%% data.diagnostics.
dumpFreq = abs(diag_frequency(1));
nDumps = round(endTime/dumpFreq);
dumpIters = round((1:nDumps)*dumpFreq/deltaT);
dumpIters = dumpIters(dumpIters > nIter0);
nDumps = length(dumpIters);

%%% Calculate time-averaged velocity
uu = readIters(exppath,'UVEL',dumpIters,deltaT,tmin,tmax,Nx,Ny,Nr);
vv = readIters(exppath,'VVEL',dumpIters,deltaT,tmin,tmax,Nx,Ny,Nr);

%%% Grid spacing matrices
DX = repmat(DXG,[1 1 Nr]);
DY = repmat(DYG,[1 1 Nr]);
DZ = repmat(reshape(DRF,[1 1 Nr]),[Nx Ny 1]);

%%% Calculate depth-averaged zonal velocity
UU = sum(uu.*DZ.*hFacW,3);

%%% Calculate barotropic streamfunction
Psi = zeros(Nx+1,Ny+1);
Psi(2:Nx+1,2:Ny+1) = -cumsum(UU.*DYG,2);
Psi = Psi(1:Nx,1:Ny);




save(fullfile('products',[expname,'_PsiBT.mat']),'Psi','uu','vv','-v7.3');

