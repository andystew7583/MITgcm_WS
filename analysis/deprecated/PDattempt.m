%%%%% Attempt to plot Potential Density ;


%%% Read experiment data
loadexp;

%%% Frequency of diagnostic output
dumpFreq = abs(diag_frequency(15));
nDumps = round(nTimeSteps*deltaT/dumpFreq);
dumpIters = round((1:nDumps)*dumpFreq/deltaT);
dumpIters = dumpIters(dumpIters > nIter0);
nDumps = length(dumpIters);

%%%%% Load Salinity and Temperature

tmin = 8*86400*365;
tmax = 9*86400*365;
theta = readIters(exppath,'THETA',dumpIters,deltaT,tmin,tmax,Nx,Ny,Nr);
salt = readIters(exppath,'SALT',dumpIters,deltaT,tmin,tmax,Nx,Ny,Nr);




tp = zeros(Nx,Ny,Nr);

%%% Calculate potential density
pd = densmdjwf(salt,theta,tp) - 1000;





[YZ,ZY]=meshgrid(yy,zz);

temp = squeeze(pd(end,:,:));

pcolor(YZ,ZY,temp')

