%%%
%%% evaluateTides.m
%%%
%%% Computes SSH variance and evaluate against obs.
%%%

%%% Options
expdir = '../experiments';
expname = 'hires_seq_onethird_RTOPO2_hifreq';
tmin = 1774980*480;
tmax = tmin+30*t1day;

%%% Load experiment
loadexp;

%%% Frequency of diagnostic output - should match that specified in
%%% data.diagnostics.
dumpFreq = abs(diag_frequency(end));
nDumps = round(endTime/dumpFreq);
dumpIters = round((1:nDumps)*dumpFreq/deltaT);
dumpIters = dumpIters(dumpIters > nIter0);
nDumps = length(dumpIters);







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PRE-DETERMINE ITERATION NUMBERS %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Determine iteration numbers to process
itersToRead = [];
times = [];
for n=1:length(dumpIters)
 
  tsecs = dumpIters(n)*deltaT;
 
  if ((tsecs >= tmin) && (tsecs < tmax))    
    itersToRead = [itersToRead dumpIters(n)];
    times = [times dumpIters(n)*deltaT];
  end
  
end
Ntime = length(itersToRead);












%%%%%%%%%%%%%%%%%%%%%%%%
%%% COMPUTE PRODUCTS %%%
%%%%%%%%%%%%%%%%%%%%%%%%

%%% Calculate time-averaged SSH and SSH variance
etan_tavg = zeros(Nx,Ny);
etansq_tavg = zeros(Nx,Ny);
for n=1:Ntime
 
  %%% Print current time to keep track of calculation
  tyears = itersToRead(n)*deltaT/86400/365; 
  [num2str(tyears) num2str(itersToRead(n))] 
  
  
  %%% Read SSH
  etan = rdmdsWrapper(fullfile(exppath,'/results/ETAN_inst'),itersToRead(n));
  if (isempty(etan))
    ['Ran out of data at n=',num2str(n),'/',num2str(nDumps),' t=',num2str(tyears),' days.']
    break;
  end
  
  %%% Add to time averages
  etan_tavg = etan_tavg + etan/Ntime;
  etansq_tavg = etansq_tavg + etan.^2/Ntime;

end

%%% Variance
etansq_fluc = etansq_tavg - etan_tavg.^2;







%%%%%%%%%%%%%%%%%%
%%% MAKE PLOTS %%%
%%%%%%%%%%%%%%%%%%


figure(48);
pcolor(XC,YC,etansq_fluc);
shading interp;
colormap haxby;