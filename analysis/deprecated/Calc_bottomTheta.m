%%%% calculate bottom temperature

%%%%%%%
%%% Plot instananeous Temp

%%% Read experiment data
setExpname
loadexp;

%%% Set true if plotting on a Mac
mac_plots = false;

% deltaT = 480;
nIter0 = 0;
%%% Frequency of diagnostic output


loadexp;
%%% Frequency of diagnostic output
dumpFreq = abs(diag_frequency(15));
deltaT = 440;
% nIter0 = 0;
nDumps = round(nTimeSteps*(deltaT/dumpFreq));
dumpIters = round((1:nDumps)*dumpFreq/deltaT);
dumpIters = dumpIters(dumpIters > nIter0);
nDumps = length(dumpIters);
 
deltaT_4 = 200;
nIter0_4 = 1;
nDumps_4 = round(nTimeSteps*(deltaT_4/dumpFreq));
dumpIters_4 = round((1:nDumps_4)*dumpFreq/deltaT_4); 
dumpIters_4 = dumpIters_4(dumpIters_4 >= nIter0_4);



deltaT_5 = 300;
nIter0_5 = 1;
nDumps_5 = round(nTimeSteps*10*(deltaT_5/dumpFreq));
dumpIters_5 = round((1:nDumps_5)*dumpFreq/deltaT_5); 
dumpIters_5 = dumpIters_5(dumpIters_5 >= nIter0_5);

deltaT2 = 400;
nDumps2 = round(nTimeSteps*10*(deltaT2/dumpFreq));
dumpIters_2 = round((1:nDumps2)*dumpFreq/deltaT2); 
dumpIters2 = dumpIters_2(dumpIters_2 >= nIter0);
nDumps2 = length(dumpIters2);



%%% Mesh grids for plotting
kmax = ones(Nx,Ny);
kmin = ones(Nx,Ny);
for i=1:Nx
  for j=1:Ny
    idx = find(squeeze(hFacC(i,j,:))>0);
    if (~isempty(idx))
      kmin(i,j) = min(idx);
      kmax(i,j) = max(idx);
    end
  end
end
kn = ones(Nx,Ny);
kp= ones(Nx,Ny);
wn = 0.5*ones(Nx,Ny);
wp = 0.5*ones(Nx,Ny);
for i=1:Nx
  for j=1:Ny
    if (sum(hFacC(i,j,:),3)==0)
      continue;
    end
    zmid = 0.5 * (SHELFICEtopo(i,j) + bathy(i,j));
    kmid = max(find(squeeze(zz)>zmid));
    if (isempty(kmid))
      continue;
    end
    kp(i,j) = kmid;
    kn(i,j) = kp(i,j) + 1;
    wp(i,j) = (zmid-zz(kn(i,j))) / (zz(kp(i,j))-zz(kn(i,j)));
    wn(i,j) = 1 - wp(i,j);
  end
end


%%% Longitudinal positions

XC;


%%% Latitudinal positions

YC;







%%% Calculate time-averaged velocity
tmin = 18*86400*360;
tmax = 26*86400*360;
tmin2 = 12*86400*360;
tmax2 = 21*86400*360;
exppath1 = '/data3/MITgcm_WS/experiments/n_34452';
exppath2 = '/data3/MITgcm_WS/experiments/a_34_20boundary';
Tc = readIters(exppath1,'THETA',dumpIters,deltaT,tmin,tmax,Nx,Ny,Nr);  %%% temperature, control cold
Tw = readIters(exppath2,'THETA',dumpIters,deltaT,tmin2,tmax2,Nx,Ny,Nr);%%%% temperature, control warm



Tc2=NaN(Nx,Ny);
Tw2=NaN(Nx,Ny);
for i = 1:Nx
    for j = 1:Ny
        Tc2(i,j)=Tc(i,j,kmax(i,j));
        Tw2(i,j)=Tw(i,j,kmax(i,j));
    end
end

Tc2(SHELFICEtopo-bathy==0)=NaN;
Tw2(SHELFICEtopo-bathy==0)=NaN;
bathy(SHELFICEtopo-bathy==0)=NaN;

topog_msk = ones(Nx,Ny);

for i=1:Nx
    for j=1:Ny        
           
           if (((SHELFICEtopo(i,j)) - (bathy(i,j)) < 0  ||  (bathy(i,j) == SHELFICEtopo(i,j)) ))
               topog_msk(i,j) = 0;
           end
    end
end
  


