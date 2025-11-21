%%%
%%% calcFRISMeltTimeSeries.m
%%%
%%% Calculates a time series of total FRIS melt rates for one of our
%%% experiments.
%%%

function [tt,SHImelt,SIprod,SHImelt_mean,XC,YC,bathy,SHELFICEtopo,SHImelt_eff,SHImelt_tend,SHImelt_diff,SHImelt_tflux,SHImelt_conv] = calcFRISMeltTimeSeries (expdir,expname,compute_eff_melt) 

  %%% This needs to be set to ensure we are using the correct output
  %%% frequency
  loadexp;
  diagfreq = diag_frequency(end);


  Cp = 3.994e3; % J/kg/K
  Lf = 3.34e5; % J/kg;
  Sref = 34.5; % g/kg

  %%% Frequency of diagnostic output
  diagnum = 1; %%% Should be arbitrary because all output is provided at the same frequency
  dumpFreq = abs(diag_frequency(diagnum));
  nDumps = round(nTimeSteps*deltaT/dumpFreq);
  dumpIters = round(nIter0+(1:nDumps)*dumpFreq/deltaT);
  dumpIters = dumpIters(dumpIters > nIter0);

  %%% To store the result
  tt = zeros(1,nDumps);
  SHImelt = NaN*ones(1,nDumps);
  SIprod = NaN*ones(1,nDumps);
  SHImelt_eff = NaN*ones(1,nDumps);
  SHImelt_conv = NaN*ones(1,nDumps);
  SHImelt_tend = NaN*ones(1,nDumps);
  SHImelt_diff = NaN*ones(1,nDumps);
  SHImelt_tflux = NaN*ones(1,nDumps);
  heat_tot = NaN*ones(1,nDumps);
  SHImelt_mean = zeros(Nx,Ny);
  tlen = 0;
  
  endTime = dumpIters(end)*deltaT
  avg_len = 1*t1year

  %%% Indices over which to integrate, i.e. defining the FRIS
  xidx = find(XC(:,1)<-29.9);
  yidx = find(YC(1,:)<-74.5);

  %%% Mask for calculating sea ice production
  ETA = defineMOCgrid(XC,YC,[],[],false,false);
  msk_SIprod = (ETA < 3.5) & (hFacC(:,:,1)>0);

  %%% TODO need to fix this
  for n=1:length(dumpIters)

    tt(n) =  dumpIters(n)*deltaT;
    tt(n)

    %%% Attempt to load melt ave per month
    SHIfwFlx=rdmdsWrapper(fullfile(exppath,'/results/SHIfwFlx'),dumpIters(n));  
    SHIhtFlx=rdmdsWrapper(fullfile(exppath,'/results/SHIhtFlx'),dumpIters(n));  
    SFLUX=rdmdsWrapper(fullfile(exppath,'/results/SFLUX'),dumpIters(n));  
    TFLUX=rdmdsWrapper(fullfile(exppath,'/results/TFLUX'),dumpIters(n));  


    if (isempty(SHIfwFlx))
      continue;
    end
    
    
    %%% Mean local melt rate
    if (tt(n) > endTime-avg_len)
      
      SHImelt_mean = SHImelt_mean + SHIfwFlx;

      %%% Increment counter
      tlen = tlen + 1;
      
    end    
    
n

    %%% Compute area-integrated freshwater flux due to melt
    SHImelt(n) =  sum(sum(SHIfwFlx(xidx,yidx).*RAC(xidx,yidx)));

    %%% Compute total sea ice production equivalent
    SIprod(n) = sum(sum(msk_SIprod.*SFLUX/Sref*rhoConst/rhoShelf.*RAC));

    if (compute_eff_melt)

      SHIhtFlx=rdmdsWrapper(fullfile(exppath,'/results/SHIhtFlx'),dumpIters(n));  
      TFLUX=rdmdsWrapper(fullfile(exppath,'/results/TFLUX'),dumpIters(n));  
  
      
  
      if (isempty(TFLUX) || isempty(SHIhtFlx))
        continue;
      end
  
      %%% Ice shelf mask
      msk = hFacC(:,:,1)==0;
  

  
      SHImelt_diff(n) =  sum(sum(SHIhtFlx(xidx,yidx).*RAC(xidx,yidx)))/Lf;
      SHImelt_tflux(n) =  sum(sum(TFLUX(xidx,yidx).*RAC(xidx,yidx).*msk(xidx,yidx)))/Lf;


      UVELTH=rdmdsWrapper(fullfile(exppath,'/results/UVELTH'),dumpIters(n)); 
      VVELTH=rdmdsWrapper(fullfile(exppath,'/results/VVELTH'),dumpIters(n)); 
      THETA=rdmdsWrapper(fullfile(exppath,'/results/THETA'),dumpIters(n)); 
      if (isempty(UVELTH) || isempty(VVELTH) || isempty(THETA))
        continue;
      end

      %%% Compute heat flux divergence
      UVELTH_zint = sum(UVELTH.*DRF.*hFacW,3).*DYG;
      VVELTH_zint = sum(VVELTH.*DRF.*hFacS,3).*DXG;
      div_htflx = UVELTH_zint([2:Nx 1],:) - UVELTH_zint(1:Nx,:) + VVELTH_zint(:,[2:Ny,1]) - VVELTH_zint(:,1:Ny);

      %%% Integrate over FRIS to compute heat flux across ice front
      
      SHImelt_conv(n) = sum(sum(div_htflx(xidx,yidx).*msk(xidx,yidx))) * rhoConst*Cp/Lf;

      %%% Separate estimate of heat content tendency
      heat_tot(n) = rho0*Cp*sum(sum(sum( repmat(msk(xidx,yidx),[1 1 Nr]).*THETA(xidx,yidx,:).*repmat(DRF,[length(xidx) length(yidx) 1]).*hFacC(xidx,yidx,:).*repmat(RAC(xidx,yidx),[1 1 Nr]) )));
     
      

    end

    




  end

  SHImelt_mean = SHImelt_mean / tlen;

  if (compute_eff_melt)    
    SHImelt_tend(2:end-1) = (heat_tot(3:end) - heat_tot(1:end-2)) / (2*t1month) / Lf;
    SHImelt_eff = SHImelt_conv + SHImelt_tend - SHImelt_diff;
  end


    
end
