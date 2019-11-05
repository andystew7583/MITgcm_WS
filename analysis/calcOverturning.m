%%%
%%% calcOverturning.m
%%%
%%% Calculates the overturning circulation, calculated using the MITgcm 
%%% 'layers' package.
%%%

%%% Load experiment
loadexp;

%%% Density bins for MOC calculation  
ptlevs = layers_bounds;
Npt = length(ptlevs)-1;

%%% Frequency of diagnostic output - should match that specified in
%%% data.diagnostics.
dumpFreq = abs(diag_frequency(1));
nDumps = round(nTimeSteps*deltaT/dumpFreq);
dumpIters = round((1:nDumps)*dumpFreq/deltaT);
dumpIters = dumpIters(dumpIters > nIter0);
nDumps = length(dumpIters);

%%% Create a finer vertical grid
ffac = 5;
Nrf = ffac*Nr;
delRf = zeros(1,Nrf); 
for n=1:Nr
  for m=1:ffac
    delRf((n-1)*ffac+m) = delR(n)/ffac;
  end
end
zz = - cumsum((delR + [0 delR(1:Nr-1)])/2);
zz_f = - cumsum((delRf + [0 delRf(1:Nrf-1)])/2);

%%% Partial cell heights on fine grid
hFacS_f = zeros(Nx,Ny,Nrf);
for k=1:Nr
  hFacS_f(:,:,ffac*(k-1)+1:ffac*k) = hFacS(:,:,k*ones(1,ffac));              
end

%%% Grid of actual vertical positions, accounting for partial cells
ZZ = zeros(Nx,Ny,Nr);
ZZ_f = zeros(Nx,Ny,Nrf);
DZ = zeros(Nx,Ny,Nr);
DZ_f = zeros(Nx,Ny,Nrf);
PP = zeros(Nx,Ny,Nr);
ZZ(:,:,1) = - delR(1)*hFacS(:,:,1)/2;
for k=2:Nr
  ZZ(:,:,k) = ZZ(:,:,k-1) - 0.5*delR(k-1)*hFacS(:,:,k-1) - 0.5*delR(k)*hFacS(:,:,k);
end       
ZZ_f(:,:,1) = - delRf(1)*hFacS_f(:,:,1)/2;
for k=2:Nrf 
  ZZ_f(:,:,k) = ZZ_f(:,:,k-1) - 0.5*delRf(k-1)*hFacS_f(:,:,k-1) - 0.5*delRf(k)*hFacS_f(:,:,k);      
end
for k=1:Nr
  DZ(:,:,k) = delR(k);
end   
for k=1:Nrf
  DZ_f(:,:,k) = delRf(k);
end   
for k=1:Nr
  PP(:,:,k) = -delR(k);
end   

%%% Matrices for vertical interpolation  
k_p = zeros(Nx,Ny,Nrf);
k_n = zeros(Nx,Ny,Nrf);
w_n = zeros(Nx,Ny,Nrf);
w_p = zeros(Nx,Ny,Nrf);
for i=1:Nx
  for j=1:Ny
  
    %%% Indices of the lowest cells
    kmax = sum(squeeze(hFacS(i,j,:))~=0);
    kmax_f = ffac*kmax;

    for k=1:Nrf

      %%% Previous and next interpolation indices
      k_p(i,j,k) = ceil(k/ffac-0.5);
      k_n(i,j,k) = k_p(i,j,k) + 1;

      %%% Fine grid cell is above highest coarse grid cell, so fine grid
      %%% gamma will just be set equal to uppermost coarse grid gamma
      if (k_p(i,j,k) <= 0)

        k_p(i,j,k) = 1;
        w_p(i,j,k) = 0;
        w_n(i,j,k) = 1;

      else

        %%% Fine grid cell is below lowest coarse grid cell, so fine grid
        %%% gamma will just be set equal to lowermost coarse grid gamma
        if (k_n(i,j,k) > kmax)

          k_n(i,j,k) = kmax;
          w_n(i,j,k) = 0;
          w_p(i,j,k) = 1;

        %%% Otherwise set weights to interpolate linearly between neighboring
        %%% coarse-grid gammas
        else

          w_p(i,j,k) = (ZZ(i,j,k_n(i,j,k))-ZZ_f(i,j,k))./(ZZ(i,j,k_n(i,j,k))-ZZ(i,j,k_p(i,j,k)));
          w_n(i,j,k) = 1 - w_p(i,j,k);

        end

      end

    end
  
  end
end

%%% Calculate time-averaged isopycnal flux, density and velocity
vflux_tavg = zeros(Nx,Ny,Npt);
h_pt_tavg = zeros(Nx,Ny,Npt);
pt_tavg = zeros(Nx,Ny,Nr);
vvel_tavg = zeros(Nx,Ny,Nr);
navg = 0;
for n=1:length(dumpIters)
 
  tyears = dumpIters(n)*deltaT/86400/365;
 
  if ((tyears >= tmin) && (tyears <= tmax))    

    [tyears dumpIters(n)]
    vflux = rdmdsWrapper(fullfile(exppath,'results/LaVH1TH'),dumpIters(n));      
    h_pt = rdmdsWrapper(fullfile(exppath,'results/LaHs1TH'),dumpIters(n));      
    pt  = rdmdsWrapper(fullfile(exppath,'/results/THETA'),dumpIters(n));                      
    vvel  = rdmdsWrapper(fullfile(exppath,'/results/VVEL'),dumpIters(n));              
    
    if (isempty(vflux) || isempty(pt) ...
        || isempty(h_pt) || isempty(vvel))
      ['Ran out of data at n=',num2str(n),'/',num2str(nDumps),' t=',num2str(tyears),' days.']
      break;
    else
      vflux_tavg = vflux_tavg + vflux;
      h_pt_tavg = h_pt_tavg + h_pt;
      pt_tavg = pt_tavg + squeeze(pt(:,:,:,1));      
      vvel_tavg = vvel_tavg + squeeze(vvel(:,:,:,1));      
      navg = navg + 1;
    end
  end
   
end

%%% Calculate the time average
if (navg == 0)
  error('No data files found');
end
vflux_tavg = vflux_tavg/navg;
h_pt_tavg = h_pt_tavg/navg;
pt_tavg = pt_tavg/navg;
vvel_tavg = vvel_tavg/navg;
pt_tavg(hFacC==0) = NaN;

%%% Interpolate potential temperature to v-gridpoints  
pt_v = NaN*pt_tavg;
pt_v(:,2:Ny,:) = 0.5* (pt_tavg(:,1:Ny-1,:) + pt_tavg(:,2:Ny,:));    

%%% Interpolate onto a finer grid         
vvel_f = zeros(Nx,Ny,Nrf);
pt_f = NaN*zeros(Nx,Ny,Nrf);
if (ffac == 1)

  %%% Shortcut if fine grid resolution = coarse grid resolution
  vvel_f = vvel_tavg;        
  pt_f = pt_v;

else   

  %%% Velocity uniform throughout each coarse grid cell to preserve
  %%% mass conservation
  for k=1:Nr
    vvel_f(:,:,ffac*(k-1)+1:ffac*k) = vvel_tavg(:,:,k*ones(1,ffac));          
  end

  %%% Linearly interpolate density
  for i=1:Nx
    for j=3:Ny-1 %%% Restrict to wet grid cells  
      pt_f(i,j,:) = w_p(i,j,:).*pt_v(i,j,squeeze(k_p(i,j,:))) + w_n(i,j,:).*pt_v(i,j,squeeze(k_n(i,j,:)));
    end
  end

end            

%%% Calculate mean fluxes within mean density surfaces
vflux_m = 0*vflux_tavg;
vdz = vvel_f.*hFacS_f.*DZ_f;
vflux_m(:,:,Npt) = vflux_m(:,:,Npt) + sum(vdz.*(pt_f>ptlevs(Npt)),3);
vflux_m(:,:,1) = vflux_m(:,:,1) + sum(vdz.*(pt_f<=ptlevs(2)),3);
for m=2:Npt-1
  vflux_m(:,:,m) = vflux_m(:,:,m) + sum(vdz.*((pt_f>ptlevs(m)) & (pt_f<=ptlevs(m+1))),3);
end   

%%% Zonally integrate meridional fluxes
vflux_xint = zeros(Ny,Npt);
vflux_m_xint = zeros(Ny,Npt);
for i=1:Nx
  vflux_xint = vflux_xint + delX(i)*squeeze(vflux_tavg(i,:,:));
  vflux_m_xint = vflux_m_xint + delX(i)*squeeze(vflux_m(i,:,:));
end

%%% Sum fluxes to obtain streamfunction
psi_pt = zeros(Ny,Npt+1);
psim_pt = zeros(Ny,Npt+1);
for m=1:Npt  
  psi_pt(:,m) = sum(vflux_xint(:,m:Npt),2);     
  psim_pt(:,m) = sum(vflux_m_xint(:,m:Npt),2);     
end
psi_pt = psi_pt/1e6;
psim_pt = psim_pt/1e6;
psie_pt = psi_pt - psim_pt;

%%% Calculate mean density surface heights
h_pt_xtavg = squeeze(nanmean(h_pt_tavg));
z_pt = 0*h_pt_xtavg;
for m=1:Npt
  z_pt(:,m) = - sum(h_pt_xtavg(:,1:m-1),2);
end

%%% Calculate zonal-mean potential temperature
pt_xtavg = squeeze(nanmean(pt_tavg(:,:,:)));
pt_f_xtavg = squeeze(nanmean(pt_f(:,:,:)));

%%% Convert to z-coordinates by mapping the streamfunction at each temp 
%%% level to the mean height of that density surface
psi_z = NaN*ones(Ny,Nrf);
psim_z = NaN*ones(Ny,Nrf);
psie_z = NaN*ones(Ny,Nrf);
for j=1:Ny  

  for k=1:Nrf

    %%% Density lies in the lowest bin
    if (pt_f_xtavg(j,k) <= ptlevs(1))
      psi_z(j,k) = psi_pt(j,1);      
      psim_z(j,k) = psim_pt(j,1);    
      psie_z(j,k) = psie_pt(j,1);    
      continue;
    end

    %%% Density lies in the highest bin
    if (pt_f_xtavg(j,k) > ptlevs(Npt))
      psi_z(j,k) = psi_pt(j,Npt);      
      psim_z(j,k) = psim_pt(j,Npt);      
      psie_z(j,k) = psie_pt(j,Npt);      
      continue;
    end    

    %%% Density lies in an intermediate bin, so find the bin and assign
    %%% the overturning streamfunction via linear interpolation
    for m=1:Npt-1
      if (pt_f_xtavg(j,k) < ptlevs(m+1))
        pt_n = ptlevs(m+1);
        pt_p = ptlevs(m);
        wp = (pt_n-pt_f_xtavg(j,k))/(pt_n-pt_p);
        wn = 1 - wp;
        psi_z(j,k) = wp*psi_pt(j,m) + wn*psi_pt(j,m+1);
        psim_z(j,k) = wp*psim_pt(j,m) + wn*psim_pt(j,m+1);
        psie_z(j,k) = wp*psie_pt(j,m) + wn*psie_pt(j,m+1);
        break;
      end
    end
    
  end
  
end

% %%% Map streamfunction back to z-coordinates
% psi_z = zeros(Ny+1,Nr+1);
% psim_z = zeros(Ny+1,Nr+1);
% psie_z = zeros(Ny+1,Nr+1);
% ZZ_psi =  zeros(Ny+1,Nr+1);
% for j=2:Ny
%   
%   %%% Index of lowest wet cell and numerical topographic depth
%   kbot = sum(hFacS(1,j,:)~=0);
%   ZZ_psi(j,:) = - cumsum([0 squeeze(hFacS(1,j,:))'.*delR]);
%   
%   %%% Skip walls
%   if (kbot == 0)
%     continue;
%   end
%   
%   %%% Loop over wet cell corners
%   for k=2:kbot
%     
%     %%% Find layer depth nearest to this grid cell corner
%     mnext = -1;
%     for m=1:Npt
%       if (ZZ_psi(j,k) > z_pt(j,m))
%         mnext = m;
%         break;
%       end
%     end
%     if (mnext == -1)
%       mnext = Npt+1;
%     end
%     
%     %%% zz_psi(k) lies above the shallowest zrho
%     if (mnext == 1)
%       zprev = 0;  
%       psi_prev = 0;
%       psim_prev = 0;
%       psie_prev = 0;
%     else
%       zprev = z_pt(j,mnext-1);
%       psi_prev = psi_pt(j,mnext-1);
%       psim_prev = psim_pt(j,mnext-1);
%       psie_prev = psie_pt(j,mnext-1);
%     end
% 
%     %%% zz_psi(k) lies deeper than the deepest zrho
%     if (mnext == Npt+1)
%       znext = ZZ_psi(j,kbot);
%       psi_next = 0;
%       psim_next = 0;
%       psie_next = 0;
%     else
%       znext = z_pt(j,mnext);
%       psi_next = psi_pt(j,mnext);
%       psim_next = psim_pt(j,mnext);
%       psie_next = psie_pt(j,mnext);
%     end
%     
%     %%% Interpolation weights
%     wprev = (znext-ZZ_psi(j,k)) / (znext-zprev);
%     wnext = 1 - wprev;
%     
%     %%% Interpolate to grid cell corner
%     psi_z(j,k) = wprev*psi_prev + wnext*psi_next;
%     psim_z(j,k) = wprev*psim_prev + wnext*psim_next;
%     psie_z(j,k) = wprev*psie_prev + wnext*psie_next;
%     
%   end
% end

%%% Store computed data for later
save([expname,'_MOC_pt.mat'],'xx','yy','zz','zz_f','hFacS_f','delRf','ptlevs', ... 
  'vvel','vvel_f','pt_tavg','pt_f', ...  
  'vflux','vflux_m', ...
  'psi_pt','psim_pt','psie_pt', ...
  'psi_z','psim_z','psie_z');