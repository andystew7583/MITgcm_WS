%%%
%%% calcCrossStreamlineFluxes.m
%%%
%%% Calculates heat and salt fluxes across barotropic streamlines
%%%

%%% Load experiment data
expdir = '../experiments';
expname = 'hires_seq_onethird_notides_RTOPO2';
% expname = 'hires_seq_onesixth_notides_RTOPO2';
% expname = 'hires_seq_onetwelfth_notides_RTOPO2';
% expname = 'hires_seq_onetwentyfourth_notides_RTOPO2';
deform_cavity = false;
loadexp;
rho0 = 1000;
Cp = 4000;

%%% 3D grid spacing matrices
DXG_3D = repmat(DXG,[1 1 Nr]);
DYG_3D = repmat(DYG,[1 1 Nr]);
DRF_3D = repmat(DRF,[Nx Ny 1]);

%%% T/S flux storage file
infname = [expname,'_TSfluxes'];
if (deform_cavity)
  infname = [infname,'_deform'];
end
infname = [infname,'.mat'];

%%% Load pre-computed horizontal flux data
load(fullfile('products',infname),'uvel_tavg','vvel_tavg','theta_tavg','salt_tavg','uvelth_tavg','uvelslt_tavg','vvelslt_tavg','vvelth_tavg','tflux_tavg','sflux_tavg');

%%% Calculate depth-averaged zonal velocity
UU = sum(uvel_tavg.*DRF_3D.*hFacW,3);

%%% Calculate barotropic streamfunction
Psi = zeros(Nx+1,Ny+1);
Psi(2:Nx+1,2:Ny+1) = -cumsum(UU.*DYG,2);
Psi = Psi(1:Nx,1:Ny);

%%% Interpolate to cell centers
PsiC = 0.25*(Psi(1:Nx,1:Ny)+Psi([2:Nx 1],1:Ny)+Psi(1:Nx,[2:Ny 1])+Psi([2:Nx 1],[2:Ny 1]))/1e6;

%%% Streamunction grid for flux calculation
pp = -2:.1:10;
Np = length(pp);

%%% Compute cross-streamline fluxes
[thflux_tot,thflux_mean,thflux_eddy] = calcMeanEddyFluxes (...
    uvel_tavg,vvel_tavg,theta_tavg,uvelth_tavg,vvelth_tavg, ...
    Nx,Ny,Np, ...  
    DXG_3D,DYG_3D,DRF_3D,hFacW,hFacS,PsiC,pp);  
[sltflux_tot,sltflux_mean,sltflux_eddy] = calcMeanEddyFluxes (...
    uvel_tavg,vvel_tavg,salt_tavg,uvelslt_tavg,vvelslt_tavg, ...
    Nx,Ny,Np, ...  
    DXG_3D,DYG_3D,DRF_3D,hFacW,hFacS,PsiC,pp);
  
%%% Heat flux across barotropic streamlines
UgradPsi = (sum(uvelth_tavg.*hFacW.*DRF_3D,3).*(PsiC(1:Nx,:)-PsiC([Nx 1:Nx-1],:))./DXC + sum(vvelth_tavg.*hFacW.*DRF_3D,3).*(PsiC(:,1:Ny)-PsiC(:,[Ny 1:Ny-1]))./DYC) ./ sqrt(((PsiC(1:Nx,:)-PsiC([Nx 1:Nx-1],:))./DXC).^2+((PsiC(:,1:Ny)-PsiC(:,[Ny 1:Ny-1]))./DYC).^2);
  
%%% Integrate surface fluxes
thflux_surf = zeros(Np,1);
for m = 1:Np
  msk = ETA<eta(m);
  eflux(m) = squeeze(sum(sum(fluxdiv.*msk,1),2));
end  

  
figure(21);
plot(pp,thflux_tot*rho0*Cp);
hold on;
plot(pp,thflux_mean*rho0*Cp);
plot(pp,thflux_eddy*rho0*Cp);
hold off

figure(22);
plot(pp,sltflux_tot*rho0);
hold on;
plot(pp,sltflux_mean*rho0);
plot(pp,sltflux_eddy*rho0);
hold off
         
Psi_plot = PsiC;
Psi_plot(sum(hFacC,3)==0) = NaN;
figure(23);
pcolor(XC,YC,Psi_plot);
shading interp
hold on
[C,h]=contour(XC,YC,Psi_plot,[-2:.2:2],'EdgeColor','k');
clabel(C,h);
hold off
shading interp;
colormap redblue(60);
colorbar;
caxis([-3 3]);
set(gca,'Color',[.8 .8 .8]);

figure(24);
pcolor(XC,YC,UgradPsi);
shading interp;
colorbar;
caxis([-20 20]);
colormap redblue

figure(25);
pcolor(XC,YC,sflux_tavg);
shading interp;
colorbar;
colormap redblue;
caxis([-1 1]*1e-3)

tflux_tavg(sum(hFacC,3)==0) = NaN;
figure(26);
pcolor(XC,YC,tflux_tavg);
shading interp;
colorbar;
colormap redblue;
caxis([-200 200])

figure(27);
contourf(XC,YC,SHELFICEtopo-bathy,[0:100:1000]);
colorbar;
colormap haxby;
caxis([0 1000])

%%%
%%% Convenience function to compute mean and eddy fluxes
%%%
function [trflux_tot,trflux_mean,trflux_eddy] = calcMeanEddyFluxes (...
  uvel,vvel,tracer,uveltr,vveltr, ...
  Nx,Ny,Neta, ...  
  DXG_3D,DYG_3D,DRF_3D,hFacW,hFacS,ETA,eta)

  %%% Calculate midpoint tracer
  tracer_u = 0.5*(tracer([1:Nx],:,:)+tracer([Nx 1:Nx-1],:,:));
  tracer_v = 0.5*(tracer(:,[1:Ny],:)+tracer(:,[Ny 1:Ny-1],:));
  
  %%% Compute eddy heat and salt fluxes on cell faces
  uveltr_mean = uvel.*tracer_u;
  vveltr_mean = vvel.*tracer_v;
  uveltr_eddy = uveltr - uveltr_mean;
  vveltr_eddy = vveltr - vveltr_mean;
  
  %%% Calculate fluxes in quasi-latitude coordinates
  trflux_tot = calcQuasiLatFluxes (...
    uveltr,vveltr, ...
    Nx,Ny,Neta, ...  
    DXG_3D,DYG_3D,DRF_3D,hFacW,hFacS,ETA,eta);
  trflux_mean = calcQuasiLatFluxes (...
    uveltr_mean,vveltr_mean, ...
    Nx,Ny,Neta, ...  
    DXG_3D,DYG_3D,DRF_3D,hFacW,hFacS,ETA,eta);
  trflux_eddy = calcQuasiLatFluxes (...
    uveltr_eddy,vveltr_eddy, ...
    Nx,Ny,Neta, ...  
    DXG_3D,DYG_3D,DRF_3D,hFacW,hFacS,ETA,eta);

end








%%%
%%% Convenience function to comute fluxes in quasi-latitude space
%%%
function eflux = calcQuasiLatFluxes (...
  uflux,vflux, ...
  Nx,Ny,Neta, ...  
  DXG_3D,DYG_3D,DRF_3D,hFacW,hFacS,ETA,eta)

  %%% Integrate fluxes verticall and horizontally over each cell face
  uflux_yzint = sum(uflux .* DYG_3D .* DRF_3D .* hFacW,3);
  vflux_xzint = sum(vflux .* DXG_3D .* DRF_3D .* hFacS,3);

  %%% Compute horizontal divergence of isopycnal fluxes
  fluxdiv = zeros(Nx,Ny);
  fluxdiv(1:Nx-1,1:Ny-1) = uflux_yzint(2:Nx,1:Ny-1) ...
                              - uflux_yzint(1:Nx-1,1:Ny-1) ...
                              + vflux_xzint(1:Nx-1,2:Ny) ...
                              - vflux_xzint(1:Nx-1,1:Ny-1);
                       
  %%% Integrate flux divergence across lines of constant eta 
  eflux = zeros(Neta,1);
  for m = 1:Neta
    msk = ETA<eta(m);
    eflux(m) = squeeze(sum(sum(fluxdiv.*msk,1),2));
  end  

end