%%%
%%% plotAABWcirculation.m
%%%
%%% Plots horizontal circulation in isopycnal layers.
%%%

%%% Options
expdir = '../experiments';
expname = 'hires_seq_onethird_RTOPO2';
loadexp;

%%% Load pre-computed data
outfname = [expname,'_AABWcirc_.mat'];
load(fullfile('products',outfname));

%%% Density bounds for water masses
dens_CDW_min = 27.7;
dens_CDW_max = 27.85;
dens_AABW = 27.85;

%%% Density grid indices for water masses
k_CDW_min = find(dens_levs==dens_CDW_min);
k_CDW_max = find(dens_levs==dens_CDW_max) - 1;
k_AASW_max = k_CDW_min - 1;
k_AABW = find(dens_levs==dens_AABW);

%%% AABW layer thickness
H_AABW_w = sum(uthic_tavg(:,:,k_AABW:end),3);
H_AABW_s = sum(vthic_tavg(:,:,k_AABW:end),3);

%%% Transports and TWA velocities in different water mass layers
hu_AASW = sum(uflux_tavg(:,:,1:k_AASW_max),3);
u_AASW = hu_AASW ./ sum(uthic_tavg(:,:,1:k_AASW_max),3);
hv_AASW = sum(vflux_tavg(:,:,1:k_AASW_max),3);
v_AASW = hv_AASW ./ sum(vthic_tavg(:,:,1:k_AASW_max),3);
hu_CDW = sum(uflux_tavg(:,:,k_CDW_min:k_CDW_max),3);
u_CDW = hu_CDW ./ sum(uthic_tavg(:,:,k_CDW_min:k_CDW_max),3);
hv_CDW = sum(vflux_tavg(:,:,k_CDW_min:k_CDW_max),3);
v_CDW = hv_CDW ./ sum(vthic_tavg(:,:,k_CDW_min:k_CDW_max),3);
hu_AABW = sum(uflux_tavg(:,:,k_AABW:end),3);
u_AABW = hu_AABW ./ H_AABW_w;
hv_AABW = sum(vflux_tavg(:,:,k_AABW:end),3);
v_AABW = hv_AABW ./ H_AABW_s;

%%% Mean/eddy decomposition
hu_AABW_mean = sum(uflux_mean(:,:,k_AABW:end),3);
u_AABW_mean = hu_AABW_mean ./ H_AABW_w;
hv_AABW_mean = sum(vflux_mean(:,:,k_AABW:end),3);
v_AABW_mean = hv_AABW_mean ./ H_AABW_s;
hu_AABW_eddy = hu_AABW - hu_AABW_mean;
hv_AABW_eddy = hv_AABW - hv_AABW_mean;
u_AABW_eddy = hu_AABW_eddy ./ H_AABW_w;
v_AABW_eddy = hv_AABW_eddy ./ H_AABW_s;
u_AABW_eddy(H_AABW_w==0) = NaN;
v_AABW_eddy(H_AABW_s==0) = NaN;

%%% Cap flow speeds
uabs_AABW_eddy = sqrt(u_AABW_eddy.^2+v_AABW_eddy.^2);
uabs_AABW_eddy_max = 0.3;
u_AABW_eddy(uabs_AABW_eddy>uabs_AABW_eddy_max) = u_AABW_eddy(uabs_AABW_eddy>uabs_AABW_eddy_max) .* uabs_AABW_eddy_max./uabs_AABW_eddy(uabs_AABW_eddy>uabs_AABW_eddy_max);
v_AABW_eddy(uabs_AABW_eddy>uabs_AABW_eddy_max) = v_AABW_eddy(uabs_AABW_eddy>uabs_AABW_eddy_max) .* uabs_AABW_eddy_max./uabs_AABW_eddy(uabs_AABW_eddy>uabs_AABW_eddy_max);


%%% Attempt at an AABW streamfunction
psi_AABW = zeros(Nx,Ny);
psi_AABW(:,2:Ny) = -cumsum(hu_AABW(:,1:Ny-1).*DYG(:,1:Ny-1),2);
psi_AABW(2:Nx,:) = cumsum(hv_AABW(1:Nx-1,:).*DYG(1:Nx-1,:),1);

%%% Diapycnal velocity
w_AABW = ( hu_AABW(1:Nx,1:Ny).*DYG(1:Nx,1:Ny) ...
         - hu_AABW([2:Nx 1],1:Ny).*DYG([2:Nx 1],1:Ny) ...
         + hv_AABW(1:Nx,1:Ny).*DXG(1:Nx,1:Ny) ...
         - hv_AABW(1:Nx,[2:Ny 1]).*DYG(1:Nx,[2:Ny 1]) ) ./ RAC;

%%% Subsampling options for quiver plots
arrowspacing = 3*round(Nx/304); %%% Scale with grid sizes
xidx = 1:arrowspacing:Nx;
yidx = 1:arrowspacing:Ny;

figure(71);
quiver(XC(xidx,yidx),YC(xidx,yidx),u_AABW(xidx,yidx),v_AABW(xidx,yidx));

figure(72);
quiver(XC(xidx,yidx),YC(xidx,yidx),hu_AABW(xidx,yidx),hv_AABW(xidx,yidx));

figure(73);
pcolor(XC,YC,w_AABW);
shading interp;
colorbar;
colormap redblue;
caxis([-2 2]*1e-4);

figure(74);
pcolor(XG,YG,psi_AABW/1e6);
shading interp;
colorbar;
colormap redblue;
caxis([-10 10]);


figure(75);
quiver(XC(xidx,yidx),YC(xidx,yidx),u_CDW(xidx,yidx),v_CDW(xidx,yidx));

figure(76);
quiver(XC(xidx,yidx),YC(xidx,yidx),hu_CDW(xidx,yidx),hv_CDW(xidx,yidx));


figure(77);
quiver(XC(xidx,yidx),YC(xidx,yidx),u_AASW(xidx,yidx),v_AASW(xidx,yidx));

figure(78);
quiver(XC(xidx,yidx),YC(xidx,yidx),hu_AASW(xidx,yidx),hv_AASW(xidx,yidx));


figure(79);
quiver(XC(xidx,yidx),YC(xidx,yidx),u_AABW_eddy(xidx,yidx),v_AABW_eddy(xidx,yidx));
hold on;
[C,h] = contour(XC,YC,bathy,[-4000 -3000 -2000 -1000 -500],'EdgeColor','k');
clabel(C,h);
hold off;

figure(80);
quiver(XC(xidx,yidx),YC(xidx,yidx),hu_AABW_eddy(xidx,yidx),hv_AABW_eddy(xidx,yidx));
hold on;
[C,h] = contour(XC,YC,bathy,[-4000 -3000 -2000 -1000 -500],'EdgeColor','k');
clabel(C,h);
hold off;

figure(81);
quiver(XC(xidx,yidx),YC(xidx,yidx),u_AABW_mean(xidx,yidx),v_AABW_mean(xidx,yidx));

figure(82);
quiver(XC(xidx,yidx),YC(xidx,yidx),hu_AABW_mean(xidx,yidx),hv_AABW_mean(xidx,yidx));

