
ETA = defineMOCgrid(XC,YC,SHELFICEtopo,bathy,false);
idx = (ETA>0 & ETA<4);idx = repmat(idx,[1 1 Nr]);

ds = 0.01;
dt = 0.025;
ss = 34.2:ds:35;
tt = -2.5:dt:0.5;
dens = zeros(length(ss),length(tt));
for i=1:Nx
for j=1:Ny
for k=1:Nr
iidx = ceil((S(i,j,k)-(34.2-ds/2))/ds);
jidx = ceil((T(i,j,k)-(-2.5-dt/2))/dt);
if (idx(i,j,k) && ~isnan(iidx) && ~isnan(jidx) && (iidx>0) && (jidx>0) && (iidx<=length(ss)) && (jidx<=length(tt)))
dens(iidx,jidx) = dens(iidx,jidx)+RAC(i,j)*DRF(k)*hFacC(i,j,k);
end
end
end
end

[TT,SS] = meshgrid(tt,ss);
figure(1);
% pcolor(SS,TT,log10(dens));
pcolor(SS,TT,dens);
shading interp;
colormap(cmocean('amp'));
colorbar;