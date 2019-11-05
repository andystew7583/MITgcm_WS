%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Plot OBCS files %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%% load OBCS files
% load OBCSfiles.mat

%Grid

defineGrid

%%%%% Input file path for OBCS SOSE generated files

input = '/data3/MITgcm_WS/newexp/differentResolutions/a_minboundary';

%%%%%meshgrid for Northern Boundary
run ../analysis/setExpname

run ../analysis/loadexp


[xz,zx] = meshgrid(xx,zz);

%%%%%meshgrid for Eastern/Western Boundary

[yz,zy] = meshgrid(yy,zz);

months = 12;

%%%%%% Read-in bathy/shelf ice files


fid = fopen(fullfile(input,bathyFile),'r','b');
h = fread(fid,[Nx Ny],'real*8')'; 
fclose(fid);

fid = fopen(fullfile(input,SHELFICEtopoFile),'r','b');
icedraft = fread(fid,[Nx Ny],'real*8')'; 
fclose(fid);


    OBNs = zeros(Nx,Nr,12);
    fid = fopen(fullfile(inputpath,OBNsFile),'r','b');
    for k=1:12
        OBNs(:,:,k) = fread(fid,[Nx Nr],'real*8');
    end
    fclose(fid);
    
    OBNs_2 = NaN(Nx,Nz,12);
    nhFAC = squeeze(hFacC(:,end,:));
    for i =1:12
        OBE = OBNs(:,:,i);
        OBE(nhFAC==0) = NaN;
        OBNs_2(:,:,i) = OBE;
    end
        



% topog_msk=NaN(Ny,Nx,Nz);
% for i=1:Nx
%     for j=1:Ny 
%         for k = 1:Nz
%            
%            if (((icedraft(j,i)) - (h(j,i)) == 0  ||  (h(j,i) == icedraft(j,i)) ))
%                topog_msk(j,i,k) = 0;
%            end
%         end
%     end
% end

% topog_msk = transpose3D(topog_msk);


%%%%%%%%%%%% Masking out respective places in OBCS files
months = 12;



%%%%%%%%%%%%%%%%%%%%%%%
%%%%% ECCO BCs %%%%%%%%


% Et_ecco = zeros(Ny,Nz,months);
% obet_ec = zeros(Ny,Nz);
% Es_ecco = zeros(Ny,Nz,months);
% obes_ec = zeros(Ny,Nz);
% Eu_ecco = zeros(Ny,Nz,months);
% obeu_ec = zeros(Ny,Nz);
% Ev_ecco = zeros(Ny,Nz,months);
% obev_ec = zeros(Ny,Nz);
% 
% 
% Nt_ecco = zeros(Nx,Nz,months);
% obnt_ec = zeros(Nx,Nz);
% Ns_ecco = zeros(Nx,Nz,months);
% obns_ec = zeros(Nx,Nz);
% Nu_ecco = zeros(Nx,Nz,months);
% obnu_ec = zeros(Nx,Nz);
% Nv_ecco = zeros(Nx,Nz,months);
% obnv_ec = zeros(Nx,Nz);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% SOSE BCS %%%%%%%%%%%
Et_sose = zeros(Ny,Nz,months);
obet_sose = zeros(Ny,Nz);
Es_sose = zeros(Ny,Nz,months);
obes_sose = zeros(Ny,Nz);
Eu_sose = zeros(Ny,Nz,months);
obeu_sose = zeros(Ny,Nz);
Ev_sose = zeros(Ny,Nz,months);
obev_sose = zeros(Ny,Nz);


Nt_sose = zeros(Nx,Nz,months);
obnt_sose = zeros(Nx,Nz);
Ns_sose = zeros(Nx,Nz,months);
obns_sose = zeros(Nx,Nz);
Nu_sose = zeros(Nx,Nz,months);
obnu_sose = zeros(Nx,Nz);
Nv_sose = zeros(Nx,Nz,months);
obnv_sose = zeros(Nx,Nz);

% topog_msk = transpose3D(topog_msk);

for i = 1:size(Nt_sose,3)
    
%     obet_ec = OBEtEc(:,:,i);
%     obet_ec(squeeze(hFacC(end,:,:))<1) = NaN;
%     obet_ec(obet_ec<-10) = NaN;
%     Et_ecco(:,:,i) = obet_ec;
%     
%     
%     obes_ec = OBEsEc(:,:,i);
%     obes_ec(squeeze(hFacC(end,:,:))<1) = NaN;
%     obes_ec(obes_ec<-10) = NaN;
%     Es_ecco(:,:,i) = obes_ec;
%     
%     obeu_ec = OBEuEc(:,:,i);
%     obeu_ec(hFacC(end,:,:)<1) = NaN;
%     obeu_ec(obeu_ec<-10) = NaN;
%     Eu_ecco(:,:,i) = obeu_ec;
%     
%     obev_ec = OBEvEc(:,:,i);
%     obev_ec(squeeze(hFacC(end,:,:))<1) = NaN;
%     obev_ec(obev_ec<-10) = NaN;
%     Ev_ecco(:,:,i) = obev_ec;
% 
%     obnt_ec = OBNtEc(:,:,i);
%     obnt_ec(squeeze(hFacC(:,end,:))<1) = NaN;
%     obnt_ec(obnt_ec<-10) = NaN;
%     Nt_ecco(:,:,i) = obnt_ec;
%     
%     obns_ec = OBNsEc(:,:,i);
%     obns_ec(squeeze(hFacC(:,end,:))<1) = NaN;
%     obns_ec(obns_ec<-10) = NaN;
%     Ns_ecco(:,:,i) = obns_ec;
%     
%     obnu_ec = OBNuEc(:,:,i);
%     obnu_ec(squeeze(hFacC(:,end,:))<1) = NaN;
%     obnu_ec(obnu_ec<-10) = NaN;
%     Nu_ecco(:,:,i) = obnu_ec;
%     
%     obnv_ec = OBNvEc(:,:,i);
%     obnv_ec(squeeze(hFacC(:,end,:))<1) = NaN;
%     obnv_ec(obnv_ec<-10) = NaN;
%     Nv_ecco(:,:,i) = obnv_ec;
    
    
    
    obet_sose = OBEt(:,:,i);
    obet_sose(squeeze(hFacC(end,:,:))<1) = NaN;
    Et_sose(:,:,i) = obet_sose;
    
    
    obes_sose = OBEs(:,:,i);
    obes_sose(squeeze(hFacC(end,:,:))<1) = NaN;
    obes_sose(obes_sose<-10) = NaN;
    Es_sose(:,:,i) = obes_sose;
    
    obeu_sose = OBEu(:,:,i);
    obeu_sose(squeeze(hFacC(end,:,:))<1) = NaN;
    Eu_sose(:,:,i) = obeu_sose;
    
    obev_sose = OBEv(:,:,i);
    obev_sose(squeeze(hFacC(end,:,:))<1) = NaN;
    Ev_sose(:,:,i) = obev_sose;

    obnt_sose = OBNt(:,:,i);
    obnt_sose(squeeze(hFacC(:,end,:))<1) = NaN;
    Nt_sose(:,:,i) = obnt_sose;
    
    obns_sose = OBNs(:,:,i);
    obns_sose(squeeze(hFacC(:,end,:))<1) = NaN;
    Ns_sose(:,:,i) = obns_sose;
    
    obnu_sose = OBNu(:,:,i);
    obnu_sose(squeeze(hFacC(:,end,:))<1) = NaN;
    Nu_sose(:,:,i) = obnu_sose;
    
    obnv_sose = OBNv(:,:,i);
    obnv_sose(squeeze(hFacC(:,end,:))<1) = NaN;
    Nv_sose(:,:,i) = obnv_sose;
    
    
    

end
%%%%%%% PloTs

%%%%%%%%%% Eastern Boundary

figure(1)
pcolor(yz,zy,(Et_sose(:,:,5)')),shading interp
colorbar
colormap jet(30)
caxis([-3 1])
axis([-72 -63 -500 0]);
title('Theta Eastern Boundary (May)');
print(figure(1),'-dpng',['/data1/MITgcm_WS/newexp/plots/obcsplots/OBEt_tdependence_may.png']);

figure(2)
pcolor(yz,zy,(Es_sose(:,:,5)')),shading interp
colorbar
colormap pmkmp(30)
caxis([34.15 35])
axis([-72 -63 -500 0]);
title('Salinity Eastern Boundary (May)');
print(figure(2),'-dpng',['/data1/MITgcm_WS/newexp/plots/obcsplots/OBEs_tdependence_may.png']);
% 
% figure(3)
% pcolor(yz,zy,nanmean(Eu_sose,3)'),shading interp
% colorbar
% colormap haxby
% title('Uvel Eastern Boundary');
% print(figure(3),'-dpng',['/data1/MITgcm_WS/newexp/plots/obcsplots/sixth/OBEu.png']);
% 
% figure(4)
% pcolor(yz,zy,nanmean(Ev_sose,3)'),shading interp
% colorbar
% colormap haxby
% title('Vvel Eastern Boundary');
% print(figure(4),'-dpng',['/data1/MITgcm_WS/newexp/plots/obcsplots/sixth/OBEv.png']);


%%%%%%%%%% Northern Boundary

figure(5)
pcolor(xx,zx,(Nt_sose(:,:,1)')),shading interp
colormap jet(100)
colorbar
caxis([-3 1])
axis([-80 20 -500 0]);
title('Theta Northern Boundary (Jan)');
print(figure(5),'-dpng',['/data1/MITgcm_WS/newexp/plots/obcsplots/OBNt_jan_tdependent.png']);

figure(6)
pcolor(xx,zx,OBNs_2(:,:,1)'),shading interp
colorbar
colormap(pmkmp(100))
axis([-80 20 -500 0]);
caxis([34.15 35])
title('Salinity Northern Boundary (Nov)');
print(figure(6),'-dpng',['/data1/MITgcm_WS/newexp/plots/obcsplots/OBNs_tdep_jan.png']);

figure(7)
pcolor(xx,zx,nanmean(Nu_sose,3)'),shading interp
colorbar
colormap jet(30)
title('Uvel Northern Boundary');
print(figure(7),'-dpng',['/data1/MITgcm_WS/newexp/plots/obcsplots/sixth/OBNu.png']);

figure(8)
pcolor(xx,zx,nanmean(Nv_sose,3)'),shading interp
colorbar
colormap jet
title('Vvel Northern Boundary');
print(figure(8),'-dpng',['/data1/MITgcm_WS/newexp/plots/obcsplots/sixth/OBNv.png']);

% figure(10)
% pcolor(yy,zy,(nanmean(Et_sose,3)-nanmean(Et_ecco,3))'),shading interp
% colorbar
% colormap jet
% caxis([-1 1])
% title('Theta Eastern Boundary, SOSE minus ECCO');
% print(figure(10),'-dpng',['/data1/MITgcm_WS/newexp/plots/obcsplots/OBEtSOSEminusECCO.png']);
% % 
% figure(11)
% pcolor(yy,zy,(nanmean(Es_sose,3)-nanmean(Es_ecco,3))'),shading interp
% colorbar
% caxis([-.5 .5])
% colormap redblue
% title('Salinity Eastern Boundary, SOSE minus ECCO');
% print(figure(11),'-dpng',['/data1/MITgcm_WS/newexp/plots/obcsplots/OBEsSOSEminusECCO.png']);
% % 
% figure(12)
% pcolor(yy,zy,(nanmean(Eu_sose,3)-nanmean(Eu_ecco,3))'),shading interp
% colorbar
% colormap jet
% title('Zonal velocity Eastern Boundary, SOSE minus ECCO');
% print(figure(12),'-dpng',['/data1/MITgcm_WS/newexp/plots/obcsplots/OBEuSOSEminusECCO.png']);
% % 
% figure(13)
% pcolor(yy,zy,(nanmean(Ev_sose,3)-nanmean(Ev_ecco,3))'),shading interp
% colorbar
% colormap jet
% title('Meridional velocity Eastern Boundary, SOSE minus ECCO');
% print(figure(13),'-dpng',['/data1/MITgcm_WS/newexp/plots/obcsplots/OBEvSOSEminusECCO.png']);
% % 
% % 
% % %%%%%%%%%% Northern Boundary
% % 
% figure(14)
% pcolor(xx,zx,(nanmean(Nt_sose,3)-nanmean(Nt_ecco,3))'),shading interp
% colorbar
% colormap jet
% caxis([-1 1])
% title('Theta Northern Boundary, SOSE minus ECCO');
% print(figure(14),'-dpng',['/data1/MITgcm_WS/newexp/plots/obcsplots/OBNtSOSEminusECCO.png']);
% % 
% figure(15)
% pcolor(xx,zx,(nanmean(Ns_sose,3)-nanmean(Ns_ecco,3))'),shading interp
% colorbar
% colormap redblue
% caxis([-.5 .5])
% title('Salinity Northern Boundary, SOSE minus ECCO');
% print(figure(15),'-dpng',['/data1/MITgcm_WS/newexp/plots/obcsplots/OBNsSOSEminusECCO.png']);
% % 
% figure(16)
% pcolor(xx,zx,(nanmean(Nu_sose,3)-nanmean(Nu_ecco,3))'),shading interp
% colorbar
% colormap jet
% title('Zonal Velocity Northern Boundary,SOSE minus ECCO');
% print(figure(16),'-dpng',['/data1/MITgcm_WS/newexp/plots/obcsplots/OBNuSOSEminusECCO.png']);
%  
% figure(17)
% pcolor(xx,zx,(nanmean(Nv_sose,3)-nanmean(Nv_ecco,3))'),shading interp
% colorbar
% colormap jet
% title('Meridional Velocity Northern Boundary,SOSE minus ECCO');
% print(figure(17),'-dpng',['/data1/MITgcm_WS/newexp/plots/obcsplots/OBNvSOSEminusECCO.png']);
%  
% % 
