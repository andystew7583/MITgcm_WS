
setExpname;
expname = 'hires_seq_onetwentyfourth_RTOPO2';
loadexp;
addpath ../utils/matlab/

iter = 46203;

Eta = rdmds(fullfile(basedir,expname,'results','Eta'),iter);
figure(1);pcolor(XC,YC,Eta);shading interp;colorbar

Heff = rdmds(fullfile(basedir,expname,'results','HEFF'),iter);
Heff(hFacC(:,:,1)==0) = NaN;
figure(2);pcolor(XC,YC,Heff);shading interp;colorbar
 

KPPdiffKzT = rdmds(fullfile(basedir,expname,'results','KPPdiffKzT'),iter);
KPPdiffKzT(hFacC==0) = NaN;
[i,j] = find(KPPdiffKzT==max(KPPdiffKzT(:)));
[j,k] = find(squeeze(KPPdiffKzT(i,:,:))==max(KPPdiffKzT(:)));
figure(3);pcolor(XC,YC,log10(KPPdiffKzT(:,:,k)));shading flat;colorbar;caxis([-2 3]);
[ZZ,YY] = meshgrid(RC,YC(i,:));
figure(4);pcolor(YY,ZZ,log10(squeeze(KPPdiffKzT(i,:,:))));shading flat;colorbar;caxis([-2 3]);


Qnet = rdmds(fullfile(basedir,expname,'results','Qnet'),iter);
Qnet(hFacC(:,:,1)==0) = NaN;
figure(5);pcolor(XC,YC,Qnet);shading flat;colorbar
[i,j]=find(abs(Qnet)==max(abs(Qnet(:))));

Vwind = rdmds(fullfile(basedir,expname,'results','VWIND'),iter);
Vwind(hFacC(:,:,1)==0) = NaN;
figure(6);pcolor(XC,YC,Vwind);shading interp;colorbar

Uwind = rdmds(fullfile(basedir,expname,'results','UWIND'),iter);
Uwind(hFacC(:,:,1)==0) = NaN;
figure(7);pcolor(XC,YC,Uwind);shading interp;colorbar

Area = rdmds(fullfile(basedir,expname,'results','AREA'),iter);
Area(hFacC(:,:,1)==0) = NaN;
figure(8);pcolor(XC,YC,Area);shading interp;colorbar


T = rdmds(fullfile(basedir,expname,'results','T'),iter);
T(hFacC==0) = NaN;
figure(9);pcolor(XC,YC,T(:,:,1));shading interp;colorbar

S = rdmds(fullfile(basedir,expname,'results','S'),iter);
S(hFacC==0) = NaN;
figure(10);pcolor(XC,YC,S(:,:,1));shading interp;colorbar

W = rdmds(fullfile(basedir,expname,'results','W'),iter);
figure(11);pcolor(XC,YC,W(:,:,25));shading interp;colorbar
colormap redblue
[iw,jw] = find(W==max(abs(W(:))));
[jw,kw] = find(squeeze(W(iw,:,:))==max(abs(W(:))));

U = rdmds(fullfile(basedir,expname,'results','U'),iter);
figure(12);pcolor(XC,YC,U(:,:,1));shading interp;colorbar
colormap redblue

V = rdmds(fullfile(basedir,expname,'results','V'),iter);
figure(13);pcolor(XC,YC,V(:,:,1));shading interp;colorbar
colormap redblue

EmPmR = rdmds(fullfile(basedir,expname,'results','EmPmR'),iter);
EmPmR(hFacC(:,:,1)==0) = NaN;
figure(14);pcolor(XC,YC,EmPmR);shading flat;colorbar


Uice = rdmds(fullfile(basedir,expname,'results','UICE'),iter);
figure(15);pcolor(XC,YC,Uice);shading flat;colorbar; colormap redblue; caxis([-4 4]);


Vice = rdmds(fullfile(basedir,expname,'results','VICE'),iter);
figure(16);pcolor(XC,YC,Vice);shading flat;colorbar; colormap redblue; caxis([-4 4]);


figure(17);pcolor(XC,YC,Uice-U(:,:,1));shading flat;colorbar; colormap redblue; caxis([-1 1]);
figure(18);pcolor(XC,YC,Vice-V(:,:,1));shading flat;colorbar; colormap redblue; caxis([-4 4]);

EXFqnet = rdmds(fullfile(basedir,expname,'results','EXFqnet'),iter);
EXFqnet(hFacC(:,:,1)==0) = NaN;
figure(19);pcolor(XC,YC,EXFqnet);shading flat;colorbar

oceaTaux = rdmds(fullfile(basedir,expname,'results','oceTAUX'),iter);
oceaTaux(hFacC(:,:,1)==0) = NaN;
figure(20);pcolor(XC,YC,oceaTaux);shading flat;colorbar

SItflux = rdmds(fullfile(basedir,expname,'results','SItflux'),iter);
SItflux(hFacC(:,:,1)==0) = NaN;
figure(21);pcolor(XC,YC,SItflux);shading flat;colorbar

TFLUX = rdmds(fullfile(basedir,expname,'results','TFLUX'),iter);
TFLUX(hFacC(:,:,1)==0) = NaN;
figure(22);pcolor(XC,YC,TFLUX);shading flat;colorbar

SIhl = rdmds(fullfile(basedir,expname,'results','SIhl'),iter);
SIhl(hFacC(:,:,1)==0) = NaN;
figure(23);pcolor(XC,YC,SIhl);shading flat;colorbar

SItices = rdmds(fullfile(basedir,expname,'results','SItices'),iter);
SItices(hFacC(:,:,1)==0) = NaN;
figure(24);pcolor(XC,YC,SItices);shading flat;colorbar

EXFhl = rdmds(fullfile(basedir,expname,'results','EXFhl'),iter);
EXFhl(hFacC(:,:,1)==0) = NaN;
figure(25);pcolor(XC,YC,EXFhl);shading flat;colorbar

EXFhs = rdmds(fullfile(basedir,expname,'results','EXFhs'),iter);
EXFhs(hFacC(:,:,1)==0) = NaN;
figure(26);pcolor(XC,YC,EXFhs);shading flat;colorbar

EXFtaux = rdmds(fullfile(basedir,expname,'results','EXFtaux'),iter);
EXFtaux(hFacC(:,:,1)==0) = NaN;
figure(27);pcolor(XC,YC,EXFtaux);shading flat;colorbar
