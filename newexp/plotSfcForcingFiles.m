%%%
%%% plotSfcForcingFiles.m
%%%
%%% Makes plots of atmospheric forcing fields for verification purposes.
%%%

defineGrid;
days_end = 3287;
days_start = 1;
fignum = 30;

uwind = zeros(EXF_Nx,EXF_Ny,length(days_start:days_end));
fid = fopen(fullfile(inputfolder,zwind),'r','b');
for k=1:days_start:days_end
    uwind(:,:,k) = fread(fid,[EXF_Nx EXF_Ny],'real*8');
end
fclose(fid);

fignum = fignum+1;
figure(fignum);
plot(squeeze(mean(mean(uwind))));
xlabel('Days');
ylabel('Mean zonal wind');

vwind = zeros(EXF_Nx,EXF_Ny,length(days_start:days_end));
fid = fopen(fullfile(inputfolder,mwind),'r','b');
for k=days_start:days_end
    vwind(:,:,k-days_start+1) = fread(fid,[EXF_Nx EXF_Ny],'real*8');
end
fclose(fid);

fignum = fignum+1;
figure(fignum);
plot(squeeze(mean(mean(vwind))));
xlabel('Days');
ylabel('Mean meridional wind');

atemp = zeros(EXF_Nx,EXF_Ny,length(days_start:days_end));
fid = fopen(fullfile(inputfolder,aTemp),'r','b');
for k=days_start:days_end
    atemp(:,:,k-days_start+1) = fread(fid,[EXF_Nx EXF_Ny],'real*8');
end
fclose(fid);

fignum = fignum+1;
figure(fignum);
plot(squeeze(mean(mean(atemp))));
xlabel('Days');
ylabel('Mean atmos. temp');

lw = zeros(EXF_Nx,EXF_Ny,length(days_start:days_end));
fid = fopen(fullfile(inputfolder,aLW),'r','b');
for k=days_start:days_end      
    lw(:,:,k-days_start+1) = fread(fid,[EXF_Nx EXF_Ny],'real*8');
end
fclose(fid);

fignum = fignum+1;
figure(fignum);
plot(squeeze(mean(mean(lw))));
xlabel('Days');
ylabel('Mean longwave');

sw = zeros(EXF_Nx,EXF_Ny,length(days_start:days_end));
fid = fopen(fullfile(inputfolder,aSW),'r','b');
for k=days_start:days_end
    sw(:,:,k-days_start+1) = fread(fid,[EXF_Nx EXF_Ny],'real*8');
end
fclose(fid);

fignum = fignum+1;
figure(fignum);
plot(squeeze(mean(mean(sw))));
xlabel('Days');
ylabel('Mean shortwave');

precip = zeros(EXF_Nx,EXF_Ny,length(days_start:days_end));
fid = fopen(fullfile(inputfolder,aPrecip),'r','b');
for k=days_start:days_end
    precip(:,:,k-days_start+1) = fread(fid,[EXF_Nx EXF_Ny],'real*8');
end
fclose(fid);

fignum = fignum+1;
figure(fignum);
plot(squeeze(mean(mean(precip))));
xlabel('Days');
ylabel('Mean precip');

apressure = zeros(EXF_Nx,EXF_Ny,length(days_start:days_end));
fid = fopen(fullfile(inputfolder,pressure),'r','b');
for k=days_start:days_end
    apressure(:,:,k-days_start+1) = fread(fid,[EXF_Nx EXF_Ny],'real*8');
end
fclose(fid);

fignum = fignum+1;
figure(fignum);
plot(squeeze(mean(mean(apressure))));
xlabel('Days');
ylabel('Mean atmos. pres.');

aq = zeros(EXF_Nx,EXF_Ny,length(days_start:days_end));
fid = fopen(fullfile(inputfolder,anewAQ),'r','b');
for k=1:days_start:days_end
    aq(:,:,k) = fread(fid,[EXF_Nx EXF_Ny],'real*8');
end
fclose(fid); 

fignum = fignum+1;
figure(fignum);
plot(squeeze(mean(mean(aq))));
xlabel('Days');
ylabel('Mean humidity');
