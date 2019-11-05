%%%%%%% extracting the last year of the surface forcing
%%%%%%%
addpath ../newexp_utils

run defineGrid.m
days_year = 365;
days_start = 1;
days_end = 3287;


%     uwind = zeros(Nx,Ny,length(days_start:days_end));
%     fid = fopen(fullfile(inputfolder,zwind),'r','b');
%     for k=days_start:days_end
%         uwind(:,:,k-days_start+1) = fread(fid,[Nx Ny],'real*8');
%     end
%     fclose(fid);
% 
% 
%     vwind = zeros(Nx,Ny,length(days_start:days_end));
%     fid = fopen(fullfile(inputfolder,mwind),'r','b');
%     for k=days_start:days_end
%         vwind(:,:,k-days_start+1) = fread(fid,[Nx Ny],'real*8');
%     end
%     fclose(fid);

    temp = zeros(Nx,Ny,length(days_start:days_end));
    fid = fopen(fullfile(inputfolder,aLW),'r','b');
    for k=days_start:days_end
        temp(:,:,k-days_start+1) = fread(fid,[Nx Ny],'real*8');
    end
    fclose(fid);
    
    temp =zeros(Nx,Ny,3287);
    for i = 1:Nx
        for j = 1:Ny
            for k = 1:days_end
                temp(i,j,k) = 0;
            end
        end
    end
    

%     lw = zeros(Nx,Ny,length(days_start:days_end));
%     fid = fopen(fullfile(inputfolder,aLW),'r','b');
%     for k=days_start:days_end
%         lw(:,:,k-days_start+1) = fread(fid,[Nx Ny],'real*8');
%     end
%     fclose(fid);
% 
%     sw = zeros(Nx,Ny,length(days_start:days_end));
%     fid = fopen(fullfile(inputfolder,aSW),'r','b');
%     for k=days_start:days_end
%         sw(:,:,k-days_start+1) = fread(fid,[Nx Ny],'real*8');
%     end
%     fclose(fid);
% 
%     precip = zeros(Nx,Ny,length(days_start:days_end));
%     fid = fopen(fullfile(inputfolder,aPrecip),'r','b');
%     for k=days_start:days_end
%         precip(:,:,k-days_start+1) = fread(fid,[Nx Ny],'real*8');
%     end
%     fclose(fid);
% 
%     apressure = zeros(Nx,Ny,length(days_start:days_end));
%     fid = fopen(fullfile(inputfolder,pressure),'r','b');
%     for k=days_start:days_end
%         apressure(:,:,k-days_start+1) = fread(fid,[Nx Ny],'real*8');
%     end
%     fclose(fid);
% 
%     aq = zeros(Nx,Ny,length(days_start:days_end));
%     fid = fopen(fullfile(inputfolder,anewAQ),'r','b');
%     for k=days_start:days_end
%         aq(:,:,k-days_start+1) = fread(fid,[Nx Ny],'real*8');
%     end
%     fclose(fid);
    

% data = uwind;
% writeDataset(data,fullfile(inputfolder,zwind),ieee,prec);
% clear data
% 
% data = vwind;
% writeDataset(data,fullfile(inputfolder,mwind),ieee,prec);
% clear data
% 
% data =aq;
% writeDataset(data,fullfile(inputfolder,anewAQ),ieee,prec);
% clear data
   
data =temp;
writeDataset(data,fullfile(inputfolder,aLW),ieee,prec);
clear data

% data =precip;
% writeDataset(data,fullfile(inputfolder,aPrecip),ieee,prec);
% clear data
% 
% data =sw;
% writeDataset(data,fullfile(inputfolder,aSW),ieee,prec);
% clear data
% 
% data =lw;
% writeDataset(data,fullfile(inputfolder,aLW),ieee,prec);
% clear data