%%%%%%% Averaging Zero Days from Sfc Forcing %%%%%%
defineGrid


%%%%%% write dataset path
addpath ../newexp_utils

% days = 26292;
days = 3287;
gendir = '/data3';

inputfr = '/data3/MITgcm_WS/newexp/differentResolutions/ardbeg_tempplus10';

forcingvar1 = zeros(Nx,Ny,days);
fid = fopen(fullfile(inputfr,aTemp),'r','b');
for k=1:days
  forcingvar1(:,:,k) = fread(fid,[Nx Ny],'real*8');
end
fclose(fid);


for i = 1:Nx
    for j = 1:Ny
        for k = 1:days
            if forcingvar1(i,j,k)<-273.15
                
                forcingvar1(i,j,k) = (forcingvar1(i,j,1));
            end
        end
    end
end



writeDataset(forcingvar1,fullfile(inputfolder,aTemp),ieee,prec);

clear forcingvar1

forcingvar1 = zeros(Nx,Ny,days);
fid = fopen(fullfile(inputfr,mwind),'r','b');
for k=1:days
  forcingvar1(:,:,k) = fread(fid,[Nx Ny],'real*8');
end
fclose(fid);


for i = 1:Nx
    for j = 1:Ny
        for k = 1:days
            if (forcingvar1(i,j,k))==0
                
                forcingvar1(i,j,k) = (forcingvar1(i,j,k+8)+forcingvar1(i,j,k-8))/2;
            end
        end
    end
end



writeDataset(forcingvar1,fullfile(inputfr,mwind),ieee,prec);

% 
clear forcingvar1




