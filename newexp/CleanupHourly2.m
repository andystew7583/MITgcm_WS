%%%%%%% Averaging Zero Days from Sfc Forcing %%%%%%
defineGrid

days = 26292;

addpath ../newexp_utils/

forcingvar1 = zeros(Nx,Ny,days);
fid = fopen(fullfile(inputfolder,'apressurefile.bin'),'r','b');
for k=1:days
  forcingvar1(:,:,k) = fread(fid,[Nx Ny],'real*8');
end
fclose(fid);
for k = 1:days
%     if k >1750 && k <= 1824 %%% 223
% 
%         forcingvar1(:,:,k) = (forcingvar1(:,:,k-8) +forcingvar1(:,:,k)+ forcingvar1(:,:,k+8))/3;
%     
%     elseif k >3872 && k <= 3824 %%% 476
% 
%         forcingvar1(:,:,k) = (forcingvar1(:,:,k-8) +forcingvar1(:,:,k)+ forcingvar1(:,:,k+8))/3;
% 
%     elseif k >10392 && k <= 10408 %%% 1300
% 
%         forcingvar1(:,:,k) = (forcingvar1(:,:,k-8) +forcingvar1(:,:,k)+ forcingvar1(:,:,k+8))/3;
% 
%     elseif k > 10496 && k<= 10510  %%%% day 1313
%         forcingvar1(:,:,k) = (forcingvar1(:,:,k+8)+ forcingvar1(:,:,k)+forcingvar1(:,:,k-8))/3;
%             
%         
%     elseif k > 11216 && k <= 11232  %%%% day 1403
% 
%         forcingvar1(:,:,k) = (forcingvar1(:,:,k+8)+ forcingvar1(:,:,k)+forcingvar1(:,:,k-8))/3;
% 
%             
%     elseif k > 12152 && k <= 12168  %%%% 1520
% 
%         forcingvar1(:,:,k) = (forcingvar1(:,:,k+8)+ forcingvar1(:,:,k)+forcingvar1(:,:,k-8))/3;
%  
%                  
% 
%     elseif k > 13848 && k<=13864 %%%% 1732
% 
%         forcingvar1(:,:,k) = (forcingvar1(:,:,k+8)+ forcingvar1(:,:,k)+forcingvar1(:,:,k-8))/3; 
        
    if k > 16776 && k<=17000 %%%% 2098

        forcingvar1(:,:,k) = (forcingvar1(:,:,k+8)+ forcingvar1(:,:,k)+forcingvar1(:,:,k-8))/3;    
            
%     elseif k > 18216 && k<= 18232 %%% 2278
% 
%         forcingvar1(:,:,k) = (forcingvar1(:,:,k+8)+ forcingvar1(:,:,k)+forcingvar1(:,:,k-8))/3;
% 
%         
%     elseif k > 18384 && k<= 18400     
%         forcingvar1(:,:,k) = (forcingvar1(:,:,k+8)+ forcingvar1(:,:,k)+forcingvar1(:,:,k-8))/3;
%   
%         
%     elseif k > 24944 && k <= 24960
%         forcingvar1(:,:,k) = (forcingvar1(:,:,k+8)+ forcingvar1(:,:,k)+forcingvar1(:,:,k-8))/3;
% 
%         
%     elseif k > 25696 && k <= 25720   
%         forcingvar1(:,:,k) = (forcingvar1(:,:,k+8)+ forcingvar1(:,:,k)+forcingvar1(:,:,k-8))/3;

    end
end
    
    
    writeDataset(forcingvar1,fullfile(inputfolder,'apressurefile.bin'),ieee,prec);


clear forcingvar1






forcingvar1 = zeros(Nx,Ny,days);
fid = fopen(fullfile(inputfolder,aTemp),'r','b');
for k=1:days
  forcingvar1(:,:,k) = fread(fid,[Nx Ny],'real*8');
end
fclose(fid);

for k = 1:days
%     if k >3792 && k <= 3816 %%% 476
% 
%         forcingvar1(:,:,k) = forcingvar1(:,:,k-8) +forcingvar1(:,:,k)+ forcingvar1(:,:,k+8)/3;
% 
% 
%     elseif k >10392 && k <= 10408 %%% 1300
% 
%         forcingvar1(:,:,k) = forcingvar1(:,:,k-8) +forcingvar1(:,:,k)+ forcingvar1(:,:,k+8)/3;
% 
%     elseif k > 10496 && k<= 10510  %%%% day 1313
%         forcingvar1(:,:,k) = (forcingvar1(:,:,k+8)+ forcingvar1(:,:,k)+forcingvar1(:,:,k-8))/3;
%             
%         
%     elseif k > 11216 && k <= 11232  %%%% day 1403
% 
%         forcingvar1(:,:,k) = (forcingvar1(:,:,k+8)+ forcingvar1(:,:,k)+forcingvar1(:,:,k-8))/3;
% 
%             
%     elseif k > 12152 && k <= 12168  %%%% 1520
% 
%         forcingvar1(:,:,k) = (forcingvar1(:,:,k+8)+ forcingvar1(:,:,k)+forcingvar1(:,:,k-8))/3;
%  
%                   
%     elseif k > 13848 && k<=13864 %%%% 1732
% 
%         forcingvar1(:,:,k) = (forcingvar1(:,:,k+8)+ forcingvar1(:,:,k)+forcingvar1(:,:,k-8))/3;
% 
%             
%             
%     elseif k > 18216 && k<= 18232 %%% 2278
% 
%         forcingvar1(:,:,k) = (forcingvar1(:,:,k+8)+ forcingvar1(:,:,k)+forcingvar1(:,:,k-8))/3;
% 
%         
%     elseif k > 18384 && k<= 18400     
%         forcingvar1(:,:,k) = (forcingvar1(:,:,k+8)+ forcingvar1(:,:,k)+forcingvar1(:,:,k-8))/3;
%   
%         
%     elseif k > 24944 && k <= 24960
%         forcingvar1(:,:,k) = (forcingvar1(:,:,k+8)+ forcingvar1(:,:,k)+forcingvar1(:,:,k-8))/3;
% 
%         
%     elseif k > 25696 && k <= 25720   
%         forcingvar1(:,:,k) = (forcingvar1(:,:,k+8)+ forcingvar1(:,:,k)+forcingvar1(:,:,k-8))/3;
    if k > 16776 && k<=17000 %%%% 2098

        forcingvar1(:,:,k) = (forcingvar1(:,:,k+8)+ forcingvar1(:,:,k)+forcingvar1(:,:,k-8))/3;  
    end
end



writeDataset(forcingvar1,fullfile(inputfolder,aTemp),ieee,prec);


clear forcingvar1





forcingvar1 = zeros(Nx,Ny,days);
fid = fopen(fullfile(inputfolder,mwind),'r','b');
for k=1:days
  forcingvar1(:,:,k) = fread(fid,[Nx Ny],'real*8');
end
fclose(fid);


for k = 1:days
%     if k >3792 && k <= 3816 %%% 476
% 
%         forcingvar1(:,:,k) = (forcingvar1(:,:,k-8) +forcingvar1(:,:,k)+ forcingvar1(:,:,k+8))/3;
% 
% 
%     elseif k >10392 && k <= 10408 %%% 1300
% 
%         forcingvar1(:,:,k) = (forcingvar1(:,:,k-8) +forcingvar1(:,:,k)+ forcingvar1(:,:,k+8))/3;
% 
%     elseif k > 10496 && k<= 10510  %%%% day 1313
%         forcingvar1(:,:,k) = (forcingvar1(:,:,k+8)+ forcingvar1(:,:,k)+forcingvar1(:,:,k-8))/3;
%             
%         
%     elseif k > 11216 && k <= 11232  %%%% day 1403
% 
%         forcingvar1(:,:,k) = (forcingvar1(:,:,k+8)+ forcingvar1(:,:,k)+forcingvar1(:,:,k-8))/3;
% 
%             
%     elseif k > 12152 && k <= 12168  %%%% 1520
% 
%         forcingvar1(:,:,k) = (forcingvar1(:,:,k+8)+ forcingvar1(:,:,k)+forcingvar1(:,:,k-8))/3;
%  
%                   
%     elseif k > 13848 && k<=13864 %%%% 1732
% 
%         forcingvar1(:,:,k) = (forcingvar1(:,:,k+8)+ forcingvar1(:,:,k)+forcingvar1(:,:,k-8))/3;
% 
%             
%             
%     elseif k > 18216 && k<= 18232 %%% 2278
% 
%         forcingvar1(:,:,k) = (forcingvar1(:,:,k+8)+ forcingvar1(:,:,k)+forcingvar1(:,:,k-8))/3;
% 
%         
%     elseif k > 18384 && k<= 18400     
%         forcingvar1(:,:,k) = (forcingvar1(:,:,k+8)+ forcingvar1(:,:,k)+forcingvar1(:,:,k-8))/3;
%   
%         
%     elseif k > 24944 && k <= 24960
%         forcingvar1(:,:,k) = (forcingvar1(:,:,k+8)+ forcingvar1(:,:,k)+forcingvar1(:,:,k-8))/3;
% 
%         
%     elseif k > 25696 && k <= 25720   
%         forcingvar1(:,:,k) = (forcingvar1(:,:,k+8)+ forcingvar1(:,:,k)+forcingvar1(:,:,k-8))/3;
    if k > 16776 && k<=17000 %%%% 2098

        forcingvar1(:,:,k) = (forcingvar1(:,:,k+8)+ forcingvar1(:,:,k)+forcingvar1(:,:,k-8))/3;  
    end
end



writeDataset(forcingvar1,fullfile(inputfolder,mwind),ieee,prec);


clear forcingvar1

%%%%%%%%%%%


forcingvar1 = zeros(Nx,Ny,days);
fid = fopen(fullfile(inputfolder,zwind),'r','b');
for k=1:days
  forcingvar1(:,:,k) = fread(fid,[Nx Ny],'real*8');
end
fclose(fid);


for k = 1:days
%     if k >3792 && k <= 3816 %%% 476
% 
%         forcingvar1(:,:,k) = (forcingvar1(:,:,k-8) +forcingvar1(:,:,k)+ forcingvar1(:,:,k+8))/3;
% 
% 
%     elseif k >10392 && k <= 10408 %%% 1300
% 
%         forcingvar1(:,:,k) = (forcingvar1(:,:,k-8) +forcingvar1(:,:,k)+ forcingvar1(:,:,k+8))/3;
% 
%     elseif k > 10496 && k<= 10510  %%%% day 1313
%         forcingvar1(:,:,k) = (forcingvar1(:,:,k+8)+ forcingvar1(:,:,k)+forcingvar1(:,:,k-8))/3;
%             
%         
%     elseif k > 11216 && k <= 11232  %%%% day 1403
% 
%         forcingvar1(:,:,k) = (forcingvar1(:,:,k+8)+ forcingvar1(:,:,k)+forcingvar1(:,:,k-8))/3;
% 
%             
%     elseif k > 12152 && k <= 12168  %%%% 1520
% 
%         forcingvar1(:,:,k) = (forcingvar1(:,:,k+8)+ forcingvar1(:,:,k)+forcingvar1(:,:,k-8))/3;
%  
%                   
%     elseif k > 13848 && k<=13864 %%%% 1732
% 
%         forcingvar1(:,:,k) = (forcingvar1(:,:,k+8)+ forcingvar1(:,:,k)+forcingvar1(:,:,k-8))/3;
% 
%             
%             
%     elseif k > 18216 && k<= 18232 %%% 2278
% 
%         forcingvar1(:,:,k) = (forcingvar1(:,:,k+8)+ forcingvar1(:,:,k)+forcingvar1(:,:,k-8))/3;
% 
%         
%     elseif k > 18384 && k<= 18400     
%         forcingvar1(:,:,k) = (forcingvar1(:,:,k+8)+ forcingvar1(:,:,k)+forcingvar1(:,:,k-8))/3;
%   
%         
%     elseif k > 24944 && k <= 24960
%         forcingvar1(:,:,k) = (forcingvar1(:,:,k+8)+ forcingvar1(:,:,k)+forcingvar1(:,:,k-8))/3;
% 
%         
%     elseif k > 25696 && k <= 25720   
%         forcingvar1(:,:,k) = (forcingvar1(:,:,k+8)+ forcingvar1(:,:,k)+forcingvar1(:,:,k-8))/3;
    if k > 16776 && k<=17000 %%%% 2098

        forcingvar1(:,:,k) = (forcingvar1(:,:,k+8)+ forcingvar1(:,:,k)+forcingvar1(:,:,k-8))/3;  
    end
end



writeDataset(forcingvar1,fullfile(inputfolder,zwind),ieee,prec);


clear forcingvar1



forcingvar1 = zeros(Nx,Ny,days);
fid = fopen(fullfile(inputfolder,aLW),'r','b');
for k=1:days
  forcingvar1(:,:,k) = fread(fid,[Nx Ny],'real*8');
end
fclose(fid);


for k = 1:days
%     if k >3792 && k <= 3816 %%% 476
% 
%         forcingvar1(:,:,k) = (forcingvar1(:,:,k-8) +forcingvar1(:,:,k)+ forcingvar1(:,:,k+8))/3;
% 
% 
%     elseif k >10392 && k <= 10408 %%% 1300
% 
%         forcingvar1(:,:,k) = (forcingvar1(:,:,k-8) +forcingvar1(:,:,k)+ forcingvar1(:,:,k+8))/3;
% 
%     elseif k > 10496 && k<= 10510  %%%% day 1313
%         forcingvar1(:,:,k) = (forcingvar1(:,:,k+8)+ forcingvar1(:,:,k)+forcingvar1(:,:,k-8))/3;
%             
%         
%     elseif k > 11216 && k <= 11232  %%%% day 1403
% 
%         forcingvar1(:,:,k) = (forcingvar1(:,:,k+8)+ forcingvar1(:,:,k)+forcingvar1(:,:,k-8))/3;
% 
%             
%     elseif k > 12152 && k <= 12168  %%%% 1520
% 
%         forcingvar1(:,:,k) = (forcingvar1(:,:,k+8)+ forcingvar1(:,:,k)+forcingvar1(:,:,k-8))/3;
%  
%                   
%     elseif k > 13848 && k<=13864 %%%% 1732
% 
%         forcingvar1(:,:,k) = (forcingvar1(:,:,k+8)+ forcingvar1(:,:,k)+forcingvar1(:,:,k-8))/3;
% 
%             
%             
%     elseif k > 18216 && k<= 18232 %%% 2278
% 
%         forcingvar1(:,:,k) = (forcingvar1(:,:,k+8)+ forcingvar1(:,:,k)+forcingvar1(:,:,k-8))/3;
% 
%         
%     elseif k > 18384 && k<= 18400     
%         forcingvar1(:,:,k) = (forcingvar1(:,:,k+8)+ forcingvar1(:,:,k)+forcingvar1(:,:,k-8))/3;
%   
%         
%     elseif k > 24944 && k <= 24960
%         forcingvar1(:,:,k) = (forcingvar1(:,:,k+8)+ forcingvar1(:,:,k)+forcingvar1(:,:,k-8))/3;
% 
%         
%     elseif k > 25696 && k <= 25720   
%         forcingvar1(:,:,k) = (forcingvar1(:,:,k+8)+ forcingvar1(:,:,k)+forcingvar1(:,:,k-8))/3;
    if k > 16776 && k<=17000 %%%% 2098

        forcingvar1(:,:,k) = (forcingvar1(:,:,k+8)+ forcingvar1(:,:,k)+forcingvar1(:,:,k-8))/3;  
    end
end

writeDataset(forcingvar1,fullfile(inputfolder,aLW),ieee,prec);


clear forcingvar1

%%%%%%%%%%%


% %%%%%%%%%%%%%%


forcingvar1 = zeros(Nx,Ny,days);
fid = fopen(fullfile(inputfolder,aSW),'r','b');
for k=1:days
  forcingvar1(:,:,k) = fread(fid,[Nx Ny],'real*8');
end
fclose(fid);


for k = 1:days
%     if k >3792 && k <= 3816 %%% 476
% 
%         forcingvar1(:,:,k) = (forcingvar1(:,:,k-8) +forcingvar1(:,:,k)+ forcingvar1(:,:,k+8))/3;
% 
% 
%     elseif k >10392 && k <= 10408 %%% 1300
% 
%         forcingvar1(:,:,k) = (forcingvar1(:,:,k-8) +forcingvar1(:,:,k)+ forcingvar1(:,:,k+8))/3;
% 
%     elseif k > 10496 && k<= 10510  %%%% day 1313
%         forcingvar1(:,:,k) = (forcingvar1(:,:,k+8)+ forcingvar1(:,:,k)+forcingvar1(:,:,k-8))/3;
%             
%         
%     elseif k > 11216 && k <= 11232  %%%% day 1403
% 
%         forcingvar1(:,:,k) = (forcingvar1(:,:,k+8)+ forcingvar1(:,:,k)+forcingvar1(:,:,k-8))/3;
% 
%             
%     elseif k > 12152 && k <= 12168  %%%% 1520
% 
%         forcingvar1(:,:,k) = (forcingvar1(:,:,k+8)+ forcingvar1(:,:,k)+forcingvar1(:,:,k-8))/3;
%  
%                   
%     elseif k > 13848 && k<=13864 %%%% 1732
% 
%         forcingvar1(:,:,k) = (forcingvar1(:,:,k+8)+ forcingvar1(:,:,k)+forcingvar1(:,:,k-8))/3;
% 
%             
%             
%     elseif k > 18216 && k<= 18232 %%% 2278
% 
%         forcingvar1(:,:,k) = (forcingvar1(:,:,k+8)+ forcingvar1(:,:,k)+forcingvar1(:,:,k-8))/3;
% 
%         
%     elseif k > 18384 && k<= 18400     
%         forcingvar1(:,:,k) = (forcingvar1(:,:,k+8)+ forcingvar1(:,:,k)+forcingvar1(:,:,k-8))/3;
%   
%         
%     elseif k > 24944 && k <= 24960
%         forcingvar1(:,:,k) = (forcingvar1(:,:,k+8)+ forcingvar1(:,:,k)+forcingvar1(:,:,k-8))/3;
% 
%         
%     elseif k > 25696 && k <= 25720   
%         forcingvar1(:,:,k) = (forcingvar1(:,:,k+8)+ forcingvar1(:,:,k)+forcingvar1(:,:,k-8))/3;
    if k > 16776 && k<=17000 %%%% 2098

        forcingvar1(:,:,k) = (forcingvar1(:,:,k+8)+ forcingvar1(:,:,k)+forcingvar1(:,:,k-8))/3;  
    end
end


writeDataset(forcingvar1,fullfile(inputfolder,aSW),ieee,prec);


clear forcingvar1

%%%%%%%%%%%

%%%%%%%%%%%%%%

forcepath = fullfile(gendir,'/MITgcm_WS/newexp/differentResolutions/3hourforcing_sixth');

forcingvar1 = zeros(Nx,Ny,days);
fid = fopen(fullfile(inputfolder,aPrecip),'r','b');
for k=1:days
  forcingvar1(:,:,k) = fread(fid,[Nx Ny],'real*8');
end
fclose(fid);


for k = 1:days
%     if k >3792 && k <= 3816 %%% 476
% 
%         forcingvar1(:,:,k) = (forcingvar1(:,:,k-8) +forcingvar1(:,:,k)+ forcingvar1(:,:,k+8))/3;
% 
% 
%     elseif k >10392 && k <= 10408 %%% 1300
% 
%         forcingvar1(:,:,k) = (forcingvar1(:,:,k-8) +forcingvar1(:,:,k)+ forcingvar1(:,:,k+8))/3;
% 
%     elseif k > 10496 && k<= 10510  %%%% day 1313
%         forcingvar1(:,:,k) = (forcingvar1(:,:,k+8)+ forcingvar1(:,:,k)+forcingvar1(:,:,k-8))/3;
%             
%         
%     elseif k > 11216 && k <= 11232  %%%% day 1403
% 
%         forcingvar1(:,:,k) = (forcingvar1(:,:,k+8)+ forcingvar1(:,:,k)+forcingvar1(:,:,k-8))/3;
% 
%             
%     elseif k > 12152 && k <= 12168  %%%% 1520
% 
%         forcingvar1(:,:,k) = (forcingvar1(:,:,k+8)+ forcingvar1(:,:,k)+forcingvar1(:,:,k-8))/3;
%  
%                   
%     elseif k > 13848 && k<=13864 %%%% 1732
% 
%         forcingvar1(:,:,k) = (forcingvar1(:,:,k+8)+ forcingvar1(:,:,k)+forcingvar1(:,:,k-8))/3;
% 
%             
%             
%     elseif k > 18216 && k<= 18232 %%% 2278
% 
%         forcingvar1(:,:,k) = (forcingvar1(:,:,k+8)+ forcingvar1(:,:,k)+forcingvar1(:,:,k-8))/3;
% 
%         
%     elseif k > 18384 && k<= 18400     
%         forcingvar1(:,:,k) = (forcingvar1(:,:,k+8)+ forcingvar1(:,:,k)+forcingvar1(:,:,k-8))/3;
%   
%         
%     elseif k > 24944 && k <= 24960
%         forcingvar1(:,:,k) = (forcingvar1(:,:,k+8)+ forcingvar1(:,:,k)+forcingvar1(:,:,k-8))/3;
% 
%         
%     elseif k > 25696 && k <= 25720   
%         forcingvar1(:,:,k) = (forcingvar1(:,:,k+8)+ forcingvar1(:,:,k)+forcingvar1(:,:,k-8))/3;
    if k > 16776 && k<=17000 %%%% 2098

        forcingvar1(:,:,k) = (forcingvar1(:,:,k+8)+ forcingvar1(:,:,k)+forcingvar1(:,:,k-8))/3;  
    end
end



writeDataset(forcingvar1,fullfile(inputfolder,aPrecip),ieee,prec);


clear forcingvar1

%%%%%%%%%%%



forcingvar1 = zeros(Nx,Ny,days);
fid = fopen(fullfile(inputfolder,anewAQ),'r','b');
for k=1:days
  forcingvar1(:,:,k) = fread(fid,[Nx Ny],'real*8');
end
fclose(fid);


for k = 1:days
%     if k >3792 && k <= 3816 %%% 476
% 
%         forcingvar1(:,:,k) = (forcingvar1(:,:,k-8) +forcingvar1(:,:,k)+ forcingvar1(:,:,k+8))/3;
% 
% 
%     elseif k >10392 && k <= 10408 %%% 1300
% 
%         forcingvar1(:,:,k) = {forcingvar1(:,:,k-8) +forcingvar1(:,:,k)+ forcingvar1(:,:,k+8))/3;
% 
%     elseif k > 10496 && k<= 10510  %%%% day 1313
%         forcingvar1(:,:,k) = (forcingvar1(:,:,k+8)+ forcingvar1(:,:,k)+forcingvar1(:,:,k-8))/3;
%             
%         
%     elseif k > 11216 && k <= 11232  %%%% day 1403
% 
%         forcingvar1(:,:,k) = (forcingvar1(:,:,k+8)+ forcingvar1(:,:,k)+forcingvar1(:,:,k-8))/3;
% 
%             
%     elseif k > 12152 && k <= 12168  %%%% 1520
% 
%         forcingvar1(:,:,k) = (forcingvar1(:,:,k+8)+ forcingvar1(:,:,k)+forcingvar1(:,:,k-8))/3;
%  
%                   
%     elseif k > 13848 && k<=13864 %%%% 1732
% 
%         forcingvar1(:,:,k) = (forcingvar1(:,:,k+8)+ forcingvar1(:,:,k)+forcingvar1(:,:,k-8))/3;
% 
%             
%             
%     elseif k > 18216 && k<= 18232 %%% 2278
% 
%         forcingvar1(:,:,k) = (forcingvar1(:,:,k+8)+ forcingvar1(:,:,k)+forcingvar1(:,:,k-8))/3;
% 
%         
%     elseif k > 18384 && k<= 18400     
%         forcingvar1(:,:,k) = (forcingvar1(:,:,k+8)+ forcingvar1(:,:,k)+forcingvar1(:,:,k-8))/3;
%   
%         
%     elseif k > 24944 && k <= 24960
%         forcingvar1(:,:,k) = (forcingvar1(:,:,k+8)+ forcingvar1(:,:,k)+forcingvar1(:,:,k-8))/3;
% 
%         
%     elseif k > 25696 && k <= 25720   
%         forcingvar1(:,:,k) = (forcingvar1(:,:,k+8)+ forcingvar1(:,:,k)+forcingvar1(:,:,k-8))/3;
    if k > 16776 && k<=17000 %%%% 2098

        forcingvar1(:,:,k) = (forcingvar1(:,:,k+8)+ forcingvar1(:,:,k)+forcingvar1(:,:,k-8))/3;  
    end
end



writeDataset(forcingvar1,fullfile(inputfolder,anewAQ),ieee,prec);


clear forcingvar1




