%%%
%%% doubleRes.m
%%%
%%% Takes an MITgcm experiment and doubles its resolution, producing input
%%% files for a new experiment at the higher resolution.
%%%
%%% The resolution-doubling process used in this file is approximate - simple
%%% nearest-neighbor interpolation that is likely to produce spurious results.
%%% The doubled-resolution solution will likely contain transient artifacts
%%% after initialization, so some period of integration will be required to
%%% remove these.
%%%

%%% For file I/O
addpath ../newexp_utils/
addpath ../utils/matlab

%%% Load experiment
expdir = '/Volumes/LaCie/UCLA/Projects/MITgcm_SM/experiments';
expname = 'SM_prod_res8km_spinup';
expiter = 1808948;
Nx = 100;
Ny = 64;
Nr = 70;

%%% Formatting
ieee='b';
prec='real*8';

%%% Pull out u,v,t,s from pickup file
A = rdmds(fullfile(expdir,expname,'results/pickup'),expiter,'siTICES');
uvtse1 = A(:,:,[1:4*Nr 8*Nr+1]);

%%% Create double-resolution arrays. Here we basically just use 
%%% nearest-neighbour interpolation because it doesn't need to be 
%%% a precise doubling of the resolution.
uvtse2 = zeros(2*Nx,2*Ny,4*Nr+1);
for i=1:Nx
  for j=2:Ny-1
    for k=1:4*Nr+1
      if (uvtse1(i,j,k)==0)
        uvtse1(i,j,k) = uvtse1(i,j,k-1);
      end
      uvtse2(2*i-1,2*j-1,k) = uvtse1(i,j,k);
      uvtse2(2*i-1,2*j,k) = uvtse1(i,j,k);
      uvtse2(2*i,2*j-1,k) = uvtse1(i,j,k);
      uvtse2(2*i,2*j,k) = uvtse1(i,j,k);      
    end
  end
end

%%% It turns out we need to fill wet cells with non-zero values 
uvtse2(:,2,:) = uvtse2(:,3,:);
uvtse2(:,2*Ny-1,:) = uvtse2(:,2*Ny-2,:); 

%%% Create input arrays
writeDataset(uvtse2(:,:,1:Nr),'./DEFAULTS/input/uVelInitFile.bin',ieee,prec);
writeDataset(uvtse2(:,:,Nr+1:2*Nr),'./DEFAULTS/input/vVelInitFile.bin',ieee,prec);
writeDataset(uvtse2(:,:,2*Nr+1:3*Nr),'./DEFAULTS/input/hydrogThetaFile.bin',ieee,prec);
writeDataset(uvtse2(:,:,3*Nr+1:4*Nr),'./DEFAULTS/input/hydrogSaltFile.bin',ieee,prec);
writeDataset(uvtse2(:,:,4*Nr+1),'./DEFAULTS/input/pSurfInitFile.bin',ieee,prec);



