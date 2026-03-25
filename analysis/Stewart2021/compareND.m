%%%
%%% compareND.m
%%%
%%% Compares different neutral densities with potential density.
%%%


%%% Load experiment data
expdir = '../experiments';
expname = 'hires_seq_onethird_RTOPO2';
loadexp;

%%% Quasi-latitude horizontal grid
ETA = defineMOCgrid (XC,YC,SHELFICEtopo,bathy,false);
ETA = repmat(ETA,[1 1 Nr]);

%%% Load pre-computed neutral density variables
load(fullfile('products',[expname,'_ND1.mat']));
ND1_ref = gg_ref;
load(fullfile('products',[expname,'_ND2.mat']));
ND2_ref = gg_ref;

%%% Compute potential density
PD0 = densjmd95(ss_ref,pt_ref,-gravity*rhoConst*RC(1)/1e4*ones(Nx,Ny,Nr))-1000;

%%% Mask for region in which JM97 neutral density variable is well defined
msk = (YC>-80) & (XC>-64);
msk = repmat(msk,[1 1 Nr]);
msk = msk*1.0;
msk(msk==0) = NaN;
ND1_ref_msk = ND1_ref.*msk;
ND2_ref_msk = ND2_ref.*msk;

%%% Compute stratifications
pp_mid = 0.5*(pp_ref(:,:,1:Nr-1)+pp_ref(:,:,2:Nr));
DRC = repmat(RC(1:Nr-1)-RC(2:Nr),[Nx Ny 1]);
Nsq_true = - gravity/rhoConst * (...
                densjmd95(ss_ref(:,:,1:Nr-1),pt_ref(:,:,1:Nr-1),pp_mid) ...
              - densjmd95(ss_ref(:,:,2:Nr),pt_ref(:,:,2:Nr),pp_mid) ) ...
            ./ DRC;
Nsq_ND1 = -gravity/rhoConst * (ND1_ref(:,:,1:Nr-1)-ND1_ref(:,:,2:Nr)) ./DRC; 
Nsq_ND2 = -gravity/rhoConst * (ND2_ref(:,:,1:Nr-1)-ND2_ref(:,:,2:Nr)) ./DRC; 

ND1min = 26.5;
ND1max = 29;
ND1step = 0.01;
PD0min = 26.9;
PD0max = 28.1;
PD0step = 0.01;
ND1PD0vol = binByVolume(ND1_ref,PD0,ones(size(ND1_ref)), ...
                  ND1min,ND1max,ND1step,PD0min,PD0max,PD0step, ...
                  RAC,DRF,hFacC);
ND1PD0vol(ND1PD0vol==0) = NaN;                
[PPP,NNN] = meshgrid(PD0min:PD0step:PD0max,ND1min:ND1step:ND1max);
figure(47);
pcolor(NNN,PPP,log10(ND1PD0vol));
shading flat;
colormap jet;
colorbar

                


figure(40);
scatter(ND1_ref(:),ND2_ref(:),10,ETA(:));
colorbar;

figure(41);
scatter(ND1_ref_msk(:),ND2_ref_msk(:),10,ETA(:));
colorbar;

figure(42);
scatter(ND1_ref_msk(:),PD0(:),10,ETA(:));
colorbar;
 
figure(43);
scatter(log10(Nsq_true(:)),log10(Nsq_ND1(:)));


figure(44);
scatter(log10(Nsq_true(:)),log10(Nsq_ND2(:)));

figure(45);
scatter(log10(Nsq_true(:)),Nsq_ND1(:)./Nsq_true(:));
