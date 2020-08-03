%%%%% Plot T/S in FRESH and REF experiments 4 panels (2 full domain, 2 fris
%%%%% cavity)

%Script to plot Temperature and salinity
setExpname
%%% Read experiment data
loadexp;

load SavedFiles/TSanalysis4panel.mat


%%% Frequency of diagnostic output
dumpFreq = abs(diag_frequency(15));
deltaT = 440;
% nIter0 = 0;
nDumps = round(nTimeSteps*(deltaT/dumpFreq));
dumpIters = round((1:nDumps)*dumpFreq/deltaT);
dumpIters = dumpIters(dumpIters > nIter0);
nDumps = length(dumpIters);
 
deltaT_4 = 200;
nIter0_4 = 1;
nDumps_4 = round(nTimeSteps*(deltaT_4/dumpFreq));
dumpIters_4 = round((1:nDumps_4)*dumpFreq/deltaT_4); 
dumpIters_4 = dumpIters_4(dumpIters_4 >= nIter0_4);



deltaT_5 = 200;
nIter0_5 = 2877120;
nDumps_5 = round(nTimeSteps*10*(deltaT_5/dumpFreq));
dumpIters_5 = round((1:nDumps_5)*dumpFreq/deltaT_5); 
dumpIters_5 = dumpIters_5(dumpIters_5 >= nIter0_5);

deltaT2 = 400;
nDumps2 = round(nTimeSteps*10*(deltaT2/dumpFreq));
dumpIters_2 = round((1:nDumps2)*dumpFreq/deltaT2); 
dumpIters2 = dumpIters_2(dumpIters_2 >= nIter0);
nDumps2 = length(dumpIters2);

%%%% Load Salinity and Temperature



%%%% load REF and FRESH data

exppath1 = '/data3/MITgcm_WS/experiments/n_34452';
exppath2 = '/data3/MITgcm_WS/experiments/n_342';

% tmin = 9*86400*360;
% tmax = 18*86400*360;
% 
% Temp_c = readIters(exppath1,'THETA',dumpIters,deltaT,tmin,tmax,Nx,Ny,Nr);
% Temp_w = readIters(exppath2,'THETA',dumpIters,deltaT,tmin,tmax,Nx,Ny,Nr);
% Salt_c = readIters(exppath1,'SALT',dumpIters,deltaT,tmin,tmax,Nx,Ny,Nr);
% Salt_w = readIters(exppath2,'SALT',dumpIters,deltaT,tmin,tmax,Nx,Ny,Nr);
% 



% %%%finding volume of all cells in grid for normalization purposes      
% volume = NaN(Nx,Ny,Nr);
% % nv=NaN(Nx,Ny,Nr);
% for i = 1:Nx
%     for j = 1:Ny
%         for k = 1:Nr
%             
%             volume(i,j,k) = (RAC(i,j)*DRF(k)*hFacC(i,j,k));
%             if volume(i,j,k)>0
%                 nv(i,j,k)=(volume(i,j,k));
%             end
%         end
%     end
% end
% 
% 
% %%%% index the FRIS area (-80 to -20) W, 74S -> !
% Temp_cn = zeros(Nx,Ny,Nr);
% Temp_wn = zeros(Nx,Ny,Nr);
% Salt_wn = zeros(Nx,Ny,Nr);
% Salt_cn = zeros(Nx,Ny,Nr);
% new_vol = zeros(Nx,Ny,Nr);
% 
% 
% 
%   for i=1:Nx
%     for j=1:Ny
%         for k = 1:Nr
%             if XC(i,j)>-80 && XC(i,j)<-30 && YC(i,j)<-75
% 
%                     if SHELFICEtopo(i,j)<0
% %                 
%                         Temp_cn(i,j,k) = Temp_c(i,j,k);
%                         Salt_cn(i,j,k) = Salt_c(i,j,k);
%                         
%                         Temp_wn(i,j,k) = Temp_w(i,j,k);
%                         Salt_wn(i,j,k) = Salt_w(i,j,k);                       
%                         new_vol(i,j,k) = volume(i,j,k);
%                     end
%                  
%              end
%                  
%              
%                   
%             
%         end
%     end
%   end
% 
% % 
% % Nx_new = size(Temp_cn,1);
% % Ny_new = size(Temp_cn,2);
% % theta_cd = reshape(Temp_cn,[1],Nx_new*Ny_new*Nr);
% % salt_cd = reshape(Salt_cn,[1],Nx_new*Ny_new*Nr);
% % 
% % theta_wd = reshape(Temp_wn,[1],Nx_new*Ny_new*Nr);
% % salt_wd = reshape(Salt_wn,[1],Nx_new*Ny_new*Nr);
% % 
% % new_vol = reshape(new_vol,[1],Nx_new*Ny_new*Nr);
% % 
% % 
% % 
% % 
% % 
% % % 
% % new_vol(new_vol ==0) = [];
% % theta_wd(theta_wd ==0) = [];
% % salt_wd(salt_wd ==0) = [];
% % 
% % theta_cd(theta_cd ==0) = [];
% % salt_cd(salt_cd ==0) = [];
% % 
% % new_vol=log(new_vol);
% % 
% salt_grid=33:.005:37;
% temp_grid=-3:.005:2;
% new_v_3445=zeros(size(salt_grid,2),size(temp_grid,2));
% new_v_34=zeros(size(salt_grid,2),size(temp_grid,2));
% 
% for i = 1:Nx
%     for j = 1:Ny
%        for k = 1:Nr
%                if Salt_cn(i,j,k)>0
%                 dist_s3445    = abs(Salt_cn(i,j,k) - salt_grid);
%                 minDist_s3445 = min(dist_s3445);
%                 idx_s3445     = find(dist_s3445 == minDist_s3445);
%             
%                 dist_s34    = abs(Salt_wn(i,j,k) - salt_grid);
%                 minDist_s34 = min(dist_s34);
%                 idx_s34     = find(dist_s34 == minDist_s34);
%                end
%                if Temp_cn(i,j,k)<0 || Temp_cn(i,j,k)>0
%            
%                 dist_t_3445    = abs(Temp_cn(i,j,k) - temp_grid);
%                 minDist_t3445 = min(dist_t_3445);
%                 idx_t3445     = find(dist_t_3445 == minDist_t3445);    
%             
%                 dist_t_34    = abs(Temp_wn(i,j,k) - temp_grid);
%                 minDist_t34 = min(dist_t_34);
%                 idx_t34     = find(dist_t_34 == minDist_t34);            
%                end
%              
%                 new_v_3445((idx_s3445),(idx_t3445))=new_v_3445((idx_s3445),(idx_t3445))+new_vol(i,j,k); 
%                 new_v_34((idx_s34),(idx_t34))=new_v_34((idx_s34),(idx_t34))+new_vol(i,j,k); 
%              
%            
%        end
%     end
% end
% 
% new_v_3445(new_v_3445==0)=NaN;
% new_v_34(new_v_34==0)=NaN;
% 
% new_v_3445=log(new_v_3445);
% new_v_34=log(new_v_34);
% 
% 
% 
% %%%meshgrid of t,s points as defined by out grids
% [s,t]=meshgrid(salt_grid,temp_grid);

%%% load CDW and HSSW values for FRESH and REF
load SavedFiles/CDWHSSWvals_343445.mat

figure(1)
clf
set(gcf,'Position',[406         168        1221         814]);

r1=subplot(2,2,1);
set(gca,'Position',[0.1013    0.6070    0.3653    0.3506]);
pcolor(s,t,(new_v_3445c'));
shading interp;
colormap jet(40);

% ylabel(m,'Log (Volume) m$^3$','interpreter','latex','fontsize',20)
xlim([33.7 35]);
xticks(33.7:.2:36);
ylim([-3 1.5]);
caxis([20 28]);
% m = colorbar;
% ylabel(m,'Log (Volume) m$^3$','interpreter','latex','fontsize',20);

title('REF, full WSRM domain', 'interpreter','latex','FontSize',16);

b = gca; legend(b,'off');
pt_maxc = 1.5;
pt_minc = -3;
ss_maxc = 37;
ss_minc = 30;

%%% Grid for contouring density
pt_stepc = (pt_maxc-pt_minc)/100;
pt_gridc = pt_minc:pt_stepc:pt_maxc;
ss_stepc = (ss_maxc-ss_minc)/100;
ss_gridc = ss_minc:ss_stepc:ss_maxc;
[PT_gridc,SS_gridc] = meshgrid(pt_gridc,ss_gridc);

%%% Calculate potential density
pdc = densmdjwf(SS_gridc,PT_gridc,0) - 1000;
hold on
[C,f] = contour(SS_gridc,PT_gridc,pdc,25.0:.1:35.0,'EdgeColor','k');
clabel(C,f)

cmap=jet;
text(33.85,-4,'\textbf{a}','fontsize',14,'fontweight','bold','interpreter','latex');
set(gca,'FontSize',16);

subplot(2,2,2);
set(gca,'Position',[0.5340    0.6070    0.3653    0.3506]);
pcolor(s,t,new_v_34w');
shading interp
% m2 = colorbar;
colormap jet(40);

xlim([33.7 35]);
xticks(33.7:.2:36);
ylim([-3 1.5]);
caxis([20 28]);
title('FRESH, full WSRM domain', 'interpreter','latex','FontSize',16);

b = gca; legend(b,'off');


pt_max = 2;
pt_min = -3;
ss_max = 37;
ss_min = 30;

%%% Grid for contouring density
pt_step = (pt_max-pt_min)/100;
pt_grid = pt_min:pt_step:pt_max;
ss_step = (ss_max-ss_min)/100;
ss_grid = ss_min:ss_step:ss_max;
[PT_grid,SS_grid] = meshgrid(pt_grid,ss_grid);

%%% Calculate potential density
pd = densmdjwf(SS_grid,PT_grid,0) - 1000;
hold on
[C,f] = contour(SS_grid,PT_grid,pd,25.0:.1:30.0,'EdgeColor','k');
clabel(C,f)
% 
% 
text(33.85,-4,'\textbf{b}','fontsize',14,'fontweight','bold','interpreter','latex');
set(gca,'FontSize',16);

subplot(2,2,3);
set(gca,'Position',[0.1013    0.135    0.3653    0.3506]);
pcolor(s,t,new_v_3445');
shading interp
% m2 = colorbar;
colormap jet(40);
% l = (min(new_vol));
% h = (max(new_vol));
% caxis([l 24]);
xlim([33.7 35]);
xticks(33.7:.2:36);
ylim([-3 1.5]);
caxis([20 28]);
title('REF, FRIS cavity only', 'interpreter','latex','FontSize',16);
% xlabel('Salinity (psu)', 'interpreter','latex','FontSize',18);
% ylabel('Temperature ($^{o}$C)', 'interpreter','latex','FontSize',18);
b = gca; legend(b,'off');
set(gca,'FontSize',16);

pt_max = 1.5;
pt_min = -3;
ss_max = 37;
ss_min = 33;

%%% Grid for contouring density
pt_step = (pt_max-pt_min)/100;
pt_grid = pt_min:pt_step:pt_max;
ss_step = (ss_max-ss_min)/100;
ss_grid = ss_min:ss_step:ss_max;
[PT_grid,SS_grid] = meshgrid(pt_grid,ss_grid);

%%% Calculate potential density
pd = densmdjwf(SS_grid,PT_grid,1000) - 1000;
hold on
[C,f] = contour(SS_grid,PT_grid,pd,25.0:.1:35.0,'EdgeColor','k');
clabel(C,f)
% 
% 
text(33.85,-4,'\textbf{c}','fontsize',14,'fontweight','bold','interpreter','latex');
plot(Salt_diff_control_hssw,Temp_diff_control_hssw,'b+','markersize',12,'LineWidth',2)
plot(Salt_diff_control_cdw,Temp_diff_control_cdw,'r+','markersize',12,'LineWidth',2)
set(gca,'FontSize',16);

subplot(2,2,4);
set(gca,'Position',[0.5340    0.135    0.3653    0.3506]);
pcolor(s,t,new_v_34');
shading interp
colormap jet(40);
% l = (min(new_vol));
% h = (max(new_vol));
% caxis([l 24]);
% ylabel(m2,'Log (Volume) m$^3$','interpreter','latex','fontsize',20)
xlim([33.7 35]);
xticks(33.7:.2:36);
ylim([-3 1.5]);
caxis([20 28]);
title('FRESH, FRIS cavity only', 'interpreter','latex','FontSize',16);

b = gca; legend(b,'off');


pt_max = 2;
pt_min = -3;
ss_max = 37;
ss_min = 30;

%%% Grid for contouring density
pt_step = (pt_max-pt_min)/100;
pt_grid = pt_min:pt_step:pt_max;
ss_step = (ss_max-ss_min)/100;
ss_grid = ss_min:ss_step:ss_max;
[PT_grid,SS_grid] = meshgrid(pt_grid,ss_grid);

%%% Calculate potential density
pd = densmdjwf(SS_grid,PT_grid,1000) - 1000;
hold on
[C,f] = contour(SS_grid,PT_grid,pd,30.0:.1:35.0,'EdgeColor','k');
clabel(C,f)
% 
% 
text(33.85,-4,'\textbf{d}','fontsize',14,'fontweight','bold','interpreter','latex');
% plot(Salt_diff_wcontrol_whssw,Temp_diff_control_whssw,'b+','markersize',12);
% plot(Salt_diff_wcontrol_wcdw,Temp_diff_control_wcdw,'r+','markersize',12);
plot(Salt_diff_control_hssw,Temp_diff_control_hssw,'b+','markersize',12,'LineWidth',2)
plot(Salt_diff_control_cdw,Temp_diff_control_cdw,'r+','markersize',12,'LineWidth',2)
set(gca,'FontSize',16);
m = colorbar;
set(m,'Position',[.92,.135,.015,.8226])
ylabel(m,'Log(m$^{3})$','interpreter','latex','fontsize',16);


[a,h]=suplabel('Temperature($^\circ$C)','y');
set(h,'interpreter','latex','fontsize',20);
set(h,'Position',[-0.0280    0.5000    0.0000]);

[a,h]=suplabel('Salinity(psu)','x');
set(h,'interpreter','latex','fontsize',20);

