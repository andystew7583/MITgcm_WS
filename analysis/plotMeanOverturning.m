%%%
%%% plotMeanOverturning.m
%%%
%%% Plots the Eulerian-Mean overturning circulation.
%%%

%%% Load experiment
expdir = '../experiments';
expname = 'hires_seq_onetwelfth_notides_RTOPO2';

%%% Options (see calcMeanOverturning)
deform_cavity = false;
% psimax = 6;
% psistep = 0.5;
psimax = 2;
psistep = 0.1;

%%% Construct output file name
outfname = [expname,'_EMOC'];
if (deform_cavity)
  outfname = [outfname,'_deform'];
end
outfname = [outfname,'.mat'];

%%% Store computed data for later
load(fullfile('products',outfname));

%%% Make the plot
psi_plot = mean(psi_EM,3)/1e6;
psi_plot = sign(psi_plot).*min(abs(psi_plot),psimax);
psi_plot(abs(psi_plot)<1e-3) = NaN;
figure(1);
clf;
set(gcf,'Position',[325         460        1023         398]);
[ZZ,EE] = meshgrid(RF,eta);
contourf(EE,-ZZ,psi_plot,[-psimax:psistep:-psistep psistep:psistep:psimax],'EdgeColor','k');
caxis([-psimax psimax]);
set(gca,'YDir','reverse');
set(gca,'YLim',[0 -RF(end)]);
colormap(cmocean('balance',round(2*psimax/psistep)));
colorbar;
xlabel('\eta');
ylabel('Potential density (kg/m^3)');
title('Eulerian-Mean Overturning streamfunction (Sv)');
set(gca,'Position',[0.0787    0.1300    0.8383    0.7950]);
set(gca,'FontSize',12);