
%%%%%%%%%%% forcing movie



run ../newexp/defineGrid.m

days =3287;
inputpath = fullfile(gendir,'/MITgcm_WS/experiments/n_342/input');
inputpath2 = fullfile(gendir,'/MITgcm_WS/experiments/n_34452/input');

forcingvar = zeros(Nx,Ny,days);
fid = fopen(fullfile(inputpath,mwind),'r','b');
for k=1:days
  forcingvar(:,:,k) = fread(fid,[Nx Ny],'real*8');
end
fclose(fid);



forcingvar2 = zeros(Nx,Ny,days);
fid = fopen(fullfile(inputpath2,mwind),'r','b');
for k=1:days
  forcingvar2(:,:,k) = fread(fid,[Nx Ny],'real*8');
end
fclose(fid);
stardt = 1100;
ender=1300;
fontsize = 12;
set(gcf,'color','w')
time = ender-stardt;

M = moviein(time);



for n= stardt:ender
% for n= 1:18
%    t= n/8;

 t = n;
  

  pcolor(XMC,YMC,(forcingvar(:,:,n)-forcingvar2(:,:,n))'),shading interp

  ['Mx value: ',num2str(max(max(max(forcingvar(:,:,n)))))]
  ['Max value: ',num2str(max(max(max(forcingvar2(:,:,n)))))]
  
  
       

%    [idx1 idx2] = find(isnan(forcingvar));
  
  
  
  
  
  
  
  
  
  
  
  
    
  %%% Finish the plot
  handle=colorbar;
  colormap jet(200);
  caxis([-20 10]);
  set(handle,'FontSize',fontsize);
  title(t);

  
  set(gca,'FontSize',fontsize);
%   annotation('textbox',[0.85 0.05 0.25 0.05],'String','$\overline{\theta}$ ($^\circ$C)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
  M(n) = getframe(gcf);
  
end




