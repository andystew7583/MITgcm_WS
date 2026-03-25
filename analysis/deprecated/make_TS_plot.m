function make_TS_plot (salt,pot_temp,depths,p_ref,pot_dens_contours)
%%%
%%% USAGE: make_TS_plot (salt, pot_temp, depths, p_ref, pot_dens_contours)
%%%
%%% Makes a temperature/salinity plot using the 2D salinity/potential 
%%% temperature data.
%%%
%%% Arguments:
%%% salt - salinity matrix
%%% pot_temp - potential temperature matrix
%%% depths - vector of depths
%%% p_ref - reference pressure for contouring density lines
%%% pot_dens_contours - vector of potential density lines to plot
%%%    
%%% salt and pot_temp must be the same size (M x N matrices). depths must 
%%% be a vector of length N. p_ref must be a positive pressure in decibars.
%%%

%%% Error checking
if (size(salt) ~= size(pot_temp))
  error('make_TS_plot: Sizes of salinity and potential temperature input data must be the same');
end
if (length(depths) ~= size(salt,2))
  error('make_TS_plot: Length of depth data must match number of columns of salinity and temperature data');
end
if (p_ref < 0)
  error('make_TS_plot: Reference pressure must be >= 0');
end

%%% Easier variable names
pt = pot_temp;
ss = salt;
Ny = size(ss,1);
Nz = size(ss,2);

%%% Get ranges of T/S
pt_max = max(max(pt)) + 1;
pt_min = min(min(pt)) - 1;
ss_max = max(max(ss)) + 0.1;
ss_min = min(min(ss)) - 0.1;

%%% Grid for contouring density
pt_step = (pt_max-pt_min)/100;
pt_grid = pt_min:pt_step:pt_max;
ss_step = (ss_max-ss_min)/100;
ss_grid = ss_min:ss_step:ss_max;
[PT_grid,SS_grid] = meshgrid(pt_grid,ss_grid);

%%% Calculate potential density
pd = densmdjwf(SS_grid,PT_grid,p_ref) - 1000;

%%% Plotting options
scrsz = get(0,'ScreenSize');
fontsize = 12;

plotloc = [0.15 0.15 0.7 0.76];
framepos = [100    500   800  800];
%%% Scatter the T/S data
handle = figure;
set(handle,'Position',framepos);
clf;
%ss
for j=1:10:Ny
  scatter(ss(j,:),pt(j,:),4,-depths);
  if (j==1)
    hold on;
  end
end
[C,h] = contour(SS_grid,PT_grid,pd,25.0:.1:30.0,'EdgeColor','k');
clabel(C,h);
hold off;
xlabel('Salinity','interpreter','latex','FontSize',fontsize+2);
ylabel('Potential temperature ($^\circ$C)','interpreter','latex','FontSize',fontsize+2);
set(gca,'Position',plotloc);
handle = colorbar;
set(handle,'Position',[0.92 0.1 0.02 .85]);
colormap(flipdim(jet(100),2));
set(handle,'FontSize',fontsize);
caxis([min(-depths) max(-depths)]);

axis([ss_min ss_max pt_min pt_max]);
annotation('textbox',[0.88 0.05 0.3 0.05],'String','Depth (m)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
set(gca,'FontSize',fontsize);

end

