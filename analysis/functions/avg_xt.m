%%%
%%% avg_xt.m
%%%
%%% Calculates the time and zonal average of the output fields from MITgcm 
%%% runs. Must be called after time-averaged fields have been calculated.
%%%
%%% NOTE: This script is really out of date now, as this kind of averaging
%%% really only makes sense in a zonally symmetric domain,
%%%

%%% Calculate zonal mean
uu_avg = zeros(Ny,Nr);
vv_avg = zeros(Ny,Nr);
ww_avg = zeros(Ny,Nr);
tt_avg = zeros(Ny,Nr);
vt_avg = zeros(Ny,Nr);
wt_avg = zeros(Ny,Nr);
usq_avg = zeros(Ny,Nr);
vsq_avg = zeros(Ny,Nr);
wsq_avg = zeros(Ny,Nr);
tsq_avg = zeros(Ny,Nr);
L_wet = zeros(Ny,Nr);
for j=1:Ny
  for k=1:Nr
    L_wet(j,k) = sum(delX'.*(hFacC(:,j,k)~=0));
    if (L_wet(j,k) == 0)
      continue;
    end
    uu_avg(j,k) = sum(uu(:,j,k).*delX')/L_wet(j,k);
    vv_avg(j,k) = sum(vv(:,j,k).*delX')/L_wet(j,k);
    ww_avg(j,k) = sum(ww(:,j,k).*delX')/L_wet(j,k);
    tt_avg(j,k) = sum(tt(:,j,k).*delX')/L_wet(j,k);
    vt_avg(j,k) = sum(vt(:,j,k).*delX')/L_wet(j,k);
    wt_avg(j,k) = sum(wt(:,j,k).*delX')/L_wet(j,k);      
    usq_avg(j,k) = sum(usq(:,j,k).*delX')/L_wet(j,k);
    vsq_avg(j,k) = sum(vsq(:,j,k).*delX')/L_wet(j,k);
    wsq_avg(j,k) = sum(wsq(:,j,k).*delX')/L_wet(j,k);
    tsq_avg(j,k) = sum(tsq(:,j,k).*delX')/L_wet(j,k);        
  end
end