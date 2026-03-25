%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Calc vFlux/uFlux
%%%%%%%%%% Trying to verify OBCS mass fluxes =0 across each boundary
%%%%%%%%%% gridpoint


%%% Load experiment
loadexp;

%%% Frequency of diagnostic output - should match that specified in
%%% data.diagnostics.
dumpFreq = abs(diag_frequency(1));
nDumps = round(nTimeSteps*deltaT/dumpFreq);
dumpIters = round((1:nDumps)*dumpFreq/deltaT);
dumpIters = dumpIters(dumpIters > nIter0);
nDumps = length(dumpIters);




%%%%%%%%% Read in Files
% vvel_tavg = zeros(Nx,Ny,Nr);
% uvel_tavg = zeros(Nx,Ny,Nr);
% navg = 0;
% for n=1:length(dumpIters)
%  
%   tyears = dumpIters(n)*deltaT/86400/365;
%  
%     [tyears  dumpIters(n)];
%     uvel  = rdmdsWrapper(fullfile(exppath,'/results/UVEL_inst'),dumpIters(n));              
%     vvel  = rdmdsWrapper(fullfile(exppath,'/results/VVEL_inst'),dumpIters(n));              
%     if (isempty(vvel) || isempty(uvel))
%       break;
%     else
%     
%     vvel_tavg = vvel_tavg + squeeze(vvel(:,:,:,1));  
%     uvel_tavg = uvel_tavg + squeeze(uvel(:,:,:,1)); 
%     navg = navg + 1;
%     end
% end

%%%%%%%%%% Eastern Boundary

EBsum = (12);
for k=1:12;
EBsum(k) = sum(sum(OBEu(:,:,k).*squeeze(hFacW(end,:,:)).*repmat(DYG(end,:)',[1 Nr]).*repmat(reshape(DRF,[1 Nr]),[Ny 1])))+sum(sum(OBNv(:,:,k).*squeeze(hFacS(:,end,:)).*repmat(DXG(:,end),[1 Nr]).*repmat(reshape(DRF,[1 Nr]),[Nx 1])));
end

NBsum = (12);
for k=1:12
NBsum(k) = sum(sum(OBNv(:,:,k).*squeeze(hFacS(:,end,:)).*repmat(DXG(:,end),[1 Nr]).*repmat(reshape(DRF,[1 Nr]),[Nx 1])))+sum(sum(OBEu(:,:,k).*squeeze(hFacW(end,:,:)).*repmat(DYG(end,:)',[1 Nr]).*repmat(reshape(DRF,[1 Nr]),[Ny 1])));
end




%%%%%%%%%% Calculate time AvErAge
% vvel_tavg = vvel_tavg/navg;
% uvel_tavg = uvel_tavg/navg;


%%%%%%%%% Calculate Fluxes

%%%%%Eastern Boundary
% uFlux_Eastern = NaN(1,Ny);
% for k = size(vvel_tavg,3)
%     uFlux_Eastern = uvel_tavg(end,:,k).*hFacW(end,:,k).*DXG(end,:).*DRF(:,:,k);
% end
% 
% %%%%%%Northern Boundary
% vFlux_Northern = NaN(Nx,1);
% for k = size(vvel_tavg,3)
%     vFlux_Northern = vvel_tavg(:,end,k).*hFacS(:,end,k).*DYG(:,end).*DRF(:,:,k);
% end



