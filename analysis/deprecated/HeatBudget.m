%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Heat Budget Analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

loadexp

%%%%%%%%%%% Read in Experiment Variables For the Heat Budget
%%% Frequency of diagnostic output - should match that specified in
%%% data.diagnostics.
dumpFreq = abs(diag_frequency(1));
nDumps = round(nTimeSteps*deltaT/dumpFreq);
dumpIters = round((1:nDumps)*dumpFreq/deltaT);
dumpIters = dumpIters(dumpIters > nIter0);
nDumps = length(dumpIters);

globalArea = 3.288196794800415e12;



CellVol=NaN(Nx,Ny,Nr);
for i = 1:Nx
    for j = 1:Ny
        for k = 1:Nr
            CellVol(i,j,k) = RAC(i,j).*DRF(k).*hFacC(i,j,k);
        end
    end
end

%%%%%%%%% Read in Files

for n=1:nDumps   
    u = rdmdsWrapper(fullfile(exppath,'results/UVEL_inst'),dumpIters(n));
    ADVr_TH  = rdmdsWrapper(fullfile(exppath,'results/ADVr_TH'),dumpIters(n));              
    ADVx_TH  = rdmdsWrapper(fullfile(exppath,'/results/ADVx_TH'),dumpIters(n));   
    ADVy_TH  = rdmdsWrapper(fullfile(exppath,'/results/ADVy_TH'),dumpIters(n));

    DFrE_TH  = rdmdsWrapper(fullfile(exppath,'/results/DFrE_TH'),dumpIters(n));              
    DFxE_TH  = rdmdsWrapper(fullfile(exppath,'/results/DFxE_TH'),dumpIters(n));   
    DFyE_TH  = rdmdsWrapper(fullfile(exppath,'/results/DFyE_TH'),dumpIters(n));

    KPPg_TH  = rdmdsWrapper(fullfile(exppath,'/results/KPPg_TH'),dumpIters(n));              
    oceQsw   = rdmdsWrapper(fullfile(exppath,'/results/oceQsw'),dumpIters(n));
    WTHMASS = rdmdsWrapper(fullfile(exppath,'/results/WTHMASS'),dumpIters(n));
    TFLUX = rdmdsWrapper(fullfile(exppath,'/results/TFLUX'),dumpIters(n));
    
  if (isempty(ADVr_TH) || isempty(ADVx_TH) || isempty(ADVy_TH) ||...
          isempty(DFrE_TH) || isempty(DFxE_TH) || isempty(DFyE_TH) ||...
          isempty(KPPg_TH) || isempty(oceQsw) || isempty(WTHMASS) ||...
          isempty(TFLUX));
    break;
  end






%%%%% Doing first five grid levels
    Adv_Tend = NaN(Nx,Ny,5);
    for i = 1:Nx
        for j = 1:Ny
            for k = 1:5
            Adv_Tend(i,j,k) = -(sum(sum(ADVr_TH(i,j,k)-sum(ADVr_TH(i,j,k+1))))/CellVol(i,j,k))+sum(sum(ADVx_TH(i+1,j,k)-sum(ADVx_TH(i,j,k)))/CellVol(i,j,k))+sum(sum(ADVy_TH(i,j+1,k)-sum(ADVy_TH(i,j,k)))/CellVol(i,j,k));
            end
        end
    end


%%%%%% Diffusion

    Dif_Tend = NaN(Nx,Ny,5);
    for i = 1:Nx
        for j = 1:Ny
            for k = 1:5
            Dif_Tend(i,j,k) = -(sum(sum(DFrE_TH(i,j,k)-sum(DFrE_TH(i,j,k+1))))/CellVol(i,j,k))+sum(sum(DFxE_TH(i+1,j,k)-sum(DFxE_TH(i,j,k)))/CellVol(i,j,k))+sum(sum(DFyE_TH(i,j+1,k)-sum(DFyE_TH(i,j,k)))/CellVol(i,j,k));
            end
        end
    end

%%%%%% KPP tendency
    KPP_tend = NaN(Nx,Ny,5);
    for i = 1:Nx
        for j = 1:Ny
             for k = 1:5 
             KPP_tend = -(sum(sum(KPPg_TH(i,j,k)-sum(KPPg_TH(i,j,k+1))))/CellVol(i,j,k));
             end
        end
    end

%%%%%%%% Depth 
    for k = 1:Nr
        depth = RF(5);
        depth1 = RF(1);
        swfrac = 0.62 * exp(depth/0.6) + (1.0 - 0.62) * exp(depth/20.0);
        swfrac1 = 0.62 * exp(depth1/0.6) + (1.0 - 0.62) * exp(depth1/20.0); 
    end



%%%%%%%% SW Heat Tendency
    Qsw_tend = NaN(Nx,Ny,5);
    for i = 1:Nx
        for j = 1:Ny
            for k = 1:5
             Qsw_tend(i,j,k) = oceQsw(i,j)/(rhoConst*Cp)/(DRF(k)* hFacC(i,j,k) )*(swfrac - swfrac1); 
            end
        end
    end



    tsurfcor = NaN(Nx,Ny);
    for i = Nx
        for j =Ny
            tsurfcor = sum(WTHMASS(i,j) * RAC(i,j)) / globalArea;
        end
    end

    Surf_corr_tend = NaN(Nx,Ny,1);
    for i = 1:Nx
        for j = 1:Ny
            Surf_corr_tend = (tsurfcor- WTHMASS(i,j)) / (DRF(1) * hFacC(ix,iy,1)); 
        end
    end

    Tflx_tend=NaN(Nx,Ny,Nr);    
    for i = 1:Nx
        for j = 1:Ny
            for k = 1:Nr
                Tflx_tend = (TFLUX(i,j) - oceQsw(i,j)) / (rhoConst * Cp * DRF(1) * hfacC(i,j,k));
            end
        end
    end

%%%%%%%%%%%%%%%%% Total Heat Budget

    surf_heat_tend = Adv_Tend(i,j,k) + Dif_Tend(i,j,k) + KPP_tend(i,j,k) + Qsw_tend(i,j,k) + Tflx_tend(i,j,k) + Surf_corr_tend(i,j,k);




end



