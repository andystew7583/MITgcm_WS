%%%%%%%% Find melt-weighted depth of ice shelf base && freezing temperature where most meltwater
%%%%%%%% forms
setExpname
run ../newexp/defineGrid.m
loadexp



tmin = 9*86400*360;
tmax = 18*86400*360;
Control_shimelt = readIters(exppath,'SHIfwFlx',dumpIters,deltaT,tmin,tmax,Nx,Ny,1);

melt=0;
botdepth = 0;
thetafrz=0;
for i=1:Nx
    for j =1:Ny
                 if XC(i,j)<-20 && XC(i,j)>-80
                     if YC(i,j)<-75
                         if SHELFICEtopo(i,j)<0  
                            if Control_shimelt(i,j)<0 
                                 Pressure = (-rho0*g*SHELFICEtopo(i,j)/Pa1dbar);
                                 thetafrz = thetafrz+(.0901 - 1.96 -(7.61e-4*Pressure))*RAC(i,j)*Control_shimelt(i,j);
                                 melt = melt+Control_shimelt(i,j)*RAC(i,j);
                                 botdepth = botdepth+SHELFICEtopo(i,j)*RAC(i,j)*Control_shimelt(i,j);

                                
                            end
                         end
                
                     end
                 end
            
            
        
    end
end





botdepth = botdepth/melt;
thetafrz=thetafrz/melt;

