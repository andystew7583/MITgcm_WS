%%%% Selecting region over which to calculate HSSW fluxes aka Ronne polynya region

setExpname
loadexp




start = NaN(Nx,Ny);
for i  = 1:Nx
    for j = 1:Ny
            if (XC(i,j)<-48) && XC(i,j)>-70&& SHELFICEtopo(i,j)==0 &&(SHELFICEtopo(i,j)-bathy(i,j)) > 0 && YC(i,j) <-71.5 
                    start(i,j) = 1;
                
            end
        

        
    end    
    
end

for i = 1:Nx
    for j=1:Ny-30
       if (XC(i,j)<-48) && XC(i,j)>-70 && SHELFICEtopo(i,j)==0 &&(SHELFICEtopo(i,j)-bathy(i,j)) > 0 && YC(i,j) <-71
        if start(i,j+30)==1
            start(i,j+30)=NaN;
        end
       
       end
      

       
       
    end
end


% %%%% finding volume of region

Vol = 0;

for i = 1:Nx
    for j = 1:Ny
            if start(i,j)==1
                Vol = Vol + RAC(i,j);
                
            end
        
    end
end

H = 300;




