%%%%Calc WM properties just outside at ice shelf face

setExpname
loadexp

yvalues = NaN(Nx,Ny,Nr);
for i = 1:Nx
    for j = 1:Ny
        if SHELFICEtopo(i,j)==0 &&  bathy(i,j)<-200
            disp(j)
            break
        end
    end
end
