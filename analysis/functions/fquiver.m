%%%
%%% fquiver.m
%%%
%%% Quiver plot with filled arrow heads.
%%%
function fquiver (x,y,u,v,headWidth,headLength,LineLength)

  %%% Make a quiver plot to generate quiver data
  hq = quiver(x,y,u,v);        
  set(hq,'Visible','off');

  %%% Get the data from regular quiver
  U = hq.UData;
  V = hq.VData;
  X = hq.XData;
  Y = hq.YData;
  
  %%% Normalize vector data
  Uabs = sqrt(U.^2+V.^2);
  Umax = nanmax(nanmax(Uabs));
  U = U./Umax;
  V = V./Umax;
 
  %%% Use annotations to re-plot quiver data
  for ii = 1:size(X,1)
    for ij = 1:size(X,2)        
      if (isnan(U(ii,ij)+V(ii,ij)))
        continue;
      end
      
      ah = annotation('arrow'); ... %[X(ii,ij) X(ii,ij)+LineLength*U(ii,ij)],[Y(ii,ij) Y(ii,ij)+LineLength*V(ii,ij)], ...          
      set(ah,'parent',gca);
      set(ah,'Position',[X(ii,ij) Y(ii,ij) LineLength*U(ii,ij) LineLength*V(ii,ij)]);
      set(ah,'headStyle','cback1','HeadLength',headLength*(Uabs(ii,ij)/Umax).^(1),'HeadWidth',headWidth*(Uabs(ii,ij)/Umax).^(1),'Color','k');
%       ah = annotation('line');
%       set(ah,'parent',gca);
%       set(ah,'Position',[X(ii,ij) Y(ii,ij) LineLength*U(ii,ij) LineLength*V(ii,ij)]);
    end
  end
  
end

