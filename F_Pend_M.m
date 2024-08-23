% Transformation of F for scalar t to F_M for vector t 
% Vector-function F at vector t is matrix-function F_M
% Algorithm from F to F_M:
%  1) y(i) ->  y(:,i)    2) ";" -> ","  in formula of F
% Note: for nonautonom ODE t is n0-dimensional column,
%       => operatios '*' -> '.*', etc
                    
function F = F_Pend_M(t,y,m,L)  
c = -1/(m*L);
g = 9.81 + 0.05*sin(2*pi.*t);
F = [y(:,3); y(:,4); c*y(:,1)*y(:,5); c*y(:,2)*y(:,5)-g; y(:,1)^2+y(:,2)^2-L^2];            
