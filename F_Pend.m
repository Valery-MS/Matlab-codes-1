% Right hand functions of Non linear ODEs y'=F(t,y)
                    
function F = F_Pend(t,y,m,L)  
c = -1/(m*L);
g = 9.81 + 0.05*sin(2*pi*t);
F = [ y(3); y(4); c*y(1)*y(5); c*y(2)*y(5)-g;  y(1)^2+y(2)^2-L^2 ];            
