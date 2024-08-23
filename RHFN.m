% Right hand functions of Non linear ODEs y'=F(t,y)
                    
function F = RHFN(t,y,A,j)  
% t - scalar, A - empty
if     j==1, F=[y(2);  y(1)*( y(2)/t-2*y(1)/t^2+2)/t^2 ]; % N 66 Romanko                     
elseif j==2, F=[2*t*y(1)*y(4);10*t*y(1)^5*y(4);2*t*y(4);-2*t*(y(3)-1)]; %Hai146 
elseif j==3                       % F with fantasy choice
   t1 = t+1;    
   t2 = 1/t1;
   st = sin(t);
   sqt = sqrt(t1);
   F=[y(4)*log(abs(y(1))) - st/sqt  + cos(t)*exp(st); ...
      sin(y(3))      - sin(atan(t)) + t2;             ...
      exp(y(2))      - t1           + 1/(t*t+1);      ...
      y(4)^2         - t2           - 0.5/(t1*sqt)   ]; 
elseif j==4              % this block = F_DNK, t & j are not used here 
   s = sin(y(1)+y(2));
   F = [ y(3); y(4); -A(1)*sin(y(1))+A(3)*s; -A(2)*sin(y(2))+A(4)*s ];
end

% function F = F_DNK(t,y,a)
% s = sin(y(1)+y(2));
% F = [ y(3); y(4); -a(1)*sin(y(1))+a(3)*s; -a(2)*sin(y(2))+a(4)*s ];

