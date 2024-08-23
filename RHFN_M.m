% Transformation of RHFN for scalar t to RHFN_M for vector t 
% Vector-function F at vector t is matrix-function F_M
% Transformation algorithm from F to F_M:
%  1) y(i) ->  y(:,i)    2) ";" -> ","  in formula of F
% Note: for nonautonom ODE t is n0-dimensional column,
%       => operatios '*' -> '.*', etc

function F = RHFN_M(t,y,A,j)   
% t - n-vector ( n = 16 )
% y - n*M-matrix
% F - n*M-matrix  [y(2);  y(1)*( y(2)/t-2*y(1)/t^2+2)/t^2 ];

t = t(:);
if     j==1, F = [ y(:,2),  y(:,1).*(y(:,2)./t-2*y(:,1)./t.^2+2)./t.^2 ]; % n66
elseif j==2,                                                 % H146
   F = [ 2*t.*y(:,1).*y(:,4),  10*t.*y(:,1).^5.*y(:,4), ...
         2*t.*y(:,4),          -2*t.*(y(:,3)-1) ];
elseif j==3                                                  % fant
   t1 = t+1;                                                 
   t2 = 1./t1;
   st = sin(t);
   sqt = sqrt(t1);
   F = [ y(:,4).*log(y(:,1)) - st./sqt      + cos(t).*exp(st), ...
         sin(y(:,3))         - sin(atan(t)) + t2,              ...
         exp(y(:,2))         - t1           + 1./(t.*t+1),     ...
         y(:,4).^2           - t2           - 0.5./(t1.*sqt)  ]; 
elseif j==4                                                   % DNK
   s = sin(y(:,1)+y(:,2));
   F = [ y(:,3), y(:,4), -A(1)*sin(y(:,1))+A(3)*s, -A(2)*sin(y(:,2))+A(4)*s ];
end
