                   % Exact solutions for non linear ODEs
% t is column of M = (1 or 16)(number of initial points)
%   M=1:  solution at next 15 points will calc by dop853, odex... in ds10
%   M=16: solition at 16 points had calculated before ds10
% es is MxN-matrix where N is number of equation = dim of y

function es = esN(t,A,j,k)   
t = t(:);
if  j == 1                              % s_n66(k,t)
   t2 = t.^2;  t3 = 5*t2.*t+1;
   es = 6*t2./t3;                       % 1st coordinate of exact solution
   if k ~= 1, es = [es, -6*t.*(t3-3)./t3.^2 ]; end  % exact solution
  
elseif j == 2                           % H146(k,t)
   t2 = t.^2; st = sin(t2); 
   es = exp(st);
   if k ~= 1, es = [es, exp(5*st), st+1, cos(t2) ]; end
   
elseif j == 3                           % fant(k,t)
   es =  exp(sin(t));
   if k ~= 1, es = [es, log(t+1), atan(t), 1./sqrt(t+1)]; end

elseif j == 4                           % DNK
   errordlg('Exact DNK-solution unknown => y0 must be given at a single point') 
end
