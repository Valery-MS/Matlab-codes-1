%                  Linear inhomogeneous  ODEs
% Calc of solution of non homogeneous Linear ODE y'=GM*y by DScheme of order10
%  - Highest Order Stable Difference Scheme(HOS_DS) 10 = 2m+4, m = 3
%    ( Mathcad-files DS10, DS12)
% alfa = ( 1 La mu 3nu 0 3nu mu La  1) 
% beta = (b4 b3 b2  b1 b0 b1 b2 b3 b4);  b = [b1,..., b5], b(5) = b0
% GM  - function-> 1*n-cell array of generating M*M-matrix G(ti) values,i=1:n,
%        adapted for vectorized computations
% t    = 1*n-double array
% y0   - M*8-double array of solution initial values
% h    - t-step size
% Lmn  = [ Lambda, Mu, Nu ]
% b    =  beta = [ b1, b2, b3, b4, b5 ],  b5 = beta0
 
function y = ds10Li( GM, t, y0, op, varargin) 

h = op{1}; ef = op{3}; Lmn = op{4}; b = op{5}; 
n = numel(t);              % n - number of t-points
[IL, M] = size(y0);        % M - number of coordinates of y = dim of G

if ef,   [G, g] = GM(t,varargin{:});
else
   G = cell(1,n);   g = zeros(M,n);
   for i = 1:n,     [ G{i}, g(:,i) ] = GM(t(i),varargin{:});  end,end

if IL == 1
   f   = op{2};  
   zam = op{6};
   setms = { @dopset @odeset};  
   mets  = { @dop853 @ode113_c};         
   warning('off','all')
   opt = setms{zam}('InitialStep',h,'RelTol',eps,'AbsTol',eps);   
   [t_, y0, nf] = mets{zam}( f, t(1:8), y0, opt, varargin{:} ); % m=3=> 2*m+2=8 
   warning('on','all');  end

y = nan(M,n);
y(:,1:8) = y0';              

hb = h*b;   La = Lmn(1);  mu = Lmn(2);  nu3 = 3*Lmn(3);

for i = 5:n-4
   S4 = -hb(4)*G{i+4};
   for j = 1:M,  S4(j,j) = 1+S4(j,j);  end
   
   r = hb(4)*g(:,i+4)+La*(y(:,i-1)-y(:,i+1)) + mu*(y(:,i-2)-y(:,i+2)) ...
     + nu3*(y(:,i-3)-y(:,i+3)) + y(:,i-4) ...
     + hb(1)*(G{i-1}*y(:,i-1)+g(:,i-1)  + G{i+1}*y(:,i+1)+g(:,i+1)) ...
     + hb(2)*(G{i-2}*y(:,i-2)+g(:,i-2)  + G{i+2}*y(:,i+2)+g(:,i+2)) ...
     + hb(3)*(G{i-3}*y(:,i-3)+g(:,i-3)  + G{i+3}*y(:,i+3)+g(:,i+3)) ...
     + hb(4)*(G{i-4}*y(:,i-4)+g(:,i-4)) + hb(5)*(G{i}*y(:,i)+g(:,i));
 
   y(:,i+4) = linsolve(S4, r); end

y = y';
