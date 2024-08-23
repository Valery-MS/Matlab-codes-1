%                Linear homogeneous ODEs, variable step size
% Calc of solution of Linear homogeneous ODE y'=G*y by Diff Scheme of 10th order - 
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
 
function y = ds10Lhv( GM, t, y0, op, varargin) 

h = op{1}; ef = op{3}; Lmn = op{4}; b = op{5}; AT = op{7}(1); RT = op{7}(2);
% Lmn: a1:a3 (ds10), a1,a2 (ds8), a1 (ds6)
% b:   b1:b5 (ds10), b1:b4 (ds8), b1:b3 (ds6) - beta_1:beta_n,beta_0
n = numel(t);             % n - number of t-points
[IL, M] = size(y0);       % M - number of coordinates of y = dim of G

if ef,   [G, g] = GM(t,varargin{:});  % g don't used
else
   G = cell(1,n);  
   for i = 1:n,     [ G{i}, g] = GM(t(i),varargin{:});  end,end
   
if IL == 1
   f   = op{2};  
   zam = op{6};
   setms = { @dopset @odeset};  
   mets  = { @dop853 @ode113_c};         
   warning('off','all')
   opt = setms{zam}('InitialStep',h,'RelTol',eps,'AbsTol',eps);    
   [t_, y0, nf] = mets{zam}( f, t(1:16), y0, opt, varargin{:} ); % m=3=> 2*m+2=8 
   warning('on','all');  end
        
y = nan(M,n);
y(:,1:16) = y0';              

L0 = Lmn(1);  m0 = Lmn(2);   n0 = 3*Lmn(3);   % La, mu, nu for ds10
L8 = Lmn(4);  m8 = 2*Lmn(5);                  % La, mu     for ds8
L6 = Lmn(6);                                  % La         for ds6

hb = h*b;      hb2 = 2*hb;
H10 = hb2(1);  H20 = hb2(2);  H30 = hb2(3);  H40 = hb2(4); H50 = hb2(5);
h10 = hb (1);  h20 = hb (2);  h30 = hb (3);  h40 = hb (4); h50 = hb (5);
H18 = hb2(6);  H28 = hb2(7);  H38 = hb2(8);  H48 = hb2(9);
h18 = hb (6);  h28 = hb (7);  h38 = hb (8);  h48 = hb (9);
H16 = hb2(10); H26 = hb2(11); H36 = hb2(12);
h16 = hb (10); h26 = hb (11)  h36 = hb (12);

p = 10; o2p = 1/2^(p-1);

for i = 5:n-4

   s8 = -h4*G{i+8};  S8 = 2*s8;
   for j = 1:M,  s8(j,j) = 1+s8(j,j); S8(j,j) = 1+S8(j,j);end
   
   r =  L0*(y(:,i+3)-y(:,i+5)) + m0*(y(:,i+2)-y(:,i+6)) ...
     +  n0*(y(:,i+1)-y(:,i+7)) + y(:,i) ...
     +    h10*(G{i+3}*y(:,i+3) + G{i+5}*y(:,i+5)) ...
     +    h20*(G{i+2}*y(:,i+2) + G{i+6}*y(:,i+6)) ...
     +    h30*(G{i+1}*y(:,i+1) + G{i+7}*y(:,i+7)) ...
     +    h40* G{i}  *y(:,i)   + G{i+4}*y(:,i+4)*h50;

   R =  L0*(y(:,i-2)-y(:,i+2)) + m0*(y(:,i-4)-y(:,i+4)) ...
     +  n0*(y(:,i-6)-y(:,i+6)) + y(:,i-8) ...
     +    H10*(G{i-2}*y(:,i-2) + G{i+2}*y(:,i+2)) ...
     +    H20*(G{i-4}*y(:,i-4) + G{i+4}*y(:,i+4)) ...
     +    H30*(G{i-6}*y(:,i-6) + G{i+6}*y(:,i+6)) ...
     +    H40* G{i-8}*y(:,i-8) + G{i}*y(:,i)*H50;
   
   y8 = linsolve(s8, r); 
   Y8 = linsolve(S8, R);   end

aby = abs(y8-Y8);
erA = o2p*;

y = y';
