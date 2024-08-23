% Numerical solution of a system of first order ODE  y'=F(x,y) on MatLab
% Hairer "Solving ODE I" (in Russian p.459)
                                          
function [t, z, varargout]  = odexO( FCN, XSF, Y, opt, varargin )

global NFCN                 % num of FCN evaluations
global DZ T NJ HH W ERR FAC A EPSD4 UROUND FAC1 FAC2 SAFE2

H    = opt{1};      % Initial stepsize
EPS  = opt{2};      % local tolerance

N = numel(Y);       % dimension of system ODE ( N <= 51 )
Farg = varargin;

DZ = zeros(N,1);  T = zeros(9,N); HH = zeros(9,1);  W = HH; 
NJ = (2:2:18)';      
A  = [ 3 7 13 21 31 43 57 73 91 ];

NMAX = 1000;        % max number of steps
z = nan(N,NMAX);    % row i = i-coordinate in NACCPT+1 points
z(:,1) = Y;        % insert y0 into matrix of numerical solution z
KM = 9;
UROUND = 1.111e-16; % smallest number satisfying 1 + UROUND > 1
                    % to be adapted by the user
FAC1 = 2E-2; FAC2 = 4; FAC3 = 0.9; FAC4 = 0.8; SAFE1 = 0.65; SAFE2 = 0.94;
                    %Initial preparations
EPSD4 = EPS*SAFE1; 
NFCN = 0; 
                    % num of computed, accepted & rejected steps
NSTEP = 0; NACCPT = 0; NREJCT = 0; 
K = max(3, min( 8, fix(-log10(EPS)*0.6+1.5)));  % ceil analog INT in fortran
X = XSF(1);  XEND = XSF(2);
HMAX = XEND-X;      % max stepsixe
t = X;
H = min([H, HMAX, (XEND-X)/2 ]);
ERR = 0; W(1) = 0;  REJECT = false;  LAST = false;
FAIL = 0;
goto = 10;
                   % Is XEND reachedin the next step ?
while true         %  goto = 10 30
   if goto == 10  
      H1 = XEND-X;
      if H1 <= UROUND, break; end
      H = min([ H, H1, HMAX ]);
      if H >= H1-UROUND, LAST = true; end
      DZ(1:N) = FCN( X, Y, Farg{:} );
      NFCN = NFCN+1;
             % 1st & last step
      if NSTEP == 0 || LAST
         NSTEP = NSTEP+1;
         for J = 1:K
            KC = J;
            MIDEX( J, X, Y, H, HMAX, N, FCN, Farg{:} )
            if J > 1 && ERR <= EPS, goto = 60; break; end,end
         if goto ~= 60, goto = 55; end,end,end 
               
   if goto <= 30     % Basic integration step   goto = 10 30 55 60          
      NSTEP = NSTEP+1;
      if NSTEP >= NMAX, FAIL = 1; break;end
      KC = K-1;
      for J = 1:KC, MIDEX( J, X, Y, H, HMAX, N, FCN, Farg{:} ); end
               % Convergence monitor
      if K == 2 || REJECT, goto = 50;
      elseif ERR <= EPS,   goto = 60;
      elseif ERR/EPS > (NJ(K+1)*NJ(K)/4)^2, goto = 100; end,end
                     %  goto = 10 30 50 55 60 100
   if goto <= 50               
      MIDEX( K, X, Y, H, HMAX, N, FCN, Farg{:} );
      KC = K;
      if ERR <= EPS, goto = 60; end,end
            % Hope for convergence in line K+1, goto = 10 30 50 55 60 100
   if goto <= 55         
      if ERR/EPS > (NJ(K+1)/2)^2, goto = 100; end,end 
  
   if goto <= 55  
      KC =K+1;
      MIDEX( KC, X, Y, H, HMAX, N, FCN, Farg{:} )
      if ERR > EPS, goto = 100; end,end 
               % Step is accepted.  Here goto = 10 30 50 55 60 100
   if goto <= 60
      X = X+H; 
      for I = 1:N, Y(I) = T(1,I); end
      NACCPT = NACCPT+1;
      t = [t; X];
      z(:,NACCPT+1) = Y;
            % Compute optimal order
      if KC == 2              
         if REJECT, KOPT = 2; 
         else       KOPT = 3; end
         goto = 80; end
                         % goto = 10 30 50 55 60 80
      if goto <= 60
         if KC <= K
            KOPT = KC;
            if W(KC-1) < W(KC)*FAC3, KOPT = KC-1; end
            if W(KC) < W(KC-1)*FAC3, KOPT = min(KC+1,KM-1); end
         else
            KOPT = KC-1;
            if KC > 3 && W(KC-2) < W(KC-1)*FAC3, KOPT = KC-2; end
            if W(KC) < W(KOPT)*FAC3, KOPT = min(KC,KM-1); end,end,end
                    % After a rejected step               
      if REJECT      % Label 80
         K = min(KOPT,KC);
         H = min(H,HH(K));
         REJECT = false;
            % Compute stepsize for next step goto = 10 30 50 55 60 80
      else
         if     KOPT <= KC,                 H=HH(KOPT);
         elseif KC<K && W(KC)<W(KC-1)*FAC4, H=HH(KC)*A(KOPT+1)/A(KC);
         else                               H=HH(KC)*A(KOPT)/A(KC); end
         K = KOPT; end
      goto = 10; 
                         % Step is rejected
   else                  % goto == 100                 
      K = min(K,KC);
      if K > 2 && W(K-1) < W(K)*FAC3, K= K-1; end
      NREJCT = NREJCT+1;
      H = HH(K);
      REJECT = true;
      goto = 30; end,end   

z = z(:,1:NACCPT+1)'; 
if nargout == 3,  varargout{:} = [ NSTEP NACCPT NFCN NREJCT FAIL ]; end
  

%****************************************************************************
%*****************************************************************************

function MIDEX( J, X, Y, H, HMAX, N, FCN, varargin)

global NFCN                    % num of FCN evaluations
global DZ T NJ HH W ERR FAC A EPSD4 UROUND FAC1 FAC2 SAFE2

Farg = varargin;
DY = nan(N,1);  YH1 = DY; YH2 = DY; 
HJ = H/NJ(J);

               % EULER starting step
for I = 1:N
   YH1(I) = Y(I);
   YH2(I) = Y(I)+HJ*DZ(I); end

              % Explicit midpoint rule
M = NJ(J)-1;
for MM = 1:M
    DY = FCN( X+HJ*MM, YH2, Farg{:} );
    for I = 1:N
        YS = YH1(I);
        YH1(I) = YH2(I);
        YH2(I) = YS+2*HJ*DY(I); end,end

               % Final smoothing step
DY = FCN( X+H, YH2, Farg{:} );             
for I = 1:N,  T(J,I) = (YH1(I) + YH2(I) + HJ*DY(I))/2; end
NFCN = NFCN+NJ(J);

               % Polinomial extrapolation
if J == 1, return,end
for L = J:-1:2
    FAC = (NJ(J)/NJ(L-1))^2 - 1;
    for I = 1:N
        T(L-1,I) = T(L,I) + (T(L,I)-T(L-1,I))/FAC; end,end

                % Scaling
ERR = 0;
for I = 1:N
    SCAL = max( [ abs(Y(I)), abs(T(1,I)), 1e-6, UROUND/EPSD4 ] );
    ERR = ERR + ( (T(1,I)-T(2,I))/SCAL )^2; end
ERR = sqrt(ERR/N);

                % Compute optimal stepsize
EXP0 = 1/(2*J-1);
FACMIN = FAC1^EXP0; 
FAC = min( FAC2/FACMIN, max( FACMIN, (ERR/EPSD4)^EXP0/SAFE2 ) );
FAC = 1/FAC;
HH(J) = min(H*FAC, HMAX );
W(J) = A(J)/HH(J);
