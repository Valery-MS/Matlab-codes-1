function [t, z, nf] = ODEX( FCN, X, Y, opt, varargin)

% odex - другой код (алгоритм тот же)
% ODEX - переписанный на МАТЛАБ код из книги Хайрера
% ODEX - Solve non-stiff differential equations. This is an extrapolation 
% algorithm (GBS), based on the explicit midpoint rule with step size 
% control, order selection and dense output. 
%
%     INPUT PARAMETERS  
%
% FCN         NAME (EXTERNAL) OF SUBROUTINE COMPUTING THE
%             VALUE OF F(X,Y):
%                  SUBROUTINE FCN(N,X,Y,F,RPAR,IPAR)
%                DOUBLE PRECISION X,Y(N),F(N)
%                F(1)=...   ETC.
%
% XSF         [ Start_X, Final_ X ] 2-array of INITIAL & FINAL X-VALUEs
%
% Y(1:N)      INITIAL VALUES FOR Y;  N - DIMENSION OF THE SYSTEM 
% 
% opt         cell-array, containing following variables;
%
% H           INITIAL STEP SIZE GUESS;
%             H=1.D0/(NORM OF F'), USUALLY 1.D-1 OR 1.D-3, IS GOOD.
%             THIS CHOICE IS NOT VERY IMPORTANT, THE CODE QUICKLY
%             ADAPTS ITS STEP SIZE. WHEN YOU ARE NOT SURE, THEN
%             STUDY THE CHOSEN VALUES FOR A FEW
%             STEPS IN SUBROUTINE "SOLOUT".
%             (IF H=0.D0, THE CODE PUTS H=1.D-4).
%
% RTOL,ATOL   RELATIVE AND ABSOLUTE ERROR TOLERANCES. THEY
%             CAN BE BOTH SCALARS OR ELSE BOTH VECTORS OF LENGTH N.
%
% ITOL        SWITCH FOR RTOL AND ATOL:
%               ITOL=0: BOTH RTOL AND ATOL ARE SCALARS.
%                 THE CODE KEEPS, ROUGHLY, THE LOCAL ERROR OF
%                 Y(I) BELOW RTOL*ABS(Y(I))+ATOL
%               ITOL=1: BOTH RTOL AND ATOL ARE VECTORS.
%                 THE CODE KEEPS THE LOCAL ERROR OF Y(I) BELOW
%                 RTOL(I)*ABS(Y(I))+ATOL(I).
%
% SOLOUT      NAME(EXTERNAL was changed to internal) OF SUBROUTINE PROVIDING 
%             THE NUMERICAL SOLUTION DURING INTEGRATION ( see lower ... (
%     
% IOUT        SWITCH FOR CALLING THE SUBROUTINE SOLOUT:
%                IOUT=0: SUBROUTINE IS NEVER CALLED
%                IOUT=1: SUBROUTINE IS USED FOR OUTPUT
%                IOUT=2: DENSE OUTPUT IS PERFORMED IN SOLOUT
%
% WORK        ARRAY OF WORKING SPACE OF LENGTH "LWORK".
%             SERVES AS WORKING SPACE FOR ALL VECTORS.
%             "LWORK" MUST BE AT LEAST
%                N*(KM+5)+5*KM+20+(2*KM*(KM+2)+5)*NRDENS
%             WHERE NRDENS=IWORK(8) (SEE BELOW) AND
%                    KM=9                IF IWORK(2)=0
%                    KM=IWORK(2)         IF IWORK(2).GT.0
%             WORK(1),...,WORK(20) SERVE AS PARAMETERS
%             FOR THE CODE. FOR STANDARD USE, SET THESE
%             PARAMETERS TO ZERO BEFORE CALLING.
%
% LWORK       DECLARED LENGTH OF ARRAY "WORK".
%
% IWORK       INTEGER WORKING SPACE OF LENGTH "LIWORK".
%             "LIWORK" MUST BE AT LEAST
%                           2*KM+21+NRDENS 
%             IWORK(1),...,IWORK(20) SERVE AS PARAMETERS
%             FOR THE CODE. FOR STANDARD USE, SET THESE
%             PARAMETERS TO ZERO BEFORE CALLING.
%
% LIWORK      DECLARED LENGTH OF ARRAY "IWORK".
%
% RPAR, IPAR  REAL AND INTEGER PARAMETERS (OR PARAMETER ARRAYS) WHICH  
%             CAN BE USED FOR COMMUNICATION BETWEEN YOUR CALLING
%             PROGRAM AND THE FCN, JAC, MAS, SOLOUT SUBROUTINES. 
%
% varargin    arguments for FCN ( if needed )
%
%-----------------------------------------------------------------------
% 
% SOPHISTICATED SETTING OF PARAMETERS
% -----------------------------------
%          SEVERAL PARAMETERS (WORK(1),...,IWORK(1),...) ALLOW
%          TO ADAPT THE CODE TO THE PROBLEM AND TO THE NEEDS OF
%          THE USER. FOR ZERO INPUT, THE CODE CHOOSES DEFAULT VALUES.
%
%    WORK(1)   UROUND, THE ROUNDING UNIT, DEFAULT 2.3D-16.
%
%    WORK(2)   MAXIMAL STEP SIZE, DEFAULT XEND-X.
%
%    WORK(3)   STEP SIZE IS REDUCED BY FACTOR WORK(3), IF THE
%          STABILITY CHECK IS NEGATIVE, DEFAULT 0.5.
%
%    WORK(4), WORK(5)   PARAMETERS FOR STEP SIZE SELECTION
%          THE NEW STEP SIZE FOR THE J-TH DIAGONAL ENTRY IS
%          CHOSEN SUBJECT TO THE RESTRICTION
%             FACMIN/WORK(5) <= HNEW(J)/HOLD <= 1/FACMIN
%          WHERE FACMIN=WORK(4)**(1/(2*J-1)) 
%          DEFAULT VALUES: WORK(4)=0.02D0, WORK(5)=4.D0
%
%    WORK(6), WORK(7)   PARAMETERS FOR THE ORDER SELECTION
%          STEP SIZE IS DECREASED IF    W(K-1) <= W(K)*WORK(6)
%          STEP SIZE IS INCREASED IF    W(K) <= W(K-1)*WORK(7)
%          DEFAULT VALUES: WORK(6)=0.8D0, WORK(7)=0.9D0
%
%    WORK(8), WORK(9)   SAFETY FACTORS FOR STEP CONTROL ALGORITHM
%         HNEW=H*WORK(9)*(WORK(8)*TOL/ERR)**(1/(J-1))
%         DEFAULT VALUES: WORK(8)=0.65D0,
%                    WORK(9)=0.94D0  IF "HOPE FOR CONVERGENCE"
%                    WORK(9)=0.90D0  IF "NO HOPE FOR CONVERGENCE"
%
%    IWORK(1)  THIS IS THE MAXIMAL NUMBER OF ALLOWED STEPS.
%          THE DEFAULT VALUE (FOR IWORK(1)=0) IS 10000.
%
%    IWORK(2)  THE MAXIMUM NUMBER OF COLUMNS IN THE EXTRAPOLATION 
%          TABLE. THE DEFAULT VALUE (FOR IWORK(2)=0) IS 9.
%          IF IWORK(2).NE.0 THEN IWORK(2) SHOULD BE .GE.3.
%
%    IWORK(3)  SWITCH FOR THE STEP SIZE SEQUENCE (EVEN NUMBERS ONLY)
%          IF IWORK(3).EQ.1 THEN 2,4,6,8,10,12,14,16,...
%          IF IWORK(3).EQ.2 THEN 2,4,8,12,16,20,24,28,...
%          IF IWORK(3).EQ.3 THEN 2,4,6,8,12,16,24,32,...
%          IF IWORK(3).EQ.4 THEN 2,6,10,14,18,22,26,30,...
%          IF IWORK(3).EQ.5 THEN 4,8,12,16,20,24,28,32,...
%          THE DEFAULT VALUE IS IWORK(3)=1 IF IOUT.LE.1;
%          THE DEFAULT VALUE IS IWORK(3)=4 IF IOUT.GE.2.
% 
%    IWORK(4)  STABILITY CHECK IS ACTIVATED AT MOST IWORK(4) TIMES IN
%          ONE LINE OF THE EXTRAP. TABLE, DEFAULT IWORK(4)=1. 
%
%    IWORK(5)  STABILITY CHECK IS ACTIVATED ONLY IN THE LINES
%          1 TO IWORK(5) OF THE EXTRAP. TABLE, DEFAULT IWORK(5)=1. 
%
%    IWORK(6)  IF  IWORK(6)=0  ERROR ESTIMATOR IN THE DENSE
%          OUTPUT FORMULA IS ACTIVATED. IT CAN BE SUPPRESSED
%          BY PUTTING IWORK(6)=1.
%          DEFAULT IWORK(6)=0  (IF IOUT.GE.2).
%
%    IWORK(7)  DETERMINES THE DEGREE OF INTERPOLATION FORMULA 
%          MU = 2 * KAPPA - IWORK(7) + 1
%          IWORK(7) SHOULD LIE BETWEEN 1 AND 6
%          DEFAULT IWORK(7)=4  (IF IWORK(7)=0).
%
%    IWORK(8)  = NRDENS = NUMBER OF COMPONENTS, FOR WHICH DENSE OUTPUT
%          IS REQUIRED
%
%    IWORK(21),...,IWORK(NRDENS+20) INDICATE THE COMPONENTS, FOR WHICH
%          DENSE OUTPUT IS REQUIRED
%
%----------------------------------------------------------------------
% OUTPUT PARAMETERS 
% ----------------- 
% X           X-VALUE FOR WHICH THE SOLUTION HAS BEEN COMPUTED
%             (AFTER SUCCESSFUL RETURN X=XEND).
%
% Y(N)        NUMERICAL SOLUTION AT X
% 
% H           PREDICTED STEP SIZE OF THE LAST ACCEPTED STEP
%
% IDID        REPORTS ON SUCCESSFULNESS UPON RETURN (if varargout is exist):
%               IDID=1  COMPUTATION SUCCESSFUL,
%               IDID=-1 COMPUTATION UNSUCCESSFUL.
%
%   IWORK(17)  NFCN    NUMBER OF FUNCTION EVALUATIONS
%   IWORK(18)  NSTEP   NUMBER OF COMPUTED STEPS
%   IWORK(19)  NACCPT  NUMBER OF ACCEPTED STEPS
%   IWORK(20)  NREJCT  NUMBER OF REJECTED STEPS (DUE TO ERROR TEST),
%                  (STEP REJECTIONS IN THE FIRST STEP ARE NOT COUNTED)  
%----------------------------------------------------------------------

H    = opt.InitialStep;      % Initial stepsize  
RTOL = opt.RelTol;      % local Relative tolerance
ATOL = opt.AbsTol;      % local Absolute tolerance
ITOL = opt.ITOL;      % ITOL=0: RTOL AND ATOL ARE SCALARS, otherwise vectors
%SOLOUT = opt{5};   % func of NUMERICAL SOLUTION DURING INTEGRATION
%IOUT   = opt{5};   % 0:SOLOUT isn't called, 1:usual output, 3:DENSE OUTPUT
WORK   = opt.WORK;    % WORKING array: NRDENS=IWORK(8), etc 
%LWORK  = opt{7};   % LENGTH OF WORK>=N*(KM+5)+5*KM+20+(2*KM*(KM+2)+5)*NRDENS 
IWORK  = opt.IWORK;    % IWORK(1:20)- PARAMETERS FOR THE CODE
%LIWORK = opt{9};   % LENGTH OF ARRAY "IWORK" 
%IPAR   = opt{12};  % array WHICH CAN BE USED FOR SOLOUT, FCN  ( not used )

Farg = varargin;
Y = Y(:);           % Y -> column-vector
N = numel(Y);       % dimension of system ODE ( N <= 51 )
t = sort(X);          % dense x-points  (RPAR = t) 
X = t(1);             % INITIAL X-VALUE

if numel(t) == 2, XEND = t(2);   IOUT = 0; t = X; % without dense calcs
else              XEND = t(end); IOUT = 3; end    % with    dense calcs
               
                   % SETTING THE PARAMETERS 
NFCN = 0;  NSTEP = 0;  NACCPT = 0;  NREJCT = 0;  ARRET = false;

                   % NMAX, THE MAXIMAL NUMBER OF STEPS
if IWORK(1)==0, NMAX = 10000;
else            NMAX = IWORK(1); 
   if NMAX <= 0
      Fail = sprintf('WRONG INPUT IWORK(1)=%d',IWORK(1));
      ARRET = true; end,end
                   
                   % KM     MAXIMUM NUMBER OF COLUMNS IN THE EXTRAPOLATION 
if IWORK(2) == 0, KM = 9;
else              KM = IWORK(2);
   if KM <= 2 
      Fail = sprintf(' CURIOUS INPUT IWORK(2)=%d',IWORK(2));
      ARRET = true; end,end
                 
                  % NSEQU     CHOICE OF STEP SIZE SEQUENCE
NSEQU = IWORK(3);
if IWORK(3) == 0 && IOUT <= 1, NSEQU = 1; end
if IWORK(3) == 0 && IOUT >= 2, NSEQU = 4; end
if NSEQU <= 0 || NSEQU >= 6
   Fail =sprintf(' CURIOUS INPUT IWORK(3)=%d',IWORK(3));
   ARRET = true; end

if NSEQU <= 3 && IOUT >= 2
   Fail = sprintf('IWORK(3) NOT COMPATIBLE WITH IOUT');
   ARRET = true; end
             
if IWORK(4) == 0, MSTAB = 1;            % MSTAB PARAMETER FOR STABILITY CHECK
else              MSTAB = IWORK(4); end 
                   
if IWORK(5) == 0, JSTAB = 2;            % JSTAB PARAMETER FOR STABILITY CHECK
else              JSTAB = IWORK(5); end
    
if IWORK(6) == 0      %IDERR  PARAMETER FOR ERROR ESTIMATION IN DENSE OUTPUT
   if IOUT <= 1, IDERR = 1; end
   if IOUT >= 2, IDERR = 0; end
else
   IDERR = IWORK(6);
   if IOUT <= 1
      Fail = sprintf('ERROR ESTIMATION IN DENSE OUTPUT NOT POSSIBLE %s %d',...
             'WRONG IWORK(6)=',IWORK(6));
      ARRET = true; end,end

if IWORK(7) == 0, MUDIF = 4;           % MUDIF
else              MUDIF = IWORK(7);
   if MUDIF <= 0 || MUDIF >= 7
      Fail =  sprintf('WRONG INPUT IWORK(7)=%d',IWORK(7));
      ARRET = true; end,end

%NRDENS = IWORK(8);       % NRDENS  NUMBER OF DENSE OUTPUT COMPONENTS
NRDENS = N;               % my change
if NRDENS < 0 || NRDENS > N
   Fail = sprintf('CURIOUS INPUT IWORK(8)=%d',IWORK(8));
   ARRET = true; end
          % IWORK(21),...,IWORK(NRDENS+20)
if NRDENS == N, IWORK(21:N+20) = 1:N; end
   % for I=1:NRDENS, IWORK(20+I) = I; end,end 

if WORK(1) == 0, UROUND = 2.3e-16; % SMALLEST NUMBER SATISFYING 1.D0+UROUND>10  
else             UROUND = WORK(1); 
   if UROUND <= 1e-35 || UROUND >= 1
      Fail = sprintf('WHICH MACHINE DO YOU HAVE? YOUR UROUND WAS:%g',WORK(1));
      ARRET = true; end,end

if WORK(2) == 0,  HMAX = XEND-X;               % MAXIMAL STEP SIZE    
else              HMAX = abs(WORK(2)); end

if WORK(3) == 0,  SAFE3 = 0.5;                 % STEP SIZE REDUCTION FACTOR
else              SAFE3 = WORK(3);
   if SAFE3 <= UROUND || SAFE3 >= 1
      Fail = sprintf('CURIOUS INPUT WORK(3)=%g',WORK(3));
      ARRET = true; end,end

if WORK(4) == 0, FAC1 = 0.02;   % FAC1,FAC2 PARAMETERS FOR STEP SIZE SELECTION
else             FAC1 = WORK(4); end

if WORK(5) == 0, FAC2 = 4;
else             FAC2 = WORK(5); end

if WORK(6) == 0, FAC3 = 0.8;  % FAC3, FAC4 PARAMETERS FOR THE ORDER SELECTION
else             FAC3 = WORK(6); end

if WORK(7) == 0, FAC4 = 0.9;
else             FAC4 = WORK(7); end


if WORK(8) == 0, SAFE1 = 0.65; % SAFE1 SAFETY FACTOR FOR STEP SIZE PREDICTION
else             SAFE1 = WORK(8); end

if WORK(9) == 0, SAFE2 = 0.94; % SAFE2 SAFETY FACTOR FOR STEP SIZE PREDICTION
else             SAFE2 = WORK(9); end

                 % PREPARE THE ENTRY-POINTS FOR THE ARRAYS IN WORK 
LFSAFE = 2*KM*KM+KM;
IEDY   = 21;
IEYH1  = IEDY+N;
IEYH2  = IEYH1+N;
IEDZ   = IEYH2+N;
IESCAL = IEDZ+N;
IET    = IESCAL+N;
IEFS   = IET+KM*N;
IEYS   = IEFS+LFSAFE*NRDENS;
IEHH   = IEYS+KM*NRDENS;
IEW    = IEHH+KM;
IEA    = IEW+KM;
IEFAC  = IEA+KM;
                 % TOTAL STORAGE REQUIREMENT
IECO   = IEFAC+2*KM;
ISTORE = IECO+(2*KM+5)*NRDENS-1;
if ISTORE > 1e6   % my change: LWORK -> 1e6
   Fail  = sprintf('INSUFFICIENT STORAGE FOR WORK, MIN. LWORK=%d',ISTORE);
   ARRET = true; end

                 % ENTRY POINTS FOR INTEGER WORKSPACE 
ICOM = 21;
IENJ = ICOM+NRDENS;    
      
                 % TOTAL REQUIREMENT 
IEIP   = IENJ+KM;
ISTORE = IEIP+KM+1-1;
if ISTORE > length(IWORK);
   Fail  = sprintf('INSUFF. STORAGE FOR IWORK, MIN. LIWORK=%d',ISTORE);
   ARRET = true; end

                 % WHEN A FAIL HAS OCCURED, WE RETURN WITH Fail (IDID=-1)
if ARRET, nf = Fail; return,end
  
                 % CALL TO CORE INTEGRATOR 
NRD   = max(1,NRDENS); 
NCOM  = max(1,(2*KM+5)*NRDENS);
ICOMP = IWORK(21:NRDENS+20);

[t,z,nf] = ODXCOR(N,FCN,  X,Y,   XEND,  HMAX,H, RTOL, ATOL,  ITOL,  KM,...
                  IOUT,   NMAX,  UROUND,...
                  NCOM,   ICOMP, NSEQU, MSTAB, JSTAB, LFSAFE,...
                  SAFE1,  SAFE2, SAFE3, FAC1,  FAC2,  FAC3,  FAC4,  IDERR,...
                  MUDIF,  NRD,   t,     NFCN,  NSTEP, NACCPT,NREJCT,Farg{:});



%****************************************************************************
%****************************************************************************
function [t,z,nf] = ...
   ODXCOR( N, FCN, X, Y,  XEND,  HMAX,H, RTOL,  ATOL,   ITOL,   KM,...
           IOUT,   NMAX,  UROUND, ...
           NCOM,   ICOMP, NSEQU, MSTAB,  JSTAB, LFSAFE,...
           SAFE1,  SAFE2, SAFE3, FAC1,   FAC2,  FAC3,   FAC4,   IDERR,  ...   
           MUDIF,  NRD,   t,     NFCN,   NSTEP, NACCPT, NREJCT, varargin )
 
global XOLDD HHH KMIT 

Farg = varargin;
                          % inner args, not used as input args
DY     = zeros(N,1);    YH1 = DY;  YH2 = DY;  DZ = DY;  SCAL = DY; 
FSAFE  = zeros(LFSAFE,NRD);
YSAFE  = zeros(KM,NRD);
T      = zeros(KM,N);
HH     = zeros(KM,1);   W = HH;  A = HH;  NJ = HH;
DENS   = zeros(NCOM,1);
IPOINT = zeros(KM+1,1);
ERRFAC = zeros(2*KM,1);

Fail = 0;
nt = numel(t);          %
z  = nan(nt,N);         % matrix of numerical solution  
z(1,:) = Y';

                 % DEFINE THE STEP SIZE SEQUENCE
if     NSEQU == 1,                  for I=1:KM, NJ(I) = 2*I;      end
elseif NSEQU == 2, NJ(1)=2;         for I=2:KM, NJ(I) = 4*I-4;    end
elseif NSEQU == 3, NJ(1:3)=[2 4 6]; for I=4:KM, NJ(I) = 2*NJ(I-2);end
elseif NSEQU == 4,                  for I=1:KM, NJ(I) = 4*I-2;    end
elseif NSEQU == 5,                  for I=1:KM, NJ(I) = 4*I;      end,end

                 % DEFINE THE A(I) FOR ORDER SELECTION
A(1) = 1+NJ(1);
for I = 2:KM,  A(I) = A(I-1)+NJ(I); end

                 % INITIAL SCALING
for I = 1:N
   if ITOL == 0, SCAL(I) = ATOL(1)+RTOL(1)*abs(Y(I));
   else          SCAL(I) = ATOL(I)+RTOL(I)*abs(Y(I)); end,end

                 % INITIAL PREPARATIONS
if XEND >= X, POSNEG =  1;
else          POSNEG = -1; end

K    = max( 2, min( KM-1, fix(-log10(RTOL(1))*0.6+1.5)));
HMAX = abs(HMAX);
H    = max(abs(H),1E-4); 
H    = POSNEG*min([H,HMAX,abs(XEND-X)/2]);
goto = 10;

if IOUT >= 1
   if IOUT >= 2
      IPOINT(1) = 0;
      for I = 1:KM  
          NJADD = 4*I-2;
          if NJ(I) > NJADD, NJADD = NJADD+1; end
          IPOINT(I+1) = IPOINT(I)+NJADD; end
      
      for MU = 1:KM*2
          ERRX = sqrt(MU/(MU+4))*0.5;
          PROD = 1/(MU+4)^2;
          for J = 1:MU, PROD = PROD*ERRX/J; end
          ERRFAC(MU) = PROD; end 
      
      IPT = 0; end
  
   %IRTRN = 0;
   %XOLD  = X;
   % [dz,IRTRN] = SOLOUT(NACCPT+1,XOLD,X,Y,N,DENS,NCOM,ICOMP,NRD, dt);                     
   % if IRTRN < 0, goto = 120; Fail = 'After SOLOUT IRTRN < 0'; end
   end

if goto == 10
   ERR    = 0; 
   ERROLD = 1e10;
   HOPTDE = POSNEG*HMAX;
   W(1)   = 0;  
   REJECT = false;
   LAST   = false; end
   
                              %  10 120
while goto ~= 120            
   ATOV = false; 
   if goto == 10
      goto = 0;
                       % IS XEND REACHED IN THE NEXT STEP?
      if 0.1*abs(XEND-X) <= abs(X)*UROUND, goto = 110; break; end
      
      H = POSNEG*min([abs(H), abs(XEND-X), HMAX, abs(HOPTDE)]);
      if (X+1.01*H-XEND)*POSNEG > 0 
         H = XEND-X;
         LAST = true; end
                                          
      if NSTEP == 0 || IOUT ~= 2, DZ = FCN( X, Y, Farg{:} ); end
      NFCN = NFCN+1;
         
                       % THE FIRST AND LAST STEP 
      if NSTEP == 0 || LAST 
         IPT = 0;
         NSTEP = NSTEP+1;
         for J = 1:K
            KC = J; 
            [H,YH2,T,HH,W,ERR,ATOV,REJECT,ERROLD,FSAFE,YSAFE,NFCN] ...
            = MIDEX(J,X,Y,H,HMAX,N,FCN,DY,YH1,YH2,DZ,T,NJ,HH,W, ...
                    ERR,A,SAFE1,UROUND,FAC1,FAC2,SAFE2,SCAL,ATOV,SAFE3,  ...
                    REJECT,KM,RTOL,ATOL,ITOL,MSTAB,JSTAB,ERROLD,FSAFE,...
                    LFSAFE,IOUT,IPT,YSAFE,ICOMP,NRD,NFCN,Farg{:});
            if     ATOV,              goto = 10; break;
            elseif J > 1 && ERR <= 1, goto = 60; break; end,end
        
         if goto == 10,  continue
         elseif goto ~= 60,      goto = 55; end,end,end

                         % BASIC INTEGRATION STEP  goto = 0 30 55 60   
   if goto <= 30
      IPT   = 0;
      NSTEP = NSTEP+1;
      if NSTEP == NMAX
         goto = 120; 
         Fail = sprintf('NSTEP = NMAX = %d, Xmax = %g, Xend = %g',NMAX,X,XEND);
         break;end
    
      KC = K-1;
      for J = 1:KC
         [H,YH2,T,HH,W,ERR,ATOV,REJECT,ERROLD,FSAFE,YSAFE,NFCN] ...
         = MIDEX(J,X,Y,H,HMAX,N,FCN,DY,YH1,YH2,DZ,T,NJ,HH,W, ...
                 ERR,A,SAFE1,UROUND,FAC1,FAC2,SAFE2,SCAL,ATOV,SAFE3,  ...
                 REJECT,KM,RTOL,ATOL,ITOL,MSTAB,JSTAB,ERROLD,FSAFE,  ...
                 LFSAFE,IOUT,IPT,YSAFE,ICOMP,NRD,NFCN,Farg{:});
         if ATOV, goto = 10; break; end,end
     
      if goto == 10, continue; end
                                                
                               % CONVERGENCE MONITOR
      if K == 2 || REJECT,                goto = 50; 
      elseif ERR <= 1,                    goto = 60;
      elseif ERR > ((NJ(K+1)*NJ(K))/4)^2, goto = 100; end,end
  
                                               
   if goto <= 50 
      [H,YH2,T,HH,W,ERR,ATOV,REJECT,ERROLD,FSAFE,YSAFE,NFCN] ...
      = MIDEX(K,X,Y,H,HMAX,N,FCN,DY,YH1,YH2,DZ,T,NJ,HH,W,     ...
              ERR,A,SAFE1,UROUND,FAC1,FAC2,SAFE2,SCAL,ATOV,SAFE3, ...
              REJECT,KM,RTOL,ATOL,ITOL,MSTAB,JSTAB,ERROLD,FSAFE, ...
              LFSAFE,IOUT,IPT,YSAFE,ICOMP,NRD,NFCN,Farg{:});
      if ATOV, goto = 10;  continue,end
      KC = K; 
      if ERR <= 1, goto = 60; end,end
         
                              % HOPE FOR CONVERGENCE IN LINE K+1
   if goto <= 55
      if ERR > (NJ(K+1)/2)^2, goto = 100; 
      else
         KC = K+1;
         [H,YH2,T,HH,W,ERR,ATOV,REJECT,ERROLD,FSAFE,YSAFE,NFCN] ...
         = MIDEX(KC,X,Y,H,HMAX,N,FCN,DY,YH1,YH2,DZ,T,NJ,HH,W,  ...
                 ERR,A,SAFE1,UROUND,FAC1,FAC2,SAFE2,SCAL,ATOV,SAFE3,...
                 REJECT,KM,RTOL,ATOL,ITOL,MSTAB,JSTAB,ERROLD,FSAFE,...
                 LFSAFE,IOUT,IPT,YSAFE,ICOMP,NRD,NFCN,Farg{:});
         if ATOV, goto = 10;  continue,end
            if ERR > 1, goto = 100; end,end,end
    
                       % STEP IS REJECTED 
   if goto == 100            
      K = min([K,KC,KM-1]);
      if K > 2 && W(K-1) < W(K)*FAC3, K = K-1; end
      NREJCT = NREJCT+1;
      H = POSNEG*HH(K);
      REJECT = true;
      goto = 30;
      continue,  end
  
                       % STEP IS ACCEPTED 
   XOLD = X;  Yold = Y;      %  goto <= 60
   X    = X+H;
   if IOUT >= 2           %  KMIT = MU OF THE PAPER                      
      KMIT = 2*KC-MUDIF+1;
      for I = 1:NRD,  DENS(I) = Y(ICOMP(I)); end 
      XOLDD = XOLD;
      HHH   = H;
      for I = 1:NRD,  DENS(NRD+I) = H*DZ(ICOMP(I)); end
      KLN = 2*NRD;
      for I = 1:NRD,  DENS(KLN+I) = T(1,ICOMP(I));  end 
         
                              % COMPUTE SOLUTION AT MID-POINT 
      for J = 2:KC,
        DBLENJ = NJ(J);
        for L = J:-1:2
          FACTOR = (DBLENJ/NJ(L-1))^2-1;
          for I = 1:NRD
            YSAFE(L-1,I)=YSAFE(L,I)+(YSAFE(L,I)-YSAFE(L-1,I))/FACTOR;end,end,end
          
      KRN = 4*NRD;
      for I = 1:NRD,  DENS(KRN+I) = YSAFE(1,I); end
            
                              % COMPUTE FIRST DERIVATIVE AT RIGHT END
      for I = 1:N,  YH1(I) = T(1,I); end
      YH2 = FCN( X, Y, Farg{:} );
      KRN = 3*NRD;
      for I = 1:NRD,  DENS(KRN+I) = YH2(ICOMP(I))*H; end     
         
                             % THE LOOP 
      for KMI = 1:KMIT       % COMPUTE KMI-TH DERIVATIVE AT MID-POINT                  
         KBEG = fix((KMI+1)/2);
         for KK = KBEG:KC
            FACNJ = (NJ(KK)/2)^(KMI-1);  
            IPT   = IPOINT(KK+1)-2*KK+KMI;
            for I = 1:NRD, YSAFE(KK,I) = FSAFE(IPT,I)*FACNJ; end,end

         for J = KBEG+1:KC
            DBLENJ = NJ(J);
            for L = J:-1:KBEG+1
               FACTOR = (DBLENJ/NJ(L-1))^2-1;
               for I = 1:NRD
                  YSAFE(L-1,I)=YSAFE(L,I)+(YSAFE(L,I)-YSAFE(L-1,I))/FACTOR;
                  end,end,end
          
         KRN = (KMI+4)*NRD;
         for I = 1:NRD, DENS(KRN+I) = YSAFE(KBEG,I)*H; end
         if KMI == KMIT, continue; end
              
                              % COMPUTE DIFFERENCES
         for KK = fix((KMI+2)/2):KC    % !!!  change to fix((KMI+2)/2)
            LBEG = IPOINT(KK+1);
            LEND = IPOINT(KK)+KMI+1;
            if KMI == 1 && NSEQU == 4, LEND = LEND+2; end
            for L = LBEG:-2:LEND
               for I = 1:NRD
                  FSAFE(L,I) = FSAFE(L,I)-FSAFE(L-2,I); end,end
              
            if KMI == 1 && NSEQU == 4 
               L = LEND-2;
               for I=1:NRD, FSAFE(L,I)=FSAFE(L,I)-DZ(ICOMP(I));end,end,end     

                                % COMPUTE DIFFERENCES
         for KK = fix((KMI+2)/2):KC            % !!!  change to fix((KMI+2)/2)
            LBEG = IPOINT(KK+1)-1;
            LEND = IPOINT(KK)+KMI+2;
            for L = LBEG:-2:LEND
               for I=1:NRD, FSAFE(L,I)=FSAFE(L,I)-FSAFE(L-2,I);end,end,end
         end % KMI

      DENS = INTERP(NRD,DENS,KMIT);
         
                               % ESTIMATION OF INTERPOLATION ERROR  
      if IDERR == 0 && KMIT >= 1
         ERRINT = 0;
         for I = 1:NRD
            ERRINT = ERRINT+(DENS((KMIT+4)*NRD+I)/SCAL(ICOMP(I)))^2;end
             
         ERRINT = sqrt(ERRINT/NRD)*ERRFAC(KMIT);
         HOPTDE = H/max((ERRINT)^(1/(KMIT+4)),0.01);
         if ERRINT > 10  
            H = HOPTDE;
            X = XOLD;
            NREJCT = NREJCT+1;
            REJECT = true;  
            goto = 10;
            continue; end,end

      for I = 1:N, DZ(I) = YH2(I); end 
      end  % if IOUT >= 2
         
      for I = 1:N, Y(I) = T(1,I); end
      NACCPT = NACCPT+1; 
      
   if IOUT < 1                      % without dens calcs                    
      t = [t; X];
      z(NACCPT+1,:) = Y'; end
      
   %else                            % only with dens calcs
   if IOUT >= 2 
      ndt = find( XOLD < t & t <= X );
      if ~isempty(ndt)
         dt = t(ndt);
         [dz,IRTRN] = SOLOUT(NACCPT+1,XOLD,X,Yold,Y,N,DENS,NCOM,ICOMP,NRD,dt);
         if IRTRN<0, goto = 120; Fail = 'After SOLOUT IRTRN < 0';break;end,
         z(ndt,:) = dz; end,end
                         % COMPUTE OPTIMAL ORDER
   if KC == 2
      KOPT = min(3,KM-1);
      if REJECT, KOPT = 2; end  
      goto = 80; end

   if goto ~= 80
      if KC <= K
         KOPT = KC; 
         if W(KC-1) < W(KC)*FAC3, KOPT = KC-1; end  
         if W(KC) < W(KC-1)*FAC4, KOPT = min(KC+1,KM-1); end
      else 
         KOPT = KC-1;
         if KC<3 && W(KC-2)<W(KC-1)*FAC3, KOPT = KC-2; end
         if W(KC) < W(KOPT)*FAC4,         KOPT = min(KC,KM-1);end,end,end

                         % AFTER A REJECTED STEP
   if REJECT               % label 80 
      K = min(KOPT,KC);
      H = POSNEG*min(abs(H),abs(HH(K)));
      REJECT = false;
   else              % COMPUTE STEPSIZE FOR NEXT STEP                    
      if     KOPT <= KC,                     H = HH(KOPT);
      elseif KC < K && W(KC) < W(KC-1)*FAC4, H = HH(KC)*A(KOPT+1)/A(KC);
      else                                   H = HH(KC)*A(KOPT)/A(KC); end
          
      K = KOPT;
      H = POSNEG*abs(H);end
     
   goto = 10;                    
end    

if ischar(Fail),  nf = Fail;    
else              nf = [NFCN, NSTEP, NACCPT, NREJCT];end


%****************************************************************************
%****************************************************************************
%          J,X,Y,H,HMAX,N,FCN,DY,YH1,YH2,DZ,T,NJ,HH,W, ...
%          ERR,A,SAFE1,UROUND,FAC1,FAC2,SAFE2,SCAL,SAFE3,  ...
%          KM,RTOL,ATOL,ITOL,MSTAB,JSTAB,ERROLD,FSAFE,...
%          LFSAFE,IOUT,IPT,YSAFE,ICOMP,NRD,NFCN,Farg{:}
function [H,YH2,T,HH,W,ERR,ATOV,REJECT,ERROLD,FSAFE,YSAFE,NFCN] = ...
   MIDEX( J, X, Y, H, HMAX, N, FCN, DY, YH1,YH2, DZ, T, NJ, HH, W, ...
          ERR,  A, SAFE1, UROUND, FAC1, FAC2, SAFE2, SCAL, ATOV,SAFE3,...
          REJECT, KM, RTOL, ATOL, ITOL, MSTAB, JSTAB, ERROLD, FSAFE, ...
          LFSAFE, IOUT, IPT, YSAFE, ICOMP, NRD, NFCN, varargin)

% variables FAC,DY,KM, LFSAVE, RPAR,IPAR don't needed, but used in Fortran-odex    

Farg = varargin;
HJ = H/NJ(J);
                   % EULER STARTING STEP
for I = 1:N,   YH1(I) = Y(I); YH2(I) = Y(I)+HJ*DZ(I);  end

                   % EXPLICIT MIDPOINT RULE  
M     = NJ(J)-1;
NJMID = NJ(J)/2;
goto  = 0;

for MM = 1:M 
   if IOUT >= 2 && MM == NJMID
      for I = 1:NRD,  YSAFE(J,I) = YH2(ICOMP(I)); end,end
                  
   DY = FCN( X+HJ*MM, YH2, Farg{:}); 
   if IOUT >= 2 && abs(MM-NJMID) <= 2*J-1 
      IPT = IPT+1;
      for I = 1:NRD, FSAFE(IPT,I) = DY(ICOMP(I)); end,end

   for I = 1:N
      YS     = YH1(I);  
      YH1(I) = YH2(I);
      YH2(I) = YS+2*HJ*DY(I); end

   if MM <= MSTAB && J <= JSTAB          % STABILITY CHECK               
      DEL1 = 0;
      for I = 1:N, DEL1 = DEL1+(DZ(I)/SCAL(I))^2; end
      DEL2 = 0;
      for I = 1:N, DEL2 = DEL2+((DY(I)-DZ(I))/SCAL(I))^2; end
      QUOT = DEL2/max(UROUND,DEL1);
      if QUOT > 4
         NFCN = NFCN+1;
         goto = 79; 
         break; end,end,end 
      
if goto ~= 79      % FINAL SMOOTHING STEP  
   DY = FCN( X+H, YH2, Farg{:});  
   if IOUT >= 2 && NJMID <= 2*J-1 
      IPT = IPT+1;
      for I = 1:NRD, FSAFE(IPT,I) = DY(ICOMP(I)); end,end

   for I = 1:N, T(J,I) = (YH1(I)+YH2(I)+HJ*DY(I))/2; end
   NFCN = NFCN+NJ(J); 
   
                  % POLYNOMIAL EXTRAPOLATION
   if J == 1, return; end
   
   DBLENJ = NJ(J);
   for L = J:-1:2
      FAC = (DBLENJ/NJ(L-1))^2-1;
      for I = 1:N,  T(L-1,I) = T(L,I)+(T(L,I)-T(L-1,I))/FAC; end,end
  
   ERR = 0;
                    % SCALING
   for I = 1:N
      T1I = max(abs(Y(I)),abs(T(1,I)));
      if ITOL == 0, SCAL(I) = ATOL(1)+RTOL(1)*T1I;
      else          SCAL(I) = ATOL(I)+RTOL(I)*T1I; end
      ERR = ERR+((T(1,I)-T(2,I))/SCAL(I))^2; end 
  
   ERR = sqrt(ERR/N);
   if     ERR*UROUND >= 1,        goto = 79; 
   elseif J > 2 && ERR >= ERROLD, goto = 79; end,end

if goto ~= 79
   ERROLD = max(4*ERR,1);
               % COMPUTE OPTIMAL STEPSIZES
   EXPO   = 1/(2*J-1);
   FACMIN = FAC1^EXPO;
   FAC    = min( FAC2/FACMIN,  max(FACMIN, (ERR/SAFE1)^EXPO/SAFE2));
   FAC    = 1/FAC;
   HH(J)  = min(abs(H)*FAC, HMAX);
   W(J)   = A(J)/HH(J);
else       
   ATOV   = true;
   H      = H*SAFE3;
   REJECT = true; end  
 
% t
%******************************************************************************
%******************************************************************************
%     COMPUTES THE COEFFICIENTS OF THE INTERPOLATION FORMULA

function   DENS = INTERP(NRD,DENS,KMIT)  

% DIMENSION DENS(NRD*(KMIT+5)),  A(0:30)       
A = zeros(1,31); 
                      % BEGIN WITH HERMITE INTERPOLATION
for I = 1:NRD
   Y0    = DENS(I);
   Y1    = DENS(2*NRD+I);
   YP0   = DENS(NRD+I);
   YP1   = DENS(3*NRD+I); 
   YDIFF = Y1-Y0;
   ASPL  = -YP1+YDIFF;
   BSPL  =  YP0-YDIFF;
   DENS(NRD+I)   = YDIFF;
   DENS(2*NRD+I) = ASPL;
   DENS(3*NRD+I) = BSPL; 
   if KMIT < 0, continue; end 
                   % COMPUTE THE DERIVATIVES OF HERMITE AT MIDPOINT
   PH0 = (Y0+Y1)*0.5+0.125*(ASPL+BSPL);
   PH1 = YDIFF+(ASPL-BSPL)*0.25;
   PH2 = -(YP0-YP1);
   PH3 = 6*(BSPL-ASPL);
                  %  COMPUTE THE FURTHER COEFFICIENTS 
   if KMIT >= 1
     A(2) = 16*(DENS(5*NRD+I)-PH1);
     if KMIT >= 3
       A(4)=16*(DENS(7*NRD+I)-PH3+3*A(2));
       if KMIT >= 5
         for IM = 5:2:KMIT
           FAC1 = IM*(IM-1)/2;
           FAC2 = FAC1*(IM-2)*(IM-3)*2;
           A(IM+1) = 16*(DENS((IM+4)*NRD+I)+FAC1*A(IM-1)-FAC2*A(IM-3));end
           end,end,end
  
   A(1) = (DENS(4*NRD+I)-PH0)*16;
   if KMIT >= 2
      A(3) = (DENS(NRD*6+I)-PH2+A(1))*16;
      if KMIT >= 4
         for IM = 4:2:KMIT
            FAC1 = IM*(IM-1)/2;
            FAC2 = IM*(IM-1)*(IM-2)*(IM-3);
            A(IM+1) = (DENS(NRD*(IM+4)+I)+A(IM-1)*FAC1-A(IM-3)*FAC2)*16;end
            end,end

   for IM = 0:KMIT, DENS(NRD*(IM+4)+I) = A(IM+1); end,end
            

%******************************************************************************
%******************************************************************************
% 
function [dz, IRTRN] = SOLOUT(NR,XOLD,X,Yold,Y,N,DENS,NCOM,ICOMP,NRD, dt)
%SOLOUT      NAME (EXTERNAL -> internal) OF SUBROUTINE PROVIDING THE
%             NUMERICAL SOLUTION DURING INTEGRATION. 
%             IF IOUT.GE.1, IT IS CALLED AFTER EVERY SUCCESSFUL STEP.
%             SUPPLY A DUMMY SUBROUTINE IF IOUT=0.   
%             SOLOUT FURNISHES THE SOLUTION "Y" AT THE NR-TH
%                GRID-POINT "X" (THEREBY THE INITIAL VALUE IS
%                THE FIRST GRID-POINT).
%             "XOLD" IS THE PRECEEDING GRID-POINT.
%             "IRTRN" SERVES TO INTERRUPT THE INTEGRATION. IF IRTRN
%                IS SET <0, ODEX WILL RETURN TO THE CALLING PROGRAM.
%       
%      -----  CONTINUOUS OUTPUT (IF IOUT=2): -----
%             DURING CALLS TO "SOLOUT", A CONTINUOUS SOLUTION
%             FOR THE INTERVAL [XOLD,X] IS AVAILABLE THROUGH
%             THE DOUBLE PRECISION FUNCTION
%                >>>   CONTEX(I,S,DENS,NCOM,ICOMP,NRD)   <<<
%             WHICH PROVIDES AN APPROXIMATION TO THE I-TH
%             COMPONENT OF THE SOLUTION AT THE POINT S. THE VALUE
%             S SHOULD LIE IN THE INTERVAL [XOLD,X].

Sdt = 1;
Ldt = length(dt);
dz  = nan(Ldt,N);

if dt(Sdt)-XOLD <= eps*abs(XOLD)
   dz(Sdt,:) = Yold;
   Sdt = Sdt+1; end

if X-dt(Ldt) <= eps*abs(X)
   dz(Ldt,:) = Y;
   Ldt = Ldt-1; end

IRTRN = 1; 
for r = Sdt:Ldt
   S = dt(r);
   for k = 1:N
      cx = CONTEX( k, S, DENS, NCOM, ICOMP, NRD); 
      if ischar(cx), break
      else           dz(r,k) = cx; end,end
   if ischar(cx), IRTRN = -1; break; end,end
   
%******************************************************************************
%******************************************************************************
%     THIS FUNCTION CAN BE USED FOR CONINUOUS OUTPUT IN CONECTION
%     WITH THE OUTPUT-SUBROUTINE FOR ODEX. IT PROVIDES AN
%     APPROXIMATION TO THE k-TH COMPONENT OF THE SOLUTION AT S.

function cx = CONTEX( k, S, DENS, NCOM, ICOMP, NRD)

global  XOLDD  HHH  KMIT          
                       % COMPUTE PLACE OF k-TH COMPONENT 
I = 0; 
for J = 1:NRD, if ICOMP(J) == k, I = J; end,end

if I == 0,  cx = sprintf('NO DENSE OUTPUT AVAILABLE FOR COMP. %d',k); 
else 
                      % COMPUTE THE INTERPOLATED VALUE 
  THETA  = (S-XOLDD)/HHH;
  THETA1 = 1-THETA;
  PHTHET = DENS(I)+ ...
        THETA*(DENS(NRD+I)+THETA1*(DENS(2*NRD+I)*THETA+DENS(3*NRD+I)*THETA1));
  if KMIT < 0, cx = PHTHET;
  else
     THETAH = THETA-0.5;
     cx = DENS(NRD*(KMIT+4)+I);
     for IM = KMIT:-1:1,    cx = DENS(NRD*(IM+3)+I)+cx*THETAH/IM; end
     cx = PHTHET+(THETA*THETA1)^2*cx; end,end
      