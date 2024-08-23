% @(t) is slowler f(t) in the body, but
% @(t) is faster if they are passed as parameters to another function.

function DSvsALi__

ticAll = tic;
                           %  PRESET CONSTANTS                         
naODEs = { 'N3_43',  'N3_46', 'N4_35', 'N8_57',  'Nord6'};
  
                % Vector-functions of ODE y'=F(t,y)
Fs={@(t,y,z) [ y(2); y(3); -(3*y(3)+4*y(1))/t-(4+1/t^2)*y(2) ],  ...    
    @(t,y,z) [ y(2); y(3); -5*y(3)/t-4*y(2)/t^2+log(t) ],     ...       
    @(t,y,z) [ y(2); y(3); y(4); -8*y(4)/t-12*y(3)/t^2-y(1)/t^4 ],... 
    @(t,y,z) N8_57( t, y, z ), ...                                  
    @(t,y,z) Nord6( t, y, z ) }; 

               % Generations matpicies of ODE y'=G*y               
Gs={ @(t,A) [0 1 0; 0 0 1; -4/t, -(4+1/t^2), -3/t ],   ...  % ds10L N3_43
     { @(t,A) [0 1 0; 0 0 1;    0,    -4/t^2,  -5/t ], ...  % ds10L N3_46
       @(t,A) [0; 0;log(t)] },                ... % ds10LI inhomogeneouse N3_46
     @(t,A) [0 1 0 0;  0 0 1 0;  0 0 0 1; -1/t^4, 0, -12/t^2, -8/t],...% N4_35
     @(t,A) N8_57_G(t,A),                                           ...% N8_57
     { @(t,A) Nord6(t,[],A), [] }};               % ds10Lc constructed  Nord6
      
  
  t0fs  = [ [1 21];   [1 21];  [1 21];   [0 20];   [0 20]];
  t0 = t0fs(1,1);
   J0 = besselj(0,t0); Y0 = bessely(0,t0); 
   J1 = besselj(1,t0); Y1 = bessely(1,t0);  JY = J0*Y1+J1*Y0;
y01 = [J0*Y0; -JY; 2*(J1*Y1-J0*Y0)+JY/t0];
y02 = [1.5; -0.25; -0.75];
   sqr5 = sqrt(5); sqr1 = (sqr5-1)*0.5; sqr2 = -(sqr5+1)*0.5;
y03 = [1;  0.5*(sqr5+1); -2*sqr5; 8*sqr5+2]; 

   %     a b  c      для этих a,b,c проверить N8_57
   % 1) -1 2 -1      sqr=5  al=2 be=-3 ga=0.5 de=-2 
   % 2) -1 2  2      sqr=4  al=3 be=-1 ga=1   de=-1
   % 3) -1 3  2      sqr=6  al=4 be=-2 ga=1   de=-1
   % 4) -1 4  2      sqr=8  al=5 be=-3 ga=1   de=-1
   % 5)  1 2  1      sqr=5  al=3 be=-2 ga=2   de=-0.5
   % 6)  1 4  4      sqr=10 al=7 be=-3 ga=2   de=-0.5
   % 7)  2 2 -1      sqr=5  al=2 be=-3 ga=2   de=-0.5
   % 8)  2 3  4      sqr=10 al=7 be=-3 ga=3   de=-0.333333
   % 9)  2 4  2      sqr=10 al=6 be=-4 ga=2   de=-0.5
   % 10) 3 3  2      sqr=10 al=6 be=-4 ga=3   de=-0.333333
   ia = 5;
   as = [ -1 -1 -1 -1 1 1  2 2 2 3 ];   a = as(ia); 
   bs = [  2  2  3  4 2 4  2 3 4 3 ];   b = bs(ia); 
   cs = [ -1  2  2  2 1 4 -1 4 2 2 ];   c = cs(ia);  abc = [a b c];
   sqr = sqrt((2*a+c)^2+4*b^2); 
   al = 0.5*(c+sqr);  ga = (a+al)/b;  
   be = 0.5*(c-sqr);  de = (a+be)/b;    gd = ga+de;
y04 = [ 2; -2; gd;-gd ];
y05 = [ 0; 1; 1; 0; 0; 1 ];
y0s = {y01,  y02,  y03, y04, y05};

Fargs = cell(size(naODEs)); Fargs{4} = [ a b c ];       
er = @(ye,y) sqrt(sum((ye(:,1)-y(:,1)).^2)/length(y)); % |ApproxSol-ExactSol| 
                            
                     % SOURCE DATA for odex                    
ITOL  = 0;
%IOUT = 3;           my cange: IOUT calcs automatically by dim(X)
WORK  = zeros(1,9);     
IWORK = zeros(1,71); % 2*KM+21+NRDENS
IWORK(1) = 15000;    % NMAX, default = 1e4
IWORK(7) = 4;        % DEGREE OF INTERPOL F-LA MUDIF=4 (DEFAULT) LIE in [1, 6]

% H    = opt{1};      % Initial stepsize  
% RTOL = opt{2};      % local Relative tolerance
% ATOL = opt{3};      % local Absolute tolerance
% ITOL = opt{4};      % ITOL=0: RTOL AND ATOL ARE SCALARS, otherwise vectors
                 % SOLOUT = opt{5};   don't used
% IOUT = opt{5};    % 0:SOLOUT isn't called, 1:usual output, 3:DENSE OUTPUT
% WORK = opt{6};    % WORKING array: NRDENS=IWORK(8), etc 
                 % LWORK    LENGTH OF WORK don't used
% IWORK = opt{7};    % IWORK(1:20)- PARAMETERS FOR THE CODE
                 % LIWORK      LENGTH OF ARRAY "IWORK" don't used

                                 
namets_c = {'dop853' 'ode113_c' 'ode45_c' 'odex'};     % Another methods
mets    = cellfun(@str2func,namets_c,'Uniform',false); 
namets  = cellfun(@(x) strrep(x,'_c',''), namets_c,'Uniform',false);  
nametsD = ['ds10' namets];                             % D = ds10+Another 
namets_ = nametsD(1:end-1);                            % _ = D-BadMet

setmsD = {@(h,RT,AT,d) DSet(h,RT,AT, d{1},d{2} ),               ...  %   ds10
          @(h,RT,AT,d) dopset('InitialSt',h,'RelT',RT,'AbsT',AT),... % 1 dop853
          @(h,RT,AT,d) odeset('InitialSt',h,'RelT',RT,'AbsT',AT),... % 2 ode113
          @(h,RT,AT,d) odeset('InitialSt',h,'RelT',RT,'AbsT',AT),... % 3 ode45 
          @(h,RT,AT,d) {h, RT, AT, ITOL,  WORK, IWORK} };            % 4 odex 
setms = setmsD(2:end);
                      % SOURCE DATA for ds10                               
m    = 3;             % number of parameters of the difference scheme
kLmn = 1;             % number of choice of coefs La,Mu,Nu

                        % SOURCE DATA
h0 = 0.01;
nODEs = [1 2 3];        % equations
nmets = 1:4;            % methods
RT = eps; AT = eps;
qrun  = 5;   fprintf('%30s qrun = %d',' ',qrun);   % q-ty of method runs
                               % Dependent CONSTANTS  
Lnmets = length(nmets);        % all methods : 1:4  
NDs= {[100 250 500 1000 2500 5000 7500 10000 12500 15000 20000]; ... % n66
      [125 250 500 1000 1200 1250 1500 2000 2500 3000 4000 8000 12000];... %H146
      [150 250 500 1000 2500 5000 7500 10000 12500 15000 20000]};    % fant
         
NAs= {[500 1000 2500 5000 7500 10000 12500 15000 20000]; ...  % n66
      [250 500 1000 1500 2000 3000 4000  8000  12000 ];  ...  % H146
      [500 1000 2500 5000 7500 10000 12500 15000 20000]};     % fant
      
N0 = [1 1 1];    % какого номера по абсциссе строить график
basa = [6:2:14 14.8 15.2 15.55]; Tail = 1;
tols = 10.^-basa;  Ltols = length(tols);
zn = '^so*d'; 
natols = cellfun(@(x) strrep(num2str(x,'e_%.3g'),'.','_'),num2cell(basa),...
         'Uniform',false);    % replace . -> _
I  = 1:Lnmets;
dJ = [2 2 1];     % how many right extra points is in ED
TOC = zeros(1,length(nODEs));

                              % LOOPS                            
for j = nODEs
   ticj = tic;
   naODE = naODEs{j}; 
   load( ['BestC_N_' naODE] );     % 'BestC','t0f','NDS','Ermin','sErmin'
   ND = NDs{j}; 
   F    = str2func( naODE );     % ODE func   
   
   G    = Gs{j};
   if length(G) == 1,     ds10 = @ds10L;
   elseif ishandle(G{2}), ds10 = @ds10LI;  
   else                   ds10 = @ds10Lc; end
   
   t0f  = t0fs(j,:); 
   Farg = Fargs{j};  
   fprintf('\n\n%30s ODEfunc = %s   t0f = [%g  %g]\n',' ',naODE,t0f);
   
   t0  = t0f(1);      tf = t0f(2);  
   y0  = y0s{j};
   ES1 = ES1s{j};   % ?
   NA = NAs{j};
   MD = length(ND);  LND = log2(ND/100);
   MA = length(NA);  LNA = log2(NA/100);  M1 = MA+1;
   ED = nan(1,MD);              TD = ED;  nD = TD;
   EA = nan(Lnmets, Ltols,M1);  TA = EA;  nA = TA;
   nB = cell(1,M1);
   hsD = (tf-t0)./ND;  
   hsA = (tf-t0)./NA;  
   hs_  = sprintf('%.4f     ',hsD);  fprintf('hsD = %3s %s\n',' ',hs_);
   naND = cellfun(@(x) num2str(x,'N_%d'),num2cell(ND),'Uniform',false);
   
                              % ds10 method
   for r = MD:-1:1
      h   = hsD(r);      Nr = ND(r);
      t   = (t0:h:tf)'; 
      Nh  = find(NDS == Nr);          % Nh = r, h's number in BestC
      Lmn = BestC{m}{Nh}(kLmn,1:m);
      op  = setmsD{1}( h, [],[],{Lmn []});
      [t1,y,nf] = ds10(G,t,y0,op,Farg);   
      tic; for l=1:qrun, [t1,y,nf] = ds10L(G,t,y0,op,Farg);end; T=toc; 
      ye=ES1(t);  TD(r)= T/qrun;   ED(r)= er(ye,y);  nD(r)= nf(1); end
  
   EDt = str2num(num2str(ED,'%.2g  '));
   TDt = str2num(num2str(TD,'%.2g  '));
   ETN = array2table([EDt; TDt; nD],'Var',naND,'Row',{'ES' 'TS' 'NS'})
   
   [mi, imi] = min(ED); 
   J   = 1:imi-dJ(j);      % exclude extra points from right hand
   LED = -log10(ED);       LnD = log10(nD); 
   
                              % All other methods
   for r = M1:-1:1          % At first: big time for big N  
      nB{r} = []; % numbers of Bad meth- that did not find the solution 
      if r == M1
         t = [t0 tf];  
         fprintf('\n\n%30s Non dense output\n',' '); 
      else
         h = hsA(r);
         t = (t0:h:tf)';
         fprintf('\n\n%30s N = %d   h = %.3g\n',' ',NA(r),h); end
      
                              % Time, Err, Nfs by Another METHODs    
      for i = 1:Lnmets        % number of method
         nm   = nmets(i);     % count  of method 
         met  = mets{nm};  
         setm = setms{nm};
         
         for s = 1:Ltols
            if strcmp('odex',namets(i))
               if r < M1,  tol = tols(s)*1e3;
               else        tol = tols(s); end
            else           tol = tols(s); end
            
            op = setm( h0, tol,tol,[]);
            T = 0;
            warning('off','all');
            [t1,y,nf] = met(F,t,y0,op,Farg);
            for l=1:qrun, tic; [t1,y,nf]= met(F,t,y0,op,Farg); T=T+toc; end
            warning('on','all');
            
            if isnumeric(t1)
               ye=ES1(t1); 
               EA(i,s,r)=er(ye,y);  TA(i,s,r)=T/qrun;  nA(i,s,r)=nf(1); 
            else 
               if s==1, BadM = {'NaN' 'NaN'};
               else     BadM = {sprintf('E^*=10^{%.2g}  T^*=%.2g',...
                               log10(EA(i,s-1,r)), TA(i,s-1,r))...
                               sprintf('nf^*=%d',nA(i,s-1,r))};end
               nB{r} = [nB{r} i]; 
               break; end,end,end

      EAt = str2num(num2str(EA(:,:,r),'%.3g  '));
      TAt = str2num(num2str(TA(:,:,r),'%.2g  '));
      ERR = array2table(EAt,'Var',natols,'Row',namets)
      TIM = array2table(TAt,'Var',natols,'Row',namets)
      NFU = array2table(nA(:,:,r),'Var',natols,'Row',namets), end
  
                      % OUTPUT  for j
                      % 1 Fig  NON DENSE  All methods                    
   LEA = -log10(EA);  
   LnA =  log10(nA);   
   figure('Name',[naODE ', ND_Tn(E)']); 
   subplot(2,1,1);          %  T(E)
   p1 = plot(LED(J),TD(J), LEA(:,:,M1)',TA(:,:,M1)');
   %title([naODE ', non dense output']);
   for p=1:Lnmets+1, p1(p).Marker = zn(p); p1(p).LineWidth = 1;end 
   ylabel('Time');   % xlabel('-log_{10}E');
   legend(nametsD,'Location','northwest'); 
   
   subplot(2,1,2);           % nf(E)
   p2 = plot(LED(J),LnD(J), LEA(:,:,M1)',LnA(:,:,M1)'); 
   for p=1:Lnmets+1, p2(p).Marker = zn(p); p2(p).LineWidth = 1;end 
   ylabel('log_{10}nf');  xlabel('-log_{10}E'); 
   legend(nametsD,'Location','northwest'); 
   
                      % 2 Fig  DENSE  exclude odex(BadMet)                      
   figure('Name',[naODE ', D_Tn(E)']); 
   k = 0;
   for r = [1 MA]
      I_ = I;   I_(nB{r}) = [];
      k = k+1;
      ps = subplot(3,1,k);      %   Time(E)
      p3 = plot(LED(J),TD(J),LEA(I_,:,r)',TA(I_,:,r)'); % exclude increasing ED
      if k == 1,  text( ps.XTick(4), ps.YTick(2), BadM{1} ); end
      for p=1:Lnmets, p3(p).Marker = zn(p); p3(p).LineWidth = 1;end 
      title( sprintf('N=%d  h=%.2g',NA(r),hsA(r)) );
      ylabel('Time');   % xlabel('-log_{10}E');
      legend(namets_,'Location','northwest');  end 
      
   ps = subplot(3,1,3);         %   nf(E) 
   r = 1;
   p4 = plot(LED(J),LnD(J), LEA(:,:,r)',LnA(:,:,r)'); 
   %text( ps.XTick(4), ps.YTick(3), BadM{2} ); 
   for p=1:Lnmets+1, p4(p).Marker = zn(p); p4(p).LineWidth = 1;end 
   title( sprintf('N=%d  h=%.2g',NA(r),hsA(r)) ); 
   ylabel('log_{10}nf');  xlabel('-log_{10}E'); 
   legend(nametsD,'Location','northwest');  

                       % 3,... Fig   
   pTA  = permute( TA(:,:,1:MA), [1 3 2]);
   pLEA = permute(LEA(:,:,1:MA), [1 3 2]);
   K = LND > LNA(N0(j))*(1-eps); % in order for the definition areas to coincide
   scalog = floor(LNA(1)) : ceil(LND(end));
   scale = num2cell(2.^scalog*100);
   L_T = Ltols-Tail;
   for s = L_T+1:Ltols
      figure('Name',sprintf('D_TE(N) %d) %s tol=%.2g',s-L_T,naODE,tols(s)));  
      ps = subplot(2,1,1);          %   Time(N)
      p5 = plot(LND(K),TD(K), LNA,pTA(I_,:,s) );  
      ps.XTickLabel = scale;
      for p=1:Lnmets, p5(p).Marker = zn(p); p5(p).LineWidth = 1;end 
      legend(namets_,'Location','northwest');
      ylabel('Time'); 
  
      ps = subplot(2,1,2);          %   E(N)
      p6 = plot(LND(K),LED(K), LNA,pLEA(I_,:,s) );
      ps.XTickLabel = scale;
      for p=1:Lnmets, p6(p).Marker = zn(p); p6(p).LineWidth = 1;end
      %title('Maximal accuracy');
      legend(namets_,'Location','southeast');
      ylabel('-log_{10}E');xlabel('N'); end
   
   TOC(j) = toc(ticj);
   end  % for j 
tocAll = toc(ticAll);   
fprintf('Total time = %s= %.3g < %.3g\n',sprintf('+ %g ',TOC),sum(TOC),tocAll)


                          % Vector-functions of ODE y'=F(t,y)
function F = N3_43( t, y, empty )  % F: eq y'=F(t,y) N3_43 of Kamke p.465
F = [ y(2); y(3); -(3*y(3)+4*y(1))/t-(4+1/t^2)*y(2) ];

function F = N3_46( t, y, empty )   % F: eq y'=F(t,y) N3_46 of Kamke p.466
F = [ y(2); y(3); -5*y(3)/t - 4*y(2)/t^2 + log(t)];

function F = N4_35( t, y, empty ) % F: y'=F Kamke p.478 homogeneous Euler eq
F = [ y(2); y(3); y(4); -8*y(4)/t - 12*y(3)/t^2 - y(1)/t^4 ];

function F = N8_57( t, y, A )  % F: y'=F  N8_57 of Kamke p.540, тOB p.16
% A = [ a b c ]
ct = A(3)*t;
bc = A(2)*cos(ct);
bs = A(2)*sin(ct);
F = [            A(1)*y(2)  +   bc*y(3) +   bs*y(4); ...
      -A(1)*y(1)            +   bs*y(3) -   bc*y(4); ...
        -bc*y(1) - bs*y(2)              + A(1)*y(4); ...
        -bs*y(1) + bc*y(2)  - A(1)*y(3) ];
    
function G = N8_57_G( t, A ) % Generetion matrix: y'=G*y Kamke p.540, тOB p.16
ct = A(3)*t;
bc = A(2)*cos(ct);
bs = A(2)*sin(ct);
n = length(t);
z = zeros(1,n);
a = A(1)*ones(1,n);
G = [ z   a  bc  bs; ...
     -a   z  bs -bc; ...
     -bc -bs z   a;...
     -bs  bc -a  z ]; 
    
function Aud = Nord6( t, y, empty )   % F: eq y'=F(y)  тOB p.17
s  = sin(t);   ss  = s.*s;      es  = exp(s);    as  = atan(s);  
               ss1 = ss+1;      os1 = 1./ss1;    Ls  = log(ss1);
c  = cos(t);   cs1 = c.*os1;    ces = c.*es;     sos = 2*s.*c.*os1;   
t1 = t+1;      ot1 = 1./t1;     r1  = sqrt(t1);  or1 = 1./r1;     
               L1  = log(t1);   cL  = cos(L1);   sL  = sin(L1).*ot1;
at = atan(t);  o21 = 1./(t.*t+1); tt = 2*t.*o21;
u = [ Ls;  at;  es;   cL; as;  o21 ];
d = [ sos; o21; ces; -sL; cs1; -tt ];      
A = [ t1  -sos   ot1 -tt   es   L1;  ...
      c   -cs1   at   cL   or1  o21; ...
     -r1  -es   -sL   ss  -as  -ces; ...
     -o21  r1    s    ces -Ls  -ot1; ...
      Ls  -at    as   t   -ss   sos; ...
     -cL   tt   -L1  -or1  cs1  sL ];
n = length(t);
if n > 1
   m  = length(u);
   AA = cell(1,n);
   for i = 1:n, AA{i} = A(:,i:n:i+(m-1)*n); end
   Aud = {AA, u, d};
else Aud = A*(y-u)+d; end    
  
                     % exact solutions
function es = ess(j,k,t,A)
if j == 1                                                            % N3_43
   J0 = besselj(0,t); Y0 = bessely(0,t); 
   es = J0.*Y0;
   if k ~= 1
      J1 = besselj(1,t); Y1 = bessely(1,t);  JY = J0.*Y1+J1.*Y0;
      es = [es; -JY; 2*(J1.*Y1-J0.*Y0)+JY./t]; end
  
elseif j == 2                                                        % N3_46
   Lt = log(t); ot = 1/t; 
   es = 1-0.5*t+ot+Lt*(0.25*t+ot);
   if k ~= 1
      ot2 = ot/t; s = 0.25-ot2;
      es = [es; -0.25+s*Lt; (s+2*Lt*ot2)*ot ]; end 
  
elseif j == 3                                                        % N4_35h
    Lt = log(t);       s5 = sqrt(5);
    s51 = (s5-1)*0.5;  t_ = t.^s51;    tp = 1./(t_*t); 
    es = t_ + Lt.*tp;
    if k ~= 1  
       ot = 1./t; ot2 = ot.*ot;  
       es = [es;       
             (  (s5- 1)*t_ +          (2-(s5+1)*Lt).*tp).*ot*0.5; ...
             ( -(s5- 2)*t_ +          (s5+2)*(Lt-1).*tp).*ot2;   ...
             ((7*s5-15)*t_ + (9*s5+19-(7*s5+15)*Lt).*tp).*ot2.*ot*0.5 ]; end
     
elseif j == 4                                                        % N8_57 
   a = A(1);  b = A(2);  c = A(3);
   sqr = sqrt((2*a+c)^2+4*b^2);        
   al = 0.5*(c+sqr);  at = al*t;  ca = cos(at);  sa = sin(at);
   be = 0.5*(c-sqr);  bt = be*t;  cb = cos(bt);  sb = sin(bt);          
   es = ca+sa+cb+sb;
   if k ~= 1
      ga = (a+al)/b;   de = (a+be)/b; 
      es = [es; sa-ca+sb-cb; ga*(sb+cb)+de*(sa+ca); ga*(-cb+sb)+de*(-ca+sa)];end
  
else                                                                 % Nord6
   s  = sin(t); 
   s2 = s.*s+1;
   es = log(s2);
   if k > 1
      ex = exp(s);
      ot2 = 1./(t.*t+1);
      t1 = t+1;
      Lt1 = log(t1);
      es = [ es; atan(t);  ex;  cos(Lt1);  atan(s); ot2 ];
      if k > 2
      c = cos(t);
      es = { es, ...
             2*s*c/s2; ot2; c*ex; -sin(Lt1)/t1; -2*t*ot2^2 }; end,end
end
