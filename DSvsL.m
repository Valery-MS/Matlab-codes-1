%  ds vs other methods for Linear ODEs
% Celecting of coeffs Arg for the best relation Eds/EO

function DSvsL
ticAll = tic;
                           %  PRESET CONSTANTS                         
naODEs = { 'N3_43',  'N4_35',  'N8_57',  'N3_46',   'Nord6'};
nHom = 3;                   % number of homogeneouse ODEs  

t0fs  = [ [1 21];   [2 22];  [0 20];    [1 21];  [0 21]];
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
   cs = [ -1  2  2  2 1 4 -1 4 2 2 ];   c = cs(ia);  

Args = cell(size(naODEs)); Args{nHom} = [ a b c ];       
er = @(ye,y) sqrt(sum((ye(:,1)-y(:,1)).^2)/length(y)); % |ApproxSol-ExactSol| 

                     % SOURCE DATA for odex    
tol_odex = 1e-9;
ITOL  = 0;
%IOUT = 3;           my change: IOUT calcs automatically by dim(X)
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
namets_c = {'odex'}; 
mets    = cellfun(@str2func,namets_c,'Uniform',false); 
namets  = cellfun(@(x) strrep(x,'_c',''), namets_c,'Uniform',false);  
nametsD = ['ds10' namets];                             % D = ds10+Another 
namets_ = nametsD(1:end-1);                            % _ = D-BadMet

setms = {@(h,R,A) dopset('InitialSt',h,'RelT',R,'AbsT',A),... % 1 dop853
         @(h,R,A) odeset('InitialSt',h,'RelT',R,'AbsT',A),... % 2 ode113
         @(h,R,A) odeset('InitialSt',h,'RelT',R,'AbsT',A),... % 3 ode45 
         @(h,R,A) struct('InitialStep',h,'RelTol',R,'AbsTol',A,...
                    'ITOL',ITOL, 'WORK',WORK, 'IWORK',IWORK)};    % 4 odex

                      % SOURCE DATA for ds10                               
m    = 3;             % number of parameters of the difference scheme
kLmn = 1;             % number of choice of coefs La,Mu,Nu

                        % SOURCE DATA
h0 = 0.01;
nODEs = [5];            % equations
nmets = 1;            % methods
%RT = eps; AT = eps;
qrun  = 0;   fprintf('%30s qrun = %d\n',' ',qrun);   % q-ty of method runs
qrun0 = max(1,qrun);
fprintf('%5s %10s %7s %4s %4s %3s %3s %12s\n',...
        'Arg','EDmi','EOmi','Rat','Nh','kLm','ND','Lmn_mi');

                        % Dependent CONSTANTS  
Lnmets = length(nmets);        % all methods : 1:4  
I   = 1:Lnmets;
%NDS = [150 250 500 1000 2500 5000 7500 10000 12500 15000 20000];    
N0 = [1 1 1 1 1];    % какого номера по абсциссе строить график
%basa = [6:2:14 14.8 15.2 15.65]; Tail = 1;
basa = [15.65];
tols = 10.^-basa;  Ltols = length(tols);
zn = '^so*d'; 
natols = cellfun(@(x) strrep(num2str(x,'e_%.3g'),'.','_'),num2cell(basa),...
         'Uniform',false);    % replace . -> _
TOC = zeros(1,length(nODEs));

                       % FLAGs
dJ  = [ 0 0 0 0 0 ] ;  % [2 2 1];       % how many right extra points is in ED
ef  = 1;               % Effective way of computing Gen.Matrix
nIL = 1;                                % 1/2 if takeoff by es/(dop853|ode113)
ILs = [8 1];           IL = ILs(nIL);   % Initial values Length of y 
pri_ds = false;        % print solution by ds10
afteCr = [ 3 0 3 0 ];  % number of out-points after Crash
expr = true;           % experiment
a_d  = [];
OUT  = false;
EO_ED = 0.5;             % = max(EO/ED)
Cra = 0;               % Crash of method
Stat = [];             % Statistic of bad parameters: iA,Nh,kL: err Ed>1
         
                              % LOOPS  
% Q1=0:2; Q2=0:2; Q3=0:2; Q4=0:2; Q5=0:2; Q6=-3:1:-1;      % for Nord6
% for a=Q1,for b=Q2,for c=Q3,for d=Q4,for e=Q5,for f=Q6,   Arg=[a b c d e f];

% fprintf('Arg = %d %d %d %d    ',Arg);
% for ia = 1:10,   Arg = [as(ia) bs(ia) cs(ia)];           % for N8_57  

%Arg=[1     1     0     2     1    -1];  Lmi=[0.2000 0.1000 0];
%LAr = length(Arg); 

for j = nODEs
   naODE = naODEs{j}; 
   %!load(naODE);
   %!for d=1:size(a_d,1)
   %!                      Arg  = a_d(d,1:LAr);  Lmi  = a_d(d,LAr+7:LAr+9); 
   %!EDmi = a_d(d,LAr+1);  EOmi_ = a_d(d,LAr+2);  Ratmi= a_d(d,LAr+3); 
   %!Nhmi = a_d(d,LAr+4);  kLmi = a_d(d,LAr+5);  NDmi = a_d(d,LAr+6); 
   %LAr = 2;  Ama = 5;                                       % Ama is Arg_max 
 %for iA=1:Ama^LAr-1, Arg = mod(fix(iA./Ama.^(0:LAr-1)),Ama);  % universal
 iA=0;for A1=[0.25 0.5 1 ], for A2=[0.1 0.25], Arg = [A1 A2 2/pi]; iA=iA+1;
 %iA=0;for A1= [0.5], for A2=[0.25], Arg = [A1 A2 2/pi]; iA=iA+1;        
   LAr = length(Arg);
   if ~any(Arg), continue,end  
   % Arg = Args{j};  
   % ticj  = tic;

   load( ['BestC_N_' naODE] );     % 'BestC','t0f','NDS','Ermin','sErmin'
   if j <= nHom, ds10L = @ds10Lh;
   else          ds10L = @ds10Li;  end

   t0f = t0fs(j,:); 
   t0  = t0f(1);      tf = t0f(2);  
   
   if expr, NA = [];
   else     NA = NDS; end 
   
   MD  = length(NDS);  LND = log2(NDS/100);
   MA  = length(NA);   LNA = log2( NA/100);  M1 = MA+1;
   ED  = nan(1,MD);              TD = ED;  nD = TD;
   EA  = nan(Lnmets, Ltols,M1);  TA = EA;  nA = TA;
   nB  = cell(1,M1);
   hsD = (tf-t0)./NDS;  
   hsA = (tf-t0)./NA;  
   naND= cellfun(@(x) num2str(x,'N_%d'),num2cell(NDS),'Uniform',false);
   
   if OUT
      hs_ = sprintf('%.4f     ',hsD); 
      fprintf('\n\n%30s func %s t0f=[%g %g]\n hsD = %s\n',' ',naODE,t0f,hs_);end
   
                              % ds10 method    
   EDmi = Inf;                         
   for Nh = 1:MD        %!Nh = Nhmi;
      h   = hsD(Nh);    
      t   = t0:h:tf;  n = length(t);
      y0I = esL(t(1:IL),Arg,j,2);      % (j,k,t,A) 
      ye  = esL(t,Arg,j,1);  % esL(t,A,j,k)
      BCm = BestC{m}{Nh};
      LL  = size(BCm,1);
      
      for kL = 1:min(50,LL)    %!  kL = kLmi;          %[1:8 10:3:min(25,LL)]
         Lmn = BCm(kL,1:m);
         op1  = DSet(h,[],ef, Lmn,[] );  
         tic; y = ds10L(@GMs,t,y0I,op1,Arg,j); T = toc;  
         tic; for l=1:qrun, y = ds10L(@GMs,t,y0I,op1,Arg,j);end; T=toc; 
         Ed = er(ye,y);
         if Ed < EDmi,  EDmi = Ed; Nhmi = Nh;  kLmi = kL; Lmi = Lmn;
         elseif Ed >1,  Stat = [Stat; [iA Nh kL]]; end,end % for kL
      
      if OUT, ED(Nh) = Ed;  TD(Nh) = T/qrun0;   nD(Nh) = n-IL;  end
         
      if pri_ds
         figure( 'Name',sprintf('%d %s %s',j,naODE, func2str(ds10L))); 
         p = Nh; if Nh > 6, p = Nh+6; end
         subplot(4,6,p); plot(t',esL(t,Arg,j,2));
         title(sprintf('%d %.3g',NDS(Nh),h))
         subplot(4,6,p+6); plot(t',y); end,end % for Nh
  
   esma = max(max(abs(esL(t,Arg,j,2))));
   
   if  OUT
      EDt = str2num(num2str(ED,'%.2g  '));
      TDt = str2num(num2str(TD,'%.2g  '));
      ETN = array2table([EDt; TDt; nD],'Var',naND,'Row',{'ES' 'TS' 'NS'})
   
      [mi, imi] = min(ED); 
      J   = 1:imi-dJ(j);      % exclude extra points from right hand
      LED = -log10(ED);       LnD = log10(nD); end
   
                      % All other methods
   for r = 1:M1        % At first: big time for big N  
      nB{r} = [];      % numbers of Bad meth - that did not find the solution 
      if r == M1
         t = [t0 tf];  
         if OUT, fprintf('\n\n%30s Non dense output\n',' '); end
      else
         h = hsA(r);
         t = (t0:h:tf)';
         fprintf('\n\n%30s N = %d   h = %.3g\n',' ',NA(r),h); end
     
      y0 = esL(t(1),Arg,j,2); 
      
                              % Time, Err, Nfs by other METHODs  
      EOmi = Inf;  
      for i = 1:Lnmets        % number of method
         nm   = nmets(i);     % count  of method 
         met_odex = strcmp(namets{i},'odex');
         met  = mets{nm};  
         setm = setms{nm};

         for s = 1:Ltols
            if met_odex
               if r < M1,  tol = tols(s)*1e5;
               else        tol = tols(s); end
            else           tol = tols(s); end
            
            op2 = setm( h0, tol,tol);
            warning('off','all');
            tic; [t1,y,nf] = met(@RHFL,t,y0,op2,Arg,j); T = toc;
            tic; for l=1:qrun, [t1,y,nf]= met(@RHFL,t,y0,op2,Arg,j); T=toc;end
            warning('on','all');

            if OUT 
               [Mr Mc] = find( y > 1.5*esma, 1 );
               if ~isempty(Mr),               Cra=1;     L=1:Mr+afteCr(i);end 
               if isnan(t1(end)), t1(end)=[]; Cra=Cra+2; L=1:length(t1);  end
               
               if Cra
                  figure('Name',sprintf('s=%d tol=%.3g %s',s,tol));
                  if Cra > 1
                     %if s==1, BadM = {'NaN' 'NaN'};
                     %else     BadM = {sprintf('E^*=10^{%.2g}  T^*=%.2g',...
                     %                log10(EA(i,s-1,r)), TA(i,s-1,r))...
                     %                sprintf('nf^*=%d',nA(i,s-1,r))};end
                     nB{r} = [nB{r} i]; 
                     break; end
              
                  plot(t1(L)',y(L,:)); title(sprintf('CRASH: %d-%s %d-%s',...
                  j,naODE,i,func2str(met))),end,end 
           
            ye = esL(t1,Arg,j,1);  
            EO = er(ye,y); 
            if EO < EOmi, EOmi = EO; mi = namets{nm};end
            
            if OUT
               EA(i,s,r) = EO;
               TA(i,s,r) = T/qrun0; 
               nA(i,s,r) = nf(1); end,end,end
   
   %fprintf('%.2g  %.2g  %s\n',EDmi,EOmi,mi);
   if EDmi<1e-5 && EDmi <= EOmi/EO_ED
      a_d  = [a_d; [ Arg EDmi EOmi EOmi/EDmi Nhmi kLmi NDS(Nhmi) Lmi ]]; 
      sArg = sprintf('%+.2g ', a_d(end,1:LAr));
      fprintf('%s   %.1g  %.1g  %.1f    %d  %d  %d  %+.2f %+.2f %+.2f\n',...
      sArg, a_d(end,LAr+1:end)); end

   if OUT
      EAt = str2num(num2str(EA(:,:,r),'%.2g  '));
      TAt = str2num(num2str(TA(:,:,r),'%.2g  '));
      ERR = array2table(EAt,'Var',natols,'Row',namets)
      TIM = array2table(TAt,'Var',natols,'Row',namets)
      NFU = array2table(nA(:,:,r),'Var',natols,'Row',namets),end,end
  
                      % OUTPUT  for j
                      % 1 Fig  NON DENSE  All methods        
   if OUT                   
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
      %if k == 1,  text( ps.XTick(4), ps.YTick(2), BadM{1} ); end
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
      legend(namets_,'Location','southwest');
      ylabel('-log_{10}E');xlabel('N'); end,end

   % TOC(j) = toc(ticj);
   
%   fprintf('%d) %g-%g  %g-%g  %g\n',d,Ed,EDmi,EOmi,EOmi_,EOmi/Ed); 
 end,end  % for iA
end  % for j

toc(ticAll) 
%fprintf('Total time = %s= %.3g < %.3g\n',sprintf('+ %g ',TOC),sum(TOC),tocAll)
a_d = sortrows( a_d, -(LAr+3)); % descend - по убыванию Rat
tip = 'Arg=[A1,...A_LAr], EDmi, EOmi, Rat=EOmi/EDmi, Nhmi, kLmi, NDmi, Lmn_mi';
save(naODE, 'a_d', 'LAr', 'tip');      

for d=1:size(a_d,1), fprintf('%.1g  %.1g  %.0f  %d %2d %d\n',a_d(d,4:9));end