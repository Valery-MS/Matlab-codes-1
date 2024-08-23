%  DSvsL demostration (ds vs other methods for Linear ODEs)
% SAVED 

function DSvsLdem_save_28_06_2019
ticAll = tic;

                           % SOURCE DATA for odex                    
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
mets    = cellfun(@str2func,namets_c,'Uniform',false); 
namets  = cellfun(@(x) strrep(x,'_c',''), namets_c,'Uniform',false);  
nametsD = ['ds10' namets];                             % D = ds10+Another 
namets_ = nametsD(1:end-1);                            % _ = D-BadMet

setms = {@(h,RT,AT,d) dopset('InitialSt',h,'RelT',RT,'AbsT',AT),... % 1 dop853
         @(h,RT,AT,d) odeset('InitialSt',h,'RelT',RT,'AbsT',AT),... % 2 ode113
         @(h,RT,AT,d) odeset('InitialSt',h,'RelT',RT,'AbsT',AT),... % 3 ode45 
         @(h,RT,AT,d) {h, RT, AT, ITOL,  WORK, IWORK} };            % 4 odex 

                           %  PRESET CONSTANTS                         
naODEs = { 'N3_43', 'N4_35', 'N8_57', 'N3_46', 'Nord6'};
nHom   = 3;                   % number of homogeneouse ODEs  
t0fs   = [ [1 21];   [2 22];  [0 20];    [1 21];  [0 21]];
er = @(ye,y) sqrt(sum((ye(:,1)-y(:,1)).^2)/length(y)); % |ApproxSol-ExactSol| 
NA  = [150 250 500 1000 2500 5000 7500 10000 12500 15000 20000];
%NDS = NA;    % NDS - from loading "BestC_N_naODE"

%*****************************************************************************
                       % SOURCE DATA                          
m     = 3;             % number of parameters of the difference scheme
h0    = 0.01;
nODEs = 3;%[1 3 4 5];  % equations
nmets = 1:4;           % methods
dBest = [1 NaN 8 9 6]; % p.47 OffBas - d according best Figs
BT    = [...       % p.48 - lengthes of Bad Tails: points with bad figs
         1 NaN 3 3 2]; 
nmBT  = 3;   % accords to ode45_c
%RT   = eps; AT = eps;
qrun  = 1;     %fprintf('%30s qrun = %d\n',' ',qrun);   % q-ty of method runs
qrun0 = max(1,qrun);
%*******************************************************************************

                       % Dependent CONSTANTS  
Lnmets = length(nmets);        % all methods : 1:4  
I      = 1:Lnmets;
%NDS   = [150 250 500 1000 2500 5000 7500 10000 12500 15000 20000];    
N0     = [1 1 1 1 1];    % какого номера по абсциссе строить график
basa   = [6:2:14 14.8 15.2 15.65]; Tail = 1;   %basa  = [15.65];
tols   = 10.^-basa;  Ltols = length(tols);
zn     = '^so*d'; 
natols = cellfun(@(x) strrep(num2str(x,'e_%.3g'),'.','_'),num2cell(basa),...
         'Uniform',false);    % replace . -> _
TOC    = zeros(1,length(nODEs));
iBT    = find(nmets==nmBT);

                       % FLAGs
dJ  = [ 0 0 0 0 0 ] ;  % [2 2 1];       % how many right extra points is in ED
ef  = 1;               % Effective way of computing Gen.Matrix
nIL = 1;                                % 1/2 if takeoff by es/(dop853|ode113)
ILs = [8 1];           IL = ILs(nIL);   % Initial values Length of y 
pri_ds = false;        % print solution by ds10
afteCr = [ 3 0 3 0 ];  % number of out-points after Crash
OUT = true;  CrashOut = false;

                              % LOOPS  
% Q1=0:2; Q2=0:2; Q3=0:2; Q4=0:2; Q5=0:2; Q6=-3:1:-1;      % for Nord6
% for a=Q1,for b=Q2,for c=Q3,for d=Q4,for e=Q5,for f=Q6,   Arg=[a b c d e f];
% LAr = 3;  Ama = 5;                                       % Ama is Arg_max 
% for j=1:Ama^LAr, Arg = mod(fix(j./Ama.^(0:LAr-1)),Ama);  % universal
% fprintf('Arg = %d %d %d\n',Arg);
% for ia = 1:10,   Arg = [as(ia) bs(ia) cs(ia)];           % for N8_57  
dmax = 3;  %  number of Difference Schemes demonstrating ds-method

fprintf('%5s %10s %7s %4s %4s %3s %3s %12s\n',...
        'Arg','EDmi','EOmi','Rat','Nh','kLm','ND','Lmn_mi');
    
for j = nODEs 
   % ticj  = tic;
   naODE = naODEs{j};
   load(['BestC_N_' naODE]);   % Load NDS
   load(naODE);    %  for N8_57/Nord5 Arg = [a b c] / [a b c d e f]
   % a_d = Arg EDmi EOmi Rat=EOmi/EDmi Nhmi kLmi NDmi Lmn_mi
   
   if j <= nHom, ds10L = @ds10Lh;
   else          ds10L = @ds10Ln;  end
   
   t0f = t0fs(j,:); 
   t0  = t0f(1);      tf = t0f(2);   
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
      fprintf('\n\n%30s func %s  t0f=[%g %g]  qrun=%d\n hsD = %s\n',...
             ' ',naODE,t0f,qrun, hs_);end

   for d = dBest(j) %1 :min(10,size(a_d,1))    
     Arg  = a_d(d,1:LAr);  Lmi  = a_d(d,LAr+7:LAr+9); 
     EDmi = a_d(d,LAr+1);  EOmi = a_d(d,LAr+2);  Ratmi= a_d(d,LAr+3); 
     Nhmi = a_d(d,LAr+4);  kLmi = a_d(d,LAr+5);  NDmi = a_d(d,LAr+6);
     
     ODEdA = sprintf('%s  d=%d Arg=%s',naODE,d,sprintf('%.2f  ',Arg)); 
     fprintf('\n %s\n ED=%.2g EO=%.2g Rat=%.2g N=%d Lmn=%.1f %.1f %.1f\n',...
     ODEdA,EDmi,EOmi,Ratmi,NDmi,Lmi);
   
     Lmn = BestC{m}{Nhmi}(kLmi,1:m);
     if any(Lmn ~= Lmi)
        fprintf('Ahtung! Lmi=%.1f %.1f %.1f Lmn=%.1f %.1f %.1f\n',Lmi,Lmn);end   
   
                              % ds10 method    
     EDmi = Inf;                         
     for Nh = MD:-1:1
        h   = hsD(Nh);    
        t   = t0:h:tf;  n = length(t);
        y0I = esL(t(1:IL),Arg,j,2);      % (j,k,t,A)          
        op  = DSet(h,[],ef, Lmi,[] );  
        tic; y = ds10L(@GMs,t,y0I,op,Arg,j); T = toc;  
        tic; for l=1:qrun, y = ds10L(@GMs,t,y0I,op,Arg,j);end; T=toc; 
        ye = esL(t,Arg,j,1);  Ed = er(ye,y);
      
        if OUT, ED(Nh) = Ed;  TD(Nh) = T/qrun0;   nD(Nh) = n-IL;  end
          
        if pri_ds
           figure( 'Name',sprintf('%d %s %s',j,naODE, func2str(ds10L))); 
           p = Nh; if Nh > 6, p = Nh+6; end
           subplot(4,6,p); plot(t',esL(t,Arg,j,2));
           title(sprintf('%d %.3g',NDS(Nh),h))
           subplot(4,6,p+6); plot(t',y); end,end
 
     esma = max(max(abs(esL(t,Arg,j,2))));
   
     if OUT
        EDt = str2num(num2str(ED,'%.2g  '));
        TDt = str2num(num2str(TD,'%.2g  '));
        ETN = array2table([EDt; TDt; nD],'Var',naND,'Row',{'ES' 'TS' 'NS'})
   
        [mi, imi] = min(ED); 
        iEL1 = find(ED<1,1);
        J   = iEL1:imi;        % exclude points: ED > 1
        LED = -log10(ED);       LnD = log10(nD); end
   
                              % All other methods
   % r=12   - indexes Non dense output
   % r=11...1 indexes N=1000*[20, 15, 12.5, 10...0.25, 0.15] 
   for r = M1:-1:1 % indexes N=20000...150,  At first: big time for big N  
      nB{r} = []; % numbers of Bad meth- that did not find the solution 
      if r == M1
         t = [t0 tf];  
         if OUT, fprintf('\n\n%30s Non dense output\n',' '); end
      else
         h = hsA(r);
         t = (t0:h:tf)';
         fprintf('\n\n%30s N = %d   h = %.3g\n',' ',NA(r),h); end
     
      y0 = esL(t(1),Arg,j,2); 
      
                              % Time, Err, Nfs by Another METHODs    
      for i = 1:Lnmets        % indexes method number i=1:4
         namet = namets{i};   met_odex = strcmp(namet,'odex');
         nm   = nmets(i);     % count  of method 
         if nm == nmBT, iBT = i; end
         met  = mets{nm};  
         setm = setms{nm};
         Cra = 0;
         EOmi  = Inf;
         
         for s = 1:Ltols  % s indexes tolerances 10^(-[6:2:14 14.8 15.2 15.65])
            if met_odex
               if r < M1,  tol = tols(s)*1e5;
               else        tol = tols(s); end
            else           tol = tols(s); end
            
            op = setm( h0, tol,tol,[]);
            warning('off','all');
            tic; [t1,y,nf] = met(@RHFL,t,y0,op,Arg,j); T = toc;
            tic; for l=1:qrun, [t1,y,nf]= met(@RHFL,t,y0,op,Arg,j);end; T=toc;
            warning('on','all');

            if OUT 
               [Mr Mc] = find( y > 1.5*esma, 1 );
               if ~isempty(Mr) || met_odex && nf(2) >= IWORK(1)
                  if CrashOut 
                     figure('Name',sprintf('s=%d tol=%.3g %s',s,tol));
                     plot(t1',y);
                     title(sprintf('Crash %s %s r=%d h=%g',naODE,namet,r,h));end

                  if met_odex || s == Ltols, nB{r}= [nB{r} i];  break;
                  else continue; end,end,end 
           
            ye = esL(t1,Arg,j,1);  
            EO = er(ye,y); 
            if EO < EOmi, EOmi = EO; end
            
            if OUT
               EA(i,s,r) = EO;
               TA(i,s,r) = T/qrun0; 
               nA(i,s,r) = nf(1); end,end,end
  
   if OUT
      EAt = str2num(num2str(EA(:,:,r),'%.2g  '));
      TAt = str2num(num2str(TA(:,:,r),'%.2g  '));
      ERR = array2table(EAt,'Var',natols,'Row',namets)
      TIM = array2table(TAt,'Var',natols,'Row',namets)
      NFU = array2table(nA(:,:,r),'Var',natols,'Row',namets),end,end
  
                      % OUTPUT  for j
                      % 1 Fig  NON DENSE  All methods        
   if OUT  
    if ~isempty(iBT),    EA(iBT,end-BT(j)+1:end,M1) = NaN; end  
            % delete Bad Tail for nmBT-meth with serial number iBT
    LEA = -log10(EA);  
    LnA =  log10(nA);   
    figure('Name',[ODEdA 'ND_Tn(E)']); 
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
   
                      % 2 Fig  DENSE  exclude Bad Methods                     
    figure('Name',[naODE ', D_Tn(E)']); 
    k = 0;
    for r = [1 MA]
      I_ = I;   I_(nB{r}) = [];
      k = k+1;
      ps = subplot(3,1,k);      %   Time(E)
      p3 = plot(LED(J),TD(J),LEA(I_,:,r)',TA(I_,:,r)'); % exclude increasing ED
      %if k == 1,  text( ps.XTick(4), ps.YTick(2), BadM{1} ); end
      for p=1:length(I_)+1,  p3(p).Marker = zn(p); p3(p).LineWidth = 1;end 
      title( sprintf('N=%d  h=%.2g',NA(r),hsA(r)) );
      ylabel('Time');   % xlabel('-log_{10}E');
      legend(['ds10' namets(I_)],'Location','northwest');  end 
      
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
    K = LND > LNA(N0(j))*(1-eps); %in order for the definition areas to coincide
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
      ylabel('-log_{10}E');xlabel('N'); end,end

   end
   end  % for j 

toc(ticAll) 
1;