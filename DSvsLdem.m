%  DSvsL demostration (ds vs other methods for Linear ODEs)

function DSvsLdem
ticAll = tic;
                           % SOURCE DATA  
                           % PARAMETERs for odex                    
ITOL     = 0;
%IOUT    = 3;           my change: IOUT calcs automatically by dim(X)
WORK     = zeros(1,9);     
IWORK    = zeros(1,71); % 2*KM+21+NRDENS
IWORK(1) = 15000;       % NMAX, default = 1e4
IWORK(7) = 4;           % DEGREE OF INTERPOL F-LA MUDIF=4 (DEFAULT) LIE in [1,6]

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

                         % PRESET CONSTANTS  
        % PARAMETRER SETting FUNCTION for methods: dop853, ode113, ode45, odex
setms = {@(h,RT,AT,d) dopset('InitialSt',h,'RelT',RT,'AbsT',AT),... % 1 dop853
         @(h,RT,AT,d) odeset('InitialSt',h,'RelT',RT,'AbsT',AT),... % 2 ode113
         @(h,RT,AT,d) odeset('InitialSt',h,'RelT',RT,'AbsT',AT),... % 3 ode45 
         @(h,RT,AT,d) {h, RT, AT, ITOL,  WORK, IWORK} };            % 4 odex

                           % FLAGs
ef  = 1; % Effective way of computing Generation Matrix
nIL = 1; % 1, takeoff by exact solution es, h=const
         % 2, takeoff by dop853 or ode113,  h=const
         % 3, takeoff by exact solution es, h=variable
         % 4, takeoff by dop853 or ode113,  h=variable
ILs      = [8 1 16 1]; % Array of Initial values Length of y for any nIL
IL       = ILs(nIL);   % Initial values Length of y 
pri_ds   = false;      % print solution by ds10
OUT      = true; 
CrashOut = false;

                        %  SELDOM                  
m      = 3;             % number of parameters of the difference scheme
h0     = 0.01;
dBest  = [1 NaN 8 9 6]; % p.47 OffBas - d according best Figs for 4 ODEs
%RT    = eps; AT = eps;
qrun   = 1;     %fprintf('%30s qrun = %d\n',' ',qrun);   % q-ty of method runs
zn     = '^so*d'; 
N0     = [1 1 1 1 1];    % какого номера по абсциссе строить график

                           % TOLERANCES
basT = [6:2:14 14.8 15.2 15.65]; % Base of Tolerances
nTT  = 1;                        % numb of elms of TolsTail
Tols = 10.^-basT;                % ...  Tail = basT(end-Tail_Tol+1:end)
nTol = length(Tols);

                          % METHODS(OTHER, not DS)
namets_c = {'dop853' 'ode113_c' 'ode45_c' 'odex'};     
mets     = cellfun(@str2func,namets_c,'Uniform',false); 
namets   = cellfun(@(x) strrep(x,'_c',''), namets_c,'Uniform',false);  
nametsD  = ['ds10' namets];                             % D = ds10+other 
namets_  = nametsD(1:end-1);                            % _ = D-BadMet

                           % ODEs             
naODEs = { 'N3_43', 'N4_35', 'N8_57', 'N3_46', 'Nord6'};
nHom   = 3;                   % number of homogeneouse ODEs  
t0fs   = [ [1 21];   [2 22];  [0 20];    [1 21];  [0 21]];
er     = @(ye,y) sqrt(sum((ye(:,1)-y(:,1)).^2)/length(y)); % |Approx-Exact|
NA     = [150 250 500 1000 2500 5000 7500 10000 12500 15000 20000];
%NDS = NA;    % NDS - from loading "BestC_N_naODE"

%*****************************************************************************
                        % OFTEN 
N_ODEs = 1;%[1 3 4 5];  % Numbers of equations
N_Mets = 1:4;           % Numbers of methods       
%******************************************************************************

                       % Dependent CONSTANTS  
qrun0  = max(1,qrun);                       
nODE   = length(N_ODEs); % quantity, number of ODEs
nMet   = length(N_Mets); % quantity, number of methods 
I      = 1:nMet;  
natols = cellfun(@(x) strrep(num2str(x,'e_%.3g'),'.','_'),num2cell(basT),...
         'Uniform',false);    % replace . -> _


                              % LOOPS  
% NDS = [150 250 500 1000 2500 5000 7500 10000 12500 15000 20000];  
% Q1=0:2; Q2=0:2; Q3=0:2; Q4=0:2; Q5=0:2; Q6=-3:1:-1;      % for Nord6
% for a=Q1,for b=Q2,for c=Q3,for d=Q4,for e=Q5,for f=Q6,   Arg=[a b c d e f];
% LAr = 3;  Ama = 5;                                       % Ama is Arg_max 
% for j=1:Ama^LAr, Arg = mod(fix(j./Ama.^(0:LAr-1)),Ama);  % universal
% fprintf('Arg = %d %d %d\n',Arg);
% for ia = 1:10,   Arg = [as(ia) bs(ia) cs(ia)];           % for N8_57  


fprintf('%5s %10s %7s %4s %4s %3s %3s %12s\n',...
        'Arg','EDmi','EOmi','Rat','Nh','kLm','ND','Lmn_mi');
    
for j = N_ODEs 
   % ticj  = tic;
   naODE = naODEs{j};
   load(['BestC_N_' naODE]);   % Load NDS
   load(naODE);    %  for N8_57/Nord5 Arg = [a b c] / [a b c d e f]
   % a_d = Arg EDmi EOmi Rat=EOmi/EDmi Nhmi kLmi NDmi Lmn_mi
   
   if j <= nHom
       if nIL <=2,  ds10L = @ds10Lh; else ds10L = @ds10Lhv; end
   elseif nIL <=2,  ds10L = @ds10Li; else ds10L = @ds10Liv; end
  
   t0f = t0fs(j,:); 
   t0  = t0f(1);      tf = t0f(2);   
   MD  = length(NDS);  LND = log2(NDS/100);
   MA  = length(NA);   LNA = log2( NA/100);  M1 = MA+1;
   ED  = nan(1,MD);              TD = ED;  nD = TD;
   EA  = nan(nMet, nTol,M1);  TA = EA;  nA = TA;
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
   % iN=12   - indexes Non dense output
   % iN=11...1 indexes N=1000*[20, 15, 12.5, 10...0.25, 0.15] 
   for iN = M1:-1:1 % indexes N=20000...150,  At first: big time for big N  
      nB{iN} = []; % numbers of Bad meth- that did not find the solution 
      if iN == M1
         t = [t0 tf];  
         if OUT, fprintf('\n\n%30s Non dense output\n',' '); end
      else
         h = hsA(iN);
         t = (t0:h:tf)';
         fprintf('\n\n%30s N = %d   h = %.3g\n',' ',NA(iN),h); end
     
      y0 = esL(t(1),Arg,j,2); 
      
                              % Time, Err, Nfs by Another METHODs    
      for im = 1:nMet        % indexes method number im=1:4
         namet = namets{im};   met_odex = strcmp(namet,'odex');
         nm    = N_Mets(im);     % count  of method 
         met   = mets{nm};  
         setm  = setms{nm};
         EOmi  = Inf;
         
         for it = 1:nTol  % it indexes tolerances 10^(-[6:2:14 14.8 15.2 15.65])
            if met_odex
               if iN < M1, tol = Tols(it)*1e5;
               else        tol = Tols(it); end
            else           tol = Tols(it); end
            
            op = setm( h0, tol,tol,[]);
            warning('off','all');
            tic; [t1,y,nf] = met(@RHFL,t,y0,op,Arg,j); T = toc;
            tic; for l=1:qrun, [t1,y,nf]= met(@RHFL,t,y0,op,Arg,j);end; T=toc;
            warning('on','all');

            if OUT 
               [Mr Mc] = find( y > 1.5*esma, 1 );
               if ~isempty(Mr) || met_odex && nf(2) >= IWORK(1)
                  if CrashOut 
                     figure('Name',sprintf('s=%d tol=%.3g %s',it,tol));
                     plot(t1',y);  title( ...
                     sprintf('Crash %s %s iN=%d h=%g',naODE,namet,iN,h));end

                  if met_odex || it == nTol, nB{iN}= [nB{iN} im];  break;
                  else continue; end,end,end 
           
            ye = esL(t1,Arg,j,1);  
            EO = er(ye,y); 
            if EO < EOmi, EOmi = EO; end
            
            if OUT
               EA(im,it,iN) = EO;
               TA(im,it,iN) = T/qrun0; 
               nA(im,it,iN) = nf(1); end,end,end
  
   if OUT
      EAt = str2num(num2str(EA(:,:,iN),'%.2g  '));
      TAt = str2num(num2str(TA(:,:,iN),'%.2g  '));
      ERR = array2table(EAt,'Var',natols,'Row',namets)
      TIM = array2table(TAt,'Var',natols,'Row',namets)
      NFU = array2table(nA(:,:,iN),'Var',natols,'Row',namets),end,end

   
                      % OUTPUT  for j
                      % 1) NON DENSE  All methods        
   if OUT  
    LEA   = -log10(EA);
    LEATT = LEA(:,nTol-nTT+1:nTol,:);  % use only in 3), but LEA changes now
    TATT  =  TA(:,nTol-nTT+1:nTol,:);  % calcs similarly
    [mi, L]  = min(EA(:,:,M1),[],2); 
    for i = 1:nMet, LEA( i, L(i)+1:end, :) = NaN; end  % exclude bad pts: Er>min
    
    LnA =  log10(nA);   
    figure('Name',[ODEdA 'ND_Tn(E)']); 
    subplot(2,1,1);          %  T(E)
    p1 = plot(LED(J),TD(J), LEA(:,:,M1)',TA(:,:,M1)');
    %title([naODE ', non dense output']);
    for p=1:nMet+1, p1(p).Marker = zn(p); p1(p).LineWidth = 1;end 
    ylabel('Time');   % xlabel('-log_{10}E');
    legend(nametsD,'Location','northwest'); 
   
    subplot(2,1,2);           % nf(E)
    p2 = plot(LED(J),LnD(J), LEA(:,:,M1)',LnA(:,:,M1)'); 
    for p=1:nMet+1, p2(p).Marker = zn(p); p2(p).LineWidth = 1;end 
    ylabel('log_{10}nf');  xlabel('-log_{10}E'); 
    legend(nametsD,'Location','northwest'); 
   
                      % 2) DENSE  exclude Bad Methods (odex)
    figure('Name',[naODE ', D_Tn(E)']); 
    k = 0;
    for iN = [1 MA]
      I_ = I;   I_(nB{iN}) = [];
      k = k+1;
      ps = subplot(3,1,k);      %   Time(E)
      p3 = plot(LED(J),TD(J),LEA(I_,:,iN)',TA(I_,:,iN)'); % excl increasing ED
      %if k == 1,  text( ps.XTick(4), ps.YTick(2), BadM{1} ); end
      for p=1:length(I_)+1,  p3(p).Marker = zn(p); p3(p).LineWidth = 1;end 
      title( sprintf('N=%d  h=%.2g',NA(iN),hsA(iN)) );
      ylabel('Time');   % xlabel('-log_{10}E');
      legend(['ds10' namets(I_)],'Location','northwest');  end 
      
    ps = subplot(3,1,3);         %   nf(E) 
    iN = 1;
    p4 = plot(LED(J),LnD(J), LEA(:,:,iN)',LnA(:,:,iN)'); 
    %text( ps.XTick(4), ps.YTick(3), BadM{2} ); 
    for p=1:nMet+1, p4(p).Marker = zn(p); p4(p).LineWidth = 1;end 
    title( sprintf('N=%d  h=%.2g',NA(iN),hsA(iN)) ); 
    ylabel('log_{10}nf');  xlabel('-log_{10}E'); 
    legend(nametsD,'Location','northwest');  

                       % 3) DENSE High Tolerances 
    pTA  = permute( TATT(:,:,1:MA), [1 3 2]);
    pLEA = permute(LEATT(:,:,1:MA), [1 3 2]);
    K    = LND > LNA(N0(j))*(1-eps);   % for the definition areas to coincide
    K(1:iEL1-1) = false;               % exclude N, where Err>1                 
    scalog      = floor(LNA(1)) : ceil(LND(end));
    scale       = num2cell(2.^scalog*100);
    
    for iT = 1:nTT % Ск высоких точностей брать от конца базы...15 15.5 16
      figure('Name',sprintf('D_TE(N) %d) %s tol=%.2g',iT,naODE,Tols(end-iT+1)));
      ps = subplot(2,1,1);          %   Time(N)
      p5 = plot(LND(K),TD(K), LNA,pTA(I_,:,iT) );  
      ps.XTickLabel = scale;
      for p=1:nMet, p5(p).Marker = zn(p); p5(p).LineWidth = 1;end 
      legend(namets_,'Location','northwest');
      ylabel('Time'); 
  
      ps = subplot(2,1,2);          %   E(N)
      p6 = plot(LND(K),LED(K), LNA,pLEA(I_,:,iT) );
      ps.XTickLabel = scale;
      for p=1:nMet, p6(p).Marker = zn(p); p6(p).LineWidth = 1;end
      %title('Maximal accuracy');
      legend(namets_,'Location','southeast');
      ylabel('-log_{10}E');xlabel('N'); end,end

   end
   end  % for j 

toc(ticAll) 
1;
%ERR 
%dop853 2.2e-06 2.1e-08   2e-10    2e-12  3.7e-14  3.8e-14  7.3e-14  5.3e-14
%ode113 2.5e-05 9.1e-08   2e-11  5.7e-12  2.3e-12  2.5e-12  9.4e-13    1e-12
%ode45  4.3e-06 4.2e-08 4.1e-10  4.1e-12  1.2e-13    2e-13  8.8e-14  3.7e-13