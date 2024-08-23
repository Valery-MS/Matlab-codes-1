%                        SAVE
% % ds vs other methods for Linear ODEs

function DSvsL_SAVE
                           %  PRESET CONSTANTS                         
naODEs = { 'N3_43',  'N4_35',  'N8_57',  'N3_46',   'Nord6'};
nHom = 3;                   % number of homogeneouse ODEs  

t0fs  = [ [1 21];   [1 21];  [0 20];    [1 21];  [0 20]];
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

setms = {@(h,RT,AT,d) dopset('InitialSt',h,'RelT',RT,'AbsT',AT),... % 1 dop853
         @(h,RT,AT,d) odeset('InitialSt',h,'RelT',RT,'AbsT',AT),... % 2 ode113
         @(h,RT,AT,d) odeset('InitialSt',h,'RelT',RT,'AbsT',AT),... % 3 ode45 
         @(h,RT,AT,d) {h, RT, AT, ITOL,  WORK, IWORK} };            % 4 odex 

                      % SOURCE DATA for ds10                               
m    = 3;             % number of parameters of the difference scheme
kLmn = 1;             % number of choice of coefs La,Mu,Nu

                        % SOURCE DATA
h0 = 0.01;
nODEs = [     5 ];        % equations
nmets = 1:4;                  % methods
%RT = eps; AT = eps;
qrun  = 5;   fprintf('%30s qrun = %d\n',' ',qrun);   % q-ty of method runs
                               % Dependent CONSTANTS  
Lnmets = length(nmets);        % all methods : 1:4  
NDs= {[150 250 500 1000 2500 5000 7500 10000 12500 15000 20000];... % N3_43
      [150 250 500 1000 2500 5000 7500 10000 12500 15000 20000];... % N4_35
      [150 250 500 1000 2500 5000 7500 10000 12500 15000 20000];... % N8_57
      [150 250 500 1000 2500 5000 7500 10000 12500 15000 20000];... % N3_46
      [150 250 500 1000 2500 5000 7500 10000 12500 15000 20000]};   % Nord6
         
NAs= {[500 1000 2500 5000 7500 10000 12500 15000 20000];... % N3_43
      [500 1000 2500 5000 7500 10000 12500 15000 20000];... % N4_35
      [500 1000 2500 5000 7500 10000 12500 15000 20000];... % N8_57
      [500 1000 2500 5000 7500 10000 12500 15000 20000];... % N3_46
      [500 1000 2500 5000 7500 10000 12500 15000 20000]};   % Nord6
      
N0 = [1 1 1 1 1];    % какого номера по абсциссе строить график
basa = [6:2:14 14.8 15.2 15.55]; Tail = 1;
tols = 10.^-basa;  Ltols = length(tols);
zn = '^so*d'; 
natols = cellfun(@(x) strrep(num2str(x,'e_%.3g'),'.','_'),num2cell(basa),...
         'Uniform',false);    % replace . -> _
I  = 1:Lnmets;
dJ = [ 0 0 0 0 0 ] ;  %[2 2 1];     % how many right extra points is in ED
TOC = zeros(1,length(nODEs));
IL = 8;
ef = 1;      % Effective way of computing Gen.Matrix

                              % LOOPS                            
for j = nODEs
   ticj  = tic;
   naODE = naODEs{j};
   load( ['BestC_N_' naODE] );     % 'BestC','t0f','NDS','Ermin','sErmin'
   ND = NDs{j}; 
   if j <= nHom, ds10L = @ds10Lh;
   else          ds10L = @ds10Ln;  end

   t0f = t0fs(j,:); 
   Arg = Args{j};  
   t0 = t0f(1);      tf = t0f(2);  
   NA = NAs{j};
   MD = length(ND);  LND = log2(ND/100);
   MA = length(NA);  LNA = log2(NA/100);  M1 = MA+1;
   ED = nan(1,MD);              TD = ED;  nD = TD;
   EA = nan(Lnmets, Ltols,M1);  TA = EA;  nA = TA;
   nB = cell(1,M1);
   hsD = (tf-t0)./ND;  
   hsA = (tf-t0)./NA; 
   hs_ = sprintf('%.4f     ',hsD); 
   naND = cellfun(@(x) num2str(x,'N_%d'),num2cell(ND),'Uniform',false);
   
   fprintf('\n%30s ODEfunc = %s   t0f = [%g  %g]\n',' ',naODE,t0f);
   fprintf('hsD = %3s %s\n',' ',hs_);
   Nh = ND(MD);
   
                              % ds10 method
   for r = MD:-1:1
      h   = hsD(r);      Nr = ND(r);
      t   = (t0:h:tf);   n = length(t);
      y0I = esL(t(1:IL),Arg,j,2);
      Nh_ = find(NDS == Nr);          % Nh = r, h's number in BestC
      if ~isempty(Nh_), Nh = Nh_; end
      Lmn = BestC{m}{Nh}(kLmn,1:m);
      op  = DSet(h,[],ef, Lmn,[] );  
      y   = ds10L(@GMs,t,y0I,op,Arg,j);   
      tic; for l=1:qrun, y = ds10L(@GMs,t,y0I,op,Arg,j);end; T=toc; 
      ye= esL(t',Arg,j,1);  TD(r)= T/qrun;   ED(r)= er(ye,y);  nD(r)= n-IL; end
  
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
     
      y0 = esL(t(1),Arg,j,2); 
      
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
            warning('off','all');
            [t1,y,nf] = met(@RHFL,t,y0,op,Arg,j);
            tic; for l=1:qrun, [t1,y,nf]= met(@RHFL,t,y0,op,Arg,j); T=toc;end
            warning('on','all');
            
            if isnumeric(t1)
               ye=esL(t1,Arg,j,1); 
               EA(i,s,r)=er(ye,y);  TA(i,s,r)=T/qrun;  nA(i,s,r)=nf(1); 
            else 
               if s==1, BadM = {'NaN' 'NaN'};
               else     BadM = {sprintf('E^*=10^{%.2g}  T^*=%.2g',...
                               log10(EA(i,s-1,r)), TA(i,s-1,r))...
                               sprintf('nf^*=%d',nA(i,s-1,r))};end
               nB{r} = [nB{r} i]; 
               break; end,end,end

      EAt = str2num(num2str(EA(:,:,r),'%.2g  '));
      TAt = str2num(num2str(TA(:,:,r),'%.2g  '));
      ERR = array2table(EAt,'Var',natols,'Row',namets)
      TIM = array2table(TAt,'Var',natols,'Row',namets)
      NFU = array2table(nA(:,:,r),'Var',natols,'Row',namets),   end
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
ED
