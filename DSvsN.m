% ds vs other methods for Non linear ODEs
% % @rkf78 is bad & deleted

function DSvsN
ticAll = tic;
                           %  PRESET CONSTANTS
naODEs = { 'n66',   'H146',   'fant' };
t0fs   = [ [1 21]; [0 4];   [0 20] ];
Args   = cell(size(naODEs));   
er = @(ye,y) sqrt(sum((ye(:,1)-y(:,1)).^2)/length(y)); % |ApproxSol-ExactSol| 
                            
                     % SOURCE DATA for ODEX                    
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
                            
namets_c = {'dop853' 'ode113_c' 'ode45_c' 'ODEX'};     % Another methods
mets    = cellfun(@str2func,namets_c,'Uniform',false); 
namets  = cellfun(@(x) strrep(x,'_c',''), namets_c,'Uniform',false);  
nametsD = ['ds10' namets];                             % D = ds10+Another 
namets_ = nametsD(1:end-1);                            % _ = D-BadMet

%setms = {@(h,RT,AT,d) dopset('InitialSt',h,'RelT',RT,'AbsT',AT),... % 1 dop853
%         @(h,RT,AT,d) odeset('InitialSt',h,'RelT',RT,'AbsT',AT),... % 2 ode113
%         @(h,RT,AT,d) odeset('InitialSt',h,'RelT',RT,'AbsT',AT),... % 3 ode45 
%         @(h,RT,AT,d) {h, RT, AT, ITOL,  WORK, IWORK} };            % 4 ODEX 
setms = {@(h,R,A) dopset('InitialSt',h,'RelT',R,'AbsT',A),...    % 1 dop853
         @(h,R,A) odeset('InitialSt',h,'RelT',R,'AbsT',A),...    % 2 ode113
         @(h,R,A) odeset('InitialSt',h,'RelT',R,'AbsT',A),...    % 3 ode45 
         @(h,R,A) struct('InitialStep',h,'RelTol',R,'AbsTol',A,...
                    'ITOL',ITOL, 'WORK',WORK, 'IWORK',IWORK)};       % 4 ODEX 
     
                      % SOURCE DATA for ds10
zam  = [2 1 2];  % zam(1:2) iteration changes; zam(3)=1/2 if dop853/ode113
m    = 3;        % number of parameters of the difference scheme
kLmn = 1;        % number of Lmnwith min Err

                        % SOURCE DATA
h0    = 0.01;
nODEs = [1 2 3];        % equations
nmets = [1 2 3 4];            %ds10 ������ �������
RT    = eps; AT = eps;
qrun  = 1;   fprintf('%30s qrun = %d',' ',qrun);   % q-ty of method runs

                               % Dependent CONSTANTS  
Lnmets = length(nmets);        % all methods : 1:4  
%NDs={[100 250 500 1000 2500 5000 7500 10000 12500 15000 20000]; ... % n66
%     [125 250 500 1000 1200 1250 1500 2000 2500 3000 4000 8000 12000];...%H146
%     [150 250 500 1000 2500 5000 7500 10000 12500 15000 20000]};    % fant
         
NAs= {[500 1000 2500 5000 7500 10000 12500 15000 20000]; ...  % n66
      [250 500 1000 1500 2000 3000 4000  8000  12000 ];  ...  % H146
      [500 1000 2500 5000 7500 10000 12500 15000 20000]};     % fant
      
N0   = [1 1 1];    % ������ ������ �� �������� ������� ������
basa = [6:2:14 14.8 15.2 15.55]; Tail = 1;
tols = 10.^-basa;  Ltols = length(tols);
zn   = '^so*d'; 
natols = cellfun(@(x) strrep(num2str(x,'e_%.3g'),'.','_'),num2cell(basa),...
         'Uniform',false);    % replace . -> _
I   = 1:Lnmets;
TOC = zeros(1,length(nODEs));
nIL = 1;                                 % 1/2 if takeoff by es/(dop853|ode113)
dJs = [2 2 1; 0 0 0];  dJ = dJs(nIL,:);  % how many right extra points is in ED
ILs = [16 1];          IL = ILs(nIL);    % Initial values Length of y 
NDS1 = [1 1 1]; % � ������ N ��������-c �������� ���������� ���. �������� 1e-8
% j=1(n66)  - � NDS(1)=100, �� �������� 1e-8
% j=2(H146) - c NDS(2)=250, �� ���� 1e-9, � ��� NDS(1)=125 ���� 1e-5
% j=3(fant) - c NDS(2)=250, �� ���� 1e-8, � ��� NDS(1)=150 ���� 1e-5

                          % LOOPS                            
for j = nODEs
   ticj = tic;
   naODE = naODEs{j}; 
   load( ['BestC_N_' naODE] );     % 'BestC','t0f','NDS','Emin','sEmin'
   % BestC{m}{Nh}-qLMN*5-matrix = [Lmn1 E1 nf1; Lmn2 E2 nf2; ... ],
   %         ranged on increasing of error
   if NDS1(j)>1, NDS(1:NDS1(j)-1) = []; end % ������� �����(�����) N, ��� 
   % ������� ��� h ����� � => �������� ���� => ������� ����������
   Arg = Args{j};                % F's argument     
   t0f = t0fs(j,:); 
   t0  = t0f(1);      tf = t0f(2);  
   NA  = NAs{j};
   MD  = length(NDS);  LND = log2(NDS/100);
   MA  = length(NA);   LNA = log2( NA/100);  M1 = MA+1;
   ED  = nan(1,MD);              TD = ED;  nD = TD;
   EA  = nan(Lnmets, Ltols,M1);  TA = EA;  nA = TA;
   nB  = cell(1,M1);
   hsD = (tf-t0)./NDS;  
   hsA = (tf-t0)./NA;  
   naND= cellfun(@(x) num2str(x,'N_%d'),num2cell(NDS),'Uniform',false);
   
   hs_ = sprintf('%.4f     ',hsD); 
   fprintf('\n\n%30s ODEfunc = %s t0f=[%g %g]\n hsD = %s\n',' ',naODE,t0f,hs_);

                              % ds10 method F,t,y0,op,varargin
   for Nh = MD:-1:1  
      h   = hsD(Nh);    
      t   = t0:h:tf; 
      y0I = esN(t(1:IL),Arg,j,2);               
      Lmn = BestC{m}{Nh}(kLmn,1:m);
      op  = DSet( h, RT,AT,Lmn, zam);
      [to,y,nf] = ds10(@RHFN,t,y0I,op,Arg,j);
      tic; for l=1:qrun, [to,y,nf] = ds10(@RHFN,t,y0I,op,Arg,j); T=toc; end
      ye=esN(t(IL:end),Arg,j,1);  yI = y(IL:end,1);  ED(Nh)= er(ye,yI);
      TD(Nh)= T/qrun;  nD(Nh)= nf(1); end
  
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
         fprintf('\n\n%30s LOOSE output\n',' '); 
      else
         h = hsA(r);
         t = (t0+(IL-1)*h : h : tf)';  
         fprintf('\n\n%30s N = %d   h = %.3g\n',' ',NA(r),h); end
     
      y0I = esN(t(1),Arg,j,2); 
                              % Time, Err, Nfs by Another METHODs    
      for i = 1:Lnmets        % number of method
         nm   = nmets(i);     % count  of method 
         met  = mets{nm};  
         setm = setms{nm};
         
         for s = 1:Ltols
            if strcmp('ODEX',namets(i))
               if r < M1,  tol = tols(s)*1e3;
               else        tol = tols(s); end
            else           tol = tols(s); end
            
            op = setm( h0, tol,tol);
            T = 0;
            warning('off','all');
            [to,y,nf] = met(@RHFN,t,y0I,op,Arg,j);
            for l=1:qrun, tic; [to,y,nf]= met(@RHFN,t,y0I,op,Arg,j);T=T+toc;end
            warning('on','all');
            
            if isnumeric(nf) && nf(1) < 50000
               ye=esN(to,Arg,j,1); 
               EA(i,s,r)=er(ye,y); TA(i,s,r)=T/qrun; nA(i,s,r)=nf(1); 
            else
               EA(i,s:Ltols,r)=NaN; TA(i,s:Ltols,r)=NaN; nA(i,s:Ltols,r)=NaN;
               break; end,end,end

      EAt = str2num(num2str(EA(:,:,r),'%.3g  '));
      TAt = str2num(num2str(TA(:,:,r),'%.2g  '));
      NM  = namets(nmets);
      ERR = array2table(EAt,'Var',natols,'Row',NM)
      TIM = array2table(TAt,'Var',natols,'Row',NM)
      NFU = array2table(nA(:,:,r),'Var',natols,'Row',NM), end
  
                      % OUTPUT  for j (number of ODEs)
                      % 1 T(p) nf(p) LOOSE   
                      % here p = -lg(Err)
   LEA = -log10(EA);  
   LnA =  log10(nA);   
   figure('Name',[naODE ', T(p) nf(p) LOOSE']); 
   subplot(2,1,1);          %  T(p)
   p1 = plot(LED(J),TD(J), LEA(:,:,M1)',TA(:,:,M1)');
   %title([naODE ', non dense output']);
   for p=1:Lnmets+1, p1(p).Marker = zn(p); p1(p).LineWidth = 1;end 
   ylabel('Time');   % xlabel('-log_{10}E');
   legend(nametsD,'Location','northwest'); 
   
   subplot(2,1,2);           % nf(p)
   p2 = plot(LED(J),LnD(J), LEA(:,:,M1)',LnA(:,:,M1)'); 
   for p=1:Lnmets+1, p2(p).Marker = zn(p); p2(p).LineWidth = 1;end 
   ylabel('log_{10}nf');  xlabel('-log_{10}E'); 
   %legend(nametsD,'Location','northwest'); 
   
                      % 2 T(p) nf(p) DENSE                   
   figure('Name',[naODE ',T(p) nf(p) DENSE']); 
   k = 0;
   if j==2, I_= 1:Lnmets-1; Lnme = Lnmets;   namet = namets_;  %H146: excl ODEX
   else     I_= I;          Lnme = Lnmets+1; namet = nametsD; end
   
   for r = [1 MA]
      k  = k+1;
      ps = subplot(3,1,k);   %   Time(E)
      p3 = plot(LED(J),TD(J),LEA(I_,:,r)',TA(I_,:,r)'); % exclude increasing ED  
      if  j==2
         if k==1, ps.YLim = [0 0.6];
         else     ps.YLim = [0 1.3]; end,end
      %if k == 1,  text( ps.XTick(4), ps.YTick(2), BadM{1} ); end
      for p=1:Lnme, p3(p).Marker = zn(p); p3(p).LineWidth = 1;end 
      title( sprintf('N=%d  h=%.2g',NA(r),hsA(r)) );
      ylabel('Time');   % xlabel('-log_{10}E');
      if r==1, legend(namet,'Location','northwest');end,end 
      
   ps = subplot(3,1,3);         %   nf(E) 
   r  = MA;
   p4 = plot(LED(J),LnD(J), LEA(I_,:,r)',LnA(I_,:,r)'); 
   %text( ps.XTick(4), ps.YTick(3), BadM{2} ); 
   for p=1:Lnme, p4(p).Marker = zn(p); p4(p).LineWidth = 1;end 
   title( sprintf('N=%d  h=%.2g',NA(r),hsA(r)) ); 
   ylabel('log_{10}nf');  xlabel('-log_{10}E'); 
   %legend(nametsD,'Location','northwest');  

                       % 3  T(N)  nf(N) DENSE  
   pTA  = permute( TA(:,:,1:MA), [1 3 2]);
   pLEA = permute(LEA(:,:,1:MA), [1 3 2]);
   K = LND > LNA(N0(j))*(1-eps); %in order for the definition areas to coincide
   scalog = floor(LNA(1)) : ceil(LND(end));
   scale  = num2cell(2.^scalog*100);

   for s = Ltols-Tail+1 : Ltols
      figure('Name',sprintf(' %s T(N) nf(N) DENSE %.2g',naODE,tols(s)));  
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
      %legend(namets_,'Location','southeast');
      ylabel('-log_{10}E');xlabel('N'); end
   
   TOC(j) = toc(ticj);
   end  % for j 
   
tocAll = toc(ticAll);   
fprintf('Total time = %s= %.3g < %.3g\n',sprintf('+ %g ',TOC),sum(TOC),tocAll)
ED 
