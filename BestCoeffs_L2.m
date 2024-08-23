% Computing stabil DS coeffs La,mu,nu for linear ODEs
% Function BestCoefs_L - only for ODEs Airy, Bessel, Rabi

function BestCoeffs_L2
Tt = tic;   fprintf('%30s\n',datetime('now'))
                           %  PRESET CONSTANTS                         
naODEs = { 'N3_43',  'N4_35',  'N8_57',  'N3_46',  'Nord6'};
nHom = 3;                                     % number of homogeneouse ODEs  

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

Args       = cell(size(naODEs));
Args{nHom} = [ a b c ];       
Args{5}    = [ 0.5  0.1   2/pi];
er = @(ye,y) sqrt(sum((ye(:,1)-y(:,1)).^2)/length(y)); % |ApproxSol-ExactSol| 

                           % SOURCE DATA
naL  = { 'La = ',   'La, mu = ',  'Lmn = ' }; 
hLa  = 0.1;         hmu  = 0.1;     hnu  = 0.1; 

NDSs={[150 250 500 1000 2500 5000 7500 10000 12500 15000 20000];... % N3_43
      [150 250 500 1000 2500 5000 7500 10000 12500 15000 20000];... % N4_35
      [150 250 500 1000 2500 5000 7500 10000 12500 15000 20000];... % N8_57
      [150 250 500 1000 2500 5000 7500 10000 12500 15000 20000];... % N3_46
      [150 250 500 1000 2500 5000 7500 10000 12500 15000 20000]};   % Nord6
                   
tolE  = [1e-5 1e-5 1e-3]; %[1e-9 1e-10];
nODEs = 5;      % numbers of considered ODEs from naODEs
zam   = [];       % zamena(change): rfn-method -> zam(rfn)-method
ms    = 3;
BestC = cell(3,1); 
nIL   = 1;                           % 1/2 if takeoff by es/(dop853|ode113)
ILs   = [8 1];      IL = ILs(nIL);    % Initial values Length of y 

                    % CALCULATED CONSTANTS                  
Las = -2+hLa : hLa : 2-hLa;   qLa = numel(Las);     % for m = 1  
mus = -2+hmu : hmu : 2-hmu;   qmu = numel(mus);     % for m = 2
nus = -2+hnu : hnu : 2-hnu;   qnu = numel(nus);     % for m = 3
qLmn = [qLa, round(qLa*qmu/4), round(qLa*qmu*qnu/8)]; 

                   % LOOP
%Q1 = 1; Q2 = [-1 1]; Q3 = [-1 1]; Q4 = 1; Q5 = 1:3; Q6 = 1;   
%for a=Q1, for b=Q2, for c=Q3,for d=Q4, for e=Q5, for f=Q6, Arg = [a b c d e f];
for j = nODEs
   naODE = naODEs{j}; fprintf('%s\n',naODE);
   file  = sprintf('BestC_N_%s', naODE);
   Arg   = Args{j};             % F's argument     
   t0f   = t0fs(j,:);  t0 = t0f(1);  tf = t0f(2);
  
   if j <= nHom, ds10L = @ds10Lh;
   else          ds10L = @ds10Ln;  end
   
   NDS  = NDSs{j};  MD = length(NDS);
   hs   = (tf-t0)./NDS;
   Ermin = Inf; 
   
   for m = ms
      fprintf('m = %d\n',m);
      BestC{m} = cell(1,MD); 
      
      for s = MD:-1:1       
         Cms = nan(qLmn(m),m+2);  
         h   = hs(s); N = NDS(s);  
         t   = t0:h:tf;
         y0I = esL(t(1:IL),Arg,j,2); 
         ye  = esL(t,Arg,j,1);
         ye2 = esL(t,Arg,j,2);
         T   = 0; E = 1;  i = 0;    
         
setm = @(h,RT,AT,d) dopset('InitialSt',h,'RelT',RT,'AbsT',AT);
met = @dop853;
h0=0.01;tol=10^(-15.65);op2=setm(h0,tol,tol,[]);warning('off','all');
y0=esL(t(1),Arg,j,2);

try tic; [t1,y2,nf]=met(@RHFL,t,y0,op2,Arg,j);warning('on','all');TD=toc;
    ErD = er(ye,y2);
catch ME,  ErD = NaN;
end
%}
         if m == 3         
            for nu = nus               
            for mu = mu3(hmu,abs(nu))   
            for La = La3(hLa,mu,nu)  
               Lmn = [La,mu,nu];
               op  = DSet(h,[],[], Lmn,zam );  
               tic; y = ds10L(@GMs,t,y0I,op,Arg,j); Ty=toc;
               Ei = er(ye,y);
               if Ei<tolE(m)
                  %fprintf('Ei=%.3g T=%.4g h=%.4g nF=%d\n',Ei,Ty,h,nF);
                  i=i+1;    T = T+Ty;
                  if Ei<E,  Lmn0 = Lmn; E = Ei; end
                  Cms(i,:) = [Lmn Ei T];end,end,end,end
         else errordlg('For m=4 program don"t ready'); end
                               %figure,plot(t,[ye,y(:,1)])
         if i
            T = T/i; 
            Cms         = Cms( ~isnan(Cms(:,m+1)), 1:m+2 );
            BestC{m}{s} = sortrows( Cms, m+1 );
            Ea  = mean(Cms(:, m+1));   
            Er1 = BestC{m}{s}(1,m+1);
            if Er1 < Ermin, Ermin = Er1; sErmin = s; end
            tLmn0 = sprintf( [naL{m} sprintf('%+.2f ',Lmn0)] );
            fprintf('N=%d h=%.4f T=%.3f %.3f  E=%.2g %.2g  Ea=%.2g %s\n',...
                     N,h,TD,T,ErD,E,Ea,tLmn0); 
          end,end,end        
      %   else BestC{m}(1:s) = []; NDS(1:s) = []; break; end,end,end
toc(Tt) 
   save(file,'BestC','NDS','t0f','Ermin','sErmin');  
end  
%end,end,end,end,end,end
