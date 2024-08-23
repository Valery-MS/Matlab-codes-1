% Computing stabil DS coeffs La,mu,nu for nonlinear ODEs

function BestCoeffs
                        %  PRESET CONSTANTS
naODEs = { 'n66',    'H146',    'fant'  };
Args   = {  [],      [],        []      };
t0fs   = [ [1 21];   [0 4];     [0 20]  ];
er     = @(ye,y) sqrt(sum((ye(:,1)-y(:,1)).^2)/length(y)); 

                           % SOURCE DATA
naL  = { 'La = ',   '(La, mu) = ',  '(La, mu, nu) = ' }; 
hLa  = 0.1;         hmu  = 0.1;     hnu  = 0.1; 

NDSs = {...
[100 250 500 1000 2500 4000 5000 6000 7500 8000 10000 12000 12500 15000 20000];...
[125 250 500 1000 1200 1250 1500 2000 2500 3000 4000  8000  12000];... % H146
[150 250 500 1000 2500 4000 5000 6000 7500 8000 10000 12000 12500 15000 20000]};
                    
tol   = [1e-5 1e-5 1e-3]; %[1e-9 1e-10];
nODEs = 1:3;        % numbers of considered ODEs from naODEs
zam   = [2 1 1];    % zamena(change): rfn-method -> zam(rfn)-method
ms    = 3;
BestC = cell(3,1);
nIL   = 1;                             % 1/2 if takeoff by es/(dop853|ode113)
ILs   = [16 1];      IL = ILs(nIL);    % Initial values Length of y 

                    % CALCULATED CONSTANTS                  
Las = -2+hLa : hLa : 2-hLa;   qLa = numel(Las);     % for m = 1  
mus = -2+hmu : hmu : 2-hmu;   qmu = numel(mus);     % for m = 2
nus = -2+hnu : hnu : 2-hnu;   qnu = numel(nus);     % for m = 3
qLmn = [qLa, round(qLa*qmu/4), round(qLa*qmu*qnu/8)]; 

                    % LOOP
for j = nODEs
   naODE = naODEs{j}; fprintf('%s\n',naODE);
   file  = sprintf('BestC_N_%s', naODE);
   Arg   = Args{j};                 % F's argument     
   t0f   = t0fs(j,:);   t0 = t0f(1);  tf = t0f(2);
   NDS   = NDSs{j};     MD = length(NDS);
   hs    = (tf-t0)./NDS;
   Emin  = Inf; 
   
   for m = ms 
      fprintf('m = %d\n',m);
      BestC{m} = cell(1,MD);   
      for s = MD:-1:1       
         Cms = nan(qLmn(m),m+2);  
         h   = hs(s); N = NDS(s);  
         t   = (t0:h:tf);
         y0I = esN(j,2,t(1:IL),Arg);
         ye  = esN(j,1,t,      Arg);             
         T   = 0; E = 1;  i = 0;        
         if m == 3         
            for nu = nus               
            for mu = mu3(hmu,abs(nu))   
            for La = La3(hLa,mu,nu)  
               Lmn = [La,mu,nu];
               op  = DSet(h,[],[], Lmn,zam );  
               tic; [t1,y,nf] = ds10(@RHFN,t,y0I,op,Arg,j); Ty=toc;
               Ei = er(ye,y);
               if Ei<tol(m)
                  %fprintf('Ei=%.3g T=%.4g h=%.4g nF=%d\n',Ei,Ty,h,nF);
                  i=i+1;    T = T+Ty;
                  if Ei<E,  Lmn0 = Lmn; E = Ei; end
                  Cms(i,:) = [Lmn Ei nf];end,end,end,end
         else errordlg('For m=4 program don"t ready'); end
                               %figure,plot(t,[ye,y(:,1)])
         if i
            T = T/i; 
            Cms         = Cms( ~isnan(Cms(:,m+1)), 1:m+2 );
            BestC{m}{s} = sortrows( Cms, m+1 );
            Er1 = BestC{m}{s}(1,m+1);
            if Er1 < Emin, Emin = Er1; sEmin = s; end
            Ea          = mean(Cms(:, m+1));   
            tLmn0 = sprintf( [naL{m} sprintf('%+.2f ',Lmn0)] );
            fprintf('N=%d h=%.4f T=%.3f  E=%.2g  Ea=%.2g %s\n',...
                     N,h,T,E,Ea,tLmn0); 
         else BestC{m}(1:s) = []; NDS(1:s) = []; break; end,end,end
 
   save(file,'BestC','NDS','Emin','sEmin');  
end