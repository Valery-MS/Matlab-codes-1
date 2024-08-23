% Самостоятельная ф-я ( без COMMODE)
% Computing stabil DS coefficients La,mu,nu for Linear ODEs Airy, Bessel, Rabi 
% & Grafics T(E) - Time of Error( or Accuracy) 

function BestCoeffs_L

FU = 1;
if FU==1                         % Airy
   file = 'Airy';
   t0 = 0;  tf = 100;                        % Start & Final time   
   omega = 1;  o2 = omega^2;                 % Angular frequency at t=1
   alfa = exp(1.0i*pi/3)*omega^(2/3);
   ugly = 3^(1/6)*gamma(2/3)/2.0;            % ugly_const 
   es1 = @(t) real(ugly*(sqrt(3)*airy(0,alfa*t)+airy(2,alfa*t))); %Ignor Im
   es  = @(t) [es1(t); ...
               real(alfa*ugly*(sqrt(3)*airy(1,alfa*t)+airy(3,alfa*t)))];
   y0  = es(t0);                             % initial state
   er = @(y,ye,n) sqrt(sum((y(1,:)-ye).^2)/n);
   GM = @(t)   [ 0  1; -o2*t   0 ];          % Generator matrix   
   f  = @(t,y) [ y(2); -o2*t*y(1)];          % func for dop853/ode113 in ds10L

elseif FU==2                     % Bessel
   file = 'Bessel';
   nBes = 1;
   t0 = 1;  tf = 100;
   es1= @(t) besselj(nBes,t);
   es = @(t) [ es1(t); 0.5*(besselj(nBes-1,t)-besselj(nBes+1,t))]; 
   y0 = es(t0);                              % initial state
   er = @(y,ye,n) sqrt(sum((y(1,:)-ye).^2)/n);
   GM = @(t)   [ 0  1; (nBes/t^2-1),        -1/t ];  % Generator matrix   
   f  = @(t,y) [ y(2); (nBes/t^2-1)*y(1) -y(2)/t ];  % func for ode113 in ds10L

elseif FU==3                     % Shrod Qubit Flip
   file = 'Qubit';  
   t0 = 0; 
   y0  = [1;0]; 
   p_amp = 0.05;   p_freq = 0.99;          % pulse ampl, pulse freq
   f_ome = sqrt(p_amp^2 + (p_freq-1)^2);   % flop omega, detuning=p_freq-1
   tf   = pi / f_ome;                      % pulse time
   es1 = @(t) [];  
   er = @(y,ye,n) abs(y(:,n)'*y(:,n)-1); 
   f_amp = p_amp^2 / f_ome^2;                   % flop ampl use for QubNrm
   Px = [0,1;1,0]; Pz = [1,0;0,-1]; % Py=[0,-1.0j;1.0j,0]; Pauli matrices 
   GM = @(t) 1.0j*(0.5*Pz+p_amp*cos(p_freq*t)*Px);    %=-i*Hamiltonian
   f  = @(t,y) GM(t)*y;
   es = @(t) solveonce(GM,t,y0);
end 
G = { GM f };
  % G_p = @(t) 1.0j*(0.5*Pz + p_amp*cos(p_freq*t)*Px); % то же, что и GM
  
                       % Difference Schemes
flag = 3; % 1: Lmn-coeffs calculation only  2: Grafics T(E) only 
          % 3: Lmn-coeffs calculation & Grafics T(E)
naL = {'La = ' '(La, mu) = ' '(La, mu, nu) = '}; 
ms = [ 3];   
hs = [0.002 0.005 0.01 0.02 0.05 0.1];
hs = 0.01;           % forced assign of hs
qh = numel(hs);

if flag == 1 || flag == 3
   met = {'m = 1:', 'm = 2:', 'm = 3:'};
   C   = cell(3,1);  Eas = nan(3,qh); 
   hLa  = 0.1;  Las = -2+hLa : hLa : 2-hLa;   qLa = numel(Las);      
   hmu  = 0.1;  mus = -2+hmu : hmu : 2-hmu;   qmu = numel(mus);
   hnu  = 0.1;  nus = -2+hnu : hnu : 2-hnu;   qnu = numel(nus); 
   qLmn = [qLa, round(qLa*qmu/4), round(qLa*qmu*qnu/8)];
   tol  = [1e-6 1e-7 1e-8]; %[1e-9 1e-10];
   
   for m = ms
      fprintf('\n'); 
      C{m} = cell(1,qh); 
      
      for s = 1:qh       
         Cms = nan(qLmn(m),m+1);   
         h = hs(s);    
         t = maketimes(t0,tf,'stepsize',h); n = numel(t);
         ye1 = es1(t);
         T = 0; E = 1;  i = 0;
         
         if m == 1
            for La = Las
               b = cDS(La);
               tic;  y = DS6L( G, t, y0, h, La, b);   T = toc+T; 
               Ei = er(y,ye1,n);
               if Ei<tol(m), i = i+1;
                  if Ei<E,  Lmn0 = La;  E = Ei; end
                  Cms(i,:) = [La Ei]; end,end
          
         elseif m == 2
            for mu = mus   
            for La = La2(hLa,abs(mu))
               Lm = [La,mu]; b = cDS(Lm);  
               tic;  y = DS8L( G, t, y0,h, Lm, b );  T = T+toc; 
               Ei = er(y,ye1,n);
               if Ei<tol(m), i = i+1;  
                  if Ei<E,  Lmn0 = Lm; E = Ei; end
                  Cms(i,:) = [Lm Ei]; end,end,end 
         else           
            for nu = nus               
            for mu = mu3(hmu,abs(nu))   
            for La = La3(hLa,mu,nu)     
               Lmn = [La,mu,nu];  b = cDS(Lmn);   
               tic; y = DS10L( G, t, y0, h, Lmn, b );  T = T+toc;
               Ei = er(y,ye1,n);  
               if Ei<tol(m), i=i+1;
                  if Ei<E,  Lmn0 = Lmn; E = Ei; end
                  Cms(i,:) = [Lmn Ei];end,end,end,end,end

         tLmn0 = sprintf([naL{m} '%s'], sprintf('%+.2f ',Lmn0));
         
         if i
            T = T/i; 
            Cms = Cms(~isnan(Cms(:,m+1)),1:m+1);
            C{m}{s} = sortrows(Cms,m+1);
            Ea = mean(Cms(:,m+1)); Eas(m,s) = Ea;
            txt = sprintf('%7s h=%.4f  T=%.2f  E=%.2g  Ea=%.2g  %s',...
                  met{m},h,T,E,Ea,tLmn0); 
            fprintf('%s\n',txt);
            
         else C{m}(s:end) = []; hs(s:end) = []; break; end,end,end
 
   save(file,'C','hs');end

if flag > 1 
   if flag == 2, load(file,'C'); end
   M = 1; 
   N = 10;  de = 0:N-1;
   hc = nan(4,N,3);
   hc(:,:,1) = ...                                                %Airy
     [0.002*1.95.^de; 0.00296*1.191.^de; 0.0057*1.33.^de; 0.0061*1.4.^de];
   hc(:,:,2) = ...                                                %Bess
     [0.007*1.62.^de; 0.022*1.4.^de; 0.044*1.3.^de; 0.012*1.6.^de];
   hc(:,:,3) = ...                                                %Qubit
     [0.025*1.4.^de; 0.045*1.4.^de; 0.045*1.25.^de; 0.1*1.3.^de];

   K = {[2 ]  [1 ] [1 ] };  %K = {[1 2] [1 150] [1 92]};  
   qK = numel(K{1}); 
   ms = 1:3; 
   fprintf('N=%d',N);  fprintf('\n');
  
                        % Метод Магнуса (m=4)
   TM = zeros(N,1);  EM = nan(N,1);
   
   if 1
      for i = 1:N, h = hc(4,i,FU);
         t=maketimes(t0,tf,'stepsize',h); n=numel(t); TS=0;
         for l=1:M,tic; F=generate(GM,t); y=evolve(y0,F); TS = toc+TS;end 
         TM(i) = TS/M;   EM(i) = er(y,es1(t),n);  end
     
      LEM = -log10(EM);  end

   Sp = {{'-.k' '-.b' ':g' '-.g'}  ...
        {'--k' '--b' ':b' '-.b'} {'-k' '-b' ':m' '-.m'}};
   DSs  = {@DS6L  @DS8L  @DS10L};  
   naLmn = {'\\lambda = ' '\\lambda, \\mu = ' '\\lambda, \\mu, \\nu = '};
   Ti = nan(N,qK,3);  Er = Ti;
   Leg = {'Magnus' '' '' ''}; W=0;
   
   for s = 1  % 1:qh        % methods 50000-5000  36000-3000   14800-2500
      if ishandle(FU), close(FU); end %0.002-0.02  0.0028-0.033 0.007-0.04
      figure(FU);hold on; plot(LEM,TM,'-r','LineWidth',1.2);
      %title(sprintf('%s eq.',file));  
      
      for m = ms 
         for i = 1:N,    h = hc(m,i,FU);
            t = maketimes(t0,tf,'stepsize',h);  n = numel(t);       
            ye1 = es1(t);
            q = 0;  TS = 0; iL = qK*(m-1)+1;
            
            for k = K{m}
               q = q+1; W=W+1;
               if mod(W,200)==0, fprintf('%.f ',W/100) ;end;  
               [ Lmn, b ] = DSet( C, m, s, k );
               Leg{iL+q} = sprintf(['m=%d %d) ' naLmn{m} '%s'],...
                           m,k,sprintf('%g  ',Lmn));
               for l=1:M
                  tic; y = DSs{m}(GM, t, n, h, Lmn, b, es); TS=TS+toc;end  
              
               Eri = er(y,ye1,n); 
               if Eri > 1, Er(i,q,m) = NaN; else Er(i,q,m) = Eri;end
               Ti(i,q,m) = TS/M+tocF; end
               end  
      
         LE = -log10(Er);
         %plot(LE(:,:,m),Ti(:,:,m));end
         for q=1:qK
            plot(LE(:,q,m),Ti(:,q,m),Sp{m}{q},'LineWidth',1.2); end
         end
 
      legend(Leg,'Location','northwest'); 
      xlabel('-log_{10}E'); ylabel('Time');
      hold off;
      end
   end 