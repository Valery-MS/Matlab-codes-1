% Самостоятельная ф-я ( без COMMODE)
% Computing best coeffs La,mu,nu for nonlinear ODEs: vdp, DNK, ...
% Analisys solution of vdp eq with mu=1.25 at best Lmn

function AccTim_n

%               %  PRESET CONSTANTS

nams = { 'm = 1:',  'm = 2:',       'm = 3:'};
naL  = { 'La = ',   '(La, mu) = ',  '(La, mu, nu) = ' }; 
hLa  = 0.1;         hmu  = 0.1;     hnu  = 0.1; 

c_a     = [344 446 92 144]*1e-4;   p1 = pi-0.4;
y0_DNK  = zeros(4,1); 
y0_fant = [1; 0; 0; 1]; % u=[exp(sin t); ln(t+1); atan(t); (t+1)^-0.5]
hs      = [ 0.002 0.005 0.01 0.02 0.05 0.1 ];

namets = { '@ode45', 'ode113', 'dop853', 'rkf78', 'DS10'};    
naODEs = {'vdpm',    'DNK',    'n58',    'n60',   'fant',   'Airy' };
                   % 1st exact solution coordinate 
ES1s = {[], ...                            % vdpm
        [], ...                            % DNK
        @(t) t./(2-t.^3).^(1/3),       ... % n58
        @(t) (1+cos(t)).*exp(tan(t/2)),... % n60
        @(t) exp(sin(t)),              ... % fant ln(t+1) atan(t) 1/sqr t+1
        [] };                              % Airy    
er = @(y,ye,n) sqrt(sum((y(1,:)-ye).^2)/n); % |ApproxSol-ExactSol| 

ESas   = { @dop853,  @ode113 };   % Approximation to unknown ES
optas  = { @dopset,  @odeset };   % data set to ESapx
naESas = { 'DP8',    '113' };

                 % SOURCE DATA
RT = eps;  AT = RT;
ms = 3; 
CN = 1;              % Coordinate Number where Initial Value iv is changed
iv = 10;             % initial values (in degrees if CN<3)
mu_vdpm  = 1;
omegAiry = 1;        % Angular frequency at t=1
tol = [1e-6 1e-7 1e-8]; %[1e-9 1e-10];

t0fs  = [ [0 20];  [0 100]; [1 1.2]; [0 p1];   [0 20];   [0 100]];
y0s   = { [2; 0],  y0_DNK,  [1; 2],  [2; 1],   y0_fant,  [1; 0] };
nODEs = [ 3 ];      % numbers of considered ODEs from naODEs



h0    = 0.001;       % initial step for ode113
hs=hs/100;


nESa  = 2;          % number of ESapx, approximation to ES

                  % CALCULATED CONSTANTS
                  
Las = -2+hLa : hLa : 2-hLa;   qLa = numel(Las);     % for m = 1  
mus = -2+hmu : hmu : 2-hmu;   qmu = numel(mus);     % for m = 2
nus = -2+hnu : hnu : 2-hnu;   qnu = numel(nus);     % for m = 3
qLmn = [qLa, round(qLa*qmu/4), round(qLa*qmu*qnu/8)]; 

Fargs = { mu_vdpm,                        c_a, [], [], [], omegAiry};                             
nargs = { sprintf(' mu_vdpm=%g',mu_vdpm), '',  '', '', '', '' };

alfa    = exp(1.0i*pi/3)*omegAiry^(2/3);
ugly    = 3^(1/6)*gamma(2/3)/2.0;        
ES1s{6} = @(t) real(ugly*(sqrt(3)*airy(0,alfa*t)+airy(2,alfa*t))); % Airy

if CN<3, iv = iv*pi/180;  end;
y0_DNK(CN)  = iv;
y0s{2}      = y0_DNK;
ESa    = ESas{nESa};
optESa = optas{nESa}('RelTol',RT, 'AbsTol',AT, 'InitialStep',h0);

                   % MAIN CYCLE of BEST COEFFS CALC
qh     = numel(hs);
BestCn = cell(3,1); 
Eas    = nan(3,qh); 

Lmns = [-2.2000    3.0000   -0.8000; ...   % Good oeffs La, mu, nu 
        -1.9000    2.7000   -0.7000; ...   % for mu_vdpm = 1.25
        -1.8000    2.5000   -0.7000; ...
        -0.7000    1.2000   -0.4000; ...
        -1.8000    2.6000   -0.7000; ...
        -1.7000    2.4000   -0.7000; ...
        -1.4000    2.1000   -0.6000; ...
        -1.7000    2.5000   -0.7000; ... 
        -1.6000    2.4000   -0.7000 ];
ks = 1:9;

for j = nODEs
   naODE = naODEs{j};
   F     = str2func(['F_' naODE ]);     % ODE func 
   F_M   = str2func(['F_' naODE '_M']); % F(y}^T in nt0 t-points for start
   file  = ['BestC_hs_' naODE];
   if j==1, file = sprintf('%s_%g_%s', file, mu_vdpm, naESas{nESa}); end
   Farg  = Fargs{j};                  % F's argument    
   narg  = nargs{j}; 
   t0f   = t0fs(j,:);  t0 = t0f(1);  tf = t0f(2);
   y0    = y0s{j};
   
   for m = ms 
      BestCn{m} = cell(1,qh); 
      
      for s = 1:qh       
         Cms = nan(qLmn(m),m+1);   
         h = hs(s);   
         t = maketimes(t0,tf,'stepsize',h); n = numel(t);
                  
         if j > 2, ye1 = ES1s{j}(t);    % ES is known
         else  end
             warning('off','all');     % ES approx by dop853 or ode113
              tic; [t_,esa] = ESa( F, t, y0, optESa, Farg); TES=toc;
              % ye1 = esa(:,1)';
              warning('on','all');  
         Eesa = er(esa',ye1,n); 
         figure, plot(t,ye1);
         title(sprintf('met=%s T=%.3g E=%g',func2str(ESa),TES,Eesa));
         
         for k = ks
            Lmn = Lmns(k,:); b = cDS(Lmn);    
            tic; y = DS10(F, t, y0, h, Lmn, b, F_M, Farg); Ty= toc;
            figure, plot(t,y(1,:));
            Ey = er(y,ye1,n);
            title(sprintf('%g) Ty=%.3g E=%g Lmn=%.2f %.2f %.f',k,Ty,Ey,Lmn));       
         end
 
                     
         T = 0; E = 1;  i = 0;
      
         if m == 3         
            for nu = nus               
            for mu = mu3(hmu,abs(nu))   
            for La = La3(hLa,mu,nu)     
               Lmn = [La,mu,nu];  b = cDS(Lmn);    
               tic; y = DS10(F, t, y0, h, Lmn, b, F_M, Farg); Ty= toc;
               Ei = er(y,ye1,n);  
               if Ei<tol(m)
                  i=i+1;    T = T+Ty;
                  if Ei<E,  Lmn0 = Lmn; E = Ei; end
                  Cms(i,:) = [Lmn Ei];end,end,end,end

         else errordlg('For m=4 program don"t ready'); end
        
         if i
            T = T/i; 
            Cms          = Cms( ~isnan(Cms(:,m+1)), 1:m+1 );
            BestCn{m}{s} = sortrows( Cms, m+1);
            Ea           = mean(Cms(:, m+1));      Eas(m,s) = Ea;
            tLmn0 = sprintf([naL{m} '%s'], sprintf('%+.2f ',Lmn0), narg);
            txt = sprintf('%7s h=%.4f  T=%.2f  E=%.2g  Ea=%.2g  %s',...
                  nams{m},h,T,E,Ea,tLmn0); 
            fprintf('%s\n',txt);
            
         else BestCn{m}(s:end) = []; hs(s:end) = []; break; end,end,end
   end
 
 save(file,'BestCn','hs');   
 