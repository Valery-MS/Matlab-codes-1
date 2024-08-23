% Computing stabil DS coefficients La,mu,nu for Linear ODEs:
%  Bessel, Airy & Shrod 

function Task_BestC(handles)

global t0 tn h es1 es er GM

ExactS_GenerF(handles);
hs  = str2num(handles.hs.String);
naL = {'La = ' '(La, mu) = ' '(La, mu, nu) = '}; 
ms  = [3 ]; 
qh  = numel(hs);

met = {'m = 1:', 'm = 2:', 'm = 3:'};
BestC= cell(3,1);  Eas = nan(3,qh);
hLa  = 0.1;  Las = -2+hLa : hLa : 2-hLa;   qLa = numel(Las);      
hmu  = 0.1;  mus = -2+hmu : hmu : 2-hmu;   qmu = numel(mus);
hnu  = 0.1;  nus = -2+hnu : hnu : 2-hnu;   qnu = numel(nus); 
qLmn = [qLa, round(qLa*qmu/4), round(qLa*qmu*qnu/8)];
tol  = [1e-6 1e-7 1e-8];    %[1e-9 1e-10];
   
for m = ms
   BestC{m} = cell(1,qh); 
      
   for s = 1:qh      
      Cms = nan(qLmn(m),m+1);   
      h = hs(s);    
      t = maketimes(t0,tn,'stepsize',h); n = numel(t);
      ye1 = es1(t);
      T = 0; E = 1;  i = 0;
                               
      if m == 1        
         for La = Las
            b = cDS(La); 
            tic;  y = DS6L( GM, t, h, La, b, es);   T = toc+T; 
            Ei = er(y,ye1,n);
            if Ei<tol(m), i = i+1;
               if Ei<E,  Lmn0 = La;  E = Ei; end
               Cms(i,:) = [La Ei]; end,end
         
      elseif m == 2          
         for mu = mus
         for La = La2(hLa,abs(mu))
            Lm = [La,mu]; b = cDS(Lm);  
            tic;  y = DS8L( GM, t, h, Lm, b, es);   T = T+toc; 
            Ei = er(y,ye1,n);
            if Ei<tol(m), i = i+1;  
               if Ei<E,  Lmn0 = Lm; E = Ei; end
               Cms(i,:) = [Lm Ei]; end,end,end 
      else                           
         for nu = nus               
         for mu = mu3(hmu,abs(nu)) 
         for La = La3(hLa,mu,nu)     
            Lmn = [La,mu,nu];  b = cDS(Lmn);  
            tic; y = DS10L( GM, t, h, Lmn, b, es);  T = T+toc;              
            Ei = er(y,ye1,n);  
            if Ei<tol(m), i=i+1;
               if Ei<E,  Lmn0 = Lmn; E = Ei; end
               Cms(i,:) = [Lmn Ei];end,end,end,end,end

      tLmn0 = sprintf([naL{m} '%s'], sprintf('%+.2f ',Lmn0));
      
      if i
         T = T/i; 
         Cms = Cms(~isnan(Cms(:,m+1)),1:m+1);
         BestC{m}{s} = sortrows(Cms,m+1);
         Ea = mean(Cms(:,m+1)); Eas(m,s) = Ea;
         txt = sprintf('%7s h=%.4f  T=%.2f  E=%.2g  Ea=%.2g  %s',...
               met{m},h,T,E,Ea,tLmn0); 
         fprintf('%s\n',txt);
         
      else BestC{m}(s:end) = []; hs(s:end) = [];  break; end,end,end

save(handles.EqNam.UserData{2},'BestC','hs');