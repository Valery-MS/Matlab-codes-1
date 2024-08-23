function Task_S_C_lin(handles)
global  t0 tn y0 es1 es er GM

handles.h1_h2.Enable = 'off';
ExactS_GenerF(handles);

M = 1; 
N = 10;  de = 0:N-1;
NEq = 2;
FBC = { 'Airy_BestC_hs_InpD', 'Bess_BestC_hs_InpD','Shrod_BestC_hs_InpD'};
load( FBC{NEq}, 'BestC','hs');

switch handles.Linear.UserData;
case 1,  hc = ...      %Airy                 
  [0.002*1.95.^de; 0.00296*1.191.^de; 0.0057*1.33.^de; 0.0061*1.4.^de];
case 2,  hc = ...      %Bess
  [0.007*1.62.^de; 0.022*1.4.^de; 0.044*1.3.^de; 0.012*1.6.^de];
case 3,  hc = ...      %Shro
  [0.025*1.4.^de; 0.045*1.4.^de; 0.045*1.25.^de; 0.1*1.3.^de];
end

K = {[2 ]  [1 ] [1 ] };   % K = {[1 2] [1 150] [1 92]};  
qK = numel(K{1});         
ms = 1:3; 
fprintf('N=%d',N);  fprintf('\n');
  
                        % Метод Магнуса (m=4)
TM = zeros(N,1);  EM = nan(N,1);  
for i = 1:N
   h = hc(4,i);
   t = maketimes(t0,tn,'stepsize',h);  n = numel(t); TS = 0;
   for l=1:M, tic; F = generate(GM,t); y = evolve(y0',F); TS = toc+TS;end 
   TM(i) = TS/M;   EM(i) = er(y,es1(t),n);  end    
LEM = -log10(EM);  
                        % LinDS
Sp = {{'-.k' '-.b' ':g' '-.g'}  ...
     {'--k' '--b' ':b' '-.b'} {'-k' '-b' ':m' '-.m'}};

DS  = {@DS11  @DS21  @DS31};  
naLmn = {'\\lambda = ' '\\lambda, \\mu = ' '\\lambda, \\mu, \\nu = '};
Ti = nan(N,qK,3);  Er = Ti;
Leg = {'Magnus' '' '' ''}; W=0;
   
for s = 1  % 1:qh       % methods   50000-5000  36000-3000   14800-2500
   if ishandle(nF), close(nF); end %0.002-0.02  0.0028-0.033 0.007-0.04 
   figure(nF);hold on; plot(LEM,TM,'-r','LineWidth',1.2);
   %title(sprintf('%s eq.',file));  
      
   for m = ms
      for i = 1:N
         h = hc(m,i);
         t = maketimes(t0,tn,'stepsize',h);  n = numel(t);    
         ye1 = es1(t);
         q = 0;  TS = 0; iL = qK*(m-1)+1;
            
         for k = K{m}
            q = q+1; W=W+1;
            if mod(W,200)==0,fprintf('%.f ',W/100);end;  
            [ Lmn, b ] = DSet( BestC, m, s, k );   
            Leg{iL+q} = sprintf...
            (['m=%d %d) ' naLmn{m} '%s'],m,k,sprintf('%g  ',Lmn));

            for l=1:M
                tic; y = DS{m}(GM, t, n, h, Lmn, b, es); TS = TS+toc; end  
            Eri = er(y,ye1,n); 
            if Eri > 1, Er(i,q,m) = NaN; else Er(i,q,m) = Eri;end
            Ti(i,q,m) = TS/M+tocF; end,end  
      
      LE = -log10(Er);
      %plot(LE(:,:,m),Ti(:,:,m));end
      for q=1:qK, plot(LE(:,q,m),Ti(:,q,m),Sp{m}{q},'LineWidth',1.2);
      end,end
 
   legend(Leg,'Location','northwest'); 
   xlabel('-log_{10}E'); ylabel('Time');
   hold off; end