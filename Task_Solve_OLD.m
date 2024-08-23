% Solving of DNK equation
% Старый работающий вариант (Old working variant of program)
function Task_Solve_OLD(handles)
global SpecInpD a NEq_  IV cT_
global GenInpD  t0 tn h1_h2 RelT AbsT nSDE y0
global Figs  yTA

y0 = y0(:);
y0i = y0;
h = h1_h2(1);   h2 = h+h;
NuMet = handles.Methods.Value; % 1-RK, 2-DP, 3-GBS, 4-nDS

                  %  Set of optional parms
   NEq = NEq_(1); m = NEq_(2); Nh = NEq_(3); ks = NEq_(4:end); kLmn = ks(1);
   switch NEq
      case 1, load(handles.Bess.UserData{2},'BestC','hs');
      case 2, load(handles.Airy.UserData{2},'BestC','hs');
      case 3, load(handles.Shrod.UserData{2},'BestC','hs');
      otherwise, errordlg('Wrong value of a file number'); end      
   zam  = [2 1 2];  % zam(1:2) iteration changes; zam(3)=1/2 if dop853/ode113
   Lmn  = BestC{m}{Nh}(kLmn,1:m);   
   op   = DSet( h, RelT,AbsT,Lmn, zam);
   meth = 'ds10';   
   
if NuMet == 3, errordlg('GBS is not yet ready');
else
   if NuMet==1, met = @ode45;  metset = @odeset;
   else         met = @dop853; metset = @dopset; end   
   opDP = metset('RelTol',RelT,'AbsTol',AbsT,'InitialStep',h);
   meth = handles.Methods.String{handles.Methods.Value}; end

                %  Calc of solution
CN = IV(1);             % Coordinate Number where Initial Value is changed
hd  = IV(4);                % h in degries 
ivd = (IV(2)-1 : IV(3))*hd; % initial fi-values(degries) 
rd  = pi/180;               % radians in degries
iv  = ivd*rd;               % initial values in radians
% NOTE: ivd(1),ivr(1) aren't used to match with TT starting with the no 2

qiv  = numel(iv);   % quantity of periods
qsBy = IV(7);   % 0-3 - q-y of subBlocks only for y plot (without Check)
PerCheck = IV(6);  % Period Checking
% % cT, IV(8) - reserve

if PerCheck
   qsBs = [7 3 2]; % Check: [no Er(NC) Er(1:2) Er] q-y of SubBlocks in qB
   qsB  = qsBs(PerCheck);
else
   qsB  = qsBy;
   qpls = [5 8 6 4 5 6];    % q-y of plots for each qsBy
   qpl  = qpls(qsBy);  end

phi = {'\phi_1', '\phi_2', 'v_1', 'v_2'};
qB = ceil((IV(3) - IV(2)+1)/qsB);    % quantity of Blocks
if CN == 1, phiIVI = '\phi_{10} = '; 
else        phiIVI = '\phi_{20} = '; end 

if     PerCheck && qsBy,  Figs = nan(2*qB,1);
elseif PerCheck || qsBy,  Figs = nan(qB,1);
else                      Figs = [];  end

p = 0;    % Figure Number
nom  = [1 2 4 5];    % graphic numbers on subplots
qiv1 = qiv+1;
TT   = cell(1,qiv1); del = TT;  
% TT{1} = [348.5 1991 2340 ];
TT{1} = [348.5 1991 2340 ];
i = 1;  del0 = 0.1;

if     tn <= 3000, tik = 100;
elseif tn <= 6000, tik = 200;
elseif tn <= 7000, tik = 500;   
else               tik = 1000; end

T1 = 349;        qT = ceil(tn/T1);      nT = round(T1/h); 
Nopp = 55;       ns = Nopp+nT*(0:qT);   t9 = ns(end)*h;
t = (t0:h:t9)';  n  = numel(t);         zer = zeros(size(t));
   
for nF = 1:qB
   if PerCheck 
      fig1 = figure('Position',[1,1,500,1000],'MenuBar','none',...
           'Name',meth,'NumberTitle','on');       
      p=p+1;  Figs(p) = fig1.Number; end        
       
   if qsBy
      fig2 = figure('Position',[1,1,500,1000], ...
           'Name',meth,'NumberTitle','on'); 
      p=p+1;  Figs(p) = fig2.Number; end
          
   for s = 1:qsB 
      T1 = TT{i}(1); 
      i = i+1;    if i > qiv1, break,end 
      r = 0;    q = 1;  
      fi0 = sprintf('%d',ivd(i));
      ivi = iv(i);  y0i(CN) = ivi;  Ey = IV(5)*ivi;
        
      tic;  
      if NuMet==4     
         y = DSm(@DNKf,@DNKv,t0,n,y0i,opDS,a,h,m)';  % According to dop853
      elseif NuMet==3, errordlg('GBS is not yet ready');
      else
         warning('off','all');   [t, y] = met(@DNKf,t,y0i,opDP,a);
         warning('on','all'); end
      Tim = toc; 
     
      yi = y(:,CN);
      Nopp = find(ivi*yi<0,1); % Number of 1st element opposed to y0
      yi(1:Nopp) = NaN;         % Turn off value y0 from yi
      neps = 4*Nopp; 
      
      for j = 1:qT-1     
         [ma, nLo] = min(abs(yi(ns(j):ns(j+1)-1)-ivi)); % Local number
         M  = ns(j)+nLo-1;      % Global number
         tM = t(M);    
         dr = tM/T1;  rdr = round(dr);
         
         M_1 = M-1;   tM_1 = t(M_1);  fta = yi(M_1);   
         M1  = M+1;   tM1  = t(M1);   ftb = yi(M1);
            
         if abs(ftb-fta)<Ey, w=0.75*tM; t1 = 0.25*tM_1+w; t2 = 0.25*tM1+w;
         elseif fta<ftb,     t1 = tM;   t2 = (tM+tM1)*0.5;
         else                t1 = (tM_1+tM)*0.5;  t2 = tM;  end
      
         % sign before dop853t = -sign(y0iC); 
         warning('off','all')
         Tj = fminbnd(@(t) -dop853t(t,CN,tM,y(M,:)',@DNKf,a),...
              t1,t2,optimset('TolX',eps));  % Period T Accuraced
         warning('on','all');         

         %[tf1, yf1] = met(@DNKf,Tj:h:Tj+(neps-1)*h,yTA',opDP,a);         
         yf = DSm(@DNKf,@DNKv,Tj,neps,yTA',opDS,a,h,m)';   % fractal
         
         nf = Nopp:neps;
         delj = max(abs(y(nf,CN)-yf(nf,CN)))/ivi;
         if     j==1,   r=r+1; TT{i}(r)=Tj;  del{i}(r)=delj;
         elseif delj<=del0
            if abs(dr-rdr)<0.05
               if j==2, r=r+1;end
               q=q+1;   TT{i}(2)=1e3*q+rdr;  del{i}(2)=delj; 
            else r=r+1; TT{i}(r)=Tj;         del{i}(r)=delj;  end
         elseif r>2 || j+1>qT,  break,end,end             
      [ i j t(M)],TT{i},del{i}
         
                  % Solution graphics
      if qsBy > 0
         figure(fig2)
         if qsBy == 3  % plot fi1, fi2 or plot fi1,...fi4)
            for k = 1:2
              ax = subplot(qpl,1,2*(s-1)+k); plot(t,y(:,k),'b',t,zer,'k'); 
              ylabel(phi{k}); 
              set(ax,'XMinorGrid','on','XMinorTick','on','XTick',0:tik:t9);
              if k==1, title([phiIVI fi0]);end,end
         elseif qsBy == 1
            for k = 1:4 
              ax = subplot(qpl,1,nom(k)); 
              plot(t,y(:,k),'b',t,zer,'k'); 
              ylabel(phi{k}); 
              set(ax,'XMinorGrid','on','XMinorTick','on','XTick',0:tik:t9);
              if k == 1
                 tTT = sprintf('%.5g  ',TT{i}); 
                 title(sprintf('%d) %s   %s',i,[phiIVI fi0],tTT)); 
                 ax = subplot(qpl,1,3); 
                 plot(t,y(:,1:2),t,zer); 
                 ylabel([phi{1} ', ' phi{2}])
                 set(ax,'XMinorGrid','on','XMinorTick','on',...
                 'XTick',0:tik:t9); end,end
         else
            ax = subplot(qsBy,1,s); plot(t,y(:,1),'b',t,zer,'k'); 
            ylabel(phi{1});
            set(ax,'XMinorGrid','on','XMinorTick','on','XTick',0:tik:t9);
            tTT = sprintf('%.5g  ',TT{i}(max(1,end-25):end)); 
            title(sprintf('%d) %s    T = %s',i-1,[phiIVI fi0],tTT));end,end
          
   if i > qiv, break,end
   %xlabel('t');
   end, end
 
if exist('fig1','var'),set(fig1,'Name',sprintf('%s  T=%.3g',meth,Tim));end
if exist('fig2','var'),set(fig2,'Name',sprintf('%s  T=%.3g',meth,Tim));end

if PerCheck 
   figure('Name',sprintf('%s  tn=%g h=%g',meth,tn,h));
   plot(ivd,TT(:,1));xlabel('\phi_0'); ylabel('T'); end

%T = [348.5  796.4  546.8  2184.2  5257  1830.4  591.5  394.3  592];
%T=[349 1842 796  547 2184 5753 5257 1830 592 247 394 880 592];
figure; plot(ivd,TT(:,1)); 
xlabel(['Initial Values ' phi{CN}])
ylabel('Period T')

save(['C:\Users\solop\Documents\MATLAB\ODE\COMMODE\' ...
handles.DNK.UserData{1}],'SpecInpD','GenInpD');