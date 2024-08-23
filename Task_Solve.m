% Solving of DNK equation
% Ќовый на фоне старого

function Task_Solve(handles)

global SpecInpD a NEq_  IV cT_
global GenInpD  t0 tn h1_h2 RelT AbsT nSDE y0
global Figs  yTA   % yTA is defined in dop853t

tic
y0 = y0(:);
h = h1_h2(1);  
NuMet = handles.Methods.Value;

if NuMet == 1
   NEq = NEq_(1); m = NEq_(2); Nh = NEq_(3); ks = NEq_(4:end); kLmn = ks(1);
   switch NEq
      case 1, load(handles.Bess.UserData{2},'BestC','hs');
      case 2, load(handles.Airy.UserData{2},'BestC','hs');
      case 3, load(handles.Shrod.UserData{2},'BestC','hs');
      otherwise, errordlg('Wrong value of a file number');
   end      
   dsm  = { @ds6  @ds8  @ds10};    
   met  = dsm{m};
   zam  = [2 1 2];  % zam(1:2) iteration changes; zam(3)=1/2 if dop853/ode113
   Lmn  = BestC{m}{Nh}(kLmn,1:m);   
   op   = DSet( h, RelT,AbsT,Lmn, zam);
   meth = func2str(met); 
else
   if     NuMet == 2,  met = @dop853; metset = @dopset;
   elseif NuMet == 3,  met = @ode113; metset = @odeset;
   elseif NuMet == 4,  met = @ode45;  metset = @odeset;  
   elseif NuMet == 5,  met = @odex;   metset = @odeset; 
   else   errordlg('Wrong value of a code number');end
    
   op = metset('RelTol',RelT,'AbsTol',AbsT,'InitialStep',h);
   meth = handles.Methods.String{handles.Methods.Value}; end 

                %  Source data
%*************************************************************************
y0i = zeros(4,1);
SZ = [459 230 926 515]; % = X Y Width Height 
RT0 = 0.1;      % max Relative Tolerance of solution y on pattern
RTp = 1e-6;     % Relative Tolerance of solution y in point
kT  = 2;        % the number of smallest periods
CN  = 1;        % Coordinate Number where Initial Value is changed - IV(1)
Ffi = 8; %Flag fi 1)CN 2)1,2 3)1,2,3,4 4)1+2,3+4 5)1,2 vs Lin 6)1 vs Lin
qBs = [5 2 1 2 2 1 4 1];  % Qty of Blocks in Figure for Ffi (array)
qGs = [1 2 4 2 2 4 1 2];  % Qty of Grafics in Block for Ffi (array)
u   = 0.1;
ivd = [3 ]; %5:u:45 48 50:u:90];  % initial values (in degrees if CN<3)
qiv = numel(ivd);  % range of fi-values(in degries)
mp_t = 1;     % multiplier of t
tf  = 6000;    tA=3200; tB=6000;  tik = 500; 
XLi = [tA tB]; t9_= tf/mp_t;      XLi_= XLi/mp_t;  tik_= tik/mp_t;  
GRAFIKI = 1;

MT = 10:17;   % мультипликаторы периодов: различие TL и T

%*************************************************************************
t = (t0:h:tf)'; t_ = t/mp_t; n  = numel(t);  zer = zeros(size(t));
t0f = [t0 tf];

if Ffi >= 5             % for Thai
   d1=a(3); d2=a(4); e1=a(1)-d1; e2=a(2)-d2; ep=e1+e2; em=e1-e2;
   ed = sqrt(em^2+4*d1*d2);
   w1 = sqrt((ep+ed)*0.5);  w2 = sqrt((ep-ed)*0.5); 
   dw = w1-w2;  TL = 2*4*pi/dw; % tz = pi/dw;
   A = [[0 0 1 0]; [0 0 0 1]; [-e1 d1 0 0]; [d2 -e2 0 0]];
   [V, D] = eig(A);
   a_= real(V(1,1)); b_= real(V(2,1)); c_= imag(V(3,1)); d_= imag(V(4,1));
   e_= real(V(1,3)); f_= real(V(2,3)); g_= imag(V(3,3)); h_= imag(V(4,3));
   cow1 = cos(w1*t); cow2 = cos(w2*t);
   af = a_*f_; be = b_*e_; bf= b_*f_;   den = af-be;
   dwt = (w1-w2)*0.5*t;
   og_ = 2/den*[af*cos(dwt) bf*sin(dwt)];   % огибающа€ ф-ции yL_
   yL_ = [(af*cow1-be*cow2),  bf*(cow1-cow2)]/den;
    %yL = [cow1+cow2, -cow1+cow2];                         - неправильно
    %cop = cos((w1+w2)*0.5*t);  - неправильно
   end
   
IB = round((tA-t0)/h)+1:round((tB-t0)/h)+1;  tIB = t(IB);
if CN < 3, iv = ivd*pi/180; end       % range of fi-values(in radians)

if qiv == 1, qB = 1;
else         qB = qBs(Ffi); end  % qty of Grafics  for y plot; - IV(7)

qG = qGs(Ffi);      % qty of Grafics  for y plot; - IV(7)
qGB = qG*qB;
qF = ceil(qiv/qB);  Figs = nan(qF,1);    % quantity of Figures
wf = round(78/h)+1;    % width of fractal
nf = round(20/h)+1:wf; % range of points on fractal
%righf = round(20/h)+1:wf; % rught hand of fractal
%nf = numel(righf);
fi    = {'$\varphi_1$',    '$\varphi_2$',    '$v_1$',    '$v_2$'};
fiCNs = {'$\varphi_{10}$', '$\varphi_{20}$', '$v_{10}$', '$v_{20}$'};
fiCN  = fiCNs{CN};
% NOTE: ivd(1),ivr(1) aren't used to match with TT starting with the no 2
% % cT, IV(6), IV(8) - reservetIB

%cA=cos(AL*t); cB=cos(BE*t); f=[cA+cB, cA-cB];
T  = nan(2*kT-1,qiv); % T(1,1) = 348.5;
Tm = nan(1,qiv);      % T for Min RTj
RT = T;               % Relative Tolerances
i = 0;  p = 0;        % Figure Number

for nF = 1:qF 
   if GRAFIKI
      fig=figure('Position',SZ,'Name',meth,'NumberTitle','on');
      p = p+1;  Figs(p) = fig.Number;  end
     
   for s = 1:qB   % Blocks
      i   = i+1;   
      fi0 = sprintf('%g',ivd(i));
      ivi = iv(i);  y0i(CN) = ivi;  Ey = RTp*ivi;  
  
      if NuMet==1, [to,y,nuf] = met(@RHFN,t,y0i,op,a,4);
      else          
         warning('off','all'); 
         [t, y] = met( @F_DNK, t, y0i, op, a);   % ранее было t0f t
         warning('on','all'); end       
      
      yi = y(:,CN);
      Nopp = find(ivi*yi<0,1); % Number of 1st element opposed to y0
      yi(1:Nopp) = NaN;        % Turn off value y0 from yi
         
      [mi0, I0] = sortrows(abs(yi(Nopp+1:end)/ivi-1));  % Local number
      mi  = mi0(mi0<RT0);    % truncated mi
      qmi = numel(mi);   
      
      q = 0;   J = Nopp+I0(1:qmi);  I = nan(1,qmi);
      while numel(J)>1 && ~isempty(J)
         q = q+1;    I(q) = J(1);               
         K = J(2:end);
         J = K( abs(K-J(1)) > 50 ); end
     
      Ms = sort(I(1:q));      
      j = 0; r = 0; 
      while r<kT && j<q    
         j = j+1;  M = Ms(j);   tM = t(M);
         if r>0
            dr = tM./T(1:r,i);
            if any(abs(dr-round(dr))<0.05), continue,end,end
           
         M_1 = M-1;   tM_1 = t(M_1);  fta = yi(M_1);   
         M1  = M+1;  
         if M == n, tM1  = tM;     ftb = yi(M); 
         else       tM1  = t(M1);  ftb = yi(M1); end
            
         if abs(ftb-fta)<Ey, w=0.75*tM; t1 = 0.25*tM_1+w; t2 = 0.25*tM1+w;
         elseif fta<ftb,     t1 = tM;   t2 = (tM+tM1)*0.5;
         else                t1 = (tM_1+tM)*0.5;  t2 = tM;  end
      
         % sign before dop853t = -sign(y0iC); ”точнение периода
         warning('off','all')
         Tj = fminbnd(@(t) -dop853t(t,CN,tM,y(M,:)',@F_DNK,a),...
              t1,t2,optimset('TolX',eps));  % Period T Accuraced
         warning('on','all');         

         %[tf1, yf1] = met(@F_DNK,Tj:h:Tj+(neps-1)*h,yTA',op,a); 
         % yf = DSm(@F_DNK, Tj,wf, yTA', h, Lmn, b, @F_DNK_M, a)'; ?? wf
         Tfj = Tj :h: Tj+(wf-1)*h;                              % T-fractal
         [Tfj, yf] = met( @F_DNK, Tfj, yTA', op, a);           % y-fractal
         %yf = DSm(@F_DNK, Tfj, yTA', h, Lmn, b, @F_DNK_M, a)'; % y-fractal
         
         RTj = max(abs(y(nf,CN)-yf(nf,CN)))/ivi;
         if RTj<RT0, r=r+1;  T(r,i)=Tj;  RT(r,i)=RTj; end,end
     
      Tm(i) = T(1,i);  if RT(2,i)<RT(1,i), Tm(i) = T(2,i); end
      T(kT+1:end,i) = T(2:kT,i)-T(1:kT-1,i);   
      
      if GRAFIKI
      figure(fig)
      if Ffi == 1         % qGB=5*1, qG=1,  plot  fi(CN)
         ax = subplot(qGB,1,s); plot(t_,y(:,CN),'k') %,t,zer,':k');
         tTT = sprintf('%.5g $\\quad$',T(~isnan(T(1:kT,i)),i));  
         %title(['\textbf{\boldmath{' sprintf('%s  %s',...
         %[fi0 '$^\circ\quad$'],tTT) '}}'],'interpreter','latex'); 
         %title(['\textbf{\boldmath{' sprintf('%s T = %s',[fiCN '=' fi0 ...
         title(['\textbf{\boldmath{' sprintf('%s  %s',[fi0 ...
         '$^\circ\quad$'],tTT) '}}'],'interpreter','latex');
         ylabel(fi{CN},'interpreter','latex');
         set(ax,'XTick',0:tik_:t9_,'XLi',XLi_);
      elseif Ffi == 2  
         for k = 1:qG    % qGB=2*2, plot fi1, fi2 
            ax = subplot(qGB,1,2*(s-1)+k); plot(t,y(:,k),'b',t,zer,'k'); 
            ylabel(fi{k},'interpreter','latex'); 
            set(ax,'XMinorGrid','on','XMinorTick','on','XTick',0:tik:tf);
            if k==1
               tTT = sprintf('%.5g $\\quad$',T(~isnan(T(1:kT,i)),i));
               title(['\textbf{\boldmath{' sprintf('%s T = %s',...
               [fiCN '=' fi0 '$^\circ\quad$'],tTT) '}}'],...
               'interpreter','latex');  end,end
      elseif Ffi == 3
         for k = 1:qG    % qGB=4*1, plot fi1,...fi4
            ax = subplot(qGB,1,k); 
            plot(t,y(:,k),'b',t,zer,'k'); 
            ylabel(fi{k},'interpreter','latex'); 
            set(ax,'XMinorGrid','on','XMinorTick','on','XTick',0:tik:tf);
            if k == 1
               tTT = sprintf('%.5g  ',T(~isnan(T(1:kT,i)),i)); 
               title(sprintf('%s (%s)',[fiCN '=' fi0 '$\quad$'],tTT),...
               'interpreter','latex');  end,end
      elseif Ffi == 4   % qG=2, plot fi1 U fi2, f3 U fi4
         for k = 1:qG  
            ax = subplot(qGB,1,2*(s-1)+k); k2 = k+k;
            plot(tIB,y(IB,k2-1),'b',tIB,y(IB,k2),'k',tIB,zer(IB),':k'); 
            ylabel([fi{k2-1} ',' fi{k2}],'interpreter','latex');
            set(ax,'XMinorGrid','on','XMinorTick','on','XTick',0:tik:tf,...
            'XLi',XLi_ );
            if k == 1
               tTT = sprintf('%.5g $\\quad$',T(~isnan(T(1:kT,i)),i));    
               title(['\textbf{\boldmath{' sprintf('%s  T = %s',...
               [fiCN '=' fi0 '$^\circ\quad$'],tTT) '}}'],...
               'interpreter','latex'); end,end       
      elseif Ffi == 5 % f1+L1, f2+L2 Comparison with linear approx (Thai)
         for k = 1:qG 
            yL = ivi*yL_; 
            ax = subplot(qGB,1,2*(s-1)+k); 
            plot(tIB,[y(IB,k),yL(IB,k)],tIB,zer(IB),':k',...
                [MT*TL MT*Tm(i)],0,'.k'); 
            ylabel([fi{k} ',' fi{k} 'Lin'],'interpreter','latex'); 
            set(ax,'XMinorGrid','on','XMinorTick','on','XTick',0:tik:tf,...
            'XLi',XLi_ ); 
            if k == 1
               tTT = sprintf('%.5g $\\quad$',T(~isnan(T(1:kT,i)),i)); 
               title(['\textbf{\boldmath{' sprintf('%s  T = %s',...
               [fiCN '=' fi0 '$^\circ\quad$'],tTT) '}}'],...
               'interpreter','latex');  end,end
      elseif Ffi == 6 % f1 vs L2 Comparison with linear approx (Thai)
         yL = ivi*yL_;
         tz=0;
         ax(1) = subplot(qGB,1,2*(s-1)+1); 
         T1 = T(1,i);  qT = round(tB/T1);
         plot(tIB,y(IB,CN),'b',tIB,zer(IB),':k',...
            tz+T1*(0:qT),0,'.k',tz+TL*(0:qT),0,'.r');
         ylabel(fi{CN},'interpreter','latex'); 
         tTT = sprintf('%.5g $\\quad$',T(~isnan(T(1:kT,i)),i));        
         title(['\textbf{\boldmath{' sprintf('%s  T = %s',...
         [fiCN '=' fi0 '$^\circ\quad$'],tTT) '}}'],'interpreter','latex'); 
        
         ax(2) = subplot(qGB,1,2*(s-1)+2); 
         plot(tIB,yL(IB,CN),'b',tIB,zer(IB),':k',...
             tz+T1*(0:qT),0,'.k',tz+TL*(0:qT),0,'.r'); 
         title(sprintf('\\omega_1=%g \\omega_2=%g \\TL=%g',w1,w2,TL))
         ylabel('$A(\cos w_1+\cos w_2$)','interpreter','latex');   
      elseif Ffi == 7 % |f1-L1| 
         yL = ivi*yL_; 
         ax = subplot(qGB,1,s); 
         plot(tIB,(y(IB,CN)-yL(IB,CN))/ivi) %,t,zer,':k');
         tTT = sprintf('%.5g $\\quad$',T(~isnan(T(1:kT,i)),i));  
         %title(['\textbf{\boldmath{' sprintf('%s  %s',...
         %[fi0 '$^\circ\quad$'],tTT) '}}'],'interpreter','latex'); 
         %title(['\textbf{\boldmath{' sprintf('%s T = %s',[fiCN '=' fi0 ...
         title(['\textbf{\boldmath{' sprintf('%s  %s',[fiCN '=' fi0 ...
         '$^\circ\quad$'],tTT) '}}'],'interpreter','latex');
         ylabel(fi{CN},'interpreter','latex');
         set(ax,'XTick',0:tik_:t9_,'XLi',XLi_);
      elseif Ffi == 8 % fi1, fi2, L1, L2, og
      for k = 1:qG 
            yL = ivi*yL_; og = ivi*og_;
            ax = subplot(qGB,1,2*(s-1)+k); 
            plot(tIB,[y(IB,k),yL(IB,k),og(IB,k)],tIB,zer(IB),':k',...
                [MT*TL MT*Tm(i)],0,'.k'); 
            ylabel([fi{k} ',' fi{k} 'Lin'],'interpreter','latex'); 
            set(ax,'XMinorGrid','on','XMinorTick','on','XTick',0:tik:tf,...
            'XLi',XLi_ ); 
            if k == 1
               tTT = sprintf('%.5g $\\quad$',T(~isnan(T(1:kT,i)),i)); 
               title(['\textbf{\boldmath{' sprintf('%s  T = %s',...
               [fiCN '=' fi0 '$^\circ\quad$'],tTT) '}}'],...
               'interpreter','latex');  end,end    
      end

      end % if GRAFIKI=1      
      if i == qiv, break,end,end
      %xlabel('\textbf{\boldmath{$t,\times 10^3$}}','interpreter','latex');
      xlabel('$t,[\times 1000],c$','interpreter','latex') %xlabel('t,');  
      end

ns = ones(14,1);  r=1;  %  bounds of periods change  
for q=1:qiv-1
   if abs(T(1,q+1)-T(1,q))>10, r=r+1; ns(r)=q+1; end,end
ns(r+1)=qiv;

%{
figure('Name',sprintf('RelTol=%g',RT0))    %,'NumberTitle','off')
%ax1=subplot(2,1,1);
plot(ivd,T(1,:));
set(gca,'XGrid','on','YGrid','on','XMinorTick','on','TickDir','out');
%title(sprintf('1-й период  RT0=%g',RT0));
ylabel('T1'); xlabel(fiCN,'interpreter','latex');
%ax2=subplot(2,1,2);
figure
plot(ivd,T(2,:));
set(gca,'XGrid','on','YGrid','on','XMinorTick','on','TickDir','out');
%title('2-й период');
ylabel('T2');xlabel(fiCN,'interpreter','latex');
%ax3=subplot(3,1,3); plot(ivd,Tm);     title('комб период');ylabel('Tm')
%}

TA = toc; TA_qiv_T1 = [TA, qiv, TA/qiv]

figure
plot(ivd,T(1,:),'*k',ivd,T(2,:),'ob');
set(gca,'XGrid','on','YGrid','on','XMinorTick','on','TickDir','in');
ylabel('T_1, T_2'); xlabel(fiCN,'interpreter','latex')

%{
i=0; 
for nF=1:4
   figure
   s=0;

   while s<3  &&  i<r
      i=i+1; ii=ns(i):ns(i+1)-1;
      if ns(i)<ns(i+1)-1
         s=s+1; ax(s)=subplot(3,1,s); 
         plot(ii,T(1,ii)); 
         ylabel('T1');  end,end
   xlabel(phiIVI(1:end-3));
   set(ax,'XGrid','on','YGrid','on','XMinorTick','on');end
%}

%{
nt1=numel(U1);     nt2=nt1+numel(U2); nt3=nt2+numel(U3); nt4=nt3+numel(U4);
nt5=nt4+numel(U5); nt6=nt5+numel(U6); nt7=nt6+numel(U7);  
N1=1:nt1; N2=nt1+1:nt2; N3=nt2+1:nt3; N4=nt3+1:nt4; N5=nt4+1:nt5;
N6=nt5+1:nt6; N7=nt6+1:nt7;
plot(U1,W(N1),U2,W(N2),U3,W(N3),U4,W(N4),U5,W(N5),U6,W(N6),U7,W(N7));
set(gca,'XGrid','on','YGrid','on','XMinorTick','on','TickDir','out');
ylabel('T1, T2'); xlabel(fiCN,'interpreter','latex')
%}
if exist('fig2','var'),set(fig,'Name',sprintf('%s  T=%.3g',meth,Tim));end

save(['C:\Users\solop\Documents\MATLAB\ODE\COMMODE\' ...
handles.DNK.UserData{1}],'SpecInpD','GenInpD');