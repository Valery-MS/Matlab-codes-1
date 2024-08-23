% 08/11/20 from ASPr comparison of search of mins
   tic
   rt = err(y(:,CN),iv);
   [m nm] = mink(rt,2);
   N1=nm(1);
   N2=nm(2);
   if abs(N2-N1) <= 2
      N2=nm(3);
      if abs(N2-N1) <= 2, N2 = nm(4); end,end 
  toc

%        Fragments of programs
% 29.10.20 from APS
   if Tol >= 0.1         
      warning('off','all');                          % for dop853, ode113
      Tfj = Tj :h: Tj+(wf-1)*h;                      % T-fractal
      [Tfj,yf] = SPmet( @F_DNK, Tfj, yTA, SPop, a); 
      warning('on','all');                           % for dop853, ode113
      ert = max(abs(y(nf,CN)-yf(nf,CN)))/abiv;  end

% 15.01.20
%{           
         warning('off','all') % met(@RHFN,t,y0I,op,Arg,j);
         Tj1 = fminbnd(@(t) -solt(@dop853,  @F_DNK,tM,t,y(M,:)',CN,op_,a),...
               t1,t2,optimset('TolX',eps)); yTA1=yTA; 
         Tj2 = fminbnd(@(t) -solt(@ode113_c,@F_DNK,tM,t,y(M,:)',CN,op_,a),...
               t1,t2,optimset('TolX',eps)); yTA2=yTA; 
         Tj3 = fminbnd(@(t) -solt(@ode45_c, @F_DNK,tM,t,y(M,:)',CN,op_,a),...
               t1,t2,optimset('TolX',eps)); yTA3=yTA;    
         warning('on','all');         
%}
         %[tf1, yf1] = met(@F_DNK,Tj:h:Tj+(neps-1)*h,',op,a); 
         % yf = DSm(@F_DNK, Tj,wf, yTA', h, Lmn, b, @F_DNK_M, a)'; ?? wf
         
% 13.01.20
   if GRAFIKI
      fig=figure('Position',SZ,'Name',meth,'NumberTitle','on');
      p = p+1;  Figs(p) = fig.Number;  end
figure(fig)

% 12/01/20
            % Исключаются соседние(если есть) с минимумами точки
            % тк они указывают на те же минимумы только с большей погрешностью  
      [a_, n_] = sortrows(am);
      
      q = 0;   J = Nopp+m_(n_);  I = J;     %  1st option
      while numel(J) > 0                    %
         q = q+1;    I(q) = J(1);           %   
         K = J(2:end);
         J = K( abs(K-J(1)) > Nopp ); end
     
     
      q = 1;   J = Nopp+m_(n_);  I = J;     %  2nd option
      while true                            %
         K = J(2:end);                      %
         J = K( abs(K-J(1)) > Nopp );
         if numel(J) == 0, break,end      
         q = q+1; I(q) = J(1);   end
      
      Ms = sort(I(1:q));         
%    The first option was chosen 


% 20.09.19
if s==1, BadM = {'NaN' 'NaN'};
               else     BadM = {sprintf('E^*=10^{%.2g}  T^*=%.2g',...
                               log10(EA(i,s-1,r)), TA(i,s-1,r))...
                               sprintf('nf^*=%d',nA(i,s-1,r))};end
               nB{r} = [nB{r} i]; 
%5 16.9.18
if k == 4
      for i = 1:n, GG{i} = G(:,i:n:i+(m-1)*n); end 
      g = d-[sum([ t1; -sos; ot1; -tt;  es;   L1 ].*u);...
         sum([ c;  -cs1; at;  cL;   or1;  o21 ].*u);...
         sum([-r1; -es; -sL;  ss;  -as;  -ces ].*u);...
         sum([-o21; r1;  s;   ces; -Ls;  -ot1 ].*u);...
         sum([ Ls; -at;  as;  t;   -ss;   sos ].*u);...
         sum([-cL;  tt; -L1; -or1;  cs1;  sL ].*u)]; 
     
   else
      for i = 1:n
         GG{i} = G(:,i:n:i+(m-1)*n); 
         g(:,i) = d(:,i)-GG{i}*u(:,i); end,end

n = 1e6;t=1:n; g=Gs{1};G=cell(1,n); H=G;
for j=1:2 
  tic; for i=1:n,G{i}=g(t(i),[]);end;t1=toc;
  tic; b=Gs1(t);for i=1:n,H{i}=b(:,i:n:i+2*n);end;t2=toc;
  fprintf(' %.3g   %.3g  ',t1,t2); end
%for j=[ 1 4 78 346],sum(sum(G{j}-H{j})),end  

n=1e4;t=1:n;y=ones(6,1);
for j=1:2
tic;for i=1:n,o1=Nord6(t(i),y,1);end;t1=toc;
tic;for i=1:n,o0=Nord6(t(i),y,0);end;t0=toc;
fprintf('  %.3g  %.3g   ',t1,t0);
end

% 4
ES1s = { @(t) besselj(0,t).*bessely(0,t),             ... % N3_43
         @(t) 1-0.5*t+1./t+log(t)*(0.25*t+1./t),      ... % N3_46
         @(t) t.^sqr1+t.^sqr2.*log(t),                    ... % N4_35h
         @(t) cos(al*t)+sin(al*t)+cos(be*t)+sin(be*t),... % N8_57
         @(t) log(1+sin(t).^2)};                          % Nord6  
     
t0 = t0fs(1,1);
   J0 = besselj(0,t0); Y0 = bessely(0,t0); 
   J1 = besselj(1,t0); Y1 = bessely(1,t0);  JY = J0*Y1+J1*Y0;
y01 = [J0*Y0; -JY; 2*(J1*Y1-J0*Y0)+JY/t0];
y02 = [1.5; -0.25; -0.75];
   sqr5 = sqrt(5); sqr1 = (sqr5-1)*0.5; sqr2 = -(sqr5+1)*0.5;
y03 = [1;  0.5*(sqr5+1); -2*sqr5; 8*sqr5+2]; 

% 3  03.12.16 Nfsk_Solve
             jm = jm+1;  Tm(jm,1) = i;
             Imu = find(mult);
             qm = numel(Imu);
             if qm==1, Tm(jm,2) = rdr(Imu);
             else      Tm(jm,2) = sum(rdr(Imu).*1000.^(qm-1:-1:0)); end        
             continue,end,end
 
 %U1=18:u:20;  U2=23:u:25; U3=28.5:u:31.5; U4=32:u:34; U5=57.5:u:58.5;
%U6=60.5:u:62.5; U7=71.5:u:78;
%ivd = [U1 U2 U3 U4 U5 U6 U7];


% 2
      % from Task_Solve
      if PerCheck == 0, TT(i,:) = TiA; break; end  % exit from wile 1
         
      TiC = TiA(1);   % True Period T Checked
                         % Ti Accuracing 
      while 1 
         %  !!!  Проверять надо не только Т/2, но и Т/3, Т/4 и т.д.
         % => надо менять весь алгоритм поиска периода:
         % начинать с малых Т и двигаться направо,
         % а не так, как сейчас - от большого Т налево
         % !!! эту прогу оставить (может, ещё пригодится)
         % переделывать копию
         %           y1_qPA,CIVI,TiC,h,F,yTiC,Farg
         erA = PerEr(y1_qPA,CIVI,TiC,h,@DNKf,yTA,a); % on qPA points(-?)
         if     erA < er0 && TiC <= TiA, TiC = TiC/2;   
         elseif erA < er0,               fnd = true;  break
         elseif TiC <= TiA,              fnd = true;  break
         else              t9 = cT*TiC; fnd = false; break; end;end 
     
      if fnd, TT(i,1) = TiC;  TT(i,2) = TiA(2);break 
      else    continue; end
      
      if PerCheck  
         figure(fig1) 
         ni = round(TiC/h_);  hi = TiC/ni;  ni1 = ni+1;  ni2 = 2*ni+1;
         op = dopset('RelTol',RelT,'AbsTol',AbsT,'InitialStep',hi);
         warning('off','all')  
         [tDP, yT] = dop853(@DNKf,[t0 TiC],y0i,op,a);
         warning('on','all');
         if NuMet==4 
            op1 = nDSset1(RelT,AbsT,hi, @DNKf,t0,y0i,a,m);  
            y = DSm(@DNKv,ni2,op1,op2)'; % Transp of y respect to dop853
         elseif NuMet==3, errordlg('GBS is not yet ready');
         else
            t2T = t0:hi:2*TT(i,1);
            op = metset('RelTol',RelT,'AbsTol',AbsT,'InitialStep',hi);
            warning('off','all') 
            [t, y] = met(@DNKf,t2T,y0i,op,a); warning('on','all');end
        
         dlt = y(1:ni1,:)-y(ni1:ni2,:);   
         er = sqrt(sum(dlt.^2)/ni1);  
         tT = t0:hi:TiC;
         
         if PerCheck == 1
            subplot(6,1,s); plot(tT,dlt(:,CIVI));
            title(sprintf('%s=%s er=%g', phiCIVI,fi0,er(CIVI)))
         elseif PerCheck == 2
            for k=1:2
               subplot(6,1,2*(s-1)+k); plot(tT,dlt(:,1:2));
               if k==1,title(sprintf('%s=%s er=%.2g %.2g',...
               phiCIVI, fi0, er(1:2)));end,end
         elseif PerCheck == 3
            for k=1:4
               subplot(6,1,4*(s-1)+k); plot(tT,dlt(:,k));
               if k==1,title(sprintf('%s=%s er=%.2g %.2g',...
               phiCIVI, fi0, er));end,end,end,end

% 1
function ARCHIVE_COMMODE
% 19.07.16 --- Executes on button press in Fig.
function Fig_Callback(hObject, eventdata, handles)
global y0 kt tt y12 
figure
den = round(pi/y0(1));
zer = zeros(size(tt));
subplot(3,1,1);  plot(tt,y12(:,1),'b',tt,zer,'k'); 
meth = handles.Methods.String{handles.Methods.Value};
dkm  = sprintf('%d, kt=%d  %.0s',den,kt,meth);
title({['\phi_{10}=\pi/' dkm] '\phi_1'})
subplot(3,1,2);  plot(tt,y12(:,2),'r',tt,zer,'k');    title('\phi_2')
subplot(3,1,3);  plot(tt,y12, tt,zer,'k'); title('\phi_1, \phi_2')
xlabel('t'); 
title('\phi_1, \phi_2');
