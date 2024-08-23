% Comparison NDS with Dormand-Prince(dop853)
function Task_NDSvsDP(handles)
global a1 y0 t0 tn h kt RelT AbsT  
global t1 t2 y12 a T2 vEmN
global Q 

y1_0 = y0(:);
a  = a1*1e-4;
h1 = h(1);  t1 = (t0:h1:tn)';                 % for DP-method
h2 = h(2);  t2 = (t0:h2:tn)'; n2 = numel(t2); % for NDS-method

                  % DP
e2 = 2*eps;
op = dopset('RelTol',e2,'AbsTol',e2,'InitialStep',h1);
warning('off','all')
tic;  [t1, y1] = dop853(@DNKf,t1,y1_0,op,a); T1 = toc;
y12 = y1(:,1:2); 

                  % NDS
cDS1 = {@cDS11   @cDS21   @cDS31};
DS1  = {@DS11na  @DS21na  @DS31na};      
vr = vEmN(1); E = vEmN(2); m = vEmN(3); Nh = vEmN(4); ks = vEmN(5:end);
switch E
   case 1, load(handles.Bess.UserData{2},'BestC','hs');
   case 2, load(handles.Airy.UserData{2},'BestC','hs');
   case 3, load(handles.Shrod.UserData{2},'BestC','hs');
   otherwise, errordlg('Wrong value of a file number'); end 

                % NDSset
Cmh  = BestC{m}{Nh};
cDSm = cDS1{m};   DSm = DS1{m};
Lmn  = Cmh(ks,1:m); b = cDSm(Lmn);
t2_0 = (t0:h2:(4*m+3)*h2)';
op = dopset('RelTol',eps,'AbsTol',eps,'InitialStep',h2);  

               % initial data computing by DP-method for NDS-method       
warning('off','all')
tic;  [t, y2_0] = dop853(@DNKf,t2_0,y1_0,op,a);  T20 = toc;
               % NDS-method
%tic
y2_0 = y2_0';
%for k = ks  
   y2 = DSm(n2,y2_0,h2,b,@DNKf,Lmn,vr,a); 
%end
%T2 = toc; 
T=T20+T2;
%y22=y2;
y22= y2(1:2,:)';

midQ = sum(Q)/(n2-16);
figure
zer1 = zeros(size(t1)); zer2 = zeros(size(t2));
subplot(5,1,1);  plot(t1,y12(:,1),'b',t1,zer1,'k');
title(sprintf('T1=%.2g kt=%d',T1,kt)) 
subplot(5,1,2);  plot(t1,y12(:,2),'r',t1,zer1,'k');    
subplot(5,1,3);  plot(t2,y22(:,1),'b',t2,zer2,'k');    
subplot(5,1,4);  plot(t2,y22(:,2),'r',t2,zer2,'k'); 
title(sprintf(' T2=%.2g+%.2g=%.2g, vr=%d, midQ=%.2g',T20,T2,T,vr,midQ)) 
subplot(5,1,5);  plot(t1,y12(:,1),t2,y22(:,1), t1,zer1,'k');
title('\phi_1, \phi_2')
xlabel('t'); 
  
warning('on','all'); 

save(handles.DNK.UserData{1},'a','y0','t0','tn','h','kt','RelT','AbsT');