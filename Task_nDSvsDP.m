% Comparison nDS with Dormand-Prince(dop853)

function Task_nDSvsDP(handles)

global t0 tn h1_h2 RelT AbsT nSDE y0
global SpecInpD NEq_   % Neq_ - Number of Equation and other
global t1 t y12 T
global nFeDS 

Fun = handles.EqNam.UserData;
F   = str2func(Fun{2}); 
F_M = str2func(Fun{3});
Farg = SpecInpD{1};
y1_0 = y0(:);
h1 = h1_h2(1);  t1 = (t0:h1:tn)';              % for DP-method
h  = h1_h2(2);  t  = (t0:h:tn)'; n = numel(t); % for nDS-method
                     % DP
e2 = 2*eps;
opDP = dopset('RelTol',e2,'AbsTol',e2,'InitialStep',h1,'Stats','on');
warning('off','all')
tic; [t1, y1, St] = dop853( F, t1, y1_0, opDP, Farg); T1 = toc;
nfDP = St.nfevals;
warning('on','all'); 
y12 = y1(:,1:2); 
                     % NDS 
E = NEq_(1); m = NEq_(2); Nh = NEq_(3); ks = NEq_(4:end); kLmn = ks(1);                     
DSs = {@DS6  @DS8  @DS10};  

switch E % save(handles.EqNam.UserData{2},'BestC','hs');
   case 1, load(handles.Bess.UserData{2},'BestC','hs');
   case 2, load(handles.Airy.UserData{2},'BestC','hs');
   case 3, load(handles.Shrod.UserData{2},'BestC','hs');
   otherwise, errordlg('Wrong value of a file number'); end 

[Lmn, b] = DSet( BestC, m, Nh, kLmn); 
tic; [y, nfDS] = DSs{m}(F, t, y0, h, Lmn, b, F_M, Farg ); T = toc; 

midQ = Qit/(n2-16);
figure
zer1 = zeros(size(t1)); zer = zeros(size(t));
subplot(5,1,1);  plot(t1,y1(:,1),'b',t1,zer1,'k');
title(sprintf('T1=%.2g',T1)) 
subplot(5,1,2);  plot(t1,y1(:,2),'r',t1,zer1,'k'); 

subplot(5,1,3);  plot(t,y(1,:),'b',t,zer,'k'); 
title(sprintf(' T=%.3g, midQ=%.2g',T,midQ))
subplot(5,1,4);  plot(t,y(2,:),'r',t,zer,'k'); 
subplot(5,1,5);  plot(t1,y1(:,1),t,y(1,:), t1,zer1,'k');
title('\phi_1, \phi_2')
xlabel('t'); 

GenInpD = {t0 tn h1_h2 RelT AbsT nSDE y0};
save(handles.DNK.UserData{1},'GenInpD','SpecInpD');