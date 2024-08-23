% Rms-error for solution at Accuraced Period
function erA = FractErr(y1_qPA,CIVI,TiC,h,F,yTiC,Farg)
warning('off','all');
op = dopset('RelTol',eps,'AbsTol',eps,'InitialStep',h);
qPA = numel(y1_qPA);
t = TiC+h*(1:qPA);
[t1, y] = dop853(F,t,yTiC,op,Farg);
warning('on','all'); 
erA = sqrt(sum( ( 1-y(:,CIVI)./y1_qPA ).^2 )/qPA);
