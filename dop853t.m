% Output f = y(n,CIVI), where n - numel of t_
% CIVI - Chaning Initial Value Index
% (tM,yM) - max point in discret solution y(:,CIVI)
function f = dop853t(t,CN,tM,yM,F,op,Farg)
global yTA  % y at Period T Accuraced
if t == tM, f = yM(CN);
else      
   op.InitialStep = (t-tM)*0.5;
   warning('off','all');
   [t_, y] = dop853(F,[tM t],yM,op,Farg);
   warning('on','all');
   yTA = y(end,:);
   f = yTA(CN);end
%{
function st = solt(met,F,tM,t,yM,CN,op,varargin)
global yTA  % y at Period T Accuraced
if t == tM, st = yM(CN);
else      
   op.InitialStep = (t-tM)*0.5;
   warning('off','all');
   [t_, y] = met(F,[tM t],yM,op,varargin{:});
   warning('on','all');
   yTA = y(end,:);
   st = yTA(CN);end
%}