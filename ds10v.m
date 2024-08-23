% DS of 10th order variable step (m=3,2m+4=10) for NonLinear ODEs y'=F(t,y)
%        see Mathcad-files DS10, DS12, formula (4.5)
% t   = [ t1; t1+h; ...; tn ] column of t-points = 0
% y0  - vector of initial values
% h   - step size of t
% Lmn - 3 base DS coeffs La,mu,nu from alfa = (1 La mu 3nu 0 3nu mu La 1) 
% b   - 5 base DS coeffs b1 b2 b3 b4 b0, beta = (b4 ... b1 b0 b1 ... b4) 
% F_M - matrix-function for F
% varargin - parameters for calc F 
% zam(1:2) - способы итераций
% zam(3) = 1/2 for dop853/ode113
%                             NOTES
% Против DS31na_1, DS31na_2 эта версия программы самая быстрая(Arh_COMMODE)
% if varargout < 0, accuracy is not achieved

% Idea - how to solve the stiff ODEs:
%   1) if error become large, remember D=|y(i)-y(i-1)|, which is too large 
%   2) go one step back, halve(уменьшать наполовину) the main step
%      h->h/2,  h/2->h/4 and repeate calcs
%   3) if solution in point i is found, we compare
%      if D-|y(i)-y(i-1)| < 0, h->2h, h/2->h
                              
function [t1,y, varargout] = ds10v(F,t,y0,op,varargin )  

h=op{1}; RT=op{2}; AT=op{3}; Lmn=op{4}; b=op{5}; zam=op{6};
t1   = NaN;              % t1 is dummy, used only for according with odex, etc
F_M  = str2func([func2str(F) '_M']);     % The function @F_M must already exist
  
t_n0 = t(1:16);   % m = 1,2,3;   4*m+4 = 8, 12, 16
M    = size(y0,2);    % number of coordinates of y
nfs  = 0;

if isvector(y0)
   nm = zam(3);
   setms = { @dopset @odeset};  
   mets  = { @dop853 @ode113_c};         
   warning('off','all')
   opt = setms{nm}('InitialStep',h,'RelTol',eps,'AbsTol',eps);   
   [t_, y0, nf] = mets{nm}( F, t_n0, y0, opt, varargin{:} ); 
   nfs = nf(1);
   warning('on','all');  end
   
n = numel(t);     % number of t-points 
y = nan(M,n);  f = y;
y(:,1:16) = y0';

hb = h*b;    hb2 = 2*hb;
H1 = hb2(1); H2  = hb2(2); H3 = hb2(3); H4 = hb2(4); H5 = hb2(5);
h1 = hb(1);  h2  = hb(2);  h3 = hb(3);  h4 = hb(4);  h5 = hb(5);
La = Lmn(1); mu = Lmn(2);  nu3 = 3*Lmn(3); 

f(:,1:16) = F_M( t_n0, y0, varargin{:})';    
qitM0 = 20;   qitMM = 200;
AT1 = 8*AT;  RT1 = 8*RT;

for i = 9:n-8 
   s2 = La *(y(:,i-2)-y(:,i+2)) + mu*(y(:,i-4)-y(:,i+4)) + ...
        nu3*(y(:,i-6)-y(:,i+6)) + y(:,i-8) + ...
        H1 *(f(:,i-2)+f(:,i+2)) + H2*(f(:,i-4)+f(:,i+4)) + ...
        H3 *(f(:,i-6)+f(:,i+6)) + H4*f(:,i-8) + H5*f(:,i); 

   s1 = La *(y(:,i+3)-y(:,i+5)) + mu*(y(:,i+2)-y(:,i+6)) + ...
        nu3*(y(:,i+1)-y(:,i+7)) + y(:,i) + ...
        h1 *(f(:,i+3)+f(:,i+5)) + h2*(f(:,i+2)+f(:,i+6)) + ...
        h3 *(f(:,i+1)+f(:,i+7)) + h4*f(:,i) + h5*f(:,i+4); 
    
   t8 = t(i+8);
   %f8  = (s1-s2)/h4;   yC = h4*fP+s1;
   %yP = 2*s1-s2;   fP = F( t8, yP, Farg{:});   % y,f in Previous step  
   %nevP = yP-h4*fP-s1;  nrmP = sum(abs(nevP)); % vector-neviazka, its norm

   yC   = 2*s1-s2;  fC = F( t8, yC, varargin{:} );  % y,f in Current step
   nevC = abs(yC-h4*fC-s1);  %nrmC = sum(abs(nevC)); % in Prev & Current steps   
   nfs  = nfs+1; 
   qit  = 0;   ay_ = 1;  ay = abs(yC);
   qitM = qitM0; 
   rfn  = 1;
   
   while any( ay_ > AT )    && any( nevC > AT )   &&  ...
         any( ay_ > RT*ay ) && any( nevC > RT*ay )  
      qit = qit+1;
      if qit == qitM
         if all( ay_ < AT1 )    || all( nevC < AT1 )   ||  ...
            all( ay_ < RT1*ay ) || all( nevC < RT1*ay ), break
         elseif  qitM < qitMM,  qitM = qitMM; 
        else      break; end
      end
      yCold = yC;     
           
      if     rfn==1,   yC = h4*fC + s1; 
      elseif rfn==2     
         nevP = yP-h4*fP-s1; 
         yP_C = yP-yC;
         nn   = nevC-nevP;       % block Conf Alushta p.18-19 
         I = nn ~= 0;
         if     all(I), yC    = yP    + yP_C   .* nevP   ./nn; 
         elseif any(I), yC(I) = yP(I) + yP_C(I).* nevP(I)./nn(I); 
                        d = 1-yP(I)./yC(I);  
                        d( d < eps ) = eps;
                        if numel(d)~=sum(~I), d = norm(d,2); end
                        yC(~I) = yC(~I).*(1+d);         
         else yC = h4*fC + s1; end,end
           
      fP   = fC;   yP = yCold;  nevP = nevC;  
      fC   = F( t8, yC, varargin{:}); 
      nevC = abs(yC-h4*fC-s1);   %nrmC = max(abs(nevC));
      nfs  = nfs+1;

      if any( nevC > nevP )
         rfn = zam(rfn);                    % return to P as nrmP < nrmC
         if rfn == 1, fC = fP; yC = yP; end,end
     
      ay_ = abs(yC-yP); ay = abs(yC);  end  % not to P
    
   if qitM >= qitMM,    break,end  % divergence

   f(:,i+8) = fC;
   y(:,i+8) = yC;  
end

y = y';   
if nargout==3, varargout{1} = sign(qitMM-qit-0.5)*nfs; end
