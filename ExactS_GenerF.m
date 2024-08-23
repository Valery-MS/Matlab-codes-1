% Building of Exact Solving es, Generation Functiom GM & error function er
function ExactS_GenerF(handles)
global t0 tn y0     % used in Bess, Shrod
global omegAiry     % Airy-SpecInpD     
global nBess        % Bess-SpecInpD   
global Pamp Pfreq   % Shrod-SpecInpD   
global es1 es er GM 

switch handles.EqNam.String
case  handles.Airy.Label          
   alfa = exp(1.0i*pi/3)*omegAiry^(2/3);
   ugly = 3^(1/6)*gamma(2/3)/2.0; % ugly_const 
   es1_ = @(t) real(ugly*(sqrt(3)*airy(0,t)+airy(2,t)));  %Ignore imag
   es2_ = @(t) real(alfa*ugly*(sqrt(3)*airy(1,t)+airy(3,t))); 
   es1  = @(t) es1_(alfa*t);
   es2  = @(t) es2_(alfa*t); 
   es   = @(t) [es1(t); es2(t)];
   er   = @(y,ye,n) sqrt(sum((y(1,:)-ye).^2)/n);
   GM   = @(t) [ 0 1; -omegAiry*omegAiry*t 0 ] ;    
           
case handles.Bess.Label
   es1 = @(t) besselj(nBess,t);
   es2 = @(t) 0.5*(besselj(nBess-1,t)-besselj(nBess+1,t)); 
   es  = @(t) [es1(t); es2(t)];
   er  = @(y,ye,n) sqrt(sum((y(1,:)-ye).^2)/n);
   GM  = @(t) [0, 1.0; (nBess/t^2-1.0), -1.0/t];   % Generator matrix  
   y0  = es(t0)';
   handles.y0.String = num2str(y0);
         
case handles.Shrod.Label
   es1 = @(t) [];  
   er  = @(y,ye,n) abs(y(:,n)'*y(:,n)-1);       
   P_x = [0,1;1,0]; P_y=[0,-1.0j;1.0j,0]; P_z=[1,0;0,-1]; % Pauli matrices 
   GM  = @(t) 1.0j*(0.5*P_z+Pamp*cos(Pfreq*t)*P_x);    %=-i*Hamiltonian
   es  = @(t) solveonce(GM,t,y0'); 
   f_omega = sqrt(Pamp^2 + (Pfreq-1)^2); % flop omega, detuning=Pfreq-1
   %f_amp = Pamp^2 / f_omega^2;           % исп в  QubNrm, QubPlot
   tn  = pi / f_omega;
   handles.tn.String = num2str(tn);
   handles.tn.Enable = 'off';
end 