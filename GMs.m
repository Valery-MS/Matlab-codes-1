% Generetion matrix G at vector t=[t1,...tn]                 
% This function is used in ds10, ds10Lh, ds10Ln                
                
function [ G, g ] = GMs( t, A, j )
n = length(t);   g = [];   
if n == 1,  O = 0;          I = 1;  
else        O = zeros(1,n); I = ones(1,n);   end

                      % Homogeneous ODEs y'=G(t)*y  ( j = 1,2,3 )
if j == 1                                     %  N3_43 of Kamke p.465
   It = -4./t;                                                %  
   Gm = [O I O; O O I; It, -(4+0.0625*It.*It), 0.75*It ];

elseif j == 2                                 % N4_35 of Kamke p.478  Euler eq
   It = 1./t; It2 = It.^2;
   Gm = [ O I O O; O O I O; O O O I; -It2.^2  O  -12*It2  -8*It ];
  % [y(2); y(3); y(4); -8*y(4)/t-12*y(3)/t^2-y(1)/t^4 ]; % N4_35 
  
elseif j == 3                                 % N8_57 Kamke p.540, ÚOB p.16
   a  = A(1)*I;
   ct = A(3)*t;
   bc = A(2)*cos(ct);
   bs = A(2)*sin(ct);
   Gm = [ O a bc bs;   -a O bs -bc;   -bc -bs O a;   -bs bc -a O ]; 
   
         % Non homogeneous ODEs y'=G(t)*y+g  ( j = 4,5 )
         % t - 1*n-array, varargout = [GG, g] for linear methods, F = G(t)*y+g
         %     Gc - 1*n-cell array of G(ti), g - m*n-double matrix
elseif j == 4                                       % N3_46 of Kamke p.466
   It = 1./t; It2 = It.*It;
   Gm = [ O I O; O O I; O, -4*It2,  -5*It];
   g  = [ O; O; log(t).*It2];
    
elseif j == 5                          % Nord6  ÚOB p.17

   o1 = 1./(t+1);      o2 = 1./(t.^2+1);  
   s  = sin(A(1)*t);   As = A(1)*s;
   c  = cos(A(1)*t);   Ac = A(1)*c;    
   e  = exp(-A(2)*t);  e1 = e.*o1;   e2 = e.*o2;   e3 = exp(-t);          
   
   u1 = (t-1).*o1;     du1 = 2*o1.^2;  
   u2 = s;             du2 = Ac;
   u3 = c;             du3 =-As;    
   u4 = e.*s;          du4 = e.*(-A(2)*s + Ac);  
   u5 = exp(s-1);      du5 = Ac.*u5;    
   u6 = A(3)*atan(t);  du6 = A(3)*o2;        
   
   Gm =[ o1   e2  u1   s  e1  o2; ...  % u1
         e2   e1  -s  o2  o1 -e3; ...  % u2
         e1   o1 -e2   c  o2  o1; ...  % u3     
         o2  -e2   s  e1  u1  e3; ...  % u4 
         u3  -u5  o1   c  u2 -o2; ...  % u5   
         u5  -o2   c -u2  u3  e3];     % u6
    
%{   
   t1 = t+1;          o1 = 1./(t+1);      o2 = 1./(t.^2+1);  
   s  = sin(A(1)*t);  e  = exp(-A(2)*t);  e1 = e.*o1;   e2 = e.*o2;  
   c  = cos(A(1)*t);  e3 = exp(-t);       
   A1s = A(1)*s;
   
   u1 = e.*c;            du1 =-e.*(    A1s + A(2)*c);  
   u2 = e.*s;            du2 = e.*(-A(2)*s + A(1)*c);  
   u3 = c./t1;           du3 =   -(    A1s + c./t1)./t1;    
   u4 = o2;              du4 =-2*t.*u4.^2;   
   u5 = A(3)-atan(t);    du5 =-u4;        
   u6 = -c;              du6 = A1s;           

   Gm =[  o1   e2  u1   I  e1  o2; ...  % u1
          e2   e1  -I  o2  o1 -e3; ...  % u2
          e1   o1 -e2   c  o2  o1; ...  % u3     
          o2  -e2   s  e1  u1  e3; ...  % u4 
          u3  -u5  o1   I  u2 -o2; ...  % u5   
          u5  -o2   c -u2  u3  e3];     % u6
    
   g = [du1-( o1.*u1 + e2.*u2 + u1.*u3 +  I.*u4 + e1.*u5 + o2.*u6 );...
        du2-( e2.*u1 + e1.*u2 -  I.*u3 + o2.*u4 + o1.*u5 - e3.*u6 );...
        du3-( e1.*u1 + o1.*u2 - e2.*u3 +  c.*u4 + o2.*u5 + o1.*u6 );...
        du4-( o2.*u1 - e2.*u2 +  s.*u3 + e1.*u4 + u1.*u5 + e3.*u6 );... 
        du5-( u3.*u1 - u5.*u2 + o1.*u3 +  I.*u4 + u2.*u5 - o2.*u6 );...
        du6-( u5.*u1 - o2.*u2 +  c.*u3 - u2.*u4 + u3.*u5 + e3.*u6 )];  


         % Formula:  g = d  - Gm*u;
%}          
        
  u = [  u1   u2   u3   u4   u5   u6 ];
  d = [ du1; du2; du3; du4; du5; du6 ];      
  g = d -[sum(reshape(Gm(1,:).*u,[n 6])');...
          sum(reshape(Gm(2,:).*u,[n 6])');...
          sum(reshape(Gm(3,:).*u,[n 6])');...
          sum(reshape(Gm(4,:).*u,[n 6])');...
          sum(reshape(Gm(5,:).*u,[n 6])');...
          sum(reshape(Gm(6,:).*u,[n 6])')];         
  

end

if n == 1, G = Gm;
else 
   G = cell(1,n);
   m = size(Gm,1);
   for i = 1:n,  G{i} = Gm(:,i:n:i+(m-1)*n); end,end
