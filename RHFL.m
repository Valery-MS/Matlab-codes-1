% Right hand functions F of Linear ODEs y'=F(y)=G(t)*y+g from Kamke       
% This function is used in dop853, ode113, ode45, odex, etc ... 

function F = RHFL(t,y,A,j)  

                      % Homogeneouse ODEs y'=G(t)*y (j = 1,2,3 )
if     j==1, F = [y(2); y(3); -(3*y(3)+4*y(1))/t-(4+1/t^2)*y(2) ];    % N3_43
elseif j==2, F = [y(2); y(3); y(4); -8*y(4)/t-12*y(3)/t^2-y(1)/t^4 ]; % N4_35
    
elseif j==3                                  % N8_57 of Kamke p.540, ÚOB p.16
   ct = A(3)*t;  % A = [ a b c ]
   bc = A(2)*cos(ct);
   bs = A(2)*sin(ct);
   F  = [           A(1)*y(2)  +   bc*y(3) +   bs*y(4); ...
         -A(1)*y(1)            +   bs*y(3) -   bc*y(4); ...
           -bc*y(1) - bs*y(2)              + A(1)*y(4); ...
           -bs*y(1) + bc*y(2)  - A(1)*y(3) ];
       
                      % Non homogeneouse ODEs y'=G(t)*y+g  (j = 4,5 )
elseif j==4,  F = [y(2); y(3); -5*y(3)/t+(log(t)-4*y(2))/t^2 ];  % N3_46
    
elseif j==5                                                  % Nord6 ÚOB p.17
   I  = 1;   % O = 0;  
   
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
   
   u = [  u1;  u2;  u3;  u4;  u5;  u6 ];
   d = [ du1; du2; du3; du4; du5; du6 ];  
   
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
         
   u = [  u1;  u2;  u3;  u4;  u5;  u6 ];
   d = [ du1; du2; du3; du4; du5; du6 ];
   
   Gm =[  o1   e2  u1   I  e1  o2; ...
          e2   e1  -I  o2  o1 -e3; ...
          e1   o1 -e2   c  o2  o1; ...      
          o2  -e2   s  e1  u1  e3; ... 
          u3  -u5  o1   I  u2 -o2; ...   
          u5  -o2   c -u2  u3  e3];
    
   s = sin(A(1)*t);   A1s = A(1)*s;   e  = exp(-A(2)*t); 
   c = cos(A(1)*t);   A4t = A(4)*t;   t1 = A(3)*t+1;      
   
   u1 = e.*c;            du1 =-e.*(    A1s + A(2)*c);  
   u2 = e.*s;            du2 = e.*(-A(2)*s + A(1)*c);  
   u3 = c./t1;           du3 =   -(    A1s + A(3)*c./t1)./t1;    
   u4 = 1./(A4t.^2+1);   du4 =-2*A(4)*A4t.*u4.^2;   
   u5 = A(5)-atan(A4t);  du5 =-A(4)*u4;        
   u6 = -c;              du6 = A1s;           
         
   u = [  u1;  u2;  u3;  u4;  u5;  u6 ];
   d = [ du1; du2; du3; du4; du5; du6 ];
   
   Gm =[  u6  -u1  u2  -u3  u4  -u5; ...
          u1  -u2  u3  -u4  u5  -u6; ...
          u2  -u3  u4  -u5  u6  -u1; ...      
          u3  -u4  u5  -u6  u1  -u2; ... 
          u4  -u5  u6  -u1  u2  -u3; ...   
          u5  -u6  u1  -u2  u3  -u4]; 
        
   s  = sin(t); ss  = s.*s;   ss1 = ss+1;     os2 = 1./ss1; ot2 = 1./(t.*t+1); 
   c  = cos(t); cc  = c.*c;   s2  = 2*s.*c;
   t1 = t+1;    ot1 = 1./t1;  st1 = sqrt(t1); or1 = 1./st1; L1  = log(t1);    
   
   u1 = A(1)*log(ss1);     du1 = A(1)*s2.*os2;  
   u2 = A(2)*cos(L1);      du2 =-A(2)*sin(L1).*ot1;
   u3 = A(3)*atan(t);      du3 = A(3)*ot2;
   u4 = A(4)*atan(s);      du4 = A(4)*c.*os2;       
   u5 = A(5)*exp(c);       du5 =-u5.*s;        
   u6 = A(6)*exp(s);       du6 = u6.*c;

   u = [  u1;  u2;  u3;  u4;  u5;  u6 ];
   d = [ du1; du2; du3; du4; du5; du6 ];

   Gm =[  u1   s  os2 cc-1   u2  du6; ... % u1
          c   u2   ss   u3  du5   -1; ... % u2
        -u3   s2   u3  du4  ot2   u6; ... % u3
         u6  ot1  du3   u4   cc   s2; ... % u4
        or1  du2   u3    1   u5 ss-1; ... % u5
        du1   u4    c -or1    s   u6 ];   % u6

    Gm =[  1  du1  ot1  du6  u3   L1; ...
         c  -L1   u2   u4 or1   u6; ...
       du2  -u3  du4   ss -u5 -du3; ...
       -u6    0    s  du3 -u1 -ot1; ...
        u1  -u2   u5    1 -ss  du1; ...
       -u4 -du6  du2 -or1   0 -du4 ];
%} 
  g = d-Gm*u;   
  F = Gm*y +g;         
end 
