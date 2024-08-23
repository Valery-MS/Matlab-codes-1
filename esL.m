                     % exact solutions for Linear ODEs
function es = esL(t,A,j,k)
t = t(:);
if j == 1                                % N3_43  homogeneous
   J0 = besselj(0,t); Y0 = bessely(0,t); 
   es = J0.*Y0;
   if k ~= 1
      J1 = besselj(1,t); Y1 = bessely(1,t);  JY = J0.*Y1+J1.*Y0;
      es = [es, -JY 2*(J1.*Y1-J0.*Y0)+JY./t]; end  

elseif j == 2                           % N4_35  homogeneous
    Lt   = log(t); 
    a    = A(1);    b = A(2);  ab = a+b*Lt;       
    c    = A(3);    d = A(4);  cd = c+d*Lt;    % d = A(4)/2;   
    s5   = sqrt(5);     
    m    = (-1-s5)*0.5; 
    n    = (-1+s5)*0.5;   % n-m = s5
    tn_m = t.^s5;  
    tm_3 = t.^(-0.5*(7+s5)); 
    tm_2 = tm_3.*t;
    tm_1 = tm_2.*t;
    es   = tm_1.*t.*(ab + tn_m.*cd);
    if k ~= 1  
       M1 = m*(m-1); M =M1*(m-2);   
       N1 = n*(n-1); N =N1*(n-2);
       es = [ es, ...       
             (tm_1.*(b            +  m*ab + tn_m.*(d            +  n*cd))), ...
             (tm_2.*(b*(2*m-1)    + M1*ab + tn_m.*(d*(2*n-1)    + N1*cd))), ...
             (tm_3.*(b*(3*m^2-6*m+2)+M*ab + tn_m.*(d*(3*n^2-6*n+2)+N*cd)))];end 
     
elseif j == 3                          % N8_57  homogeneous 
   a = A(1);  b = A(2);  c = A(3);
   sqr = sqrt((2*a+c)^2+4*b^2);        
   al = 0.5*(c+sqr);  at = al*t;  ca = cos(at);  sa = sin(at);
   be = 0.5*(c-sqr);  bt = be*t;  cb = cos(bt);  sb = sin(bt);          
   es = ca+sa+cb+sb;
   if k ~= 1
      ga = (a+al)/b;   de = (a+be)/b; 
      es = [es, sa-ca+sb-cb, ga*(sb+cb)+de*(sa+ca), ga*(-cb+sb)+de*(-ca+sa)];end
  
elseif j == 4                         % N3_46     non homogeneous
   Lt = log(t); ot = 1./t;  
   b  = A(2);   ab = A(1)+b*Lt; 
   es = ab.*ot + 0.25*t.*(Lt-2) + A(3);
   if k ~= 1
      ot2 = ot./t; 
      es = [es, (b-ab).*ot2+0.25*(Lt-1), (-3*b+2*ab).*(ot2.*ot)+0.25*ot ]; end 
  
elseif j == 5                         % Nord6     non homogeneous
   es = (t-1)./(t+1); 
   if k > 1 
      At = A(1)*t; s = sin(At);    
      es = [ es, s, cos(At), exp(-A(2)*t).*s, exp(s-1), A(3)*atan(t) ]; end
end
%{
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

   u1 = e.*c;           
   u2 = e.*s;           
   u3 = c./t1;          
   u4 = o2;             
   u5 = A(3)-atan(t);    
   u6 = -c;     

   c = cos(A(1)*t);   A1s = A(1)*s;   e  = exp(-A(2)*t); 
   s = sin(A(1)*t);   A4t = A(4)*t;   t1 = A(3)*t+1;      
   
   u1 = e.*c;            du1 =-e.*(    A1s + A(2)*c);  
   u2 = e.*s;            du2 = e.*(-A(2)*s + A(1)*c);  
   u3 = c./t1;           du3 =   -(    A1s + A(3)*c./t1)./t1;    
   u4 = 1./(A4t.^2+1);   du4 =-2*A(4)*A4t.*u4.^2;   
   u5 = pi/2-atan(A4t);  du5 =-A(4)*u4;        
   u6 = -c;              du6 = A1s;        

   s = sin(t);  c = cos(t);  e = exp(t);  t1 = t+1;  ot = 1./(t.*t+1);    
   
   u1 = A(1)*s.*e;     du1 = A(1)*e.*(c-s);  
   u2 = A(2)*c.*e;     du2 =-A(2)*e.*(c+s);
   u3 = A(3)*s./t1;    du3 = A(3)*(c-s./t1)./t1;  
   u4 = A(4)*c./t1;    du4 =-A(4)*(s+c./t1)./t1;       
   u5 = A(5)*ot;       du5 =-2*t.*u5.*ot;        
   u6 = A(6)-atan(t);  du6 = ot;             % A(6) = pi/2

   s  = sin(t); ss  = s.*s;   ss1 = ss+1;     os2 = 1./ss1; ot2 = 1./(t.*t+1); 
   c  = cos(t); s2  = s-2;
   t1 = t+1;    ot1 = 1./t1;  L1  = log(t1);    
   
   u1 = A(1)*log(s2);     du1 = A(1)*c./s2;  
   u2 = A(2)*cos(L1);     du2 =-A(2)*sin(L1).*ot1;
   u3 = A(3)*atan(t);     du3 = A(3)*ot2;
   u4 = A(4)*atan(s);     du4 = A(4)*c.*os2;       
   u5 = A(5)*exp(c);      du5 =-u5.*s;        
   u6 = A(6)*ot2;         du6 =-2*t.*ot2.^2;

   u1 = A(1)*log(ss1);     du1 = A(1)*s2.*os2;  
   u2 = A(2)*cos(L1);      du2 =-A(2)*sin(L1).*ot1;
   u3 = A(3)*atan(t);      du3 = A(3)*ot2;
   u4 = A(4)*atan(s);      du4 = A(4)*c.*os2;       
   u5 = A(5)*exp(c);       du5 =-u5.*s;        
   u6 = A(6)*exp(s);       du6 = u6.*c;
%} 