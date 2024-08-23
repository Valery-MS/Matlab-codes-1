% К-ты отстоят от границы на величину d: hLa < d < 2hLa
% A=Ldm=Lq-L0 <La< Ldp=Lq L0=B
% массив Las содержит числа из интервала [A,B], 
% 1) либо кратные hLa,
% 2) либо кратные одному из элементов массива hs, меньшему hLa,
% 3) либо Las = (A B)/2, если B-A< min(hs)

function Las = La3(hLa,mu,nu)
p3 = (2-mu)/3; 
if nu==0,  B = p3^1.5;  A = -B;
else
   nu2 = nu*nu;  L0 = (p3+nu2)^1.5;
   Lq  = abs(nu)*(1+mu-2*nu2);
   B = Lq+L0;  A = Lq-L0;  end

A0 = floor(A);  
hs = [0.5 0.2 0.1];   % 0.05 0.002 0.001 0.0005];
hs = hs(hs<=hLa);
for h = hs
   M = A0:h:B;
   Las = M(A<M & M<B);
   if ~isempty(Las)
      if numel(Las) > 2
         Las = Las(2:end-1);
         if nu<0, Las = -flip(Las); end
         return,end,end,end
Las = [];