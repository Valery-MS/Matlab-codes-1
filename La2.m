% К-ты отстоят от границы на величину d: hLa < d < 2hLa
% A=4*amu-3 < La< mu^2 1=B
% массив Las содержит числа из интервала [A,B], 
% 1) либо кратные hLa,
% 2) либо кратные одному из элементов массива hs, меньшему hLa,
% 3) либо Las = (A B)/2, если B-A< min(hs)
% amu = abs(mu)


function Las = La2(hLa,amu)
A = 4*amu-3;  B = amu*amu+1;
A0 = floor(A);  
hs = [0.5 0.2 0.1];     % 0.05 0.002 0.001 0.0005];
hs = hs(hs<=hLa);
for h = hs
   M = A0:h:B;
   Las = M(A<M & M<B);
   if ~isempty(Las)
      if numel(Las) > 2, Las = Las(2:end-1); return,end,end,end
Las = [];       

