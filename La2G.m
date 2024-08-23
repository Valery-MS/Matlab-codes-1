% выч к-тов Las по mu включая к-ты, близкие к границе(G)

function Las = La2G(hLa,amu)
A = 4*amu-3;  B = amu*amu+1;
A0 = floor(A);  
hs = [0.5 0.2 0.1 0.05 0.002 0.001 0.0005];
hs = hs(hs<=hLa);
for h = hs
   M = A0:h:B;
   Las = M(A<M & M<B);
   if ~isempty(Las), return; end,end
Las = (A+B)/2;