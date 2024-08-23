                        %  Out of BestCoeffs
function BesTab(BCm,NDS,t0f,k)

L   = length(BCm);
hs  = (t0f(2)-t0f(1))./NDS;
Tab = zeros(L*k+L-1,7);

for s = 1:L
   nr = size(BCm{s},1);
   if nr < k, k = nr; end,end

k1  = k+1;

for s = 1:L
   for j = 1:k
       hB = str2num(num2str([hs(s) BCm{s}(j,4)],'%.4f %.2g'));
       Tab( (s-1)*k1+j,: ) = [ NDS(s) hB BCm{s}(j,[5 1:3]) ];end,end

array2table(Tab,'Var',{'N' 'h' 'Err' 'La','mu' 'nu'})
