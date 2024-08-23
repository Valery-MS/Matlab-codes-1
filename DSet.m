% DSet - parameter Set for ds10 for solving  ODE y'=F(y)

% h    - t-step
% RT_f - Relative Tolerance / func f for dop853/ode113 for linear ODE (ds10L)
% AT   - Absolute  Tolerance 
% Lmn  - coefs Lambda, mu, nu        for which BestC are calculated
% zam  -  zam(1:2) iteration changes; zam(3)=1/2 if dop853/ode113
% cDS(Lmn) - coeffs beta

function op = DSet( h, RT_f, AT, Lmn,zam)
                         %  b    
op = { h,RT_f, AT, Lmn, cDS(Lmn), zam};             