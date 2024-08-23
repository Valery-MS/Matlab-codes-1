% !! It's change to according with other programms

function [t1,y,nf] = ode113_c(F,t,y0,op,varargin)
[t1,y,nf] = ode113(F,t,y0,op,varargin{:});
nf = nf([ 3 1 2 ]);