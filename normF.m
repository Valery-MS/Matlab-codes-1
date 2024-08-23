function nF = normF(F,t,y,varargin)
f  = F(t,y,varargin{:});
nF = f*f';