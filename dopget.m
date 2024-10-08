function o = dopget(options,name,default)
% ---------------------------
%  DOPGET Get GNI OPTIONS parameters.
%  VAL = DOPGET(OPTIONS,'NAME') extracts the value of the named property
%  from the options structure OPTIONS. The syntax is the same as for ODEGET.
%
%   See also DOPSET DOP54 DOP54D DOP853 DOP853D.
% ---------------------------

if nargin < 2
  error('Not enough input arguments.');
end

if ~isempty(options) & ~isa(options,'struct')
  error('First argument must be an options structure created with GNISET.');
end

if isempty(options)
  if nargin == 3
    o = default;
  else
    o = [];
  end
  return;
end

Names = ['RelTol            '; ...
         'AbsTol            '; ...
         'NormControl       '; ...
         'NonNegative       '; ...
         'OutputFcn         '; ...
         'OutputSel         '; ...
         'DenseOutputFcn    '; ...
         'DenseOutputSel    '; ...
         'Stats             '; ...
         'InitialStep       '; ...
         'MaxStep           '; ...
         'Events            '; ...
         'Mass              '; ...
         'MStateDependence  '; ...
         'MaxIter           '];

[m,n] = size(Names);
names = lower(Names);

lowName = lower(name);
j = strmatch(lowName,names);
if isempty(j)               % if no matches
  error(sprintf(['Unrecognized property name ''%s''.  ' ...
                 'See GNISET for possibilities.'], name));
elseif length(j) > 1            % if more than one match
  % Check for any exact matches (in case any names are subsets of others)
  k = strmatch(lowName,names,'exact');
  if length(k) == 1
    j = k;
  else
    msg = sprintf('Ambiguous property name ''%s'' ', name);
    msg = [msg '(' deblank(Names(j(1),:))];
    for k = j(2:length(j))'
      msg = [msg ', ' deblank(Names(k,:))];
    end
    msg = sprintf('%s).', msg);
    error(msg);
  end
end

eval(['o = options.' Names(j,:) ';']);

if (nargin == 3) & isempty(o)
  o = default;
end
