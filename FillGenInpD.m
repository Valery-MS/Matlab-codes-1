% Filling of panel GenInpD: t t0 h RelT AbsT nSDE y0 y0'
function  f = FillGenInpD(handles,GenInpD)
global t0 tn h1_h2 RelT AbsT nSDE y0

t0 = GenInpD{1}; handles.t0.String = num2str(t0); handles.t0.Enable = 'on';
tn = GenInpD{2}; handles.tn.String = num2str(tn); handles.tn.Enable = 'on';
h1_h2 = GenInpD{3}; handles.h1_h2.String = num2str(h1_h2);
                    handles.h1_h2.Enable = 'on';
RelT  = GenInpD{4}; handles.RelT.String  = num2str(RelT,'%.1g');
                    handles.RelT.Enable  = 'on';
AbsT  = GenInpD{5}; handles.AbsT.String  = num2str(AbsT,'%.1g');
                    handles.AbsT.Enable  = 'on';
nSDE  = GenInpD{6}; handles.nSDE.String  = num2str(nSDE);
                    handles.nSDE.Enable  = 'on';
y0    = GenInpD{7}; handles.y0.String    = num2str(y0,'%.3g  ');
                    handles.tn.Enable    = 'on';