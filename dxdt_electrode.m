function f = dXdT_electrode(t,Jel,tsim,Jox)
%% This function describes the delay in  electrode measurements of VO2
tau = 4; % seconds
f = (interp1(tsim,Jox,t)- Jel)/tau;