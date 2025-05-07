function [xdot] = InversePODE_L(t,x,A,B,pend1,i)

amp = pend1.amp(i);
t_pulse = pend1.t_pulse;

u = input_f(t,amp,t_pulse);
xdot= A*x+ B*u(t); 

end
