function [xdot] = InversePODE_R(t,x,A,B,t_pulse)

u = input_d(t,t_pulse);
xdot= A*x+ B*u(t); 

end
