function [xdot]= pendulumODE1(t,x, gamma, w_n,l_n, input_pendulum)

x1= x(1); 
x2 = x(2);

u = input_pendulum(t, 1, 0);

xdot1 = x2;
xdot2 = -gamma(1)*x2 - w_n(1) * sin(x1) - (1/l_n(1))*cos(x1) * u;

xdot = [xdot1; xdot2]; 
% riflettere se ha senso mettere tempo di 10s 

end

