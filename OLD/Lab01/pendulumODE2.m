function [xdot]= pendulumODE2(t,x, gamma, w_n,l_n)

x1= x(1); 
x2 = x(2);

u = input_pendulum_impulse(t,1,-1)
disp(['t = ', num2str(t), ', u = ', num2str(u)]); % Stampa del valore di impulso
xdot1 = x2;
xdot2 = -gamma(1)*x2 - w_n(1) * sin(x1) - (1/l_n(1))*cos(x1) * u;

xdot = [xdot1; xdot2]; 
% riflettere se ha senso mettere tempo di 10s 

end