function [xdot] = InversePODE(t,x,pend1,i)

M = pend1.M;
m = pend1.m;
g = pend1.g;
b = pend1.b;
l = pend1.l;
c = pend1.c;
amp(i) = pend1.amp(i);
t_pulse = pend1.t_pulse;

x1 = x(1);
x2 = x(2);
x3 = x(3);
x4 = x(4);

u = input_f(t,amp(i),t_pulse);

xdot1 = x2;
xdot2 = ((b/l)*x4*cos(x3) - m*g*sin(x3)*cos(x3) + m*l*(x4^2)*sin(x3) -c*x2 +u(t))/(M+m-(m*cos(x3)^2));
xdot3 = x4;
xdot4 = ((-m)/(M+m)*(x4^2)*cos(x3)*sin(x3) +c/((M+m)*l)*x2*cos(x3) - u(t)/((M+m)*l)*cos(x3) -b/(m*l^2)*x4 +g/l*sin(x3))/(1-m/(M+m)*(cos(x3)^2));

xdot= [xdot1; xdot2; xdot3; xdot4];
end

