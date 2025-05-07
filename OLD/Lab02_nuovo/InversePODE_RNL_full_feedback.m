function [xdot] = InversePODE_RNL_full_feedback(t,x,l,m,g,M,K,A2,B2,t_pulse)

x1 = x(1);
x2 = x(2);
x3 = x(3);
x4 = x(4);

k1 = K(1);
k2 = K(2);
k3 = K(3);
k4 = K(4);

d = input_d(t,t_pulse);


xdot1 = x2;
xdot2 = (+ m*g*sin(x3)*cos(x3) - m*l*(x4^2)*sin(x3)  +...
    ((-k1*x1-k2*x2-k3*x3-k4*x4) + d(t)))/(M+m-(m*cos(x3)^2));
xdot3 = x4;
xdot4 = ((-m)/(M+m)*(x4^2)*cos(x3)*sin(x3)  + ((-k1*x1-k2*x2-k3*x3-k4*x4) + d(t))...
    /((M+m)*l)*cos(x3) +g/l*sin(x3))/(1-m/(M+m)*(cos(x3)^2));


xdot= [xdot1; xdot2; xdot3; xdot4];
 % xdot = (A2-B2*K)x - B2(K*x - d(t));

end
