function [A,B,C,D] = state_space1(pend1)

M = pend1.M;
m = pend1.m;
g = pend1.g;
b = pend1.b;
l = pend1.l;
c = pend1.c;

A = [0 1 0 0; 0 -c/M (-m*g)/M -b/(M*l); 0 0 0 1; 0 -c/(M*l) (-g*(M+m))/(M*l) (-b*(M+m))/(M*m*(l^2))];
B = [0; 1/M; 0; 1/(M*l)];
C = [1 0 0 0; 0 0 1 0];
D = 0;
end

