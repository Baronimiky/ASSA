function [A2,B2,C2,D2] = state_space2(pend1)

M = pend1.M;
m = pend1.m;
g = pend1.g;
l = pend1.l;


A2 = [0 1 0 0; 0 0 (m*g)/M 0; 0 0 0 1; 0 0 (g*(M+m))/(M*l) 0];
B2 = [0; 1/M; 0; 1/(M*l)];
C2 = [1 0 0 0; 0 0 1 0];
D2 = 0;

end

