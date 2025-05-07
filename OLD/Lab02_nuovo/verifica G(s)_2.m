clear all
close all
clc
syms c M b m l g s
A = [0 1 0 0; 0 0 (m*g)/M 0; 0 0 0 1; 0 0 (g*(M+m))/(M*l) 0];
B = [0; 1/M; 0; 1/(M*l)];
C = [1 0 0 0; 0 0 1 0];
D = 0;

G = C*inv(s*eye(4)-A)*B;
G1 = simplify(G(1))
