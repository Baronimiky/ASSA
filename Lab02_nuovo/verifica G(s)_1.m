
syms c M b l m g s
A = [0 1 0 0; 0 -c/M (-m*g)/M -b/(M*l); 0 0 0 1; 0 -c/(M*l) (-g*(M+m))/(M*l) (-b*(M+m))/(M*m*(l^2))];
B = [0; 1/M; 0; 1/(M*l)];
C = [1 0 0 0; 0 0 1 0];
D = 0;

G = C*inv(s*eye(4)-A)*B
G1 = simplify(G(1))
G2 = simplify(G(2))
