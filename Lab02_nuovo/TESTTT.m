clear all
close all
clc
% parametri

M = 3; % massa carrello [Kg]
m = 0.1; % massa pendolo [Kg]
l = 0.75; % lunghezza pendolo [Kg]
g = 9.81; % acc. di gravità [m/s^2]
c = 0.1; % attrito viscoso carrello [N/m/s]
b = 0.001; % attrito viscoso pendolo [Nm/rad/s]

x0 = [0; 0; pi; 0]; % cond. iniziali
t_i = 0; % tempo iniziale [s]
t_f = 20; % tempo finale [s]
h= 0.001; % discrezione
deltat = t_i:h:t_f; % intervallo di tempo
amp = [1 5 25]; % amplificazione forza pulsante
t_pulse = [1 2]; % tempo di applicazione dello scalino

%% sistema non lineare

deltat2 = tiledlayout(3,2);
for i=1:3
    odefun = @(t,x) InversePODE(t,x,b,l,m,g,M,c,amp(i),t_pulse);
    [t_n, x_n] = ode45(odefun,t_i:0.1:t_f,x0);
    nexttile; plot(t_n,x_n(:,1), 'LineWidth', 2); % spostamento
    ylabel('x [m]');
    title(sprintf('spostamento del Carrello = %g N',amp(i)));
    grid on;
    nexttile; plot(t_n,x_n(:,3)*(180/pi), 'LineWidth', 2); % angolo in gradi
    xlabel('Tempo (s)');
    ylabel('theta (°)');
    title(sprintf('Angolo del pendolo = %g N',amp(i)));
    grid on;
end

% Plot della forza applicata con un amplificazione di 1 N
F = (deltat > 1 & deltat < 2) * 1; % 1 N tra 1 s e 2 s, 0 altrove
figure
plot(deltat, F, 'LineWidth', 2);
xlabel('Tempo (s)');
ylabel('Forza (N)');
title('Forza sul Carrello amp = 1 N');
grid on;

%% sistema linearizzato
clc
A = [0 1 0 0; 0 -c/M (-m*g)/M -b/(M*l); 0 0 0 1; 0 -c/(M*l) (-g*(M+m))/(M*l) (-b*(M+m))/(M*m*(l^2))];
B = [0; 1/M; 0; 1/(M*l)];
C = [1 0 0 0; 0 0 1 0];
D = [0];

x0 = [0;0;0;0]


figure
deltat2 = tiledlayout(3,2);

for i=1:3
    
    odefun_l = @(t,x) InversePODE_L(t,x,A,B,amp(i),t_pulse);
    [t_n, x_n_l] = ode45(odefun_l,deltat,x0);

    nexttile;
    plot(t_n,x_n_l(:,1),'LineWidth', 2); % spostamento
    xlabel('Tempo (s)');
    ylabel('x [m]');
    title(sprintf('spostamento del Carrello linearizzato = %g N',amp(i)));
    grid on;
    nexttile;
    plot(t_n,x_n_l(:,3)*180/pi,'LineWidth', 2); % angolo
    xlabel('Tempo (s)');
    ylabel('theta (°)');
    title(sprintf('Angolo del pendolo linearizzato = %g N',amp(i)));
    grid on;
end

%% funzione di trasferimento

states = {'x' 'x_dot' 'theta' 'theta_dot'};
inputs = {'forza f'};
outputs = {'x' 'theta'};

sys_mimo = ss(A,B,C,D,'statename',states,'inputname',inputs,'outputname',outputs);

G = tf(sys_mimo);

for i=1:2
    figure
    pzmap(G(i));
    yline(0);
    xline(0);
    grid on
    figure
    bode(G(i));
    grid on
end



%% task 2

A2 = [0 1 0 0; 0 0 (m*g)/M 0; 0 0 0 1; 0 0 (g*(M+m))/(M*l) 0];
B2 = [0; 1/M; 0; 1/(M*l)];
C2 = C;
D2 = D;
x02 = [0;0;0;0];

inputs2 = {'input u di controllo'};

sys_mimo2 = ss(A2,B2,C2,D2,'statename',states,'inputname',inputs2,'outputname',outputs);

G2 = tf(sys_mimo2);
kp = (m+M)*g; % guadagno regolatore P
figure
pzmap(G2(2));
yline(0);
xline(0);
grid on

p = pole(G2(2));
z = zero(G2(2));

Ft_P = G2(2)/(1+kp*G2(2)); % funzione di trasferimento da d a theta con regolatore p
[A3,B3,C3,D3] = tf2ss(cell2mat(Ft_P.Numerator),cell2mat(Ft_P.Denominator)); % spazio agli stati regolato

%% controllo stabilità del sistema regolato ( kp deve tendere all'infinito per renderlo stabile)

t_pulse2 = 0.1; % % tempo di applicazione dello scalino

odefun_l = @(t,x) InversePODE_R(t,x,A3,B3,t_pulse2);
[t_n, x_n_l] = ode45(odefun_l,deltat,x02);
figure
plot(t_n,x_n_l(:,3)*180/pi,'LineWidth', 2); % angolo
xlabel('Tempo (s)');
ylabel('theta (°)');
title(sprintf('Angolo del pendolo con regolatore P = %g N',amp(1)));
grid on;

% Plot della forza applicata con un amplificazione di 1 N
F = (deltat > 0 & deltat < t_pulse2) * 1; % 1 N tra 1 s e 2 s, 0 altrove
figure
plot(deltat, F, 'LineWidth', 2);
xlabel('Tempo (s)');
ylabel('Forza (N)');
title('disturbo sul carrello amp = 1 N');
grid on;
xlim([0 5]);

%% Task 2.3 D

% Parametri
a = 1/(M*l); 
b = (m+M)*g/(M*l);
k_min = (M+m)*g+1;
omega = sqrt(k_min*a-b);    % Frequenza angolare 
deltat2 = -0.05:0.001:100; % Intervallo di tempo

% Inizializza theta(t)
theta = zeros(size(deltat2));

% Calcola theta(t) per i diversi intervalli
for i = 1:length(deltat2)
    if deltat2(i) < 0
        theta(i) = 0;
    elseif deltat2(i) < 0.1
        theta(i) = (a / omega^2) * (1 - cos(omega * deltat2(i)));
    else
        theta(i) = (a / omega^2) * ((1 - cos(omega * deltat2(i))) - (1 - cos(omega * (deltat2(i) - 0.1))));
    end
end

% Plot
figure;
plot(deltat2, theta, 'LineWidth', 2);
xlabel('t [s]');
ylabel('rad');
title('Andamento di \theta(t) definita a tratti');
grid on;

%% Task 2.4 B

psi = 0.456;
wn = 3.63;  % [rad/s]

Kd = 2*psi*wn / a;
Kp = (wn^2+b) /a;

wd = wn*sqrt(1-psi^2);
phi = acos(psi);
deltat3 = 0:0.001:5;

theta_step_v = zeros(size(deltat3));

theta_step = @(t) a/wn^2 * (1-1/sqrt(1-psi^2)*exp(-psi*wn*t)*sin(wd*t+phi));

for i = 1:length(deltat3)
    theta_step_v(i) = theta_step(deltat3(i));
end

figure
plot(deltat3,theta_step_v)
grid on 


%% Task 2.4 C

% plot forza di controllo

u = @(t) dirac(t) - exp(-psi*wn*t)*((wn^2*(1-2*psi^2) - b)* sin(wd*t)/wd + 2*psi*wn*cos(wd*t));
u_v = zeros(size(deltat3));

for i = 1:length(deltat3)
    u_v(i) = u(deltat3(i));
end

figure
plot(deltat3,u_v);
grid on 
title('Forza di controllo lineare')


x = @(t) a*l*dirac(t) - 2*a*l*psi*wn*exp(-psi*wn*t)*cos(wd*t) - ...
    exp(-psi*wn*t)*sin(wd*t)*(a*l*wd + (a*g)/wd);

x_v = zeros(size(deltat3));

for i = 1:length(deltat3)
    x_v(i) = x(deltat3(i));
end

figure
plot(deltat3,x_v);
grid on
title('cart position lineare') % plot cart position 

%%  NON LINEARE

s = tf('s') % definisce la s di Laplace
R = Kp + Kd*s;
t_pulse2 = 0.1;

odefun_NL = @(t,x) InversePODE_RNL(t,x,b,l,m,g,M,c,Kp,Kd,t_pulse2);


x02 = [0;0;0;0]
[t_n_NL, x_n_l_NL] = ode45(odefun_NL,deltat3,x02);


figure
plot(deltat3,x_n_l_NL(:,1))
title('cart position')
grid on 
xlabel('t[s]')
ylabel('x [m]')
F = -Kp*x_n_l_NL(:,3) - Kd * x_n_l_NL(:,4);

figure
plot(deltat3,F)
title('Forza')
grid on 
xlabel('t[s]')
ylabel('F [N]')

figure
plot(deltat3,x_n_l_NL(:,3) .* 180/pi)
title('angolo')
grid on 
xlabel('t[s]')
ylabel('\theta [°]')
%% 
clc

H_s = (s^2-b) / (s^2+2*psi*wn*s + wn^2);
H_t = @(t) ilaplace(H_s,s,t);

H_t_delay =@(t) ilaplace(H_s*exp(-0.1*s),s,t);

u_t = @(t) heaviside(t-0.1);
G  =@(t) H_t(t)-H_t_delay(t) * u_t(t)

G(1)



