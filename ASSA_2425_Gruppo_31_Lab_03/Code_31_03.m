%% CODE

% Task 1 (risoluzione equazioni)
clear 
close all
clc

% dati 
h_trim = 0;
rho_h_trim = 1.225;

S = 32.8; % [m^2]
m = 7484.4; % [kg]
I_yy = 84309; % [kgm^2]
c = 2.29; % [m]
alpha_t = 0; % [°]
delta_z_t = 0.378; % [m]
V_kts = 120; % [kts]
V = 61.74; % [m/s] 
g = 9.81; % [m/s^2]


x0 = [0.5, -0.5, 10000];  % vettore di guess iniziale, fsolve ha bisogno di un vettore di valori 
% di stima iniziali per poi effettuare delle iterazioni fino a convergenza.
% In teoria andrebbero messi valori che potrebbero essere reali (ho preso
% spunto dalle soluzioni). Mentre una buona guess per T potrebbe essere il
% drag, T=D. 
% il vettore è x0 = [alpha delta T]

% By default fsolve chooses the trust-region dogleg algorithm. pagina
% mathworks di fsolve


F_trim = @(x) trim_equations(x, m, I_yy, rho_h_trim, V, S,delta_z_t,alpha_t,c,g);
[x_trim] = fsolve(F_trim, x0);

alpha_trim = x_trim(1) * 180/pi; % trasformo in °
theta_trim = alpha_trim; % essendo  gamma nullo (volo livellato) 
delta_trim = x_trim(2) * 180/pi; % trasformo in °


% Da x_trim ottengo alpha_trim delta_trim e T_trim, per ottenere u e w
% sfrutto la conoscenze della velocità V: V = sqrt(u^2 + w^2) e l'ultima
% equazione della dinamica longitudinale: usin(alpha) = w cos(alpha),

u = V / (sqrt(1+tan(x_trim(1))^2 )); 
w = u* tan(x_trim(1));

% Definisco le equazioni della dinamica longitudinale: 

% Task 2 

coeff = struct();  
coeff.CL_0   = 0.895;   % termine costante di portanza
coeff.CL_alpha = 5.01;  % CL derivato rispetto all'angolo d'attacco
coeff.CL_e   = 0.722;   % CL derivato rispetto alla deflessione dell'equilibratore
    
coeff.CD_0   = 0.177;   % drag parassita
coeff.CD_alpha = 0.232; % drag lineare con alpha
coeff.CD_alpha2 = 1.393; % drag quadratico con alpha
    
coeff.CM_0   = -0.046;   % momento a zero alpha
coeff.CM_alpha = -1.087; % CM derivato del momento rispetto a alpha
coeff.CM_q   = -7.055;   % CM derivato rispetto al tasso di variazione di alpha
coeff.CM_e   = -1.88; 

tempo = 0:0.01:5; 
x0_ver = [u; w; 0; theta_trim*pi/180; h_trim];
T_trim = x_trim(3);

Odefun = @(t,x) long_dyn(t,x, m, I_yy, rho_h_trim, V, S,delta_z_t, alpha_t,c,g,coeff,delta_trim*pi/180,T_trim);

[~, x_out] = ode45(Odefun, tempo, x0_ver);

for i = 1:5
    figure
    plot(tempo, x_out(:,i), 'LineWidth', 2);  % Usa direttamente x_out
    grid on
    xlabel('Tempo(s)')

    switch i
        case 1
            ylabel('Velocità in avanti "u" (m/s)');
        case 2
            ylabel('Velocità verticale "w" (m/s)');
        case 3
            ylabel('Velocità angolare di beccheggio "q" (rad/s)');
        case 4
            ylabel('Angolo di beccheggio "theta (°)');
        case 5
            ylabel('Quota "h" (m)');
        end
    % Calcola la media dei dati per centrare l'asse Y attorno ai dati
    y_center = mean(x_out(:,i));        % Media dei dati
    range_fixed = 0.1;                   % Imposta il range fisso dell'asse Y (±0.05)
    step = 0.05;                         % Step tra le tacche (0.05)
    
    % Calcola i limiti dell'asse Y
    y_min = y_center - range_fixed / 2;
    y_max = y_center + range_fixed / 2;
    
    % Imposta le tacche ogni 0.05
    yt = y_min:step:y_max;
    
    % Modifica l'asse Y
    ylim([y_min y_max]);  % Limiti fissi
    yticks(yt);            % Tacche ogni 0.05
    yline(0, '--k');       % Aggiungi linea orizzontale a y=0 (facoltativo)
end

alpha = atan(x_out(:,2)/x_out(:,1));
alpha_deg = alpha.* 180/pi;
figure
plot(tempo, alpha_deg, 'LineWidth', 2);
xlabel('tempo(s)'); 
ylabel('Angolo di attacco "alpha" (°)')
 % Calcola la media dei dati per centrare l'asse Y attorno ai dati
    y_center = mean(alpha_deg(:,1));        % Media dei dati
    range_fixed = 0.1;                   % Imposta il range fisso dell'asse Y (±0.05)
    step = 0.05;                         % Step tra le tacche (0.05)
    
    % Calcola i limiti dell'asse Y
    y_min = y_center - range_fixed / 2;
    y_max = y_center + range_fixed / 2;
    
    % Imposta le tacche ogni 0.05
    yt = y_min:step:y_max;
    
    % Modifica l'asse Y
    ylim([y_min y_max]);  % Limiti fissi
    yticks(yt);            % Tacche ogni 0.05




% Task 3 

tempo2 = 0:0.1:100; 

diff_delta_e = -1; % [°]


Odefun3 = @(t,x) long_dyn(t, x, m, I_yy, rho_h_trim, V, S, delta_z_t, alpha_t, c, g, coeff, ...
                          pos_equilibratore(t, delta_trim*pi/180, diff_delta_e*pi/180), ...
                          T_trim);

[t_out3, x_out3] = ode45(Odefun3, tempo2, x0_ver);


x_out3(:,4) = [x_out3(:,4)*180/pi]; % per ottenere theta in gradi 

for i = 1:5
   figure(i)
   plot(tempo2,x_out3(:,i)); 
   grid on
    xlabel('Tempo(s)')

    switch i
        case 1
            ylabel('Velocità in avanti "u" (m/s)');
        case 2
            ylabel('Velocità verticale "w" (m/s)');
        case 3
            ylabel('Velocità angolare di beccheggio "q" (rad/s)');
        case 4
            ylabel('Angolo di beccheggio "theta (°)');
        case 5
            ylabel('Quota "h" (m)');
        end
end

alpha = atan(x_out3(:,2)./x_out3(:,1));
alpha_deg = alpha.* 180/pi;
figure
plot(tempo2, alpha_deg);
xlabel('tempo(s)'); 
ylabel('Angolo di attacco "alpha" (°)')
grid on


% Task 4

sim_Lab3_4 = sim('Model_31_03_1.slx');

figure();
tiledlayout(3,2);

nexttile;
plot(sim_Lab3_4.tout,sim_Lab3_4.x1);
xlabel('Tempo (s)');
ylabel('Velocità in avanti "u" (m/s)');

nexttile;
plot(sim_Lab3_4.tout,sim_Lab3_4.x2);
xlabel('Tempo (s)');
ylabel('Velocità verticale "w" (m/s)');

nexttile;
plot(sim_Lab3_4.tout,sim_Lab3_4.x3);
xlabel('Tempo (s)');
ylabel('Velocità angolare di beccheggio "q" (rad/s)');

nexttile;
plot(sim_Lab3_4.tout,sim_Lab3_4.x4);
xlabel('Tempo (s)');
ylabel('Angolo di beccheggio "theta (°)');

nexttile;
plot(sim_Lab3_4.tout,sim_Lab3_4.x5);
xlabel('Tempo (s)');
ylabel('Quota "h" (m)');

% task 5
% linearizzazione

% a mano



% task 6
clc

    Xu     = -0.057076461;
    Xw     =  0.125051144;
    Xq     =  0.00;
    Xwdot  =  0.00;

    Zu     = -0.30507998;
    Zw     = -0.86334716;
    Zq     =  0.00;
    Zwdot  =  0.00;

    Mu     = -0.0014736;
    Mw     = -0.036668742;
    Mq     = -0.5442436;
    Mwdot  =  0.00;

    Xde    =  0.00;
    Zde    = -7.38503134;
    Mde    = -3.90965145;

    theta_trim = 0.8426;
    q_trim = 0;
    w_trim = 0.9051;
    u_trim = 61.7334;

    % Forze aerodinamiche in asse body 
    Xt = cos(alpha_t);
    Zt = sin(alpha_t);
    Mt = delta_z_t;

 A = [Xu Xw-q_trim -w_trim -g*cos(theta_trim*pi/180) 0;
      (q_trim+Zu) (Zw) (Zq+u_trim) (-g*sin(theta_trim*pi/180)) 0 ; 
        Mu Mw Mq 0 0;
        0 0 1 0 0;
        sin(theta_trim*pi/180) -cos(theta_trim*pi/180) 0 u_trim*cos(theta_trim*pi/180)+w_trim*sin(theta_trim*pi/180) 0];
 B = [Xde Xt; Zde Zt; Mde Mt; 0 0; 0 0]; 
 C = [1 0 0 0 0; 0 1 0 0 0; 0 0 1 0 0; 0 0 0 1 0; 0 0 0 0 1];
 D = 0;

states = {'u' 'w' 'q' 'theta' 'h'};
inputs = {'delta_e' 'T'};
outputs = {'u' 'w' 'q' 'theta' 'h'};

sys = ss(A,B,C,D,'statename',states,'inputname',inputs,'outputname',outputs);

[wn, psi, autoval] = damp(sys); % frequenza naturale, smorzamento e autovalori

G = tf(sys);

for i=1:5
figure;

% 🔧 Allarga la figura per spazio laterale
set(gcf, 'Units', 'normalized', 'Position', [0.1 0.1 0.6 0.7]);

% Calcola poli e zeri
[p, z] = pzmap(G(i,1));

% Traccia poli e zeri
h1 = plot(real(p), imag(p), 'xr', 'MarkerSize', 12, 'LineWidth', 2); hold on;
h2 = plot(real(z), imag(z), 'ob', 'MarkerSize', 12, 'LineWidth', 2);

% Imposta limiti PRIMA
xlim([-15 15]);
ylim([-10 10]);
drawnow;  % forza il rendering per aggiornare i limiti

% Salva limiti reali
xl = xlim;
yl = ylim;

% Usa plot per disegnare gli assi (più robusto)
plot(xl, [0 0], 'k-', 'LineWidth', 1);  % asse x
plot([0 0], yl, 'k-', 'LineWidth', 1);  % asse y

% Griglia & legenda
axis equal;
grid on;
xlabel('Re');
ylabel('Im');
ingresso = G.InputName{1,1};
uscita = G.OutputName{i,1};
title(sprintf('Poli e Zeri della G(s) da "%s" a "%s"', ingresso, uscita));
legend([h1, h2], {'Poli', 'Zeri'}, 'Location', 'northeast');

% Crea testo con poli/zeri
poli_str = arrayfun(@(x) sprintf('%.2f %+.2fi', real(x), imag(x)), p, 'UniformOutput', false);
zeri_str = arrayfun(@(x) sprintf('%.2f %+.2fi', real(x), imag(x)), z, 'UniformOutput', false);
testo = ['Poli:'; poli_str; {''}; 'Zeri:'; zeri_str];

% Annotation ben spostata a destra, ma più stretta
annotation('textbox', [0.7 0.15 0.25 0.7], ...  % Ridotto a 0.7
    'String', testo, ...
    'FitBoxToText', 'on', ...
    'BackgroundColor', 'white', ...
    'EdgeColor', 'black', ...
    'FontSize', 10);
hold off;
    figure
    bode(G(i,1));
    h = findall(gcf,'Type','line');
    set(h, 'LineWidth', 2);
    grid on 
end




% task 7

x0_ver_lin = [61.7334;0.9051;0;0.84*pi/180; 0]; % [u,w,q,theta,h]


tempo2 = 0:0.1:300;


Odefun4 = @(t,x) long_dyn_lin(t,x,delta_z_t, alpha_t,g, ...
                          diff_delta_e*pi/180);
x0 = [0;0;0;0;0];
[t_out5, x_out5] = ode45(Odefun4, tempo2, x0);


for i = 1:4
   figure
   plot(t_out5,x0_ver_lin(i)+x_out5(:,i)); 
end

%
% 
% grafico il confronto molto rapidamente e molto male
figure
i = 5;

plot(t_out5,x_out5(:,i)); 
hold on
plot(tempo2,x_out3(:,i)-x0_ver_lin(i))





% task 8
% simulink linearizzato

sim_Lab3_8 = sim('Model_31_03_2.slx');


figure
plot(sim_Lab3_8.tout,sim_Lab3_8.x1);

figure
plot(sim_Lab3_8.tout,sim_Lab3_8.x2);

figure
plot(sim_Lab3_8.tout,sim_Lab3_8.x3);

figure
plot(sim_Lab3_8.tout,sim_Lab3_8.x4);

figure
plot(sim_Lab3_8.tout,sim_Lab3_8.x5);





%% FUNCTIONS

% Task 1

function F = trim_equations(x, m, I_yy, rho, V, S,delta_z_t, alpha_t,c,g)
    % x = [alpha; delta_e; T]
    alpha = x(1);
    delta_e = x(2);
    T = x(3);

    % dalla quarta equazione:
    q = 0; % theta_dot = q

    % dalla definizione degli angoli: theta = gamma + alpha, ma gamma = 0
    % (volo livellato, traiettoria // all'orizzonte)

    theta = alpha;

    % Coefficienti aerodinamici
    CL_0   = 0.895;   % termine costante di portanza
    CL_alpha = 5.01;  % CL derivato rispetto all'angolo d'attacco
    CL_e   = 0.722;   % CL derivato rispetto alla deflessione dell'equilibratore
    
    CD_0   = 0.177;   % drag parassita
    CD_alpha = 0.232; % drag lineare con alpha
    CD_alpha2 = 1.393; % drag quadratico con alpha
    
    CM_0   = -0.046;   % momento a zero alpha
    CM_alpha = -1.087; % CM derivato del momento rispetto a alpha
    CM_q   = -7.055;   % CM derivato rispetto al tasso di variazione di alpha
    CM_e   = -1.88;    % CM derivato del momento rispetto a delta_e
    
    % Coefficienti
    CL = CL_0 + CL_alpha * alpha + CL_e * delta_e;
    CD = CD_0 + CD_alpha * alpha + CD_alpha2 * alpha^2;
    CM = CM_0 + CM_alpha * alpha + CM_e * delta_e + CM_q * c/V * q;

    % Forze
    D = 0.5 * rho * V^2 * S * CD; % Resistenza
    M = 0.5 * rho * V^2 * S * c * CM; % Momento rispetto al CG
    L = 0.5 * rho * V^2 * S * CL; % Portanza 
    S = 0; % Devianza

    % D S L sono validi eslcusivamente in asse wind, è quindi necessario
    % trasportarli, attraverso una matrice di rotazione, in asse body,
    % poiche le equazioni della dinamica longitudinale sono in tale sistema
    % di riferimento

    f_A = [-D; S; -L]; % Forze aerodinamiche

    % Matrice di passaggio asse body-wind
    C_bw = [cos(alpha) 0 -sin(alpha); 0 1 0; sin(alpha) 0 cos(alpha)];

    % Forze aerodinamiche in asse body 
    X_A = [1 0 0] * (C_bw*f_A);
    X_T = T*cos(alpha_t);
    Z_A = [0 0 1] * (C_bw*f_A);
    Z_T = T*sin(alpha_t);
    M_A = M;
    M_T = delta_z_t * T;

    % Equazioni dinamica longitudinale
    
    eq1 = -g*sin(theta) + 1/m * X_A + 1/m*X_T;
    eq2 = g*cos(theta) + 1/m * Z_A + 1/m *Z_T;
    eq3 = 1/I_yy * (M_A+M_T);

   F = [eq1, eq2, eq3];

end


% Task 2
function [xdot] = long_dyn(t,x, m, I_yy, rho, V, S,delta_z_t, alpha_t,c,g,coeff,delta_e,T)

    u = x(1); 
    w = x(2);
    q = x(3);
    theta = x(4);
    h = x(5);

    alpha = atan(x(2)/x(1)); 
    V = sqrt((x(1))^2 + (x(2))^2);

    % coefficienti 
    CL = coeff.CL_0 + coeff.CL_alpha * alpha + coeff.CL_e * delta_e;
    CD = coeff.CD_0 + coeff.CD_alpha * alpha + coeff.CD_alpha2 * alpha^2;
    CM = coeff.CM_0 + coeff.CM_alpha * alpha + coeff.CM_e * delta_e + coeff.CM_q * c/V * q;


    % Forze
    D = 0.5 * rho * V^2 * S * CD; % Resistenza
    M = 0.5 * rho * V^2 * S * c * CM; % Momento rispetto al CG
    L = 0.5 * rho * V^2 * S * CL; % Portanza 
    S = 0; % Devianza

    % D S L sono validi eslcusivamente in asse wind, è quindi necessario
    % trasportarli, attraverso una matrice di rotazione, in asse body,
    % poiche le equazioni della dinamica longitudinale sono in tale sistema
    % di riferimento

    f_A = [-D; S; -L]; % Forze aerodinamiche

    % Matrice di passaggio asse body-wind
    C_bw = [cos(alpha) 0 -sin(alpha); 0 1 0; sin(alpha) 0 cos(alpha)];

    % Forze aerodinamiche in asse body 
    X_A = [1 0 0] * (C_bw*f_A);
    X_T = T*cos(alpha_t);
    Z_A = [0 0 1] * (C_bw*f_A);
    Z_T = T*sin(alpha_t);
    M_A = M;
    M_T = delta_z_t * T;

    
    xdot1 = -q*w -g*sin(theta) + 1/m * X_A + 1/m*X_T; 
    xdot2 = q*u + g*cos(theta) + 1/m * Z_A + 1/m *Z_T;
    xdot3 = 1/I_yy * (M_A + M_T);
    xdot4 = q; 
    xdot5 = u * sin(theta) - w* cos(theta); 

    xdot= [xdot1; xdot2; xdot3; xdot4;xdot5];


end


function [delta_e3] = pos_equilibratore(t,delta_trim, diff_delta_e)

    if t>=5 && t<=10 
        delta_e3 =  delta_trim + diff_delta_e;
    else 
        delta_e3 =  delta_trim;
    end
end


% Task 5

function [xdot] = long_dyn_lin(t,x,delta_z_t, alpha_t,g,diff_delta_e)
    
    Xu     = -0.057076461;
    Xw     =  0.125051144;
    Xq     =  0.00;
    Xwdot  =  0.00;

    Zu     = -0.30507998;
    Zw     = -0.86334716;
    Zq     =  0.00;
    Zwdot  =  0.00;

    Mu     = -0.0014736;
    Mw     = -0.036668742;
    Mq     = -0.5442436;
    Mwdot  =  0.00;

    Xde    =  0.00;
    Zde    = -7.38503134;
    Mde    = -3.90965145;

    theta_trim = 0.8426;
    q_trim = 0;
    w_trim = 0.9051;
    u_trim = 61.7334;

    % Forze aerodinamiche in asse body 
    Xt = cos(alpha_t);
    Zt = sin(alpha_t);
    Mt = delta_z_t;

    % calcolo il \Delta delta_e 
   
     if t>=5 && t<=10 
        delta_e =  ( diff_delta_e);
    else 
        delta_e = 0;
     end


  A = [Xu Xw-q_trim -w_trim -g*cos(theta_trim*pi/180) 0;
      (q_trim+Zu) (Zw) (Zq+u_trim) (-g*sin(theta_trim*pi/180)) 0 ; 
        Mu Mw Mq 0 0;
        0 0 1 0 0;
        sin(theta_trim*pi/180) -cos(theta_trim*pi/180) 0 u_trim*cos(theta_trim*pi/180)+w_trim*sin(theta_trim*pi/180) 0];

   B = [Xde Xt; Zde Zt; Mde Mt; 0 0; 0 0]; 

   xdot = A*x + B*[delta_e; 0];

end

