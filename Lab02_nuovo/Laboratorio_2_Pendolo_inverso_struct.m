%% parametri e struct
clear all
close all
clc

M = 3; % massa carrello [Kg]
m = 0.1; % massa pendolo [Kg]
l = 0.75; % lunghezza pendolo [Kg]
g = 9.81; % acc. di gravità [m/s^2]
c = 0.1; % attrito viscoso carrello [N/m/s]
b = 0.001; % attrito viscoso pendolo [Nm/rad/s]

x0 = [0; 0; pi; 0]; % cond. iniziali
x02 = [0; 0; 0; 0]; % cond. iniziali linearizzate (delta_x)
t_i = 0; % tempo iniziale [s]
t_f = 20; % tempo finale [s]
h= 0.001; % discrezione
deltat = t_i:h:t_f; % intervallo di tempo
amp = [1 5 25]; % amplificazione forza pulsante
t_pulse = [1 2]; % tempo di applicazione dello scalino

a_c = 1/(M*l); 
b_c = (m+M)*g/(M*l);
k_min = (M+m)*g-1;
omega = sqrt(k_min*a_c-b_c);    % Frequenza angolare 
deltat2 = 0:0.001:100; % Intervallo di tempo

kp_m = (m+M)*g; % guadagno regolatore P marginalmente stabile

t_pulse2 = 0.1;  % % tempo di applicazione dello scalino

psi = 0.456;
wn = 3.63;  % [rad/s]
Kd = 2*psi*wn / a_c;
Kp = (wn^2+b_c) /a_c;
wd = wn*sqrt(1-psi^2);
phi = acos(psi);
deltat3 = 0:0.001:5;

pend1 = struct('M',M,'m',m,'l',l,'g',g,'c',c,'b',b,'x0',x0,'x02',x02,'deltat',deltat,'amp',amp,'t_pulse',t_pulse,...
     'a_c',a_c,'b_c',b_c,'k_min',k_min,'omega',omega,'deltat2',deltat2,'kp_m',kp_m,'t_pulse2',t_pulse2,...
     'psi',psi,'wn',wn,'Kd',Kd,'Kp',Kp,'wd',wd,'phi',phi,'deltat3',deltat3);
%%
% Task 1.1/1.2 (sistema non lineare)
layout = tiledlayout(3,2);
% Ingrandisci subito la figura
set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
for i=1:3
    odefun = @(t,x) InversePODE(t,x,pend1,i);
    [t_n, x_n] = ode45(odefun,deltat,x0);
    nexttile; plot(t_n,x_n(:,1), 'LineWidth', 2); % spostamento
    ylabel('Spostamento [m]');
    title(sprintf('Spostamento del Carrello = %g N',pend1.amp(i)));
    grid on;
    nexttile; plot(t_n,x_n(:,3)*(180/pi), 'LineWidth', 2); % angolo in gradi
    xlabel('Tempo (s)');
    ylabel('Theta (°)');
    title(sprintf('Angolo del pendolo = %g N',pend1.amp(i)));
    grid on;
end
sgtitle('Risposta non lineare all impulso', 'FontSize', 16, 'FontWeight', 'bold', 'Color', 'red');
set(gcf, 'PaperPositionMode', 'auto');

% Plot della forza applicata con un amplificazione di 1 N
F = (pend1.deltat > 1 & pend1.deltat < 2) * 1; % 1 N tra 1 s e 2 s, 0 altrove
figure
plot(pend1.deltat, F, 'LineWidth', 2);
xlabel('Tempo (s)');
ylabel('Forza (N)');
title('Forza sul Carrello amp = 1 N');
grid on;
% Task 1.1/1.2 (scelta metodo di integrazione)
layout = tiledlayout(3,2);

t_int_45 = zeros(1,3);
t_int_15s = zeros(1,3);
t_int_78 = zeros(1,3); 
t_int_45m = zeros(1,3);
options1 = odeset('AbsTol', 1e-15, 'RelTol', 1e-12, 'MaxStep', 0.001);
options2 = odeset('AbsTol', 1e-12, 'RelTol', 1e-9, 'MaxStep', 0.01);
x_n_113_n = [];
x_n_45_n = []; 
x_n_45m_n = []; 
x_n_78_n = [];

for i=1:3
    
    odefun = @(t,x) InversePODE(t,x,pend1,i);
    tic
    [t_n, x_n] = ode45(odefun,deltat,x0);
    x_n_45_n = [x_n_45_n, x_n]; 
    t_int_45(i) = toc;
    tic
    [t_n_stiff, x_n_es_stiff] = ode15s(odefun,deltat,x0);
    t_int_15s(i) = toc;
    
    % Verifichiamo che il problema sia non stiff anche se si capisce già da
    % un breve osservazione iniziale. 
    % Ora, una volta verificato che il problema sia non stiff...
    % possiamo verificare le performance del metodo attraverso un
    % metodo molto preciso che approssima molto bene la soluzione analitica
    
    [t_n_113, x_n_113] = ode113(odefun,deltat,x0,options1);
     x_n_113_n = [x_n_113_n,x_n_113];

    % metodo 45 modificato 
    tic
    [t_n, x_n_45_m] = ode45(odefun,deltat,x0, options2);
    t_int_45m(i) = toc;
    x_n_45m_n = [x_n_45m_n, x_n_45_m];

    % metodo 78 
    tic
    [t_n, x_n_78] = ode78(odefun,deltat,x0);
    t_int_78(i) = toc;
    x_n_78_n = [x_n_78_n, x_n_78];


    nexttile; plot(t_n,x_n(:,1), 'LineWidth', 2); % spostamento
    ylabel('Spostamento [m]');
    title(sprintf('Spostamento del Carrello = %g N',pend1.amp(i)));
    grid on;
    nexttile; plot(t_n,x_n(:,3)*(180/pi), 'LineWidth', 2); % angolo in gradi
    xlabel('Tempo (s)');
    ylabel('Theta (°)');
    title(sprintf('Angolo del pendolo = %g N',pend1.amp(i)));
    grid on;

end

% Plot della forza applicata con un amplificazione di 1 N

F = (pend1.deltat > 1 & pend1.deltat < 2) * 1; % 1 N tra 1 s e 2 s, 0 altrove
figure
plot(pend1.deltat, F, 'LineWidth', 2);
xlabel('Tempo (s)');
ylabel('Forza (N)');
title('Forza sul Carrello amp = 1 N');
grid on;


% Calcolo degli errori assoluti

metodi = {'ode45', 'ode45m', 'ode78'};
x_all = {x_n_45_n, x_n_45m_n, x_n_78_n};
err = struct();

for k = 1:numel(metodi)
    err.(metodi{k}).spostamento = zeros(length(t_n), 3);
    err.(metodi{k}).theta = zeros(length(t_n), 3);
end

% Loop su ogni forzante
for i = 1:3
    offset = 4*(i-1); % offset colonne: 0, 4, 8
    for k = 1:numel(metodi)
        xk = x_all{k};
        err.(metodi{k}).spostamento(:,i) = abs(x_n_113_n(:,1+offset) - xk(:,1+offset));
        err.(metodi{k}).theta(:,i)       = abs(x_n_113_n(:,3+offset) - xk(:,3+offset));
    end
end

% Calcolo errore massimi

for k = 1:numel(metodi)
    fprintf('\nErrori massimi con %s:\n', metodi{k});
    for i = 1:3
        max_sp = max(err.(metodi{k}).spostamento(:,i));
        max_th = max(err.(metodi{k}).theta(:,i));
        fprintf('Forzante %d N: Spostamento = %.2e, Theta = %.2e\n', pend1.amp(i), max_sp, max_th);
    end
end

% Visualizzazione separata per ogni metodo con tutte le forzanti

for k = 1:numel(metodi)
    figure('Name', sprintf('Errori - %s', metodi{k}), 'NumberTitle', 'off');
    tiledlayout(3,2);

    for i = 1:3
        nexttile;
        plot(t_n, err.(metodi{k}).spostamento(:,i), 'LineWidth', 2);
        ylabel('Spostamento [m]');
        title(sprintf('Errore Spostamento = %g N (%s)', pend1.amp(i), metodi{k}));
        grid on;

        nexttile;
        plot(t_n, err.(metodi{k}).theta(:,i), 'LineWidth', 2);
        xlabel('Tempo (s)');
        ylabel('Theta (°)');
        title(sprintf('Errore Angolo = %g N (%s)', pend1.amp(i), metodi{k}));
        grid on;
    end
end
%%
% Task 1.3 (Simulink non lineare)
figure
layout = tiledlayout(3,2);

for F=pend1.amp

    % con l'utilizzo dell'approccio con blocchi di integrazione
    sim_Lab2_1_3_1 = sim('Model_31_02_1.slx');

    nexttile;
    plot(sim_Lab2_1_3_1.tout,sim_Lab2_1_3_1.x1,'LineWidth', 2); % spostamento
    xlabel('Tempo (s)');
    ylabel('Spostamento [m]');
    title(sprintf('Spostamento del Carrello con F = %g N',F));
    grid on;
    nexttile;
    plot(sim_Lab2_1_3_1.tout,sim_Lab2_1_3_1.theta,'LineWidth', 2); % angolo
    xlabel('Tempo (s)');
    ylabel('Theta (°)');
    title(sprintf('Angolo del pendolo con F = %g N',F));
    grid on;

end

figure
layout = tiledlayout(3,2);

for F=pend1.amp

    % con l'utilizzo di M-function
    sim_Lab2_1_3_2 = sim('Model_31_02_2.slx');
    
    nexttile;
    plot(sim_Lab2_1_3_2.tout,sim_Lab2_1_3_2.x1,'LineWidth', 2); % spostamento
    xlabel('Tempo (s)');
    ylabel('Spostamento [m]');
    title(sprintf('Spostamento del Carrello con F = %g N',F));
    grid on;
    nexttile;
    plot(sim_Lab2_1_3_2.tout,sim_Lab2_1_3_2.theta,'LineWidth', 2); % angolo
    xlabel('Tempo (s)');
    ylabel('Theta (°)');
    title(sprintf('Angolo del pendolo con F = %g N',F));
    grid on;

end


%% Task 1.4/1.5 (sistema linearizzato)
[A,B,C,D] = state_space1(pend1);

figure
layout = tiledlayout(3,2);
% Ingrandisci subito la figura
set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
for i=1:3
    odefun_l = @(t,x) InversePODE_L(t,x,A,B,pend1,i);
    [t_n, x_n_l] = ode45(odefun_l,pend1.deltat,pend1.x02);
    nexttile;
    plot(t_n,x_n_l(:,1),'LineWidth', 2); % spostamento
    xlabel('Tempo (s)');
    ylabel('Spostamento [m]');
    title(sprintf('Spostamento del Carrello linearizzato con F = %g N',pend1.amp(i)));
    grid on;
    nexttile;
    plot(t_n,(x_n_l(:,3)+pi)*180/pi,'LineWidth', 2); % angolo, sommo pi per ritrovare theta(t) e non delta_theta(t)
    xlabel('Tempo (s)');
    ylabel('Theta (°)');
    title(sprintf('Angolo del pendolo linearizzato con F = %g N',pend1.amp(i)));
    grid on;
end
sgtitle('Risposta lineare all impulso', 'FontSize', 16, 'FontWeight', 'bold', 'Color', 'red');
set(gcf, 'PaperPositionMode', 'auto');

% Task 1.4 (differenze tra lin. e non lin.)
[A,B,C,D] = state_space1(pend1);

figure
layout = tiledlayout(3,2);
% Ingrandisci subito la figura
set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
for i=1:3
    odefun_l = @(t,x) InversePODE_L(t,x,A,B,pend1,i);
    [t_n, x_n_l] = ode45(odefun_l,pend1.deltat,pend1.x02); %lineare
    nexttile;
    plot(t_n,x_n_l(:,1),'LineWidth', 2); % spostamento
    xlabel('Tempo (s)');
    ylabel('Spostamento [m]');
    title(sprintf('Spostamento del Carrello linearizzato con F = %g N',pend1.amp(i)));
    grid on;
    hold on;
    odefun = @(t,x) InversePODE(t,x,pend1,i);
    [t_n, x_n] = ode45(odefun,deltat,x0); % non lineare
    plot(t_n,x_n(:,1), 'LineWidth', 2); % spostamento
    ylabel('Spostamento [m]');
    title(sprintf('Spostamento del Carrello = %g N',pend1.amp(i)));
    grid on;
    legend('Lineare','Non lineare');
    nexttile;
    plot(t_n,(x_n_l(:,3)+pi)*180/pi,'LineWidth', 2); % angolo, sommo pi per ritrovare theta(t) e non delta_theta(t)
    xlabel('Tempo (s)');
    ylabel('Theta (°)');
    title(sprintf('Angolo del pendolo linearizzato con F = %g N',pend1.amp(i)));
    grid on;
    hold on;
    plot(t_n,x_n(:,3)*(180/pi), 'LineWidth', 2); % angolo in gradi
    xlabel('Tempo (s)');
    ylabel('Theta (°)');
    title(sprintf('Angolo del pendolo = %g N',pend1.amp(i)));
    grid on;
    legend('Lineare','Non lineare');
end
sgtitle('differenze tra risposta lineare e non lineare', 'FontSize', 16, 'FontWeight', 'bold', 'Color', 'red');
set(gcf, 'PaperPositionMode', 'auto');


%% Task 1.4 (funzione di trasferimento,pzmap)
states = {'x' 'x_dot' 'theta' 'theta_dot'};
inputs = {'forza f'};
outputs = {'x' 'theta'};

sys_mimo = ss(A,B,C,D,'statename',states,'inputname',inputs,'outputname',outputs);

G = tf(sys_mimo);

for i=1:2
figure;

% 🔧 Allarga la figura per spazio laterale
set(gcf, 'Units', 'normalized', 'Position', [0.1 0.1 0.6 0.7]);

% Calcola poli e zeri
[p, z] = pzmap(G(1));

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
title('Poli e Zeri della G(s) da u a \theta');
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
    bode(G(i));
    h = findall(gcf,'Type','line');
    set(h, 'LineWidth', 2);
    grid on 
end

% zoom mappa poli e zeri G_1(s)
figure
pzmap(G(1))
yline(0);
xline(0);
grid on

% zoom mappa poli e zeri G_2(s)
figure
pzmap(G(2))
yline(0);
xline(0);
grid on

% G_1(s) asintotico 
clc
asymp(sys_mimo(1,:),0.01,1000);
legend("Reale","Asintotico");

% % G_2(s) asintotico 
clc
asymp(sys_mimo(2,:),0.01,1000);
legend("Reale","Asintotico");


% Task 1.6 (Simulink linearizzato) 
% comparare su latex
figure
layout = tiledlayout(3,2);

for F=pend1.amp

    % con l'utilizzo di M-function
    sim_Lab2_1_6 = sim('Model_31_02_3.slx');
    
    nexttile;
    plot(sim_Lab2_1_6.tout,sim_Lab2_1_6.x1,'LineWidth', 2); % spostamento
    xlabel('Tempo (s)');
    ylabel('Spostamento [m]');
    title(sprintf('Spostamento del Carrello linearizzato con F = %g N',F));
    grid on;
    nexttile;
    plot(sim_Lab2_1_6.tout,sim_Lab2_1_6.theta,'LineWidth', 2); % angolo
    xlabel('Tempo (s)');
    ylabel('Theta (°)');
    title(sprintf('Angolo del pendolo linearizzato con F = %g N',F));
    grid on;

end

% task 2.1/2.2 (non lineare, spazio stati non lineare, lineare, spazio stati lineare,tf,pzmap)
[A2,B2,C2,D2] = state_space2(pend1);

states = {'x' 'x_dot' 'theta' 'theta_dot'};
inputs2 = {'input u di controllo'};
outputs = {'x' 'theta'};

sys_mimo2 = ss(A2,B2,C2,D2,'statename',states,'inputname',inputs2,'outputname',outputs);

G2 = tf(sys_mimo2);
figure
pzmap(G2(2));
yline(0);
xline(0);
grid on
p = pole(G2(2));
z = zero(G2(2));

Ft_P = G2(2)/(1+pend1.kp_m*G2(2)); % funzione di trasferimento da d a theta con regolatore p

figure
rlocus(sys_mimo2(2,:))
title('Luogo delle radici diretto');
ax = gca;             % Ottieni l'asse corrente
lines = findall(ax, 'Type', 'line');  % Trova tutte le linee nel plot
set(lines, 'LineWidth', 2);           % Modifica lo spessore

figure
rlocus(-sys_mimo2(2,:))
title('Luogo delle radici inverso');
ax = gca;             % Ottieni l'asse corrente
lines = findall(ax, 'Type', 'line');  % Trova tutte le linee nel plot
set(lines, 'LineWidth', 2);           % Modifica lo spessore


% Task 2.3 (controllo stabilità del sistema regolato) 
% Plot del  applicata con un amplificazione di 1 N
F = (pend1.deltat > 0 & pend1.deltat < pend1.t_pulse2) * 1; % 1 N tra 0 s e 0.1 s, 0 altrove
figure
plot(pend1.deltat, F, 'LineWidth', 2);
xlabel('Tempo (s)');
ylabel('Forza (N)');
title('Disturbo sul carrello amp = 1 N');
grid on;
xlim([0 5]);


% Task 2.3 D (dimostrazione valore minimo per stabilità marginale)
% Inizializza theta(t)
theta = zeros(size(pend1.deltat2));

% Calcola theta(t) per i diversi intervalli
for i = 1:length(pend1.deltat2)
    if pend1.deltat2(i) < 0
        theta(i) = 0;
    elseif deltat2(i) < 0.1
        theta(i) = (pend1.a_c / pend1.omega^2) * (1 - cos(pend1.omega * pend1.deltat2(i)));
    else
        theta(i) = (pend1.a_c / pend1.omega^2) * ((1 - cos(pend1.omega * pend1.deltat2(i))) - (1 - cos(pend1.omega * (pend1.deltat2(i) - 0.1))));
    end
end

% Plot
figure;
plot(pend1.deltat2, theta, 'LineWidth', 2);
xlabel('t [s]');
ylabel('rad');
title('Andamento di \theta(t) definita a tratti con K_p = (M+m)g-1');
grid on;


% Task 2.4 B ( determinazione kp e kd a mano)
% fare diagramma blocchi closed-loop system 
theta_step_v = zeros(size(pend1.deltat3));

theta_step = @(t) pend1.a_c/pend1.wn^2 * (1-1/sqrt(1-pend1.psi^2)*exp(-pend1.psi*pend1.wn*t)*sin(pend1.wd*t+pend1.phi));

for i = 1:length(pend1.deltat3)
    theta_step_v(i) = theta_step(pend1.deltat3(i));
end

figure
plot(pend1.deltat3,theta_step_v)
grid on 

%% Task 2.4 C (plot della risposta non lineare del sistema controllato dal PD) 
s = tf('s'); % definisce la s di Laplace
R = Kp + Kd*s;

odefun_NL = @(t,x) InversePODE_RNL(t,x,pend1);

[t_n_NL, x_n_l_NL] = ode45(odefun_NL,pend1.deltat3,pend1.x02);

figure
plot(pend1.deltat3,x_n_l_NL(:,1),'LineWidth',2);
title('Posizione del carrello')
grid on
exportgraphics(gcf, 'pos_car_2.pdf', 'ContentType', 'vector', 'BackgroundColor','none')
figure
plot(pend1.deltat3,x_n_l_NL(:,3)*180/pi,'LineWidth',2)
title('Angolo del pendolo')
grid on
exportgraphics(gcf, 'ang_pen_2.pdf', 'ContentType', 'vector', 'BackgroundColor','none')

F = -pend1.Kp*x_n_l_NL(:,3) - pend1.Kd * x_n_l_NL(:,4);

figure
plot(pend1.deltat3,F,'LineWidth',2)
title('Forza')
grid on
exportgraphics(gcf, 'F2.pdf', 'ContentType', 'vector', 'BackgroundColor','none')
%% Task 2.4 D (Punto C con utilizzo di Simulink)
d = 1;  % [N] Disturbo

figure
layout = tiledlayout(1,2);

% con l'utilizzo di M-function
sim_Lab2_4d = sim('Model_31_02_4.slx');

nexttile;
plot(sim_Lab2_4d.tout,sim_Lab2_4d.x1,'LineWidth', 2); % spostamento
xlabel('Tempo (s)');
ylabel('Spostamento [m]');
title(sprintf('Spostamento con d = %g N',d));
grid on;
nexttile;
plot(sim_Lab2_4d.tout,sim_Lab2_4d.control_force,'LineWidth', 2); % forza
xlabel('Tempo (s)');
ylabel('Forza di controllo (N)');
title(sprintf('Forza di controllo [N]'));
grid on;

%% Task 2.5

Co = ctrb(A2, B2);
rank(Co); 

p_des = [-4,-5, -1.6553 + 3.2306i , -1.6553 - 3.2306i];
K= place(A2, B2, p_des);

% K = place(A,B,p) pone i poli desiderati a loop chiuso p calcolando una 
% matrice di guadagno del feedback di stato K
odefun_NL_full_feedback = @(t,x) InversePODE_RNL_full_feedback(t,x,l,m,g,M,K,A2,B2,t_pulse2);
[t_n_NL_full , x_n_l_NL_full] = ode45(odefun_NL_full_feedback,deltat3,x02);

figure
plot(t_n_NL_full,x_n_l_NL_full(:,1),'LineWidth',2)
grid on
hold on
plot(t_n_NL_full,x_n_l_NL(:,1),'LineWidth',2)
ylabel('Posizione del carrello');
title('Posizione del carrello full feedback');
legend('Posizione del carrello full','Posizione del carrello partial','location','southwest');

figure
plot(t_n_NL_full,x_n_l_NL_full(:,3).*180/pi,'LineWidth',2)
grid on
hold on
plot(t_n_NL_full,x_n_l_NL(:,3).* 180/pi,'LineWidth',2)

ylabel('Angolo del pendolo');
title('Angolo del pendolo full feedback');
legend('Angolo del pendolo full','Angolo del pendolo partial');


k1 = K(1);
k2 = K(2);
k3 = K(3);
k4 = K(4);

F_full_feedback = -k1*x_n_l_NL_full(:,1)-k2*x_n_l_NL_full(:,2)...
    -k3*x_n_l_NL_full(:,3)-k4*x_n_l_NL_full(:,4); 
figure
plot(t_n_NL_full,F_full_feedback,'LineWidth',2);
ylabel('Forza ');
title('Forza full feedback');
grid on