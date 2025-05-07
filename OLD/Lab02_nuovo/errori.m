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
t_f = 2000; % tempo finale [s]
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

% Task 1.1/1.2 (sistema non lineare)
layout = tiledlayout(3,2);
% Ingrandisci subito la figura
set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
for i=1:3
    odefun = @(t,x) InversePODE(t,x,pend1,i);
    [t_n1, x_n1] = ode113(odefun,deltat,x0);
    [t_n, x_n] = ode78(odefun,deltat,x0);
    nexttile; plot(t_n,abs(x_n(:,1)-x_n1(:,1)), 'LineWidth', 2); % spostamento
    ylabel('Spostamento [m]');
    title(sprintf('Spostamento del Carrello = %g N',pend1.amp(i)));
    grid on;
    nexttile; plot(t_n,abs(x_n(:,3)*(180/pi)-x_n1(:,3)*(180/pi)), 'LineWidth', 2); % angolo in gradi
    xlabel('Tempo (s)');
    ylabel('Theta (°)');
    title(sprintf('Angolo del pendolo = %g N',pend1.amp(i)));
    grid on;
end
sgtitle('Errore tra i diversi metodi all impulso', 'FontSize', 16, 'FontWeight', 'bold', 'Color', 'red');
set(gcf, 'PaperPositionMode', 'auto');

