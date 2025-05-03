clear all
close all
clc

rho = 1000; % density
d = 1; 
h = 1;
g = 9.81; % acceleration
damping_ratio = 0.003;

n = 10;
eps = [1.841, 5.329, 8.531];
eps_old = eps(end);
for i = 4:n
    eps_new = eps_old + pi;
    eps =[eps, eps_new]; % cercare cos'è epsilon!!!
    eps_old = eps_new;
end

m_liq = (pi/4)*rho*d^2 *h; % liquid mass

m_n = []; % n-pendulum mass
l_n = []; % n-lenght of the pendulum
h_n = []; % n-distance from the CoG
w_n = []; % natural frequency of the n-pendulum in rad/s
fn_n = []; % natural frequency of the n-pendulum in Hz
m0 = 0; % static mass
h0 = 0; % distance from the CoG of the static mass

for i = 1:n
    m_n = [m_n, m_liq*((d*tanh(2*eps(i)*(h/d))/(eps(i)*(eps(i)^2 -1)*h)))];
    l_n = [l_n, d/(2*eps(i)*tanh(2*eps(i)*h/d))];
    h_n = [h_n, (h/2 - (d/(2*eps(i)))*(tanh(eps(i)*h/d) - 1/(sinh(2*eps(i)*h/d))))];
    w_n = [w_n, sqrt(g/l_n(i))];
    fn_n =[fn_n, w_n(i)/(2*pi)];
    m0 = m0 + m_n(i);
end

m0 = m_liq - m0;

for i = 1:n
      h0 = h0 + ((m_n(i)*(h_n(i)-l_n(i)))/m0);
end
  
    
% Task 2a

A_non_smorzata = zeros(2*n,2*n);
B_non_smorzata = zeros(2*n,1);
C_non_smorzata = zeros(n,2*n);
D_non_smorzata = zeros(n,1);

A_smorzata = zeros(2*n,2*n);
C_smorzata = zeros(n,2*n);

gamma = zeros(1,n);

for i = 2:2:2*n
    gamma(1,i/2) = 2*damping_ratio*m_n(i/2)*l_n(i/2)^2*w_n(i/2);

    A_non_smorzata(i-1,i)=1;
    A_non_smorzata(i,i-1)=-g/l_n(i/2);
    B_non_smorzata(i) = -1/l_n(i/2);
    C_non_smorzata (i/2,i-1) = m_n(i/2)*g;
    D_non_smorzata(i/2,1) = - (m0) ;

    A_smorzata(i-1,i)=1;
    A_smorzata(i,i-1)=-g/l_n(i/2);
    A_smorzata(i,i) = -gamma(1,i/2)/(m_n(i/2)*l_n(i/2)^2);

    C_smorzata (i/2,i-1) = m_n(i/2)*g; % segno invertito rispetto equazioni analitiche 
    C_smorzata(i/2,i) = gamma(1,i/2)/l_n(i/2); 
end 


B_smorzata = B_non_smorzata;
D_smorzata = D_non_smorzata;

% Task 3

sys = ss(A_smorzata(1:2,1:2),B_smorzata(1:2,1),C_smorzata(1,1:2),D_smorzata(1));
tempo  = 0:0.01:10;

% soluzione con IMPULSE di matlab

impulse(sys, tempo);
[y_impulse, t_out] = impulse(sys, tempo);


-gamma(1)/l_n(1)^2 - gamma(1)/(m_n(1)*l_n(1)^2)*m0
(m_n(1)+ m0)*g/l_n(1)
G1 = tf(sys)

% SOLUZIONE ESATTA IMPULSE


y = @(t) C_smorzata(1,1:2)*expm(A_smorzata(1:2,1:2)*t)*B_smorzata(1:2,1) % + D_smorzata(1)*Dirac;

y_vals_impulse = arrayfun(y, tempo); % Calcola y per ogni valore di tempo
%y_vals_impulse(1) = y_vals_impulse(1) + D_smorzata(1); % per simulare il delta di Dirac ma non so quanto sia corretto

figure(2)
plot(tempo, y_vals_impulse);

% METODO NUMERICO IMPULSE 

u_impulse = input_pendulum_impulse(tempo, 1, 0.01);
odefun_impulse = @(t,x) A_smorzata(1:2,1:2)*x; % B_smorzata(1:2,1)*u_impulse(t);
% [t_n_impulse, theta_n_impulse] = ode45(odefun_impulse,tempo,
% B_smorzata(1:2,1)); 
% ode45 standard
tic
[t_n_impulse_45, theta_n_impulse_45] = ode45(odefun_impulse,tempo,B_smorzata(1:2,1));
tempo_esecuzione_45 = toc

for ii = 1:length(tempo)
    F_impulse_45(ii) = C_smorzata(1,1)*theta_n_impulse_45(ii,1)+C_smorzata(1,2)...
        *theta_n_impulse_45(ii,2); %+ D_smorzata(1)* u_impulse(ii/100);
end

% ode45 modificato
options = odeset(MaxStep = 0.01, RelTol = 1e-6, AbsTol = 1e-9);
tic;
[t_n_impulse_45m, theta_n_impulse_45m] = ode45(odefun_impulse,tempo,B_smorzata(1:2,1),options);
tempo_esecuzione_45m = toc

for ii = 1:length(tempo)
    F_impulse_45m(ii) = C_smorzata(1,1)*theta_n_impulse_45m(ii,1)+C_smorzata(1,2)...
        *theta_n_impulse_45m(ii,2); %+ D_smorzata(1)* u_impulse(ii/100);
end

% ode78 standard

tic
[t_n_impulse_78, theta_n_impulse_78] = ode78(odefun_impulse,tempo,B_smorzata(1:2,1));
tempo_esecuzione_78 = toc

for ii = 1:length(tempo)
    F_impulse_78(ii) = C_smorzata(1,1)*theta_n_impulse_78(ii,1)+C_smorzata(1,2)...
        *theta_n_impulse_78(ii,2); %+ D_smorzata(1)* u_impulse(ii/100);
end


figure
plot(tempo,F_impulse_45,'LineWidth',1.5) 
hold on
plot(tempo, F_impulse_45m,'LineWidth',1.5)
hold on
plot(tempo, F_impulse_78,'LineWidth',1.5)
legend('ode45', 'ode45 modificato', 'ode78')
xlabel('tempo (s)')
ylabel('F_x (N)')
title('Confronto output diversi metodi')


error_abs_45 = abs(y_vals_impulse - F_impulse_45);
error_abs_45m = abs(y_vals_impulse - F_impulse_45m);
error_abs_78 = abs(y_vals_impulse - F_impulse_78);

max_err_45 = max(error_abs_45)
max_err_45m = max(error_abs_45m)
max_err_78 = max(error_abs_78)


figure
plot(tempo, error_abs_45)
xlabel('tempo(s)')
ylabel('errore')
title('Errore assoluto impulso ode45')
figure
plot(tempo, error_abs_45m)
xlabel('tempo(s)')
ylabel('errore')
title('Errore assoluto impulso ode45m')
figure
plot(tempo, error_abs_78)
xlabel('tempo(s)')
ylabel('errore')
title('Errore assoluto impulso ode78')


% Forma modale di impulse 

[V,Lambda] = eig(A_smorzata(1:2,1:2)); % V matrice di autovettori di A 
   
lambda = diag(Lambda) % lambda matrice diag di autovalori
lambda_inv = diag(1./(lambda));
V_inv = V\B_smorzata(1:2,1);
z = @(t) lambda_inv*((expm(lambda*t)-eye(2)) * V_inv);


u_impulse = input_pendulum_impulse(tempo, 1, 0.01);

for i = 1:length(tempo)
    eLambda = diag(exp(lambda*tempo(i)));
    y_modal_impulse(i) = C_smorzata(1,1:2)*V*(eLambda)* V_inv; %+ D_smorzata(1)*u_impulse(tempo(i));
end

figure
plot(tempo,y_modal_impulse(1,:))

% CICLO ITERAZIONE DELL'ERRORE     
% (differenza tra due configurazioni successive normalizzate)

clc
tempo = 0:0.01:10;
u_impulse = input_pendulum_impulse(tempo, 1, 0.01);
j_ottimo = 0;
toll = 0.01; % 1% di errore

freq = logspace(-1, 2, 1000); % Frequenza da 0.1 a 100 rad/s
w = 2*pi*freq; % Conversione in rad/s se freq era in Hz


y_vec = zeros(n, length(tempo));

    for j = 1:n
        A_non_smorzata = zeros(2*j,2*j);
        B_non_smorzata = zeros(2*j,1);
        C_non_smorzata = zeros(j,2*j);
        D_non_smorzata = zeros(j,1);
        
        A_smorzata = zeros(2*j,2*j);
        C_smorzata = zeros(j,2*j);
        
        gamma = zeros(1,j);
 
        for i = 2:2:2*j
            gamma(1,i/2) = 2*damping_ratio*m_n(i/2)*l_n(i/2)^2*w_n(i/2);
        
            A_non_smorzata(i-1,i)=1;
            A_non_smorzata(i,i-1)=-g/l_n(i/2);
            B_non_smorzata(i,1) = -1/l_n(i/2);
            C_non_smorzata (i/2,i-1) = -m_n(i/2)*g;
            D_non_smorzata(i/2,1) = - (m0) ;
        
            A_smorzata(i-1,i)=1;
            A_smorzata(i,i-1)=-g/l_n(i/2);
            A_smorzata(i,i) = -gamma(1,i/2)/(m_n(i/2)*l_n(i/2)^2);
        
            C_smorzata (i/2,i-1) = m_n(i/2)*g; 
            C_smorzata(i/2,i) = gamma(1,i/2)/l_n(i/2); 
            B_smorzata = B_non_smorzata;
            D_smorzata = D_non_smorzata;
        end 

     
            y = @(t) C_smorzata*expm(A_smorzata*t)*B_smorzata;
            
          

        for i = 1:length(tempo)
            y_vec(j,i) = sum (y(tempo(i)));
        end

      

        if j > 1
        % Calcola l'errore relativo usando la norma L2
                err = norm(y_vec(j,:) - y_vec(j-1,:)) / norm(y_vec(j-1,:));
            if err < toll
                j_ottimo = j;
                break;
            end
        end

            if j_ottimo > 0
                break;
            end

    end

        if (j_ottimo == 0)
            disp('tolleranza non rispettata');
        end

% soluzione con STEP di matlab 

clc
j_ottimo = 4; 
tempo = 0:0.01:100;
A_smorzata_step = A_smorzata(1:2*j_ottimo,1:2*j_ottimo);
B_smorzata_step = B_smorzata(1:2*j_ottimo,1);
C_smorzata_step = C_smorzata(1:j_ottimo,1:2*j_ottimo);
D_smorzata_step = D_smorzata(1:j_ottimo,1);
C_smorzata_step = ones(1,4) * C_smorzata_step; % sommo tutti i contributi di F
% relativi ai singoli pendoli 
D_smorzata_step = D_smorzata_step(1); % D (m0) la sommo solo una volta come
% da definizione di F nel libro sullo sloshing
sys_step = ss(A_smorzata_step,B_smorzata_step,C_smorzata_step,D_smorzata_step);
y_step = step(sys_step,tempo);
step(sys_step,tempo);
exportgraphics(gcf, 'Risposta_step.pdf', 'ContentType', 'vector', 'BackgroundColor','none')
iiiii = step(sys_step,tempo);

% SOLUZIONE HEAVISIDE STEP

y_vals_heaviside_step = zeros(1,length(tempo));
for i = 1:j_ottimo
    c = gamma(i) / (m_n(i)*l_n(i)^2); 
    temp = sqrt(w_n(i)^2-0.25*c^2);
    y_heaviside_step = @(t) m_n(i)*(exp(-c/2*t)*cos(temp*t)+...
        (-0.5*c)/(temp)*exp(-c/2*t)* sin(temp*t)) 
    y_vals_heaviside_step = y_vals_heaviside_step + arrayfun(y_heaviside_step, tempo);
end
y_vals_heaviside_step = y_vals_heaviside_step + (-sum(m_n(1:j_ottimo)) - m0); % aggiungo il contributo
% m0 una sola volta 
figure
plot(tempo,y_vals_heaviside_step)
% soluzione esatta con expm STEP
y_vals_step_esatta = zeros(1,length(tempo));
y_esatta_step = @(t) C_smorzata(1:j_ottimo,1:2*j_ottimo)*inv(A_smorzata(1:2*j_ottimo,1:2*j_ottimo))...
    *(expm(A_smorzata(1:2*j_ottimo,1:2*j_ottimo)*t)-eye(2*j_ottimo))...
    *B_smorzata(1:2*j_ottimo,1);
for i = 1:length(tempo)
    y_vals_step_esatta(i) = sum(y_esatta_step(tempo(i)));
end
y_vals_step_esatta  = y_vals_step_esatta + D_smorzata(1);

figure
plot(tempo, y_vals_step_esatta)
xlabel('Tempo(s)')
ylabel('F_x (N)')
grid on
exportgraphics(gcf, 'Risposta_esatta_step.pdf', 'ContentType', 'vector', 'BackgroundColor','none')
figure
plot(tempo, y_vals_step_esatta-iiiii')
grid on
exportgraphics(gcf, 'Confronto_step_esatta.pdf', 'ContentType', 'vector', 'BackgroundColor','none')

% Forma modale STEP
[V,Lambda] = eig(A_smorzata(1:2*j_ottimo,1:2*j_ottimo)); % V matrice di autovettori di A 
   
lambda = diag(Lambda) % lambda matrice diag di autovalori
lambda_inv = diag(1./(lambda));
V_inv = V\B_smorzata(1:2*j_ottimo,1);
z = @(t) lambda_inv*((expm(lambda*t)-eye(2*j_ottimo)) * V_inv);
y_modal_step= zeros(1,length(tempo));
u_step = input_pendulum(tempo, 1, 0);
for i = 1:length(tempo)
    eLambda = diag(exp(lambda*tempo(i)));
    y_modal_step(i) = sum(C_smorzata(1:j_ottimo,1:2*j_ottimo)*V*lambda_inv*...
        (eLambda-eye(2*j_ottimo)) * V_inv,1)+ D_smorzata(1)*u_step(tempo(i));
end
% per y_modal_impulse la formula non è esplicitamente sulle slide ma l'ho
% estrapolata per comparazione dalla forma generale: 
% risposta imp = C*expm(A*t)*B
% risposta step = C*A^-1 * (expm(A*t)-I)*B
% dove expm(A*t) = V*expm(lambda*t)*V^(-1)
plot(tempo,y_modal_step(1,:))

% Metodo integrazione STEP
x0 = zeros(2*j_ottimo, 1);  % considerando 8 pendoli devo assegnare le 8 
%codnizioni ini
u_step = input_pendulum(tempo, 1, 0)

odefun_step = @(t,x) A_smorzata(1:2*j_ottimo,1:2*j_ottimo)*x +...
    B_smorzata(1:2*j_ottimo,1)*u_step(t);

% ode 45
tic 
[t_n_step, theta_n_step] = ode45(odefun_step,tempo, x0);
t_45 = toc
F_step = @(t) (sum(C_smorzata(1:j_ottimo,1:2*j_ottimo)*theta_n_step',1))  +...
     D_smorzata(1)*u_step(t);

 % ode 45 modificato 
tic 
options = odeset(MaxStep = 0.01, RelTol = 1e-6, AbsTol = 1e-9); 
[t_n_step_m, theta_n_step_m] = ode45(odefun_step,tempo, x0, options);
t45m = toc
F_step_m = @(t) (sum(C_smorzata(1:j_ottimo,1:2*j_ottimo)*theta_n_step_m',1))  +...
     D_smorzata(1)*u_step(t);

% ode 78
tic
[t_n_step_78, theta_n_step_78] = ode78(odefun_step,tempo, x0);
t_78 = toc
F_step_78 = @(t) (sum(C_smorzata(1:j_ottimo,1:2*j_ottimo)*theta_n_step_78',1))  +...
     D_smorzata(1)*u_step(t);

figure
plot(tempo,F_step(tempo),'LineWidth',1.5) 
hold on
plot(tempo, F_step_m(tempo),'LineWidth',1.5)
hold on
plot(tempo, F_step_78(tempo),'LineWidth',1.5)
legend('ode45', 'ode45 modificato', 'ode78')
xlabel('tempo (s)')
ylabel('F_x (N)')


err_abs_45 = abs(y_vals_step_esatta - F_step(tempo)); 
max_err_45 = max(err_abs_45)


err_abs_45m = abs(y_vals_step_esatta - F_step_m(tempo)); 
max_err_45m = max(err_abs_45m) 

err_abs_78 = abs(y_vals_step_esatta - F_step_78(tempo)); 
max_err_78 = max(err_abs_78)

figure
plot(tempo, err_abs_45)
xlabel('tempo (s)')
ylabel('Errore assoluto')
figure
plot(tempo, err_abs_45m)
xlabel('tempo (s)')
ylabel('Errore assoluto')
figure
plot(tempo, err_abs_78)
xlabel('tempo (s)')
ylabel('Errore assoluto')

%  TASK 6 

j = j_ottimo;
            A_non_smorzata = zeros(2*j,2*j);
            B_non_smorzata = zeros(2*j,1);
            C_non_smorzata = zeros(j,2*j);
            D_non_smorzata = zeros(j,1);
            
            A_smorzata = zeros(2*j,2*j);
            C_smorzata = zeros(j,2*j);
        
            gamma = zeros(1,j);

        for i = 2:2:2*j

            gamma(1,i/2) = 2*damping_ratio*m_n(i/2)*l_n(i/2)^2*w_n(i/2);
        
            A_non_smorzata(i-1,i)=1;
            A_non_smorzata(i,i-1)=-g/l_n(i/2);
            B_non_smorzata(i,1) = -1/l_n(i/2);
            C_non_smorzata (i/2,i-1) = -m_n(i/2)*g;
            D_non_smorzata(i/2,1) = - (m0) ;
        
            A_smorzata(i-1,i)=1;
            A_smorzata(i,i-1)=-g/l_n(i/2);
            A_smorzata(i,i) = -gamma(1,i/2)/(m_n(i/2)*l_n(i/2)^2);
        
            C_smorzata (i/2,i-1) = m_n(i/2)*g; 
            C_smorzata(i/2,i) = gamma(1,i/2)/l_n(i/2); 
            B_smorzata = B_non_smorzata;
            D_smorzata = D_non_smorzata;
        end 
        
        C_smorzata_tot = ones(1,j_ottimo) * C_smorzata; % sommo tutti i contributi di F
        % relativi ai singoli pendoli 
        D_smorzata_tot = D_smorzata(1); % D (m0) la sommo solo una volta come
        % da definizione di F nel libro sullo sloshing

        s = tf('s') % definisce la s di Laplace
        sys_freq = ss(A_smorzata,B_smorzata, C_smorzata_tot,D_smorzata_tot)
        G_freq = tf(sys_freq)
        G_freq_spostamento = G_freq*s^2
        
        % va moltiplicato per s^2 perche ora l'ingresso 
        % non è piu l'accelerazione bensi lo spostamento, dividendo per
        % tale termine, la pendenza dell'ampiezza nel grafico di bode
        % risulterà -40db/dec mentre si avra un ritardo di fase iniziale
        % -180° 
        % G = ouscita/ingresso = F / x  
        % con F = C*x(t) + D*u e x(t) = (sI-A)^(-1) * B * u

        % quindi F = C*(sI-A)^(-1)*B*u + D*u --> in Laplace ottengo

        % G = (C*(sI-A)^(-1)*B + D)*x_dot_dot / (x_dot_dot/s^2)
        
        % G = (C*(sI-A)^(-1)*B + D)*s^2 che corrisponde alla tf con
        % ingresso accelerazione * s^2

       
        bode(G_freq_spostamento)
       hold on
      
      %  Modello Frozen Fluid (strano, magari ho interpretato male)

      M_tot = sum(m_n) + m0; % massa totale del sistema che NON partecipa allo sloshing

      G_no_slosh = M_tot*s^2 % funzione di trasferimento è M*x_dot_dot/x 
      % in Laplace diventa M*s^2x/x quindi M*s^2

      bode(G_no_slosh) 


% TASK 7 (Modello Simulink)
x0 = [0 0 0 0 0 0 0 0]';
figure
layout = tiledlayout(2,1);

% con l'utilizzo di spazio agli stati
Model_1 = sim('Model_31_01_1.slx');

nexttile;
plot(Model_1.tout,Model_1.u,'LineWidth', 2); % accelerazione
xlabel('Tempo [s]');
ylabel('Accelerazione [m/s^2]');
title(sprintf('Accelerazione applicata sul tank'));
grid on;
nexttile;
plot(Model_1.tout,Model_1.force,'LineWidth', 2); % forza
xlabel('Tempo [s]');
ylabel('Forza [N]');
title(sprintf('Forza dovuta all''accelerazione sul tank'));
grid on;


% Confronto tra i vari metodi risolutivi
figure(20)
plot(Model_1.tout,Model_1.force-iiiii,'LineWidth', 2);

figure(30)
plot(Model_1.tout,Model_1.force-y_vals_step_esatta','LineWidth', 2);

figure(40)
plot(Model_1.tout,Model_1.force-y_vals_heaviside_step','LineWidth', 2);


% function utilizzate:
% function input impulso 
function  u = input_pendulum_impulse(t, amp, t0)
    % 
    % if  t==t0
    %      force = amp;
    %  else 
    %      force = 0; 
    %  end
    %  u = force; 

    u = @(t) (t < t0 ) *amp; % Impulso breve di durata ~1ms
end

% function input step
function  u = input_pendulum(t, amp, t0)
 
if t>= t0
     force = amp;
 else 
     force = 0; 
 end
 u = @(t) force; 
end


% function ODE 1 pendolo

function [xdot]= pendulumODE1(t,x, gamma, w_n,l_n, input_pendulum)

x1= x(1); 
x2 = x(2);

u = input_pendulum(t, 1, 0);

xdot1 = x2;
xdot2 = -gamma(1)*x2 - w_n(1) * sin(x1) - (1/l_n(1))*cos(x1) * u;

xdot = [xdot1; xdot2]; 
% riflettere se ha senso mettere tempo di 10s 

end

% function ODE 2 pendolo

function [xdot]= pendulumODE2(t,x, gamma, w_n,l_n)

x1= x(1); 
x2 = x(2);

u = input_pendulum_impulse(t,1,-1)
disp(['t = ', num2str(t), ', u = ', num2str(u)]); % Stampa del valore di impulso
xdot1 = x2;
xdot2 = -gamma(1)*x2 - w_n(1) * sin(x1) - (1/l_n(1))*cos(x1) * u;

xdot = [xdot1; xdot2]; 

end
