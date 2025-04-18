clear all
close all
clc

% da pag 43 (fr. nat.) a pag. 48 (equazioni)
% oss1 = la prima massa del primo pendolo con theta=0 è sotto il CoG
% oss2= sistema di riferimento differente tra slide e doc; le formule fanno
% riferimento al DOGE
% guardare introduzione!!! todos


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

% ATTENZIONE: il modal damping ratio non è direttamente il gamma, ma è
% l'epsilon di mecc aerospaziale 2 quindi --->

% gamma = 2*damping_ratio*m(i)*l(i)^2*w_n(i)

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

% creo modello spazio degli stati attraverso sys e valuto nell'intervallo
% la risposta all'impulso con la funzione integrata di matlab (considero
% solo il primo pendolo come richiesto)
% 

sys = ss(A_smorzata(1:2,1:2),B_smorzata(1:2,1),C_smorzata(1,1:2),D_smorzata(1));
tempo  = 0:0.01:10;

figure
[k, t_out] = impulse(sys, tempo);%% MATLAB esclude il contributo di D (legato al delta di
% Dirac) perche infintesimale , t_out serve per non plottare 

%
clc
-gamma(1)/l_n(1)^2 - gamma(1)/(m_n(1)*l_n(1)^2)*m0
(m_n(1)+ m0)*g/l_n(1)
G1 = tf(sys)

% confronto con soluzione analitica (intesa come esatta, non come a mano)
 
y = @(t) C_smorzata(1,1:2)*expm(A_smorzata(1:2,1:2)*t)*B_smorzata(1:2,1) % + D_smorzata(1);
%se sommassi in questo modo, starei considerando la risposta ad uno scalino
%non impulso; D va moltiplicato per il delta di Dirac, ma non è
%visualizzabile in matlab

y_vals = arrayfun(y, tempo); % Calcola y per ogni valore di tempo
y_vals(1) = y_vals(1) + D_smorzata(1); % per simulare il delta di Dirac ma non so quanto sia corretto

figure(2)
plot(tempo, y_vals);





%% risoluzione analitica con heaviside


clear
clc


tempo  = 0:0.01:10;

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

c = gamma(1) / (m_n(1)*l_n(1)^2);

temp = sqrt(w_n(1)^2-0.25*c^2)

y_heaviside = @(t) -c*m_n(1)*exp(-c/2*t)*cos(temp*t)+(0.5*c^2*m_n(1)-w_n(1)^2*m_n(1))/(temp)*exp(-c/2*t)* sin(temp*t)

y_vals_heaviside = arrayfun(y_heaviside, tempo);

% plot(tempo, y_vals_heaviside);

j_ottimo = 0;
toll = 1e-16;

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
    
        C_smorzata (i/2,i-1) = m_n(i/2)*g; % segno invertito rispetto equazioni analitiche 
        C_smorzata(i/2,i) = gamma(1,i/2)/l_n(i/2); 
        B_smorzata = B_non_smorzata;
        D_smorzata = D_non_smorzata;
    end 
    

    sys = ss(A_smorzata,B_smorzata,C_smorzata,D_smorzata);
    

    [k, t_out] = impulse(sys,tempo);

    y = @(t) C_smorzata*expm(A_smorzata*t)*B_smorzata;
    y_vals = arrayfun(y, tempo,'UniformOutput', false); % Calcola y per ogni valore di tempo
    % crea un array di celle che a loro volta contengono un vettore 10x1
    % (ogni cella contiene 
    
        for h = 2:size(k,1)
            if(abs(y_vals{h}(j)-k(h,j))<toll && (y_vals{h}(j)-k(h,j)~= 0))
                j_ottimo = j;
                break
            end
        end
    
        if(j_ottimo==j)
            break
        end
    
end

    if (j_ottimo == 0)
        disp('tolleranza non rispettata');
    end

%% 


y_vals{10}
size(k,1)

y_vals_2 = arrayfun(y, tempo,'UniformOutput', false)

y_vals_2{1}(3)
%% 
%% 
%% 
%% 


clear
clc





tempo  = 0:0.01:10;

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


y_vals_heaviside=zeros(1001,10);
% plot(tempo, y_vals_heaviside);

j_ottimo = 0;
toll = 1e-1; %non va bene, è irrealistica e probabilmente qualsiasi metodo
%numerico ha un errore intrinseco maggiore --> trovare un altro criterio

for j = 1:n
    
        A_non_smorzata = zeros(2*j,2*j);
        B_non_smorzata = zeros(2*j,1);
        C_non_smorzata = zeros(j,2*j);
        D_non_smorzata = zeros(j,1);
        
        A_smorzata = zeros(2*j,2*j);
        C_smorzata = zeros(j,2*j);
        
        gamma = zeros(1,j);
    for i = 2:2:2*j

        c = gamma(j) / (m_n(j)*l_n(j)^2);
        temp = sqrt(w_n(j)^2-0.25*c^2);  
        y_heaviside = @(t) -c*m_n(1)*exp(-c/2*t)*cos(temp*t)+(0.5*c^2*m_n(1)-w_n(1)^2*m_n(1))/(temp)*exp(-c/2*t)* sin(temp*t);
 
        y_vals_heaviside(:,j) = arrayfun(y_heaviside, tempo);

        gamma(1,i/2) = 2*damping_ratio*m_n(i/2)*l_n(i/2)^2*w_n(i/2);
    
        A_non_smorzata(i-1,i)=1;
        A_non_smorzata(i,i-1)=-g/l_n(i/2);
        B_non_smorzata(i,1) = -1/l_n(i/2);
        C_non_smorzata (i/2,i-1) = -m_n(i/2)*g;
        D_non_smorzata(i/2,1) = - (m0) ;
    
        A_smorzata(i-1,i)=1;
        A_smorzata(i,i-1)=-g/l_n(i/2);
        A_smorzata(i,i) = -gamma(1,i/2)/(m_n(i/2)*l_n(i/2)^2);
    
        C_smorzata (i/2,i-1) = m_n(i/2)*g; % segno invertito rispetto equazioni analitiche 
        C_smorzata(i/2,i) = gamma(1,i/2)/l_n(i/2); 
        B_smorzata = B_non_smorzata;
        D_smorzata = D_non_smorzata;
    end 
    

    sys = ss(A_smorzata,B_smorzata,C_smorzata,D_smorzata);
    

    [y_impulse, t_out] = impulse(sys,tempo);

    y = @(t) C_smorzata*expm(A_smorzata*t)*B_smorzata;
    y_vals = arrayfun(y, tempo,'UniformOutput', false); % Calcola y per ogni valore di tempo
    % crea un array di celle che a loro volta contengono un vettore jx1
    % (ogni cella contiene i valori di F per ogni pendolo all'istante t)
    
        % for h = 2:size(k,1)
        %     if(abs(y_vals{h}(j)^2-k(h,j)^2)<toll && (y_vals{h}(j)-k(h,j)~= 0))
        %         j_ottimo = j;
        %         break
        %     end
        % end

    rmse = sqrt(mean((y_vals_heaviside(:,j) - y_impulse(:,j)).^2))
        if rmse < toll
            j_ottimo = j;
            break;
        end
    
        if(j_ottimo==j)
            break
        end
    
end

    if (j_ottimo == 0)
        disp('tolleranza non rispettata');
    end

%% 


y_vals{10}
size(y_impulse,1)

y_vals_2 = arrayfun(y, tempo,'UniformOutput', false)
%% 

y_vals_heaviside

%% risposta STEP

odefun = @(x,t) pendulumODE1(x,t, gamma, w_n,l_n, @input_pendulum)

[t_n, theta_n] = ode45(odefun,tempo, [0 0]);

figure
plot(t_n,theta_n(:,1)); 
hold on 
plot(t_n,theta_n(:,2))

%% risposta IMPULSE


odefun = @(x,t) pendulumODE2(x,t, gamma, w_n,l_n)

[t_n, theta_n] = ode45(odefun,tempo, [0 0]);

figure
plot(t_n,theta_n(:,1)); 
hold on 
plot(t_n,theta_n(:,2))
