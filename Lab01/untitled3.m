


%% 





% function [xdot]= pendulumODE1(t,x, gamma, w_n,l_n, input_pendulum)
% 
% x1= x(1); 
% x2 = x(2);
% 
% u = input_pendulum(t, 1, 0);
% 
% xdot1 = x2;
% xdot2 = -gamma(1)*x2 - w_n(1) * sin(x1) - (1/l_n(1))*cos(x1) * u;
% 
% xdot = [xdot1; xdot2]; 
% 
% 
% % riflettere se ha senso mettere tempo di 10s 
% 
% end
% 


%% 
clc
close all
gamma(1,1) = 2*damping_ratio*m_n(1)*l_n(1)^2*w_n(1);
odefun = @(t,x) pendulumODE1(t,x, gamma, w_n,l_n, @input_pendulum);

[t_n, theta_n] = ode45(odefun,tempo, [0 0]);

figure
plot(t_n,theta_n(:,1)); 
hold on 
plot(t_n,theta_n(:,2));

F = m_n(1)*g.*sin(theta_n(:,1)).*cos(theta_n(:,1)) + m_n(1)*l_n(1).*sin(theta_n(:,1)).*theta_n(:,2).^2 ...
    +(gamma(1)/l_n(1)).*theta_n(:,2).*cos(theta_n(:,1))- m0-m_n(1) + m_n(1)*(cos(theta_n(:,1))).^2;

figure
plot(tempo(1,1:(end)),F(:,1),'k')

%% 
clc

u_step = input_pendulum(tempo, 1, 0)

odefun_step = @(t,x) A_smorzata(1:2,1:2)*x + B_smorzata(1:2,1)*u_step

[t_n_step, theta_n_step] = ode45(odefun_step,tempo, [0 0]);

F_step = C_smorzata(1,1:2)*theta_n_step' + D_smorzata(1)*u_step

plot(tempo,F_step)

%% impulse

clc
close all

u_impulse = input_pendulum_impulse(tempo, 1, 0.01);

odefun_impulse = @(t,x) A_smorzata(1:2,1:2)*x + B_smorzata(1:2,1)*u_impulse(t)

% [t_n_impulse, theta_n_impulse] = ode45(odefun_impulse,tempo,
% B_smorzata(1:2,1)); % scegliendo B come condizioni iniziali, si avvia il
% sistema nella direzione dell'impulso, con 0,0 invece l'impulso (unitario)
% non è sufficiente per inenscare il sistema oppure lo smorzamento è
% eccessivo 

[t_n_impulse, theta_n_impulse] = ode45(odefun_impulse,tempo,[0,0]);


for ii = 1:length(tempo)
    F_impulse(ii) = C_smorzata(1,1)*theta_n_impulse(ii,1)+C_smorzata(1,2)...
        *theta_n_impulse(ii,2) + D_smorzata(1)* u_impulse(ii/100);
end

figure(7)
plot(tempo,F_impulse)




% %% risposta IMPULSE
% clc
% close all
% 
% odefun_step = @(x,t) pendulumODE2(x,t, gamma(1), w_n(1),l_n(1))
% 
% [t_n_impulse, theta_n_impulse] = ode45(odefun_step,tempo, [0 0]);
% 
% figure
% plot(t_n_impulse,theta_n_impulse(:,1)); 
% plot(t_n_impulse,theta_n_impulse(:,2)); 

 %% 
