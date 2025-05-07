%% Soluzione integrazione con diversi metodi da ode45
%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% STEP %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%


clc
j_ottimo = 4; 
x0 = zeros(2*j_ottimo, 1);  % considerando 8 pendoli devo assegnare le 8 
%codnizioni ini
u_step = input_pendulum(tempo, 1, 0)

odefun_step = @(t,x) A_smorzata(1:2*j_ottimo,1:2*j_ottimo)*x +...
    B_smorzata(1:2*j_ottimo,1)*u_step(t)
tic 
[t_n_step, theta_n_step] = ode45(odefun_step,tempo, x0);
t_45 = toc
F_step = @(t) (sum(C_smorzata(1:j_ottimo,1:2*j_ottimo)*theta_n_step',1))  +...
     D_smorzata(1)*u_step(t);
 % anche in questo caso sommo tutti i contributi dei singoli pendoli stando
 % attento a non sommare molteplici volte il contributo costante della
 % massa m0 che chiaramente deve essere considerato una sola volta

tic 
options = odeset(MaxStep = 0.1, RelTol = 1e-6, AbsTol = 1e-9); 
[t_n_step_m, theta_n_step_m] = ode45(odefun_step,tempo, x0, options);
t45m = toc
F_step_m = @(t) (sum(C_smorzata(1:j_ottimo,1:2*j_ottimo)*theta_n_step_m',1))  +...
     D_smorzata(1)*u_step(t);

plot(tempo,F_step(tempo), tempo,F_step_m(tempo));


err_abs_45 = abs(y_vals_step_esatta - F_step(tempo)); 
max_err_45 = max(err_abs_45) 

err_rel_45 = abs(err_abs_45 / y_vals_step_esatta);
max_err_rel_45 = max(err_rel_45)

err_abs_45m = abs(y_vals_step_esatta - F_step_m(tempo)); 
max_err_45m = max(err_abs_45m) 

err_rel_45m = abs(err_abs_45 / y_vals_step_esatta);
max_err_rel_45m = max(err_rel_45m)

plot(tempo, err_abs_45, tempo, err_abs_45m)





%% 

%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% IMPULSE %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%


u_impulse = input_pendulum_impulse(tempo, 1, 0.01);
odefun_impulse = @(t,x) A_smorzata(1:2,1:2)*x; % B_smorzata(1:2,1)*u_impulse(t);
% [t_n_impulse, theta_n_impulse] = ode45(odefun_impulse,tempo,
% B_smorzata(1:2,1)); % scegliendo B come condizioni iniziali, si avvia il
% sistema nella direzione dell'impulso, con 0,0 invece l'impulso (unitario)
% non è sufficiente per inenscare il sistema oppure lo smorzamento è
% eccessivo 
tic;


[t_n_impulse_45m, theta_n_impulse_45m] = ode45(odefun_impulse,tempo,B_smorzata(1:2,1),options);
tempo_esecuzione_45m = toc

for ii = 1:length(tempo)
    F_impulse_45m(ii) = C_smorzata(1,1)*theta_n_impulse_45m(ii,1)+C_smorzata(1,2)...
        *theta_n_impulse_45m(ii,2); %+ D_smorzata(1) u_impulse(ii/100);
end


tic
[t_n_impulse_45, theta_n_impulse_45] = ode45(odefun_impulse,tempo,B_smorzata(1:2,1));
tempo_esecuzione_45 = toc

for ii = 1:length(tempo)
    F_impulse_45(ii) = C_smorzata(1,1)*theta_n_impulse_45(ii,1)+C_smorzata(1,2)...
        *theta_n_impulse_45(ii,2); %+ D_smorzata(1) u_impulse(ii/100);
end

figure(7)
plot(tempo,F_impulse_45, tempo, F_impulse_45m)
legend('45', '45m')


error_abs_45 = abs(y_vals_impulse - F_impulse_45);
error_abs_45m = abs(y_vals_impulse - F_impulse_45m);

max_err_45 = max(error_abs_45)
max_err_45m = max(error_abs_45m)