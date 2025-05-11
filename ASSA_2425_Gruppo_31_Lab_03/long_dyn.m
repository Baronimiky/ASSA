function [xdot] = long_dyn(x,t, m, I_yy, rho, V, S,delta_z_t, alpha_t,c,g,coeff,delta_e,T)

    u = x(1); 
    w = x(2);
    q = x(3);
    theta = x(4);
    h = x(5);

    alpha = theta; %theta = gamma + alpha ma gamma = 0
    
    
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

    xdot = zeros(5,1);
    xdot(1) = -q*w -g*sin(theta) + 1/m * X_A + 1/m*X_T; 
    xdot(2) = q*u + g*cos(theta) + 1/m * Z_A + 1/m *Z_T;
    xdot(3) = 1/I_yy * (M_A+M_T);
    xdot(4) = q; 
    xdot(5) = u * sin(theta) - w* cos(theta); 

end
