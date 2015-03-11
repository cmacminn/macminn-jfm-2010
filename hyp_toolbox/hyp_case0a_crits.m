function [xCs,tCs,gCs] = hyp_case0a_crits(M,gamma,Nf,Ns,gmin)
    % Calculate key values for case 0
    %   Note that Ns/Nf --> infty (positive slope only) is contained here
    
    % 0a:  a shock forms on the left and hits the bottom before a peak forms (no peak)
    % LL_LS, L_LR, L_RL, L_RS, L_RR
    
    % ---------------------------------------------------------
    % Key Values
    % ---------------------------------------------------------
    g_s = sqrt(M*(M-1)*Nf/Ns+M);
    xLS = -1/((M-1)*Nf/Ns+1);
    xRS =  1/((M-1)*Nf/Ns+1);
    
    % Injection
    t1 = 1;
    x1 = M;
    g1 = 1;
    
    % Collision times
    tLL_LS = 1 + (1-gamma)*(M-1)/((M-1)*Nf+Ns);
    xLL_LS = xLS;
    tRL_LR = 1 - 2*((1-gamma)/gamma)/(Nf-Ns);
    xRL_LR = -(2-gamma)/(M*gamma);
    
    % Retreat, part a
    t2a = tLL_LS;
    x2a = xLL_LS;
    g2a = g_s;
    
    % Retreat, part b
    c1 = (M*Nf+(M/(M-1))*Ns);
    c2 = ((1/(M-1))*Ns);
    c3 = -M;
    c4 = 1/(1-gamma);
    gm0 = g_s;
    tm0 = tLL_LS;
    gmF = M;
    t2b = hyp_shock(gamma,c1,c2,c3,c4,gm0,tm0,gmF);
    g2b = M;
    x2b = -(M/(g2b^2)) + (-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)/(g2b^2))*(t2b-1);
    
    % Chase
    t3 = (-(Nf-Ns)+2*(1-gamma)-(Nf-Ns)*(1-gamma)*(t2b-1)+M*Nf*t2b)/((M-1)*Nf+Ns);
    x3 = x2b + (Nf/(1-gamma))*(t3-t2b);
    g3 = M;
    
    % Sweep, part a
    c1 = (1/(1-gamma))*(M*Nf+(M/(M-1))*Ns);
    c2 = (1/(1-gamma))*((1/(M-1))*Ns);
    c3 = M;
    c4 = 1;
    gm0 = g3;
    tm0 = t3;
    gmF = g_s;
    t4a = hyp_shock(gamma,c1,c2,c3,c4,gm0,tm0,gmF);
    x4a = xRS;
    g4a = g_s;    
    
    % Sweep, part b
    g4b = min(gmin,g4a); % Cut-off shock height at which we say the plume is completely trapped
    
    c1 = (M*Nf+(M/(M-1))*Ns);
    c2 = ((1/(M-1))*Ns);
    c3 = M;
    c4 = 1/(1-gamma);
    gm0 = g4a;
    tm0 = t4a;
    gmF = g4b;
    t4b = hyp_shock(gamma,c1,c2,c3,c4,gm0,tm0,gmF);
    x4b = (M/(gmF^2)) + (-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)/(gmF^2)).*(t4b-1.);
    
    xCs = [x1,x2a,x2b,x3,x4a,x4b];
    tCs = [t1,t2a,t2b,t3,t4a,t4b];
    gCs = [g1,g2a,g2b,g3,g4a,g4b];
    
end