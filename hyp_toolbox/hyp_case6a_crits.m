function [xCs,tCs,gCs] = hyp_case6a_crits(M,gamma,Nf,Ns,gmin)
    % Calculate key values for case 6a
    %   Note that Ns/Nf-->-infty (negative slope only) is contained here
    
    % 6a:  a shock forms on the right and hits the bottom before a peak forms (no peak)
    % RR_RS, R_R0, R_RL, R_LR, R_L0, R_LS, R_LL
    
    % ---------------------------------------------------------
    % Key Values
    % ---------------------------------------------------------
    g_s = sqrt(M*(M-1)*Nf/Ns+M);
    g_0 = M*(M-1)*Nf/Ns + M; % zero of the flux function
    xLS = -1/((M-1)*Nf/Ns+1);
    xRS =  1/((M-1)*Nf/Ns+1);
    
    % Injection
    t1 = 1;
    x1 = M;
    g1 = 1;
    
    % Collision times
    tRR_RS = 1 - (1-gamma)*(M-1)/((M-1)*Nf+Ns);
    xRR_RS = xRS;
    tLR_RL = 1 + 2*((1-gamma)/gamma)/(Nf-Ns);
    xLR_RL = (2-gamma)/(M*gamma);
    
    % Event 2a:  shock forms at xRS
    t2a = tRR_RS;
    x2a = xRR_RS;
    g2a = g_s;
    
    % Event 2b:  shock changes direction at xR0
    g2b = g_0;
    c1 = (M*Nf+(M/(M-1))*Ns);
    c2 = ((1/(M-1))*Ns);
    c3 = M;
    c4 = 1/(1-gamma);
    gm0 = g2a;
    tm0 = t2a;
    gmF = g2b;
    t2b = hyp_shock(gamma,c1,c2,c3,c4,gm0,tm0,gmF);
    x2b = (M./(g2b^2)) + (-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(g2b^2))*(t2b-1);
    
    % Event 2c:  shock collides with xRL
    g2c = M;
    c1 = (M*Nf+(M/(M-1))*Ns);
    c2 = ((1/(M-1))*Ns);
    c3 = M;
    c4 = 1;
    gm0 = g2b;
    tm0 = t2b;
    gmF = g2c;
    t2c = hyp_shock(gamma,c1,c2,c3,c4,gm0,tm0,gmF);
    x2c = (M./(g2c^2)) + (-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(g2c^2))*(t2c-1);
    
    % Event 3:  shock collides with xLR
    t3 = (1/M+(1/(1-gamma))*(Nf-Ns)/M + x2c - Nf*t2c)/(-Nf+(1/(1-gamma))*(Nf-Ns)/M);
    x3 = x2c + Nf*(t3-t2c);
    g3 = M;
    
    % Event 4a:  shock changes direction at xL0
    g4a = g_0;
    c1 = (1/(1-gamma))*(M*Nf+(M/(M-1))*Ns);
    c2 = (1/(1-gamma))*((1/(M-1))*Ns);
    c3 = -M;
    c4 = (1-gamma);
    gm0 = g3;
    tm0 = t3;
    gmF = g4a;
    t4a = hyp_shock(gamma,c1,c2,c3,c4,gm0,tm0,gmF);
    x4a = -(M./(g4a.^2)) + (1/(1-gamma))*(-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(g4a.^2)).*(t4a-1);
    
    % Event 4b:  shock collides with xRS
    g4b = g_s;
    c1 = (1/(1-gamma))*(M*Nf+(M/(M-1))*Ns);
    c2 = (1/(1-gamma))*((1/(M-1))*Ns);
    c3 = -M;
    c4 = 1;
    gm0 = g4a;
    tm0 = t4a;
    gmF = g4b;
    t4b = hyp_shock(gamma,c1,c2,c3,c4,gm0,tm0,gmF);
    x4b = xLS;
    
    % Event 4c:  shock collides with xLL
    g4c = min(g4b,gmin); % Set a cut-off shock height at which we say the plume is completely trapped
    c1 = (M*Nf+(M/(M-1))*Ns);
    c2 = (1/(M-1))*Ns;
    c3 = -M;
    c4 = 1/(1-gamma);
    gm0 = g4b;
    tm0 = t4b;
    gmF = g4c;
    t4c = hyp_shock(gamma,c1,c2,c3,c4,gm0,tm0,gmF);
    x4c = -(M./(g4c.^2)) + (-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(g4c.^2)).*(t4c-1.);
    
    xCs = [x1,x2a,x2b,x2c,x3,x4a,x4b,x4c];
    tCs = [t1,t2a,t2b,t2c,t3,t4a,t4b,t4c];
    gCs = [g1,g2a,g2b,g2c,g3,g4a,g4b,g4c];
    
end