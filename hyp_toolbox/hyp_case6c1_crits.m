function [xCs,tCs,gCs] = hyp_case6c1_crits(M,gamma,Nf,Ns,gmin)
    % Calculate key values for case 6c1
    %   Note that Ns/Nf-->-infty (negative slope only) is contained here
    
    % 6c1:  peak forms first, then shock forms, hits R0, and collides with peak
    % LR_RL, RR_RS, R_R0, R_p, R_L0, R_LS, R_LL
    
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
    
    % Event 2a:  peak forms at xLR_RL
    g2a = M;
    t2a = tLR_RL;
    x2a = xLR_RL;
    
    % Event 2b:  shock formation at xRS
    t2b = tRR_RS;
    x2b = xRR_RS;
    g2b = g_s;
    
    % Event 2c:  shock changes direction at xR0
    g2c = g_0;
    c1 = (M*Nf+(M/(M-1))*Ns);
    c2 = ((1/(M-1))*Ns);
    c3 = M;
    c4 = 1/(1-gamma);
    gm0 = g2b;
    tm0 = t2b;
    gmF = g2c;
    t2c = hyp_shock(gamma,c1,c2,c3,c4,gm0,tm0,gmF);
    x2c = (M/(g2c^2)) + (-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(g2c.^2))*(t2c-1);
    
    % Event 2d:  shock collides with peak
    g = sym('g','positive');
    t2d_p = 1 + (2*(1-gamma)*M/(gamma*g^2))*((-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)/g^2)^-1);
    dt2d_p = diff(t2d_p);
    t2d_p = inline(vectorize(t2d_p));
    dt2d_p = inline(vectorize(dt2d_p));
    
    c1 = (M*Nf+(M/(M-1))*Ns);
    c2 = ((1/(M-1))*Ns);
    c3 = M;
    c4 = 1;
    gm0 = g2c;
    tm0 = t2c;
    t_err = 1;
    g2d_next = (g_0+M)/2;
    while abs(t_err) > 1E-10
        g2d = g2d_next;
        gmF = g2d;
        R = t2d_p(g2d) - hyp_shock(gamma,c1,c2,c3,c4,gm0,tm0,gmF);
        t_err = R/t2d_p(g2d);
        g2d_next = g2d_next - 0.1*(1/dt2d_p(g2d))*R;
    end
    
    t2d = t2d_p(g2d);
    x2d = (M./(g2d.^2)) + (-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(g2d.^2))*(t2d-1);
    g2d = g2d;
    
    % Event 4a:  shock changes direction at xR_0
    g4a = max(g_0,gmin); % Set a cut-off shock height at which we say the plume is completely trapped
    c1 = (1/(1-gamma))*(M*Nf+(M/(M-1))*Ns);
    c2 = (1/(1-gamma))*((1/(M-1))*Ns);
    c3 = -M;
    c4 = (1-gamma);
    gm0 = g2d;
    tm0 = t2d;
    gmF = g4a;
    t4a = hyp_shock(gamma,c1,c2,c3,c4,gm0,tm0,gmF);
    x4a = -(M./(g4a.^2)) + (1/(1-gamma))*(-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(g4a.^2)).*(t4a-1.);
    
    if gmin<g4a
        % Event 4b:  shock collides with xRS
        g4b = max(g_s,gmin);
        c1 = (1/(1-gamma))*(M*Nf+(M/(M-1))*Ns);
        c2 = (1/(1-gamma))*((1/(M-1))*Ns);
        c3 = -M;
        c4 = 1;
        gm0 = g4a;
        tm0 = t4a;
        gmF = g4b;
        t4b = hyp_shock(gamma,c1,c2,c3,c4,gm0,tm0,gmF);
        x4b = -(M./(g4b.^2)) + (1/(1-gamma))*(-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(g4b.^2)).*(t4b-1.);
    else
        g4b = g4a;
        t4b = t4a;
        x4b = x4a;
    end
    
    if gmin<g4b
        % Event 4c:  shock collides with xRR
        g4c = gmin; % Set a cut-off shock height at which we say the plume is completely trapped
        c1 = (M*Nf+(M/(M-1))*Ns);
        c2 = (1/(M-1))*Ns;
        c3 = -M;
        c4 = 1/(1-gamma);
        gm0 = g4b;
        tm0 = t4b;
        gmF = g4c;
        t4c = hyp_shock(gamma,c1,c2,c3,c4,gm0,tm0,gmF);
        x4c = -(M./(g4c.^2)) + (-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(g4c.^2)).*(t4c-1.);
    else
        g4c = g4b;
        t4c = t4b;
        x4c = x4b;
    end
    
    xCs = [x1,x2a,x2b,x2c,x2d,x4a,x4b,x4c];
    tCs = [t1,t2a,t2b,t2c,t2d,t4a,t4b,t4c];
    gCs = [g1,g2a,g2b,g2c,g2d,g4a,g4b,g4c];
    
end