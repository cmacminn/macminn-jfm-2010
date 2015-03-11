function [xCs,tCs,gCs] = hyp_case6c2_crits(M,gamma,Nf,Ns,gmin)
    % Calculate key values for case 6c2
    %   Note that Ns/Nf-->-infty (negative slope only) is contained here
    
    % 6c2:  peak forms first and eats R0 before shock gets there, then the two collide
    % LR_RL, RR_RS, R_p, R_L0, R_LS, R_LL
    % p_R0 happens either before or after RR_RS, but is a non-event
    
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
    
    % Event 2a: peak forms at xLR_RL
    g2a = M;
    t2a = tLR_RL;
    x2a = xLR_RL;
    
    % Event 2b: shock forms at xRS
    t2b = tRR_RS;
    x2b = xRR_RS;
    g2b = g_s;
    
    % Event 2c: shock collides with peak
    g = sym('g','positive');
    t2c_p = 1 + (2*(1-gamma)*M/(gamma*g^2))*((-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)/g^2)^-1);
    dt2c_p = diff(t2c_p);
    t2c_p = inline(vectorize(t2c_p));
    dt2c_p = inline(vectorize(dt2c_p));
    
    c1 = (M*Nf+(M/(M-1))*Ns);
    c2 = ((1/(M-1))*Ns);
    c3 = M;
    c4 = (1/(1-gamma));
    gm0 = g2b;
    tm0 = t2b;
    t_err = 1;
    g2c_next = g_0+0.01*(M-g_0);
    while abs(t_err) > 1E-10
        g2c = g2c_next;
        gmF = g2c;
        R = t2c_p(g2c) - hyp_shock(gamma,c1,c2,c3,c4,gm0,tm0,gmF);
        t_err = R/t2c_p(g2c);
        g2c_next = g2c_next - 0.5*(1/dt2c_p(g2c))*R;
    end
    
    t2c = t2c_p(g2c);
    x2c = (M./(g2c.^2)) + (-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(g2c.^2))*(t2c-1);
    g2c = g2c;
    
    % Event 4a:  shock collides with xRS
    g4a = max(g_s,gmin);
    c1 = (1/(1-gamma))*(M*Nf+(M/(M-1))*Ns);
    c2 = (1/(1-gamma))*((1/(M-1))*Ns);
    c3 = -M;
    c4 = 1;
    gm0 = g2c;
    tm0 = t2c;
    gmF = g4a;
    t4a = hyp_shock(gamma,c1,c2,c3,c4,gm0,tm0,gmF);
    x4a = -(M./(g4a.^2)) + (1/(1-gamma))*(-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(g4a.^2)).*(t4a-1.);
    
    if gmin<g4a
        % Event 4b:  shock collides with xRR
        g4b = gmin; % Set a cut-off shock height at which we say the plume is completely trapped
        c1 = (M*Nf+(M/(M-1))*Ns);
        c2 = (1/(M-1))*Ns;
        c3 = -M;
        c4 = 1/(1-gamma);
        gm0 = g4a;
        tm0 = t4a;
        gmF = g4b;
        t4b = hyp_shock(gamma,c1,c2,c3,c4,gm0,tm0,gmF);
        x4b = -(M./(g4b.^2)) + (-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(g4b.^2)).*(t4b-1.);
    else
        g4b = g4a;
        t4b = t4a;
        x4b = x4a;
    end
    
    xCs = [x1,x2a,x2b,x2c,x4a,x4b];
    tCs = [t1,t2a,t2b,t2c,t4a,t4b];
    gCs = [g1,g2a,g2b,g2c,g4a,g4b];
    
end