function [xCs,tCs,gCs] = hyp_case0b_crits(M,gamma,Nf,Ns,gmin)
    % Calculate key values for case 0
    %   Note that Ns/Nf --> infty (positive slope only) is contained here
    
    % 0b:  a shock forms on the left, then a peak
    % LL_LS, RL_LR, L_p, L_RS, L_RR
    
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
    x2a = xLS;
    g2a = g_s;
    
    % Retreat, part b
    t2b = tRL_LR;
    x2b = xRL_LR;
    g2b = M;
    
    % Retreat, part c
    c1 = (M*Nf+(M/(M-1))*Ns);
    c2 = ((1/(M-1))*Ns);
    c3 = -M;
    c4 = 1/(1-gamma);
    gm0 = g_s;
    tm0 = tLL_LS;
    
    g = sym('g','positive');
    t2c_p = 1 - (2*M*(1-gamma)/(gamma*g^2))*((-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)/g^2)^-1);
    dt2c_p = diff(t2c_p);
    t2c_p = inline(vectorize(t2c_p));
    dt2c_p = inline(vectorize(dt2c_p));
    
    t_err = 1;
    g2c_next = (g_s+M)/2;
    
    limiter = 0.5;
    while abs(t_err) > 1E-10
        g2c = g2c_next;
        R = t2c_p(g2c) - hyp_shock(gamma,c1,c2,c3,c4,gm0,tm0,g2c);
        g2c_next = g2c_next - limiter*(1/dt2c_p(g2c))*R;
        if ( imag(g2c_next)==0 && g2c_next<=M )
            t_err = R/t2c_p(g2c);
        else
            t_err = 1;
            g2c_next = (g_s+M)/2;
            limiter = limiter/2;
        end
    end
    
    t2c = t2c_p(g2c);
    x2c = -(M/(g2c^2))*t1 + (-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)/(g2c^2))*(t2c-t1);    
    
    
    % Sweep, part a
    c1 = (1/(1-gamma))*(M*Nf+(M/(M-1))*Ns);
    c2 = (1/(1-gamma))*((1/(M-1))*Ns);
    c3 = M;
    c4 = 1;
    gm0 = g2c;
    tm0 = t2c;
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
    
    xCs = [x1,x2a,x2b,x2c,x4a,x4b];
    tCs = [t1,t2a,t2b,t2c,t4a,t4b];
    gCs = [g1,g2a,g2b,g2c,g4a,g4b];
    
end