function [xCs,tCs,gCs] = hyp_case4_crits(M,gamma,Nf,Ns,gmin)
    % Calculate key values for case 4
    
    % ---------------------------------------------------------
    % Key Values
    % ---------------------------------------------------------
    % Injection
    t1 = 1;
    x1 = M;
    g1 = 1;
    
    % Retreat, part a
    t2a = 1 + 2*((1-gamma)/gamma)*(1/(Nf-Ns));
    x2a = (1/M)*(2/gamma-1);
    g2a = M;
    
    % Retreat, part b
    t2b = 1 - (M-1)/((M-1)*Nf+Ns);
    x2b = Ns/((M-1)*Nf+Ns);
    g2b = sqrt( M*(M-1)*Nf/Ns + M - (2*M*(M-1)*(1-gamma)/(gamma*(t2b-1)))/Ns );
    
    % Sweep
    g4 = min(g2b,gmin); % Set a cut-off shock height at which we say the plume is completely trapped
    c1 = (1/(1-gamma))*(M*Nf+(M/(M-1))*Ns);
    c2 = (1/(1-gamma))*((1/(M-1))*Ns);
    c3 = -M;
    c4 = (1-gamma);
    gm0 = g2b;
    tm0 = t2b;
    gmF = g4;
    
    t4 = hyp_shock(gamma,c1,c2,c3,c4,gm0,tm0,gmF);
    x4 = -M/(g4^2) + (1/(1-gamma))*(-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)/(g4^2)).*(t4-1.);
    
    xCs = [x1,x2a,x2b,x4];
    tCs = [t1,t2a,t2b,t4];
    gCs = [g1,g2a,g2b,g4];
    
    xCs = [x1,x2a,x2b,x4];
    tCs = [t1,t2a,t2b,t4];
    gCs = [g1,g2a,g2b,g4];
    
end