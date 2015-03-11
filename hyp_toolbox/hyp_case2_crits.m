function [xCs,tCs,gCs] = hyp_case2_crits(M,gamma,Nf,Ns,gmin)
    % Calculate key values for case 2
    
    % ---------------------------------------------------------
    % Key Values
    % ---------------------------------------------------------
    % Injection
    t1 = 1;
    x1 = M;
    g1 = 1;
    % Retreat, part a
    t2a = 1. + 2.*(1.-gamma)/(gamma*(Nf-Ns));
    x2a = (2.-gamma)/(gamma*M);
    g2a = M;
    % Retreat, part b
    t2b = 1 + (1-gamma)*(M-1)/((M-1)*Nf+Ns);
    x2b = -Ns/((M-1)*Nf+Ns);
    g2b = sqrt(M*((M-1)*Nf/Ns+1)*(1-2/gamma));
    % Sweep
    % See below
    
    % Calculate x4, t4 from the analytical solution
    g4 = min(g2b,gmin); % Set a cut-off shock height at which we say the plume is completely trapped
    c1 = (M*Nf+(M/(M-1))*Ns);
    c2 = ((1/(M-1))*Ns);
    c3 = M;
    c4 = 1/(1-gamma);
    gm0 = g2b;
    tm0 = t2b;
    gmF = g4;
    
    t4 = hyp_shock(gamma,c1,c2,c3,c4,gm0,tm0,gmF);
    x4 = (M/(gmF.^2)) + (-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gmF.^2)).*(t4-1.);
    
    xCs = [x1,x2a,x2b,x4];
    tCs = [t1,t2a,t2b,t4];
    gCs = [g1,g2a,g2b,g4];
    
end