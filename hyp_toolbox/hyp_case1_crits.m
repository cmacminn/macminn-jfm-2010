function [xCs,tCs,gCs] = hyp_case1_crits(M,gamma,Nf,Ns,gmin)
    % Calculate key values for case 1
    %   Note that Ns/Nf = 0 (flow only, Ruben's solution) is contained here
    
    % ---------------------------------------------------------
    % Key Values
    % ---------------------------------------------------------
    % Injection
    t1 = 1;
    x1 = M;
    g1 = 1;
    % Retreat
    t2 = 1. + (1.-gamma)*(M-1.)/((M-1.)*Nf+Ns);
    x2 = -Ns/((M-1.)*Nf+Ns);
    g2 = M;
    % Chase
    t3 = 1.+(1.-gamma)*(M+1.)/((M-(1.-gamma))*Nf+(1.-gamma)*Ns);
    x3 = ((2.-gamma)*Nf-(1.-gamma)*Ns)/((M-(1.-gamma))*Nf+(1.-gamma)*Ns);
    g3 = M;
    % Sweep
    % See below
    
    % Calculate x4, t4 from the analytical solution
    g4 = gmin; % Cut-off shock height at which we say the plume is completely trapped
    
    c1 = (M*Nf+(M/(M-1))*Ns);
    c2 = ((1/(M-1))*Ns);
    c3 = M;
    c4 = 1/(1-gamma);
    gm0 = g3;
    tm0 = t3;
    gmF = g4;
    
    t4 = hyp_shock(gamma,c1,c2,c3,c4,gm0,tm0,gmF);
    x4 = (M./(gmF).^2) + (-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gmF.^2)).*(t4-1.);
    
    xCs = [x1,x2,x3,x4];
    tCs = [t1,t2,t3,t4];
    gCs = [g1,g2,g3,g4];
    
end