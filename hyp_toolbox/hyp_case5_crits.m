function [xCs,tCs,gCs] = hyp_case5_crits(M,gamma,Nf,Ns,gmin)
    % Calculate key values for case 5
    %   Note that this case only exists if M<sqrt(2/gamma-1)
    
    % ---------------------------------------------------------
    % Key Values
    % ---------------------------------------------------------
    % Injection
    t1 = 1;
    x1 = M;
    g1 = 1;
    
    % Retreat
    t2 = 1-(M-1)/((M-1)*Nf+Ns);
    x2 = M+(M*Nf+Ns)*(t2-1);
    g2 = M;
    
    % Chase
    t3 = (x2+1/M+(1/(1-gamma))*((Nf-Ns)/M)-Nf*t2)/((1/(1-gamma))*((Nf-Ns)/M)-Nf);
    x3 = x2+Nf*(t3-t2);
    g3 = M;
    
    % Sweep
    % Calculate x4, t4 from the analytical solution
    g4 = gmin; % Set a cut-off shock height at which we say the plume is completely trapped
    c1 = (1/(1-gamma))*(M*Nf+(M/(M-1))*Ns);
    c2 = (1/(1-gamma))*((1/(M-1))*Ns);
    c3 = -M;
    c4 = (1-gamma);
    gm0 = g3;
    tm0 = t3;
    gmF = g4;
    
    t4 = hyp_shock(gamma,c1,c2,c3,c4,gm0,tm0,gmF);
    x4 = -(M/(g4^2)) + (1/(1-gamma))*(-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)/(g4^2)).*(t4-1.);
    
    xCs = [x1,x2,x3,x4];
    tCs = [t1,t2,t3,t4];
    gCs = [g1,g2,g3,g4];
    
end