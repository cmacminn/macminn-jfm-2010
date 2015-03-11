function [xCs,tCs,gCs] = hyp_case3_crits(M,gamma,Nf,Ns,gmin)
    % Calculate key values for case 3
    %   Note that Ns/Nf = -(M-1) is contained here
    
    % ---------------------------------------------------------
    % Key Values
    % ---------------------------------------------------------
    % Injection
    t1 = 1;
    x1 = M;
    g1 = 1;
    
    % Retreat
    t2 = 1. + 2.*(1.-gamma)/(gamma*(Nf-Ns));
    x2 = (2.-gamma)/(gamma*M);
    g2 = M;
    
    % ...
    g4 = gmin; % Set a cut-off shock height at which we say the plume is completely trapped
    t4 = 1 + 2*M/((gamma/(1-gamma))*(M*Nf+(M/(M-1))*Ns-(1/(M-1))*Ns*g4^2));
    x4 = -M/g4^2 + (1/(1-gamma))*(-(1/(M-1))*Ns + (M*Nf+(M/(M-1))*Ns)/g4^2)*(t4-1);
    % t4 = 1 + 2*M*(1-gamma)/(gamma*(M*Nf+Ns));
    % x4 = M*(2/gamma-1);
    
    xCs = [x1,x2,x4];
    tCs = [t1,t2,t4];
    gCs = [g1,g2,g4];
    
end