function [tms] = hyp_shock(gamma,c1,c2,c3,c4,g0,t0,gms)
    % Shock height as a function of time, general case (based on an implicit analytical expression)
    % Specifically: For a shock that has initial height g0 at time t0, 
    % calculate and return the times tms at which the shock has heights gms.
    % 
    % After much algebra, the constants c1,c2,c3,c4 are defined as follows:
    %
    % c1 = M*Nf + (M/(M-1))*Ns if the shock is colliding with a drainage front
    %   OR
    %    = (above)/(1-gamma) if colliding with an imbibition front
    %
    % c2 = (1/(M-1))*Ns if the shock is colliding with a drainage front
    %   OR
    %    = (above)/(1-gamma) if colliding with an imbibition front
    %
    % c3 = -M if the shock is colliding with a LEFT front
    %   OR
    %    = M if colliding with a RIGHT front
    %
    % c4 = 1 if the shock is a drainage front colliding with a drainage front,
    %        or if imibition front colliding with imbibition front
    %   OR
    %    = (1-gamma) if drainage front colliding with imbibition front
    %   OR
    %    = 1/(1-gamma) if imbibition front colliding with drainage front
    %
    % CWM, 2009
    
    c5 = c2/c1;
    
    if c4==1
        I_gms = 2*log(gms./(gms-1));
        I_g0  = 2*log(g0 ./(g0 -1));
    elseif c5==0
        I_gms = log((gms./(gms-(1-gamma))).^2);
        I_g0 = log((g0./(g0-(1-gamma))).^2);
    else
        A = sqrt(-4*c5*(1-c4)-c4^2);
        I_gms = log((gms.^2)./(1-c4*gms-c5*(1-c4)*gms.^2)) + (2*c4/A)*atan((-2*c5*(1-c4)*gms-c4)/A);
        I_g0  = log((g0^2)/(1-c4*g0-c5*(1-c4)*g0^2)) + (2*c4/A)*atan((-2*c5*(1-c4)*g0-c4)/A);
    end
    
    tms = (1-c3/c1) + (t0-(1-c3/c1))*exp(I_gms-I_g0);
    
end