function [Lc,Tc,M,gamma,Nf,Ns,Ng] = hyp_params(Sgr,Swc,krg,H,phi,k,theta,Un,rho_w,rho_g,mu_w,mu_g,Vi,which_Tc)
    % Calculate and return characteristic scales and dimensionless parameters for dimensional input parameters.
    %
    % INPUTS (Careful about the units!)
    %   ** Multiphase flow **
    %   Sgr: Residual gas saturation [-]
    %   Swc: Connate (residual) water saturation [-]
    %   krg: End-point relative permeability to CO2 (for Sw=Swc) [-] 
    %
    %   ** Hydrogeology **
    %   H: Aquifer thickness [m]
    %   phi: Aquifer porosity [-]
    %   k: Aquifer permeability [m^2]
    %   theta: Aquifer slope [degrees]
    %   Un: Characteristic velocity of regional groundwater flow [m/s]
    %   
    %   ** Fluids **
    %   rho_w: Groundwater (brine) density [kg/m^3]
    %   rho_g: CO2 density [kg/m^3]
    %   mu_w: Groundwater (brine) viscosity [Pa-s]
    %   mu_g: CO2 viscosity [Pa-s]
    %
    %   ** CO2 Injection **
    %   Vi = Qi*Ti = volume of CO2 injected *per unit depth into the page* [m^3/m]
    %   (See the paper, Figure 1 and first paragraph of S2)
    %
    %   ** Time scale **
    %   which_Tc = 'f' for flow-based time scale or 's' for slope-based time scale (see below)
    %
    % OUTPUTS
    %   Lc: Characteristic length scale [m]
    %   Tc: Characteristic time scale [s]
    %   M: Mobility ratio (below Eq. 2.9) [-]
    %   gamma: Capillary trapping number (below Eq. 2.2) [-]
    %   Nf: Flow number (Eqs. 2.10) [-]
    %   Ns: Slope number (Eqs. 2.10) [-]
    %   Ng: Gravity number (Eqs. 2.10) [-]
    %
    % To make a dimensionless length dimensional, multiply by Lc: l_dim = Lc*l_nondim
    % To make a dimensionless time dimensional (post-injection),
    %   !! subtract one !! and then multiply by Tc: t_dim = Tc(t_nondim-1)
    
    theta = (pi/180)*2.5;
    
    M = (krg/mu_g)/(1/mu_w);
    gamma = Sgr/(1-Swc);
    
    g = 9.81; % Gravity [m/s^2]
    kappa = (rho_w-rho_g)*g*k*(krg/mu_g)/((1-Swc)*phi); % Characteristic buoyancy velocity [m/s]
    
    % Characteristic length scale
    Lc = (Vi/2)/((1-Swc)*phi*H); % [m]
    
    % Characteristic time scale [s]
    if strcmp(which_Tc,'f')
        Tc = (Vi/2)/(Un*H); % Based on groundwater flow (Nf==1)
    elseif strcmp(which_Tc,'s')
        Tc = Lc/(kappa*sin(theta)); % Based on slope (Ns==1)
    else
        error(['Input parameter which_Tc must be either ''f'' for flow-based time scale or ''s'' for slope-based time scale. You chose ' which_Tc ', which is not a valid option.'])
    end
    
    Nf = (Tc*Un*H)/(Vi/2);
    Ns = (Tc/Lc)*kappa*sin(theta);
    Ng = (Tc/Lc)*kappa*cos(theta)*(H/Lc);
    
end