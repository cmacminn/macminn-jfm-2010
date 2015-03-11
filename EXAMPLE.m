clear all; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Paramters for Region a of the Mt. Simon Sandstone from Szulczewski et al., PNAS 2012 (Table S2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ** Multiphase flow **
Sgr = 0.3; % Residual gas saturation [-]
Swc = 0.4; % Connate (residual) water saturation [-]
krg = 0.6; % End-point relative permeability to CO2 (for Sw=Swc) [-] 

% ** Hydrogeology **
H = 400; % Aquifer thickness [m]
phi = 0.2; % Aquifer porosity [-]
k = (1E-15)*100; % Aquifer permeability [100 mD --> m^2]
theta = 0.5; % Aquifer slope [degrees]
Un = (1E-2/(365.25*24*60*60))*1; % Characteristic velocity of regional groundwater flow [1 cm/yr --> m/s]

% ** Fluids **
rho_w = 1000; % Groundwater (brine) density [kg/m^3]
rho_g = 700; % CO2 density [kg/m^3]
mu_w = 0.8E-3; % Groundwater (brine) viscosity [Pa-s]
mu_g = 0.06E-3; % CO2 viscosity [Pa-s]

%   ** CO2 Injection **
Mi = 10*1E9*1E3; % Total mass of CO2 injected [10 GT --> kg]
W = 200E3; % Length of the well array [m]
Vi = (Mi/rho_g)/W; % volume of CO2 injected *per unit depth into the page* [kg --> m^3/m]

% which_Tc = 'f'; % Flow-based time scale
which_Tc = 's'; % Slope-based time scale

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get the characteristic scales and the dimensionless parameters
[Lc,Tc,M,gamma,Nf,Ns,Ng] = hyp_params(Sgr,Swc,krg,H,phi,k,theta,Un,rho_w,rho_g,mu_w,mu_g,Vi,which_Tc);

% Plot the characteristics diagram and get the critical positions and times
[xCs,tCs,eff] = hyp_chars(M,gamma,Nf,Ns);

% Get the plume shape at a selection of times over the whole lifetime, and make a movie
N = 10000;
hs = [0.,logspace(-6.,0.,N)];
ts = linspace(tCs(1),tCs(end),100);
[xss,hss,xmss,hmss] = hyp_plume(M,gamma,Nf,Ns,ts,hs);

xLmost = min([-xCs(1),xCs]);
xRmost = max(xCs);
for ti = 1:length(ts)
    
    figure(2);
    hold off
    xlabel('\xi')
    ylabel('1-\eta')
    fill(1.1*[xLmost,xLmost,xRmost,xRmost],[1,0,0,1],[1,1,1],'linewidth',1.5);
    axis([1.1*xLmost,1.1*xRmost,0,1])
    hold on
    
    xs = xss(ti,:);
    hs = hss(ti,:);
    xms = xmss(ti,:);
    hms = hmss(ti,:);
    
    fill([xms(1),xms,xms(end),xms(end),xms(1)],[1-hms(1),1-hms,1-hms(end),1,1],[200,200,200]/255,'linewidth',1.5)
    plot(xms,1-hms,'k-','linewidth',1.5)
    fill([xs(1),xs,xs(end),xs(end),xs(1)],[1-hs(1),1-hs,1-hs(end),1,1],[121,121,121]/255,'linewidth',1.5)
    plot([0,0],[0,1],'k-')
    box on;
    drawnow()
    
end

% Display the key numbers, converted to meaningful units
disp(['The plume is completely trapped after ' num2str(Tc*(tCs(end)-1)/(365*24*60*60)) ' years, having migrated ' num2str(Lc*xCs(end)/1000) ' km.'])