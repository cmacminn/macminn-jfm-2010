function [xs,gs,xms,gms] = hyp_case5_plume(xCs,tCs,gCs,M,gamma,Nf,Ns,t,gs,gmin)
    % Calculate plume shape for case 5
    %   Note that this case only exists if M<sqrt(2/gamma-1)
    
    % Be quiet?  (Suppress nonessential output?)
    switches_quiet = 1;
    
    % ---------------------------------------------------------
    % Key Values
    % ---------------------------------------------------------
    x1 = xCs(1); x2 = xCs(2); x3 = xCs(3); x4 = xCs(4);
    t1 = tCs(1); t2 = tCs(2); t3 = tCs(3); t4 = tCs(4);
    g1 = gCs(1); g2 = gCs(2); g3 = gCs(3); g4 = gCs(4);
    
    % ---------------------------------------------------------
    % Calculate the plume shape at time t
    % ---------------------------------------------------------
    
    gLs = gs;
    gRs = gs;
    
    gLms = gs;
    gRms = gs;
    
    if t<0
        
        % ?
        error('You asked for the plume shape for some t<0... ?')
        xs = zeros(1,2*length(gs));
        gs = [gLs,fliplr(gRs)];
        xms = xs;
        gms = gs;
        
    elseif t==0
        
        xs = zeros(1,2*length(gs));
        gs = [gLs,fliplr(gRs)];
        xms = xs;
        gms = gs;
        
    elseif 0<t && t<=t1
        
        % Injection
        if switches_quiet==0
            disp(['Injection:  ' num2str(t) ' < t1 = ' num2str(t1)])
        end
        xLs = -(M./(gs.^2))*t;
        xRs =  (M./(gs.^2))*t;
        xs = [xLs,fliplr(xRs)];
        gs = [gLs,fliplr(gRs)];
        
        xms = xs;
        gms = gs;
        
    elseif t1<t && t<=t2
        
        % Retreat
        if switches_quiet==0
            disp(['Retreat:  t1 = ' num2str(t1) ' < t = ' num2str(t) ' < t2 = ' num2str(t2)])
        end
        xLs = -(M./(gs.^2)) + (1/(1-gamma))*(-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gs.^2))*(t-1);
        xRs =  (M./(gs.^2)) +               (-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gs.^2))*(t-1);
        xs = [xLs,fliplr(xRs)];
        gs = [gLs,fliplr(gRs)];
        
        % End-of-injection location of the left front
        xLms = -(M./(gLms.^2));
        % Current location of the right front
        xRms = xRs;
        
        xms = [xLms,fliplr(xRms)];
        gms = [gLms,fliplr(gRms)];
        
    elseif t2<t && t<=t3
        
        % Chase
        if switches_quiet==0
            disp(['Chase:  t2 = ' num2str(t2) ' < t = ' num2str(t) ' < t3 = ' num2str(t3)])
        end
        % gLs(:) = g3;
        xLs = -(M./(gLs.^2)) + (1/(1-gamma))*(-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gLs.^2))*(t-1);
        xRs = (x2 + (Nf/(1-gamma))*(t-t2))*ones(1,length(gs));
        
        xs = [xLs,fliplr(xRs)];
        gs = [gLs,fliplr(gRs)];
        
        % End-of-injection location of the left front
        xLms = -(M./(gLms.^2));
        % Current location of the right front
        xRms = xRs;
        
        xms = [xLms,fliplr(xRms)];
        gms = [gLms,fliplr(gRms)];
        
    elseif t3<t && t<t4
        
        % Sweep
        if switches_quiet==0
            disp(['Sweep:  t3 = ' num2str(t3) ' < t = ' num2str(t) ' < t4 = ' num2str(t4)])
        end
        c1 = (1/(1-gamma))*(M*Nf+(M/(M-1))*Ns);
        c2 = (1/(1-gamma))*((1/(M-1))*Ns);
        c3 = -M;
        c4 = (1-gamma);
        gm0 = g3;
        tm0 = t3;
        
        tm = t;
        g_lim = g4;
        t_lim = t4;
        gm_t = hyp_shockHeight(gamma,c1,c2,c3,c4,gm0,tm0,tm,g_lim,t_lim);
        xm_t = -(M./(gm_t.^2)) + (1/(1-gamma))*(-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gm_t.^2))*(t-1);
        
        gLs(gLs>gm_t) = gm_t;
        % gRs(:) = gm_t;
        gRs(gRs>gm_t) = gm_t;
        xLs = -(M./(gLs.^2)) + (1/(1-gamma))*(-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gLs.^2))*(t-1);
        xRs = xm_t*ones(1,length(gs));
        
        xs = [xLs,fliplr(xRs)];
        gs = [gLs,fliplr(gRs)];
        
        % End-of-injection location of the left front
        xLms = -(M./(gLms.^2));
        % Current location of the right front ...
        xRms = xm_t*ones(1,length(gRms));
        % ... except for waves that have not yet been eaten by the shock
        gm0 = g3;
        tm0 = t3;
        gmF = gm_t;
        gm_ts = gLms(gLms>gm_t);
        tms = hyp_shock(gamma,c1,c2,c3,c4,gm0,tm0,gm_ts);
        xm_ts = -(M./(gm_ts.^2)) + (1/(1-gamma))*(-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gm_ts.^2)).*(tms-1);
        % xRms = -(M./(gRms.^2)) + (1/(1-gamma))*(-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gLms.^2)).*(t-1);
        xRms(gLms>gm_t) = xm_ts;
        
        xms = [xLms,fliplr(xRms)];
        gms = [gLms,fliplr(gRms)];
        
    elseif t4<=t
        
        xs = zeros(1,2*length(gs));
        gs = [gLs,fliplr(gRs)];
        
        % End-of-injection location of the left front
        xLms = -(M./(gLms.^2));
        % Current location of the right front ...
        xRms = xm_t*ones(1,length(gRms));
        % ... except for waves that have not yet been eaten by the shock
        gm_t = 1;
        xm_t = -(M./(gm_t.^2)) + (1/(1-gamma))*(-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gm_t.^2))*(t4-1);
        c1 = (1/(1-gamma))*(M*Nf+(M/(M-1))*Ns);
        c2 = (1/(1-gamma))*((1/(M-1))*Ns);
        c3 = -M;
        c4 = (1-gamma);
        gm0 = g3;
        tm0 = t3;
        gmF = gm_t;
        gm_ts = gLms(gLms>gm_t);
        tms = hyp_shock(gamma,c1,c2,c3,c4,gm0,tm0,gm_ts);
        xm_ts = -(M./(gm_ts.^2)) + (1/(1-gamma))*(-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gm_ts.^2)).*(tms-1);
        % xRms = -(M./(gRms.^2)) + (1/(1-gamma))*(-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gLms.^2)).*(t-1);
        xRms(gLms>gm_t) = xm_ts;
    end
    
end