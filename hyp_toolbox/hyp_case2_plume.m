function [xs,gs,xms,gms] = hyp_case2_plume(xCs,tCs,gCs,M,gamma,Nf,Ns,t,gs,gmin)
    % Calculate plume shape for case 2
    
    % Be quiet?  (Suppress nonessential output?)
    switches_quiet = 1;
    
    % ---------------------------------------------------------
    % Key Values
    % ---------------------------------------------------------
    x1 = xCs(1); x2a = xCs(2); x2b = xCs(3); x4 = xCs(4);
    t1 = tCs(1); t2a = tCs(2); t2b = tCs(3); t4 = tCs(4);
    g1 = gCs(1); g2a = gCs(2); g2b = gCs(3); g4 = gCs(4);
    
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
        
    elseif t1<t && t<=t2a
        
        % Retreat, part a
        if switches_quiet==0
            disp(['Retreat, part a:  t1 = ' num2str(t1) ' < t = ' num2str(t) ' < t2a = ' num2str(t2a)])
        end
        xLs = -(M./(gs.^2)) + (1/(1-gamma))*(-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gs.^2))*(t-1);
        xRs =  (M./(gs.^2)) +               (-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gs.^2))*(t-1);
        xs = [xLs,fliplr(xRs)];
        gs = [gLs,fliplr(gRs)];
        
        xLms = -(M./(gLms.^2));
        xRms = xRs;
        xms = [xLms,fliplr(xRms)];
        gms = [gLms,fliplr(gRms)];
        
    elseif t2a<t && t<=t2b
        
        % Retreat, part b
        if switches_quiet==0
            disp(['Retreat, part b:  t2a = ' num2str(t2a) ' < t = ' num2str(t) ' < t2b = ' num2str(t2b)])
        end
        gp_t = sqrt(((M-1)/Ns)*((M*Nf+(M/(M-1))*Ns)-2*M*(1-gamma)/(gamma*(t-1))));
        
        gLs(gLs>gp_t) = gp_t;
        gRs(gRs>gp_t) = gp_t;
        xLs = -(M./(gLs.^2)) + (1/(1-gamma))*(-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gLs.^2))*(t-1);
        xRs =  (M./(gRs.^2)) +               (-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gRs.^2))*(t-1);
        
        xs = [xLs,fliplr(xRs)];
        gs = [gLs,fliplr(gRs)];
        
        gp_ts = gLms(gLms>gp_t);
        tps = 1 + 2*M./((gamma/(1-gamma))*(M*Nf+(M/(M-1))*Ns-(1/(M-1))*Ns*gp_ts.^2));
        xLms = -(M./(gLms.^2));
        xRms = xRs;
        xRms(gRms>gp_t) =  (M./(gp_ts.^2)) +               (-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gp_ts.^2)).*(tps-1);
        
        xms = [xLms,fliplr(xRms)];
        gms = [gLms,fliplr(gRms)];
        
    elseif t2b<t && t<=t4
        
        % Sweep
        if switches_quiet==0
            disp(['Sweep:  t2b = ' num2str(t2b) ' < t = ' num2str(t) ' < t4 = ' num2str(t4)])
        end
        c1 = (M*Nf+(M/(M-1))*Ns);
        c2 = ((1/(M-1))*Ns);
        c3 = M;
        c4 = 1/(1-gamma);
        gm0 = g2b;
        tm0 = t2b;
        
        tm = t;
        g_lim = g4;
        t_lim = t4;
        gm_t = hyp_shockHeight(gamma,c1,c2,c3,c4,gm0,tm0,tm,g_lim,t_lim);
        xm_t = (M./(gm_t.^2)) + (-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gm_t.^2))*(t-1);
        
        gLs(:) = gm_t;
        gRs(gRs>gm_t) = gm_t;
        xLs = xm_t*ones(1,length(gs));
        xRs = (M./(gRs.^2)) + (-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gRs.^2))*(t-1);
        
        xs = [xLs,fliplr(xRs)];
        gs = [gLs,fliplr(gRs)];
        
        % End-of-injection location of the left front
        xLms = -(M./(gLms.^2));
        % Current location of the right front ...
        xRms = (M./(gRms.^2)) + (-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gRms.^2)).*(t-1);
        % ... except for waves that were eaten by the peak ...
        gp_ts = gRms(gRms>=g2b);
        tps = 1 + 2*M./((gamma/(1-gamma))*(M*Nf+(M/(M-1))*Ns-(1/(M-1))*Ns*gp_ts.^2));
        xp_ts = (M./(gp_ts.^2)) + (-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gp_ts.^2)).*(tps-1);
        xRms(gRms>=g2b) = xp_ts;
        % ... and except for waves that were eaten by the shock
        gmF = gm_t;
        gm_ts = gRms(logical((gRms>gm_t).*(gRms<g2b)));
        tms = hyp_shock(gamma,c1,c2,c3,c4,gm0,tm0,gm_ts);
        xm_ts = (M./(gm_ts.^2)) + (-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gm_ts.^2)).*(tms-1);
        xRms(logical((gRms>gm_t).*(gRms<g2b))) = xm_ts;
        
        xms = [xLms,fliplr(xRms)];
        gms = [gLms,fliplr(gRms)];
        
    elseif t4<t
        
        xs = zeros(1,2*length(gs));
        gs = [gLs,fliplr(gRs)];
        
        % End-of-injection location of the left front
        xLms = -(M./(gLms.^2));
        % Current location of the right front ...
        xRms = (M./(gRms.^2)) + (-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gRms.^2)).*(t-1);
        % ... except for waves that were eaten by the peak ...
        gp_ts = gRms(gRms>=g2b);
        tps = 1 + 2*M./((gamma/(1-gamma))*(M*Nf+(M/(M-1))*Ns-(1/(M-1))*Ns*gp_ts.^2));
        xp_ts = (M./(gp_ts.^2)) + (-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gp_ts.^2)).*(tps-1);
        xRms(gRms>=g2b) = xp_ts;
        % ... and except for waves that were eaten by the shock
        gm_t = 1;
        c1 = (M*Nf+(M/(M-1))*Ns);
        c2 = ((1/(M-1))*Ns);
        c3 = M;
        c4 = 1/(1-gamma);
        tm0 = t2b;
        gm0 = g2b;
        gmF = gm_t;
        gm_ts = gRms(logical((gRms>gm_t).*(gRms<g2b)));
        tms = hyp_shock(gamma,c1,c2,c3,c4,gm0,tm0,gm_ts);
        xm_ts = (M./(gm_ts.^2)) + (-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gm_ts.^2)).*(tms-1);
        xRms(logical((gRms>gm_t).*(gRms<g2b))) = xm_ts;
        
        xms = [xLms,fliplr(xRms)];
        gms = [gLms,fliplr(gRms)];
        
    end
    
end