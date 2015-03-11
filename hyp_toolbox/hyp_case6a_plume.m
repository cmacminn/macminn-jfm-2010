function [xs,gs,xms,gms] = hyp_case6a_plume(xCs,tCs,gCs,M,gamma,Nf,Ns,t,gs,gmin)
    % Calculate plume shape for case 6a
    %   Note that Ns/Nf-->-infty (negative slope only) is contained here
    
    % 6a:  a shock forms on the right and hits the bottom before a peak forms (no peak)
    % RR_RS, R_R0, R_RL, R_LR, R_L0, R_LS, R_LL
    
    % Be quiet?  (Suppress nonessential output?)
    switches_quiet = 0;
    
    % ---------------------------------------------------------
    % Key Values
    % ---------------------------------------------------------
    x1 = xCs(1); x2a = xCs(2); x2b = xCs(3); x2c = xCs(4); x3 = xCs(5); x4a = xCs(6); x4b = xCs(7); x4c = xCs(8);
    t1 = tCs(1); t2a = tCs(2); t2b = tCs(3); t2c = tCs(4); t3 = tCs(5); t4a = tCs(6); t4b = tCs(7); t4c = tCs(8);
    g1 = gCs(1); g2a = gCs(2); g2b = gCs(3); g2c = gCs(4); g3 = gCs(5); g4a = gCs(6); g4b = gCs(7); g4c = gCs(8);
    
    g_s = sqrt(M*(M-1)*Nf/Ns+M);
    g_0 = M*(M-1)*Nf/Ns + M; % zero of the flux function
    xLS = -1/((M-1)*Nf/Ns+1);
    xRS =  1/((M-1)*Nf/Ns+1);
    
    % Collision times
    tRR_RS = 1 - (1-gamma)*(M-1)/((M-1)*Nf+Ns);
    xRR_RS = xRS;
    tLR_RL = 1 + 2*((1-gamma)/gamma)/(Nf-Ns);
    xLR_RL = (2-gamma)/(M*gamma);
    
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
        
    elseif t==0
        
        xs = zeros(1,2*length(gs));
        gs = [gLs, fliplr(gRs)];
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
        gs = [gLs, fliplr(gRs)];
        xms = xs;
        gms = gs;
        
    elseif t1<t && t<=t2a
        
        % Retreat, part a
        if switches_quiet==0
            disp(['Retreat, part a:  t1 = ' num2str(t1) ' < t = ' num2str(t) ' < t2a = ' num2str(t2a)])
        end
        xLs = -(M./(gs.^2)) + (1*(gs<=g_s) + (1/(1-gamma))*(gs> g_s)).*(-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gs.^2))*(t-1);
        xRs =  (M./(gs.^2)) + ((1/(1-gamma))*(gs<=g_s) + 1*(gs> g_s)).*(-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gs.^2))*(t-1);
        xs = [xLs,fliplr(xRs)];
        gs = [gLs,fliplr(gRs)];
        
        % End-of-injection location of the left front ...
        xLms = -(M./(gLms.^2));
        % ... except left-moving waves
        xLms(gLms<g_s) = xLs(gLms<g_s);
        
        % Current location of the right front ...
        xRms = xRs;
        % ... except left-moving waves
        xRms(gRms<g_s) = (M./(gRms(gRms<g_s).^2));
        
        xms = [xLms,fliplr(xRms)];
        gms = [gLms,fliplr(gRms)];
        
    elseif t2a<t && t<=t2b
        
        % Retreat, part b
        if switches_quiet==0
            disp(['Retreat, part b:  t2a = ' num2str(t2a) ' < t = ' num2str(t) ' < t2b = ' num2str(t2b)])
        end
        
        c1 = (M*Nf+(M/(M-1))*Ns);
        c2 = ((1/(M-1))*Ns);
        c3 = M;
        c4 = 1/(1-gamma);
        gm0 = g2a;
        tm0 = t2a;
        
        tm = t;
        g_lim = g2b;
        t_lim = t2b;
        gm = hyp_shockHeight(gamma,c1,c2,c3,c4,gm0,tm0,tm,g_lim,t_lim);
        
        gRs(gRs<gm) = gm;
        xLs = -(M./(gLs.^2)) + (1*(gLs <=g_s) + (1/(1-gamma))*(gLs> g_s)).*(-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gLs.^2))*(t-1);
        xRs =  (M./(gRs.^2)) + ((1/(1-gamma))*(gRs<=g_s) + 1*(gRs> g_s)).*(-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gRs.^2))*(t-1);
        
        xs = [xLs,fliplr(xRs)];
        gs = [gLs,fliplr(gRs)];
        
        % End-of-injection location of the left front ...
        xLms = -(M./(gLms.^2));
        % ... except left-moving waves ...
        xLms(gLms<g_s) = xLs(gLms<g_s);
        
        % Current location of the right front ...
        xRms = xRs;
        % ... except left-moving waves ...
        xRms(gRms<g_s) = (M./(gRms(gRms<g_s).^2));
        % ... and right-moving waves that have been eaten by the shock
        gm_ts = gRms(logical((g_s<gRms).*(gRms<=gm)));
        tms = hyp_shock(gamma,c1,c2,c3,c4,gm0,tm0,gm_ts);
        xm_ts = (M./(gm_ts.^2)) + (-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gm_ts.^2)).*(tms-1);
        xRms(logical((g_s<gRms).*(gRms<=gm))) = xm_ts;
        
        xms = [xLms,fliplr(xRms)];
        gms = [gLms,fliplr(gRms)];
        
    elseif t2b<t && t<=t2c
        
        % Retreat, part c
        if switches_quiet==0
            disp(['Retreat, part c:  t2b = ' num2str(t2b) ' < t = ' num2str(t) ' < t2c = ' num2str(t2c)])
        end
        
        % Shock travels to the right
        c1 = (M*Nf+(M/(M-1))*Ns);
        c2 = ((1/(M-1))*Ns);
        c3 = M;
        c4 = 1;
        gm0 = g2b;
        tm0 = t2b;
        
        tm = t;
        g_lim = g2c;
        t_lim = t2c;
        gm = hyp_shockHeight(gamma,c1,c2,c3,c4,gm0,tm0,tm,g_lim,t_lim);
        
        gRs(gRs<gm) = gm;
        xLs = -(M./(gLs.^2)) + (1*(gLs <=g_s) + (1/(1-gamma))*(gLs> g_s)).*(-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gLs.^2))*(t-1);
        xRs =  (M./(gRs.^2)) + ((1/(1-gamma))*(gRs<=g_s) + 1*(gRs> g_s)).*(-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gRs.^2))*(t-1);
        
        xs = [xLs,fliplr(xRs)];
        gs = [gLs,fliplr(gRs)];
        
        % End-of-injection location of the left front ...
        xLms = -(M./(gLms.^2));
        % ... except left-moving waves
        xLms(gLms<g_s) = -(M./(gLms(gLms<g_s).^2)) + (-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gLms(gLms<g_s).^2))*(t-1);
        
        % Current location of the right front ...
        xRms = xRs;
        % ... except left-moving waves ...
        xRms(gRms<g_s) = (M./(gRms(gRms<g_s).^2));
        % ... except right-moving waves that were eaten by the shock before it changed direction...
        c1 = (M*Nf+(M/(M-1))*Ns);
        c2 = ((1/(M-1))*Ns);
        c3 = M;
        c4 = (1/(1-gamma));
        gm0 = g2a;
        tm0 = t2a;
        gm_ts = gRms(logical((g_s<gRms).*(gRms<=g_0)));
        tms = hyp_shock(gamma,c1,c2,c3,c4,gm0,tm0,gm_ts);
        xm_ts = (M./(gm_ts.^2)) + (-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gm_ts.^2)).*(tms-1);
        xRms(logical((g_s<gRms).*(gRms<=g_0))) = xm_ts;
        
        xms = [xLms,fliplr(xRms)];
        gms = [gLms,fliplr(gRms)];
        
    elseif t2c<t && t<=t3
        
        % Chase
        if switches_quiet==0
            disp(['Chase:  t2c = ' num2str(t2c) ' < t = ' num2str(t) ' < t3 = ' num2str(t3)])
        end
        
        xLs = -(M./(gLs.^2)) + (1*(gLs <=g_s) + (1/(1-gamma))*(gLs> g_s)).*(-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gLs.^2))*(t-1);
        gRs(:) = M;
        xRs = (x2b + Nf*(t-t2b))*ones(size(gs));
        
        xs = [xLs,fliplr(xRs)];
        gs = [gLs,fliplr(gRs)];
        
        % End-of-injection location of the left front ...
        xLms = -(M./(gLms.^2));
        % ... except left-moving waves
        xLms(gLms<g_s) = xLs(gLms<g_s);
        
        % Current location of the right front ...
        xRms = xRs;
        % ... except left-moving waves ...
        xRms(gRms<g_s) = (M./(gRms(gRms<g_s).^2));
        % ... except right-moving waves that were eaten by the shock before it changed direction ...
        c1 = (M*Nf+(M/(M-1))*Ns);
        c2 = ((1/(M-1))*Ns);
        c3 = M;
        c4 = (1/(1-gamma));
        gm0 = g2a;
        tm0 = t2a;
        gm_ts = gRms(logical((g_s<gRms).*(gRms<=g_0)));
        tms = hyp_shock(gamma,c1,c2,c3,c4,gm0,tm0,gm_ts);
        xm_ts = (M./(gm_ts.^2)) + (-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gm_ts.^2)).*(tms-1);
        xRms(logical((g_s<gRms).*(gRms<=g_0))) = xm_ts;
        % All taller waves than g_s are now part of the chase shock
        xRms = max([xRms;max(xRs)*ones(size(gRms))]);
        
        xms = [xLms,fliplr(xRms)];
        gms = [gLms,fliplr(gRms)];
        
    elseif t3<t && t<=t4a
        
        % Sweep, part a
        if switches_quiet==0
            disp(['Sweep, part a:  t3 = ' num2str(t3) ' < t = ' num2str(t) ' < t4a = ' num2str(t4a)])
        end
        
        %   Evolve shock through left thicknesses from g3 to g_0 (shock travels to the right)
        % This portion only happens if g3<g_0 (i.e., if t3<t4a)
        c1 = (1/(1-gamma))*(M*Nf+(M/(M-1))*Ns);
        c2 = (1/(1-gamma))*((1/(M-1))*Ns);
        c3 = -M;
        c4 = (1-gamma);
        gm0 = g3;
        tm0 = t3;
        
        tm = t;
        g_lim = g4a;
        t_lim = t4a;
        gm = hyp_shockHeight(gamma,c1,c2,c3,c4,gm0,tm0,tm,g_lim,t_lim);
        
        gLs(gLs>gm) = gm;
        xLs = -(M./(gLs.^2)) + (1*(gLs<=g_s) + (1/(1-gamma))*(gLs> g_s)).*(-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gLs.^2))*(t-1);
        gRs(:) = gm;
        xRs = -(M./(gRs.^2)) + (1/(1-gamma))*(-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gRs.^2))*(t-1);
        
        xs = [xLs,fliplr(xRs)];
        gs = [gLs,fliplr(gRs)];
        
        % End-of-injection location of the left front ...
        xLms = -(M./(gLms.^2));
        % ... except left-moving waves
        xLms(gLms<g_s) = xLs(gLms<g_s);
        
        % Current location of the right front ...
        xRms = xRs;
        % ... except left-moving waves ...
        xRms(gRms<g_s) = (M./(gRms(gRms<g_s).^2));
        % ... except right-moving waves that were eaten by the shock before it changed direction ...
        c1 = (M*Nf+(M/(M-1))*Ns);
        c2 = ((1/(M-1))*Ns);
        c3 = M;
        c4 = (1/(1-gamma));
        gm0 = g2a;
        tm0 = t2a;
        gm_ts = gRms(logical((g_s<gRms).*(gRms<=g_0)));
        tms = hyp_shock(gamma,c1,c2,c3,c4,gm0,tm0,gm_ts);
        xm_ts = (M./(gm_ts.^2)) + (-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gm_ts.^2)).*(tms-1);
        xRms(logical((g_s<gRms).*(gRms<=g_0))) = xm_ts;
        % ... or after ...
        c1 = (M*Nf+(M/(M-1))*Ns);
        c2 = ((1/(M-1))*Ns);
        c3 = M;
        c4 = 1;
        gm0 = g2b;
        tm0 = t2b;
        gm_ts = gRms(logical((g_0<gRms).*(gRms<=M)));
        tms = hyp_shock(gamma,c1,c2,c3,c4,gm0,tm0,gm_ts);
        xm_ts = (M./(gm_ts.^2)) + (-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gm_ts.^2)).*(tms-1);
        xRms(logical((g_0<gRms).*(gRms<=M))) = xm_ts;
        % ... and right-moving waves of the left front that travel past the right front
        c1 = (1/(1-gamma))*(M*Nf+(M/(M-1))*Ns);
        c2 = (1/(1-gamma))*((1/(M-1))*Ns);
        c3 = -M;
        c4 = (1-gamma);
        gm0 = g3;
        tm0 = t3;
        gm_ts = gLms(gm<=gLms);
        tms = hyp_shock(gamma,c1,c2,c3,c4,gm0,tm0,gm_ts);
        xm_ts = -(M./(gm_ts.^2)) + (1/(1-gamma))*(-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gm_ts.^2)).*(tms-1);
        xRms(gm<=gRms) = max([xRms(gm<=gRms);xm_ts]);
        xRms(gRms<gm) = max([max(xm_ts)*ones(size(gRms(gRms<gm)));xRms(gRms<gm)]);
        
        xms = [xLms,fliplr(xRms)];
        gms = [gLms,fliplr(gRms)];
        
    elseif t4a<t && t<=t4b
        
        % Sweep, part b
        if switches_quiet==0
            disp(['Sweep, part b:  t4a = ' num2str(t4a) ' < t = ' num2str(t) ' < t4b = ' num2str(t4b)])
        end
        
        %   Evolve shock through left thicknesses from g4a to g_s (shock travels to the left)
        c1 = (1/(1-gamma))*(M*Nf+(M/(M-1))*Ns);
        c2 = (1/(1-gamma))*((1/(M-1))*Ns);
        c3 = -M;
        c4 = 1;
        gm0 = g4a;
        tm0 = t4a;
        
        tm = t;
        g_lim = g4b;
        t_lim = t4b;
        gm = hyp_shockHeight(gamma,c1,c2,c3,c4,gm0,tm0,tm,g_lim,t_lim);
        
        gLs(gLs>gm) = gm;
        gRs(:) = gm;
        xLs = -(M./(gLs.^2)) + (1*(gLs<=g_s) + (1/(1-gamma))*(gLs> g_s)).*(-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gLs.^2))*(t-1);
        xRs = -(M./(gRs.^2)) + (1/(1-gamma))*(-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gRs.^2))*(t-1);
        
        xs = [xLs,fliplr(xRs)];
        gs = [gLs,fliplr(gRs)];
        
        % End-of-injection location of the left front ...
        xLms = -(M./(gLms.^2));
        % ... except left-moving waves
        xLms(gLms<g_s) = xLs(gLms<g_s);
        
        % Current location of the right front ...
        xRms = xRs;
        % ... except left-moving waves ...
        xRms(gRms<g_s) = (M./(gRms(gRms<g_s).^2));
        % ... except right-moving waves that were eaten by the shock before it changed direction ...
        c1 = (M*Nf+(M/(M-1))*Ns);
        c2 = ((1/(M-1))*Ns);
        c3 = M;
        c4 = (1/(1-gamma));
        gm0 = g2a;
        tm0 = t2a;
        gm_ts = gRms(logical((g_s<gRms).*(gRms<=g_0)));
        tms = hyp_shock(gamma,c1,c2,c3,c4,gm0,tm0,gm_ts);
        xm_ts = (M./(gm_ts.^2)) + (-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gm_ts.^2)).*(tms-1);
        xRms(logical((g_s<gRms).*(gRms<=g_0))) = xm_ts;
        % ... or after ...
        c1 = (M*Nf+(M/(M-1))*Ns);
        c2 = ((1/(M-1))*Ns);
        c3 = M;
        c4 = 1;
        gm0 = g2b;
        tm0 = t2b;
        gm = M;
        gm_ts = gRms(logical((g_0<gRms).*(gRms<=gm)));
        tms = hyp_shock(gamma,c1,c2,c3,c4,gm0,tm0,gm_ts);
        xm_ts = (M./(gm_ts.^2)) + (-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gm_ts.^2)).*(tms-1);
        xRms(logical((g_0<gRms).*(gRms<=gm))) = xm_ts;
        % ... and right-moving waves of the left front that traveled past the right front
        c1 = (1/(1-gamma))*(M*Nf+(M/(M-1))*Ns);
        c2 = (1/(1-gamma))*((1/(M-1))*Ns);
        c3 = -M;
        c4 = (1-gamma);
        gm0 = g3;
        tm0 = t3;
        gm_ts = gLms(g_0<=gLms);
        tms = hyp_shock(gamma,c1,c2,c3,c4,gm0,tm0,gm_ts);
        xm_ts = -(M./(gm_ts.^2)) + (1/(1-gamma))*(-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gm_ts.^2)).*(tms-1);
        xRms(g_0<=gRms) = max([xRms(g_0<=gRms);xm_ts]);
        xRms(gRms<g_0) = max([max(xm_ts)*ones(size(gRms(gRms<g_0)));xRms(gRms<g_0)]);
        
        xms = [xLms,fliplr(xRms)];
        gms = [gLms,fliplr(gRms)];
        
    elseif t4b<t && t<=t4c
        
        % Sweep, part c
        if switches_quiet==0
            disp(['Sweep, part c:  t4b = ' num2str(t4b) ' < t = ' num2str(t) ' < t4c = ' num2str(t4c)])
        end
        %   Evolve shock through left thicknesses from g_s (g4b) to 1 (shock travels to the left)
        c1 = (M*Nf+(M/(M-1))*Ns);
        c2 = ((1/(M-1))*Ns);
        c3 = -M;
        c4 = 1/(1-gamma);
        gm0 = g4b;
        tm0 = t4b;
        
        tm = t;
        g_lim = g4c;
        t_lim = t4c;
        gm = hyp_shockHeight(gamma,c1,c2,c3,c4,gm0,tm0,tm,g_lim,t_lim);
        
        gLs(gLs>gm) = gm;
        gRs(:) = gm;
        xLs = -(M./(gLs.^2)) + (1*(gLs<=g_s) + (1/(1-gamma))*(gLs> g_s)).*(-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gLs.^2))*(t-1);
        xRs = -(M./(gRs.^2)) + (1*(gRs<=g_s) + (1/(1-gamma))*(gRs> g_s)).*(-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gRs.^2))*(t-1);
        
        xs = [xLs,fliplr(xRs)];
        gs = [gLs,fliplr(gRs)];
        
        % End-of-injection location of the left front ...
        xLms = -(M./(gLms.^2));
        % ... except left-moving waves
        xLms(gLms<g_s) = xLs(gLms<g_s);
        % ... except left-moving waves that have been eaten by the shock
        gm_ts = gLms(logical((gm<gLms).*(gLms<g_s)));
        tms = hyp_shock(gamma,c1,c2,c3,c4,gm0,tm0,gm_ts);
        xm_ts = -(M./(gm_ts.^2)) + (-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gm_ts.^2)).*(tms-1);
        xLms(logical((gm<gLms).*(gLms<g_s))) = xm_ts;
        
        % Current location of the right front ...
        xRms = xRs;
        % ... except left-moving waves ...
        xRms(gRms<g_s) = (M./(gRms(gRms<g_s).^2));
        % ... except right-moving waves that were eaten by the shock before it changed direction ...
        c1 = (M*Nf+(M/(M-1))*Ns);
        c2 = ((1/(M-1))*Ns);
        c3 = M;
        c4 = (1/(1-gamma));
        gm0 = g2a;
        tm0 = t2a;
        gm_ts = gRms(logical((g_s<gRms).*(gRms<=g_0)));
        tms = hyp_shock(gamma,c1,c2,c3,c4,gm0,tm0,gm_ts);
        xm_ts = (M./(gm_ts.^2)) + (-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gm_ts.^2)).*(tms-1);
        xRms(logical((g_s<gRms).*(gRms<=g_0))) = xm_ts;
        % ... or after it changed direction ...
        c1 = (M*Nf+(M/(M-1))*Ns);
        c2 = ((1/(M-1))*Ns);
        c3 = M;
        c4 = 1;
        gm0 = g2b;
        tm0 = t2b;
        gm_ts = gRms(logical((g_0<gRms).*(gRms<=M)));
        tms = hyp_shock(gamma,c1,c2,c3,c4,gm0,tm0,gm_ts);
        xm_ts = (M./(gm_ts.^2)) + (-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gm_ts.^2)).*(tms-1);
        xRms(logical((g_0<gRms).*(gRms<=M))) = xm_ts;
        % ... and right-moving waves of the left front that traveled past the right front
        c1 = (1/(1-gamma))*(M*Nf+(M/(M-1))*Ns);
        c2 = (1/(1-gamma))*((1/(M-1))*Ns);
        c3 = -M;
        c4 = (1-gamma);
        gm0 = g3;
        tm0 = t3;
        gm_ts = gLms(g_0<=gLms);
        tms = hyp_shock(gamma,c1,c2,c3,c4,gm0,tm0,gm_ts);
        xm_ts = -(M./(gm_ts.^2)) + (1/(1-gamma))*(-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gm_ts.^2)).*(tms-1);
        xRms(g_0<=gRms) = max([xRms(g_0<=gRms);xm_ts]);
        xRms(gRms<g_0) = max([max(xm_ts)*ones(size(gRms(gRms<g_0)));xRms(gRms<g_0)]);
        
        xms = [xLms,fliplr(xRms)];
        gms = [gLms,fliplr(gRms)];
        
    else
        xs = zeros(1,2*length(gs));
        gs = [gLs,fliplr(gRs)];
        
        % End-of-injection location of the left front ...
        xLms = -(M./(gLms.^2));
        % ... except left-moving waves
        xLms(gLms<g_s) = -(M./(gLms(gLms<g_s).^2)) + (1*(gLms(gLms<g_s)<=g_s) + (1/(1-gamma))*(gLms(gLms<g_s)> g_s)).*(-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gLms(gLms<g_s).^2))*(t4c-1);
        % ... except right-moving waves that have been eaten by the shock
        c1 = (M*Nf+(M/(M-1))*Ns);
        c2 = ((1/(M-1))*Ns);
        c3 = -M;
        c4 = 1/(1-gamma);
        gm0 = g4b;
        tm0 = t4b;
        gm = 1;
        gm_ts = gLms(logical((gm<gLms).*(gLms<g_s)));
        tms = hyp_shock(gamma,c1,c2,c3,c4,gm0,tm0,gm_ts);
        xm_ts = -(M./(gm_ts.^2)) + (-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gm_ts.^2)).*(tms-1);
        xLms(logical((gm<gLms).*(gLms<g_s))) = xm_ts;
        
        % Current location of the right front ...
        xRms =  (M./(gRms.^2)) + (1*(gRms> g_s) + (1/(1-gamma))*(gRms<=g_s)).*(-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gRms.^2))*(t4a-1);
        % ... except left-moving waves ...
        xRms(gRms<g_s) = (M./(gRms(gRms<g_s).^2));
        % ... except right-moving waves that were eaten by the shock before it changed direction ...
        c1 = (M*Nf+(M/(M-1))*Ns);
        c2 = ((1/(M-1))*Ns);
        c3 = M;
        c4 = (1/(1-gamma));
        gm0 = g2a;
        tm0 = t2a;
        gm_ts = gRms(logical((g_s<gRms).*(gRms<=g_0)));
        tms = hyp_shock(gamma,c1,c2,c3,c4,gm0,tm0,gm_ts);
        xm_ts = (M./(gm_ts.^2)) + (-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gm_ts.^2)).*(tms-1);
        xRms(logical((g_s<gRms).*(gRms<=g_0))) = xm_ts;
        % ... or after it changed direction ...
        c1 = (M*Nf+(M/(M-1))*Ns);
        c2 = ((1/(M-1))*Ns);
        c3 = M;
        c4 = 1;
        gm0 = g2b;
        tm0 = t2b;
        gm = M;
        gm_ts = gRms(logical((g_0<gRms).*(gRms<=gm)));
        tms = hyp_shock(gamma,c1,c2,c3,c4,gm0,tm0,gm_ts);
        xm_ts = (M./(gm_ts.^2)) + (-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gm_ts.^2)).*(tms-1);
        xRms(logical((g_0<gRms).*(gRms<=gm))) = xm_ts;
        % ... and right-moving waves of the left front that traveled past the right front
        c1 = (1/(1-gamma))*(M*Nf+(M/(M-1))*Ns);
        c2 = (1/(1-gamma))*((1/(M-1))*Ns);
        c3 = -M;
        c4 = (1-gamma);
        gm0 = g3;
        tm0 = t3;
        gm_ts = gLms(g_0<=gLms);
        tms = hyp_shock(gamma,c1,c2,c3,c4,gm0,tm0,gm_ts);
        xm_ts = -(M./(gm_ts.^2)) + (1/(1-gamma))*(-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gm_ts.^2)).*(tms-1);
        xRms(g_0<=gRms) = max([xRms(g_0<=gRms);xm_ts]);
        xRms(gRms<g_0) = max([max(xm_ts)*ones(size(gRms(gRms<g_0)));xRms(gRms<g_0)]);
        
        xms = [xLms,fliplr(xRms)];
        gms = [gLms,fliplr(gRms)];
        
    end
    
end