function [xs,gs,xms,gms] = hyp_case6c2_plume(xCs,tCs,gCs,M,gamma,Nf,Ns,t,gs,gmin)
    % Calculate plume shape for case 6c2
    %   Note that Ns/Nf-->-infty (negative slope only) is contained here
    
    % 6c2:  peak forms first and eats R0 before shock gets there, then the two collide
    % LR_RL, RR_RS, R_p, R_L0, R_LS, R_LL
    % p_R0 happens either before or after RR_RS, but is a non-event
    
    % Be quiet?  (Suppress nonessential output?)
    switches_quiet = 1;
    
    % ---------------------------------------------------------
    % Key Values
    % ---------------------------------------------------------
    x1 = xCs(1); x2a = xCs(2); x2b = xCs(3); x2c = xCs(4); x4a = xCs(5); x4b = xCs(6);
    t1 = tCs(1); t2a = tCs(2); t2b = tCs(3); t2c = tCs(4); t4a = tCs(5); t4b = tCs(6);
    g1 = gCs(1); g2a = gCs(2); g2b = gCs(3); g2c = gCs(4); g4a = gCs(5); g4b = gCs(6);
    
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
        
        gp = sqrt( M*(M-1)*Nf/Ns + M - (2*(1-gamma)*M*(M-1)/(gamma*(t-1)))/Ns );
        
        gLs(gp<gLs) = gp;
        gRs(gp<gRs) = gp;
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
        % ... except right-moving waves that were eaten by the peak ...
        gp_ts = gRms(gp<gRms);
        tps = 1 + (2*(1-gamma)*M./(gamma*gp_ts.^2)).*((-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./gp_ts.^2).^-1);
        xp_ts = (M./(gp_ts.^2)) + (-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gp_ts.^2)).*(tps-1);
        xRms(gp<gRms) = xp_ts;
        
        xms = [xLms,fliplr(xRms)];
        gms = [gLms,fliplr(gRms)];
        
    elseif t2b<t && t<=t2c
        
        % Retreat, part c
        if switches_quiet==0
            disp(['Retreat, part c:  t2b = ' num2str(t2b) ' < t = ' num2str(t) ' < t2c = ' num2str(t2c)])
        end
        
        c1 = (M*Nf+(M/(M-1))*Ns);
        c2 = ((1/(M-1))*Ns);
        c3 = M;
        c4 = 1/(1-gamma);
        gm0 = g2b;
        tm0 = t2b;
        
        tm = t;
        g_lim = g2c;
        t_lim = t2c;
        gm = hyp_shockHeight(gamma,c1,c2,c3,c4,gm0,tm0,tm,g_lim,t_lim);
        
        gp = sqrt( M*(M-1)*Nf/Ns + M - (2*(1-gamma)*M*(M-1)/(gamma*(t-1)))/Ns );
        
        gLs(gp<gLs) = gp;
        gRs(gRs<gm) = gm;
        gRs(gp<gRs) = gp;
        xLs = -(M./(gLs.^2)) + (1*(gLs <=g_s) + (1/(1-gamma))*(gLs> g_s)).*(-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gLs.^2))*(t-1);
        xRs =  (M./(gRs.^2)) + ((1/(1-gamma))*(gRs<=g_s) + 1*(gRs> g_s)).*(-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gRs.^2))*(t-1);
        
        xs = [xLs,fliplr(xRs)];
        gs = [gLs,fliplr(gRs)];
        
        % End-of-injection location of the left front ...
        xLms = -(M./(gLms.^2));
        % ... except left-moving waves
        xLms(gLms<g_s) = -(M./(gLms(gLms<g_s).^2)) + (-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gLms(gLms<g_s).^2))*(t-1);
        
        % Current location of the right front ...
        xRms =  (M./(gRms.^2)) + ((1/(1-gamma))*(gRms<=g_s) + 1*(gRms> g_s)).*(-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gRms.^2))*(t-1);
        % ... except left-moving waves ...
        xRms(gRms<g_s) = (M./(gRms(gRms<g_s).^2));
        % ... except right-moving waves that were eaten by the peak ...
        gp_ts = gRms(gp<gRms);
        tps = 1 + (2*(1-gamma)*M./(gamma*gp_ts.^2)).*((-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./gp_ts.^2).^-1);
        xp_ts = (M./(gp_ts.^2)) + (-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gp_ts.^2)).*(tps-1);
        xRms(gp<gRms) = xp_ts;
        % ... except right-moving waves that were eaten by the shock before it changed direction...
        c1 = (M*Nf+(M/(M-1))*Ns);
        c2 = ((1/(M-1))*Ns);
        c3 = M;
        c4 = (1/(1-gamma));
        gm0 = g2b;
        tm0 = t2b;
        gm_ts = gRms(logical((g_s<gRms).*(gRms<=gm)));
        tms = hyp_shock(gamma,c1,c2,c3,c4,gm0,tm0,gm_ts);
        xm_ts = (M./(gm_ts.^2)) + (-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gm_ts.^2)).*(tms-1);
        xRms(logical((g_s<gRms).*(gRms<=gm))) = xm_ts;
        
        xms = [xLms,fliplr(xRms)];
        gms = [gLms,fliplr(gRms)];
        
    elseif t2c<t && t<=t4a
        
        % Sweep, part a
        if switches_quiet==0
            disp(['Sweep, part a:  t2c = ' num2str(t2c) ' < t = ' num2str(t) ' < t4a = ' num2str(t4a)])
        end
        
        %   Evolve shock through left thicknesses from g2c to g4a (shock travels to the left)
        c1 = (1/(1-gamma))*(M*Nf+(M/(M-1))*Ns);
        c2 = (1/(1-gamma))*((1/(M-1))*Ns);
        c3 = -M;
        c4 = 1;
        gm0 = g2c;
        tm0 = t2c;
        
        tm = t;
        g_lim = g4a;
        t_lim = t4a;
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
        % ... except right-moving waves that were eaten by the peak ...
        gp_ts = gRms(g2c<gRms);
        tps = 1 + (2*(1-gamma)*M./(gamma*gp_ts.^2)).*((-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./gp_ts.^2).^-1);
        xp_ts = (M./(gp_ts.^2)) + (-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gp_ts.^2)).*(tps-1);
        xRms(g2c<gRms) = xp_ts;
        % ... except right-moving waves that were eaten by the shock before it changed direction ...
        c1 = (M*Nf+(M/(M-1))*Ns);
        c2 = ((1/(M-1))*Ns);
        c3 = M;
        c4 = (1/(1-gamma));
        gm0 = g2b;
        tm0 = t2b;
        gm_ts = gRms(logical((g2b<gRms).*(gRms<=g2c)));
        tms = hyp_shock(gamma,c1,c2,c3,c4,gm0,tm0,gm_ts);
        xm_ts = (M./(gm_ts.^2)) + (-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gm_ts.^2)).*(tms-1);
        xRms(logical((g2b<gRms).*(gRms<=g2c))) = xm_ts;
        
        xms = [xLms,fliplr(xRms)];
        gms = [gLms,fliplr(gRms)];
        
    elseif t4a<t
        t = min(t,t4b);
        
        % Sweep, part b
        if switches_quiet==0
            disp(['Sweep, part b:  t4a = ' num2str(t4a) ' < t = ' num2str(t)])
        end
        %   Evolve shock through left thicknesses from g_s (g4a) to 1 (shock travels to the left)
        c1 = (M*Nf+(M/(M-1))*Ns);
        c2 = ((1/(M-1))*Ns);
        c3 = -M;
        c4 = 1/(1-gamma);
        gm0 = g4a;
        tm0 = t4a;
        
        tm = t;
        g_lim = g4b;
        t_lim = t4b;
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
        % ... except right-moving waves that were eaten by the peak ...
        gp_ts = gRms(g2c<gRms);
        tps = 1 + (2*(1-gamma)*M./(gamma*gp_ts.^2)).*((-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./gp_ts.^2).^-1);
        xp_ts = (M./(gp_ts.^2)) + (-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gp_ts.^2)).*(tps-1);
        xRms(g2c<gRms) = xp_ts;
        % ... except right-moving waves that were eaten by the shock before it changed direction ...
        c1 = (M*Nf+(M/(M-1))*Ns);
        c2 = ((1/(M-1))*Ns);
        c3 = M;
        c4 = (1/(1-gamma));
        gm0 = g2b;
        tm0 = t2b;
        gm_ts = gRms(logical((g4a<gRms).*(gRms<=g2c)));
        tms = hyp_shock(gamma,c1,c2,c3,c4,gm0,tm0,gm_ts);
        xm_ts = (M./(gm_ts.^2)) + (-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gm_ts.^2)).*(tms-1);
        xRms(logical((g4a<gRms).*(gRms<=g2c))) = xm_ts;
        
        xms = [xLms,fliplr(xRms)];
        gms = [gLms,fliplr(gRms)];
        
    end
    
end