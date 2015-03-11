function [xs,gs,xms,gms] = hyp_case0a_plume(xCs,tCs,gCs,M,gamma,Nf,Ns,t,gs,gmin)
    % Calculate plume shape for case 0a
    %   Note that Ns/Nf --> infty (positive slope only) is contained here
    
    % 0a:  a shock forms on the left and hits the bottom before a peak forms (no peak)
    % LL_LS, L_LR, L_RL, L_RS, L_RR
    
    % Be quiet?  (Suppress nonessential output?)
    switches_quiet = 1;
    
    % ---------------------------------------------------------
    % Key Values
    % ---------------------------------------------------------
    x1 = xCs(1); x2a = xCs(2); x2b = xCs(3); x3 = xCs(4); x4a = xCs(5); x4b = xCs(6);
    t1 = tCs(1); t2a = tCs(2); t2b = tCs(3); t3 = tCs(4); t4a = tCs(5); t4b = tCs(6);
    g1 = gCs(1); g2a = gCs(2); g2b = gCs(3); g3 = gCs(4); g4a = gCs(5); g4b = gCs(6);
    
    g_s = sqrt(M*(M-1)*Nf/Ns+M);
    xLS = -1/((M-1)*Nf/Ns+1);
    xRS =  1/((M-1)*Nf/Ns+1);
    tLL_LS = 1 + (1-gamma)*(M-1)/((M-1)*Nf+Ns);
    xLL_LS = xLS;
    tRL_LR = 1 - 2*((1-gamma)/gamma)/(Nf-Ns);
    xRL_LR = -(2-gamma)/(M*gamma);
    
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
        gLms = gs;
        gRms = gs;
        
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
        xLs = -(M./(gs.^2)) + ((1/(1-gamma))*(gs<=g_s)+1*(gs> g_s)).*(-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gs.^2))*(t-1);
        xRs =  (M./(gs.^2)) + ((1/(1-gamma))*(gs> g_s)+1*(gs<=g_s)).*(-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gs.^2))*(t-1);
        xs = [xLs,fliplr(xRs)];
        gs = [gLs,fliplr(gRs)];
        
        % End-of-injection location of the left front ...
        xLms = -(M./(gLms.^2));
        % ... except left-moving waves
        xLms(gLms>g_s) = xLs(gLms>g_s);
        
        % Current location of the right front ...
        xRms = xRs;
        % ... except left-moving waves
        xRms(gRms>g_s) = (M./(gRms(gRms>g_s).^2));
        
        xms = [xLms,fliplr(xRms)];
        gms = [gLms,fliplr(gRms)];
        
    elseif t2a<t && t<t2b
        
        % Retreat, part b
        if switches_quiet==0
            disp(['Retreat, part b:  t2a = ' num2str(t2a) ' < t = ' num2str(t) ' < t2b = ' num2str(t2b)])
        end
        c1 = (M*Nf+(M/(M-1))*Ns);
        c2 = ((1/(M-1))*Ns);
        c3 = -M;
        c4 = 1/(1-gamma);
        gm0 = g_s;
        tm0 = tLL_LS;
        
        tm = t;
        g_lim = g2b;
        t_lim = t2b;
        gm_t = hyp_shockHeight(gamma,c1,c2,c3,c4,gm0,tm0,tm,g_lim,t_lim);
        
        gLs(gLs<gm_t) = gm_t;
        xLs = -(M./(gLs.^2)) + ((1/(1-gamma))*(gLs<=g_s)+1*(gLs> g_s)).*(-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gLs.^2))*(t-1);
        xRs =  (M./(gRs.^2)) + ((1/(1-gamma))*(gRs> g_s)+1*(gRs<=g_s)).*(-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gRs.^2))*(t-1);
        
        xs = [xLs,fliplr(xRs)];
        gs = [gLs,fliplr(gRs)];
        
        % End-of-injection location of the left front ...
        xLms = -(M./(gLms.^2));
        % ... except left-moving waves ...
        xLms(gLms>g_s) = xLs(gLms>g_s);
        % ... except left-moving waves that have been eaten by the shock
        gm_ts = gLms(logical((gLms>g_s).*(gLms<=gm_t)));
        tms = hyp_shock(gamma,c1,c2,c3,c4,gm0,tm0,gm_ts);
        xm_ts = -(M./(gm_ts.^2)) + (-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gm_ts.^2)).*(tms-1);
        xLms(logical((gLms>g_s).*(gLms<=gm_t))) = xm_ts;
        
        % Current location of the right front ...
        xRms = xRs;
        % ... except left-moving waves
        xRms(gRms>g_s) = (M./(gRms(gRms>g_s).^2));
        
        xms = [xLms,fliplr(xRms)];
        gms = [gLms,fliplr(gRms)];
        
    elseif t2b<=t && t<=t3
        
        % Chase
        if switches_quiet==0
            disp(['Chase:  t2b = ' num2str(t2b) ' < t = ' num2str(t) ' < t3 = ' num2str(t3)])
        end
        gLs(gLs<=g3) = g3;
        xLs = (x2b + (Nf/(1-gamma))*(t-t2b))*ones(1,length(gs));
        xRs =  (M./(gRs.^2)) + ((1/(1-gamma))*(gRs> g_s)+1*(gRs<=g_s)).*(-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gRs.^2)).*(t-1);
        
        xs = [xLs,fliplr(xRs)];
        gs = [gLs,fliplr(gRs)];
        
        % End-of-injection location of the left front ...
        xLms = -(M./(gLms.^2));
        % ... except left-moving waves that have been eaten by the shock (all of them)
        c1 = (M*Nf+(M/(M-1))*Ns);
        c2 = ((1/(M-1))*Ns);
        c3 = -M;
        c4 = 1/(1-gamma);
        gm0 = g_s;
        tm0 = tLL_LS;
        gm_t = M;
        gm_ts = gLms(logical((gLms>g_s).*(gLms<=gm_t)));
        tms = hyp_shock(gamma,c1,c2,c3,c4,gm0,tm0,gm_ts);
        xm_ts = -(M./(gm_ts.^2)) + (-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gm_ts.^2)).*(tms-1);
        xLms(logical((gLms>g_s).*(gLms<=gm_t))) = xm_ts;
        
        % Current location of the right front ...
        xRms = xRs;
        % ... except left-moving waves
        xRms(gRms>g_s) = (M./(gRms(gRms>g_s).^2));
        
        xms = [xLms,fliplr(xRms)];
        gms = [gLms,fliplr(gRms)];
        
    elseif t3<t && t<=t4a
        
        % Sweep, part a
        if switches_quiet==0
            disp(['Sweep, part a:  t3 = ' num2str(t3) ' < t = ' num2str(t) ' < t4a = ' num2str(t4a)])
        end
        c1 = (1/(1-gamma))*(M*Nf+(M/(M-1))*Ns);
        c2 = (1/(1-gamma))*((1/(M-1))*Ns);
        c3 = M;
        c4 = 1;
        gm0 = g3;
        tm0 = t3;
        
        tm = t;
        g_lim = g4a;
        t_lim = t4a;
        gm_t = hyp_shockHeight(gamma,c1,c2,c3,c4,gm0,tm0,tm,g_lim,t_lim);
        xm_t = (M./(gm_t.^2)) + ((1/(1-gamma))*(gm_t> g_s)+1*(gm_t<=g_s)).*(-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gm_t.^2))*(t-1);
        
        gLs(:) = gm_t;
        gRs(gRs>gm_t)=gm_t;
        xLs = xm_t*ones(1,length(gs));
        xRs = (M./(gRs.^2)) + ((1/(1-gamma))*(gRs> g_s)+1*(gRs<=g_s)).*(-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gRs.^2)).*(t-1);
        
        xs = [xLs,fliplr(xRs)];
        gs = [gLs,fliplr(gRs)];
        
        % End-of-injection location of the left front ...
        xLms = -(M./(gLms.^2));
        % ... except left-moving waves that have been eaten by the shock (all of them)
        c1 = (M*Nf+(M/(M-1))*Ns);
        c2 = ((1/(M-1))*Ns);
        c3 = -M;
        c4 = 1/(1-gamma);
        gm0 = g_s;
        tm0 = tLL_LS;
        gm_tL = M;
        gm_ts = gLms(logical((gLms>g_s).*(gLms<=gm_tL)));
        tms = hyp_shock(gamma,c1,c2,c3,c4,gm0,tm0,gm_ts);
        xm_ts = -(M./(gm_ts.^2)) + (-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gm_ts.^2)).*(tms-1);
        xLms(logical((gLms>g_s).*(gLms<=gm_tL))) = xm_ts;
        
        % Current location of the right front ...
        xRms = xRs;
        % ... except left-moving waves
        xRms(gRms>g_s) = (M./(gRms(gRms>g_s).^2));
        
        xms = [xLms,fliplr(xRms)];
        gms = [gLms,fliplr(gRms)];
        
    elseif t4a<t && t<=t4b
        % Sweep, part b
        if switches_quiet==0
            disp(['Sweep, part b:  t4a = ' num2str(t4a) ' < t = ' num2str(t) ' < t4b = ' num2str(t4b)])
        end
        c1 = (M*Nf+(M/(M-1))*Ns);
        c2 = ((1/(M-1))*Ns);
        c3 = M;
        c4 = 1/(1-gamma);
        gm0 = g4a;
        tm0 = t4a;
        
        tm = t;
        g_lim = 1;
        t_lim = t4b;
        gm_t = hyp_shockHeight(gamma,c1,c2,c3,c4,gm0,tm0,tm,g_lim,t_lim);
        xm_t = (M./(gm_t.^2)) + ((1/(1-gamma))*(gm_t> g_s)+1*(gm_t<=g_s)).*(-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gm_t.^2))*(t-1);
        
        gLs(:) = gm_t;
        gRs(gRs>gm_t)=gm_t;
        xLs = xm_t*ones(1,length(gs));
        xRs = (M./(gRs.^2)) + ((1/(1-gamma))*(gRs> g_s)+1*(gRs<=g_s)).*(-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gRs.^2)).*(t-1);
        
        xs = [xLs,fliplr(xRs)];
        gs = [gLs,fliplr(gRs)];
        
        % End-of-injection location of the left front ...
        xLms = -(M./(gLms.^2));
        % ... except left-moving waves that have been eaten by the shock (all of them)
        c1 = (M*Nf+(M/(M-1))*Ns);
        c2 = ((1/(M-1))*Ns);
        c3 = -M;
        c4 = 1/(1-gamma);
        gm0 = g_s;
        tm0 = tLL_LS;
        gm_tL = M;
        gm_ts = gLms(logical((gLms>g_s).*(gLms<=gm_tL)));
        tms = hyp_shock(gamma,c1,c2,c3,c4,gm0,tm0,gm_ts);
        xm_ts = -(M./(gm_ts.^2)) + (-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gm_ts.^2)).*(tms-1);
        xLms(logical((gLms>g_s).*(gLms<=gm_tL))) = xm_ts;
        
        % Current location of the right front ...
        xRms = xRs;
        % ... except left-moving waves ...
        xRms(gRms>g_s) = (M./(gRms(gRms>g_s).^2));
        % ... except right-moving waves that have been eaten by the shock
        c1 = (M*Nf+(M/(M-1))*Ns);
        c2 = ((1/(M-1))*Ns);
        c3 = M;
        c4 = 1/(1-gamma);
        gm0 = g4a;
        tm0 = t4a;
        gm_tR = gm_t;
        gm_ts = gRms(logical((gRms<g_s).*(gRms>=gm_tR)));
        tms = hyp_shock(gamma,c1,c2,c3,c4,gm0,tm0,gm_ts);
        xm_ts = (M./(gm_ts.^2)) + (-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gm_ts.^2)).*(tms-1);
        xRms(logical((gRms<g_s).*(gRms>=gm_tR))) = xm_ts;
        
        xms = [xLms,fliplr(xRms)];
        gms = [gLms,fliplr(gRms)];
        
    elseif t4b<t
        
        xs = zeros(1,2*length(gs));
        gs = [gLs,fliplr(gRs)];
        xms = gs;
        gms = gs;
        
        % End-of-injection location of the left front ...
        xLms = -(M./(gLms.^2));
        % ... except left-moving waves that have been eaten by the shock (all of them)
        c1 = (M*Nf+(M/(M-1))*Ns);
        c2 = ((1/(M-1))*Ns);
        c3 = -M;
        c4 = 1/(1-gamma);
        gm0 = g_s;
        tm0 = tLL_LS;
        gm_tL = M;
        gm_ts = gLms(logical((gLms>g_s).*(gLms<=gm_tL)));
        tms = hyp_shock(gamma,c1,c2,c3,c4,gm0,tm0,gm_ts);
        xm_ts = -(M./(gm_ts.^2)) + (-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gm_ts.^2)).*(tms-1);
        xLms(logical((gLms>g_s).*(gLms<=gm_tL))) = xm_ts;
        
        % End-of-injection location of the right front ...
        xRms = (M./(gRms.^2));
        % ... except right-moving waves that have been eaten by the shock (all of them)
        c1 = (M*Nf+(M/(M-1))*Ns);
        c2 = ((1/(M-1))*Ns);
        c3 = M;
        c4 = 1/(1-gamma);
        gm0 = g4a;
        tm0 = t4a;
        gm_tR = 1;
        gm_ts = gRms(logical((gRms<g_s).*(gRms>=gm_tR)));
        tms = hyp_shock(gamma,c1,c2,c3,c4,gm0,tm0,gm_ts);
        xm_ts = (M./(gm_ts.^2)) + (-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gm_ts.^2)).*(tms-1);
        xRms(logical((gRms<g_s).*(gRms>=gm_tR))) = xm_ts;
        
        xms = [xLms,fliplr(xRms)];
        gms = [gLms,fliplr(gRms)];
        
    end
    
end