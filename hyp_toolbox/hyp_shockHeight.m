function gm = hyp_shockHeight(gamma,c1,c2,c3,c4,g0,t0,tm,g_lim,t_lim)
    % Shock height as a function of time, general case (based on an implicit analytical expression)
    % Specifically: For a shock that has initial height g0 at time t0, 
    % calculate and return the shock heights gms at times tms.
    % 
    % This uses numerical iteration to invert the implicit expressions implemented in hyp_shock.
    % See hyp_shock for defintions for c1,c2,c3,c4.
    %
    % g_lim and t_lim are limits on the height and time (ie, we know the shock won't be shorter than gm at time tm).
    %
    % CWM, 2009
    
    gm_n = g0;
    tm_n = tm;
    gm_next = (g0+g_lim)/2;
    
    g_err = 1;
    while abs(g_err) > 1E-10
        if isnan(gm_next)
            gm = g_lim;
        else
            gm = gm_next;
        end
        tm_next = hyp_shock(gamma,c1,c2,c3,c4,g0,t0,gm_next);
        if tm_next>t_lim
            tm_next = t_lim;
        end
        R = tm_next - tm;
        g_err = R/g0;
        
        dtmdgm = (tm_next-tm_n)/(gm_next-gm_n);
        
        gm_n = gm_next;
        tm_n = tm_next;
        
        limiter = 0.25;
        gm_next = gm_n - limiter*(1/dtmdgm)*R;
        if g_lim>g0
            while (gm_next>g_lim || gm_next<g0)
                limiter = limiter/2;
                gm_next = gm_n - limiter*(1/dtmdgm)*R;
            end
        else
            while (gm_next<g_lim || gm_next>g0)
                limiter = limiter/2;
                gm_next = gm_n - limiter*(1/dtmdgm)*R;
            end
        end
    end
    
end