function hyp_case6c2_chars(xCs,tCs,gCs,M,gamma,Nf,Ns,gs,gmin)
    % Draw characteristics for case 6c2
    %   Note that Ns/Nf-->-infty (negative slope only) is contained here
    
    % 6c2:  peak forms first and eats R0 before shock gets there, then the two collide
    % LR_RL, RR_RS, R_p, R_L0, R_LS, R_LL
    % p_R0 happens either before or after RR_RS, but is a non-event
    
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
    % Draw the characteristics
    % ---------------------------------------------------------
    
    % Injection -------------------------------------------------
    hyp_inj_chars(M,gamma,t1,gs);
    x1Ls = -(M./(gs.^2));
    x1Rs =  (M./(gs.^2));
    
    % Retreat, part a -------------------------------------------------
    x2aLs = -(M./(gs.^2)) + (1*(gs<=g_s) + (1/(1-gamma))*(gs> g_s)).*(-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gs.^2))*(t2a-1);
    x2aRs =  (M./(gs.^2)) + ((1/(1-gamma))*(gs<=g_s) + 1*(gs> g_s)).*(-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gs.^2))*(t2a-1);
    
    for i=1:length(gs)
        plot([x1Ls(i),x2aLs(i)],[t1,t2a],'r-','linewidth',1.5)
        plot([x1Rs(i),x2aRs(i)],[t1,t2a],'g-','linewidth',1.5)
    end
    plot([x2a],[t2a],'b*')
    
    % Retreat, part b -------------------------------------------------
    gp2b = sqrt( M*(M-1)*Nf/Ns + M - (2*(1-gamma)*M*(M-1)/(gamma*(t2b-1)))/Ns );
    
    t2bLs = t2b*ones(size(gs));
    t2bLs(gp2b<gs) = 1 + (2*(1-gamma)*M./(gamma*gs(gp2b<gs).^2)).*((-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./gs(gp2b<gs).^2).^-1);
    x2bLs = -(M./(gs.^2)) + (1*(gs<=g_s) + (1/(1-gamma))*(gs> g_s)).*(-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gs.^2)).*(t2bLs-1);
    
    t2bRs = t2b*ones(size(gs));
    t2bRs(gp2b<gs) = 1 + (2*(1-gamma)*M./(gamma*gs(gp2b<gs).^2)).*((-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./gs(gp2b<gs).^2).^-1);
    x2bRs =  (M./(gs.^2)) + ((1/(1-gamma))*(gs<=g_s) + 1*(gs> g_s)).*(-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gs.^2)).*(t2bRs-1);
    
    for i=1:1:length(gs)
        plot([x2aRs(i),x2bRs(i)],[t2a,t2bRs(i)],'g-','linewidth',1.5)
        plot([x2aLs(i),x2bLs(i)],[t2a,t2bLs(i)],'r-','linewidth',1.5)
        if gp2b<gs(i)
            plot(x2bRs(i),t2bRs(i),'bo')
        end
    end
    
    gps = linspace(g2a,gp2b,1000);
    tps = 1 + (2*(1-gamma)*M./(gamma*gps.^2)).*((-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./gps.^2).^-1);
    xps = (M./(gps.^2)) + (-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gps.^2)).*(tps-1);
    plot(xps,tps,'b-','linewidth',1.5)
    
    plot([x2b],[t2b],'b*')
    
    % Retreat, part c -------------------------------------------------
    c1 = (M*Nf+(M/(M-1))*Ns);
    c2 = ((1/(M-1))*Ns);
    c3 = M;
    c4 = (1/(1-gamma));
    gm0 = g2b;
    tm0 = t2b;
    gmF = g2c;
    
    t2cLs = t2c*ones(size(gs));
    t2cLs(g2c<=gs) = 1 + (2*(1-gamma)*M./(gamma*gs(g2c<gs).^2)).*((-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./gs(g2c<gs).^2).^-1);
    t2cLs(gp2b<=gs) = t2bLs(gp2b<=gs);
    x2cLs = -(M./(gs.^2)) + (1*(gs<=g_s) + (1/(1-gamma))*(gs> g_s)).*(-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gs.^2)).*(t2cLs-1);
    
    t2cRs = zeros(size(gs));
    t2cRs(g2c<=gs) = 1 + (2*(1-gamma)*M./(gamma*gs(g2c<gs).^2)).*((-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./gs(g2c<gs).^2).^-1);
    t2cRs(gp2b<gs) = t2bRs(gp2b<gs);
    t2cRs(gs<g2c) = hyp_shock(gamma,c1,c2,c3,c4,gm0,tm0,gs(gs<g2c));
    t2cRs(gs<g2b) = t2bRs(gs<g2b);
    x2cRs =  (M./(gs.^2)) + ((1/(1-gamma))*(gs<=g_s) + 1*(gs> g_s)).*(-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gs.^2)).*(t2cRs-1);
    
    for i=1:length(gs)
        plot([x2bLs(i),x2cLs(i)],[t2bLs(i),t2cLs(i)],'r-','linewidth',1.5)
        plot([x2bRs(i),x2cRs(i)],[t2bRs(i),t2cRs(i)],'g-','linewidth',1.5)
        if g2b<gs(i) && gs(i)<gp2b<gs(i)
            plot(x2cRs(i),t2cRs(i),'bo')
        end
    end
    
    gms = linspace(gm0,gmF,1000);
    tms = hyp_shock(gamma,c1,c2,c3,c4,gm0,tm0,gms);
    xms = (M./(gms.^2)) + (-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gms.^2)).*(tms-1);
    plot(xms,tms,'b-','linewidth',1.5)
    
    gps = linspace(gp2b,g2c,1000);
    tps = 1 + (2*(1-gamma)*M./(gamma*gps.^2)).*((-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./gps.^2).^-1);
    xps = (M./(gps.^2)) + (-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gps.^2)).*(tps-1);
    plot(xps,tps,'b-','linewidth',1.5)
    
    plot([x2c],[t2c],'b*')
    
    % Sweep, part a ---------------------------------------------------
    c1 = (1/(1-gamma))*(M*Nf+(M/(M-1))*Ns);
    c2 = (1/(1-gamma))*((1/(M-1))*Ns);
    c3 = -M;
    c4 = 1;
    gm0 = g2c;
    tm0 = t2c;
    
    t4aLs = t4a*ones(size(gs));
    t4aLs(g4a<gs) = hyp_shock(gamma,c1,c2,c3,c4,gm0,tm0,gs(g4a<gs));
    t4aLs(g2c<gs) = t2cLs(g2c<gs);
    x4aLs = -(M./(gs.^2)) + (1*(gs<=g_s) + (1/(1-gamma))*(gs> g_s)).*(-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gs.^2)).*(t4aLs-1);
    
    for i=1:length(gs)
        plot([x2cLs(i),x4aLs(i)],[t2cLs(i),t4aLs(i)],'r-','linewidth',1.5);
        if g4a<gs(i) && gs(i)<g2c
            plot(x4aLs(i),t4aLs(i),'bo')
        end
    end
    
    gmF = g4a;
    gms = linspace(gm0,gmF,1000);
    tms = hyp_shock(gamma,c1,c2,c3,c4,gm0,tm0,gms);
    xms =  -(M./(gms.^2)) + (1/(1-gamma))*(-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gms.^2)).*(tms-1.);
    
    plot(xms,tms,'b-','linewidth',1.5)
    plot([x4a],[t4a],'b*')
    
    if gmin<g4a
        % Sweep, part b ---------------------------------------------------
        c1 = (M*Nf+(M/(M-1))*Ns);
        c2 = ((1/(M-1))*Ns);
        c3 = -M;
        c4 = 1/(1-gamma);
        gm0 = g4a;
        tm0 = t4a;
        
        t4bLs = zeros(size(gs));
        t4bLs(g4a<=gs) = t4aLs(g4a<=gs);
        t4bLs(gs<g4a) = hyp_shock(gamma,c1,c2,c3,c4,gm0,tm0,gs(gs<g4a));
        x4bLs = -(M./(gs.^2)) + (1*(gs<=g_s) + (1/(1-gamma))*(gs> g_s)).*(-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gs.^2)).*(t4bLs-1);
        
        for i=1:1:length(gs)
            plot([x4aLs(i),x4bLs(i)],[t4aLs(i),t4bLs(i)],'r-','linewidth',1.5);
            if gs(i)<g_s
                plot(x4bLs(i),t4bLs(i),'bo')
            end
        end
        
        gmF = g4b;
        gms = linspace(gm0,gmF,1000);
        tms = hyp_shock(gamma,c1,c2,c3,c4,gm0,tm0,gms);
        xms =  -(M./(gms.^2)) + (-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gms.^2)).*(tms-1.);
        
        plot(xms,tms,'b-','linewidth',1.5)
        plot([x4b],[t4b],'b*')
    end
    
    axis([1.1*min([-x1,xCs]),1.1*max(xCs),0.,1.1*t4b])
    
end