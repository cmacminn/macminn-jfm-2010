function hyp_case6a_chars(xCs,tCs,gCs,M,gamma,Nf,Ns,gs,gmin)
    % Draw characteristics for case 6a
    %   Note that Ns/Nf-->-infty (negative slope only) is contained here
    
    % 6a:  a shock forms on the right and hits the bottom before a peak forms (no peak)
    % RR_RS, R_R0, R_RL, R_LR, R_L0, R_LS, R_LL
    
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
    c1 = (M*Nf+(M/(M-1))*Ns);
    c2 = ((1/(M-1))*Ns);
    c3 = M;
    c4 = 1/(1-gamma);
    gm0 = g2a;
    tm0 = t2a;
    gmF = g2b;
    
    t2bLs = t2b*ones(size(gs));
    x2bLs = -(M./(gs.^2)) + (1*(gs<=g_s) + (1/(1-gamma))*(gs> g_s)).*(-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gs.^2)).*(t2bLs-1);
    
    t2bRs = t2b*ones(size(gs));
    t2bRs(gs<g_0) = hyp_shock(gamma,c1,c2,c3,c4,gm0,tm0,gs(gs<g_0));
    t2bRs(gs<g_s) = t2a;
    x2bRs =  (M./(gs.^2)) + ((1/(1-gamma))*(gs<=g_s) + 1*(gs> g_s)).*(-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gs.^2)).*(t2bRs-1);
    
    for i=1:1:length(gs)
        plot([x2aRs(i),x2bRs(i)],[t2a,t2bRs(i)],'g-','linewidth',1.5)
        plot([x2aLs(i),x2bLs(i)],[t2a,t2bLs(i)],'r-','linewidth',1.5)
        if g_s<gs(i) && gs(i)<g_0
            plot(x2bRs(i),t2bRs(i),'bo')
        end
    end
    
    gms = linspace(gm0,gmF,1000);
    tms = hyp_shock(gamma,c1,c2,c3,c4,gm0,tm0,gms);
    xms = (M./(gms.^2)) + (-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gms.^2)).*(tms-1);
    plot(xms,tms,'b-','linewidth',1.5)    
    
    plot([x2b],[t2b],'b*')
    
    % Retreat, part c -------------------------------------------------
    c1 = (M*Nf+(M/(M-1))*Ns);
    c2 = ((1/(M-1))*Ns);
    c3 = M;
    c4 = 1;
    gm0 = g2b;
    tm0 = tms(end);
    gmF = g2c;
    
    t2cLs = t2c*ones(size(gs));
    x2cLs = -(M./(gs.^2)) + (1*(gs<=g_s) + (1/(1-gamma))*(gs> g_s)).*(-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gs.^2)).*(t2cLs-1);
    
    t2cRs = zeros(size(gs));
    t2cRs(g_0<gs) = hyp_shock(gamma,c1,c2,c3,c4,gm0,tm0,gs(g_0<gs));
    t2cRs(gs<=g_0) = t2bRs(gs<=g_0);
    x2cRs =  (M./(gs.^2)) + ((1/(1-gamma))*(gs<=g_s) + 1*(gs> g_s)).*(-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gs.^2)).*(t2cRs-1);
    
    for i=1:1:length(gs)
        plot([x2bLs(i),x2cLs(i)],[t2bLs(i),t2cLs(i)],'r-','linewidth',1.5)
        plot([x2bRs(i),x2cRs(i)],[t2bRs(i),t2cRs(i)],'g-','linewidth',1.5)
        if g_0<gs(i)
            plot(x2cRs(i),t2cRs(i),'bo')
        end
    end
    
    gms = linspace(gm0,gmF,1000);
    tms = hyp_shock(gamma,c1,c2,c3,c4,gm0,tm0,gms);
    xms = (M./(gms.^2)) + (-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gms.^2)).*(tms-1);
    plot(xms,tms,'b-','linewidth',1.5)
    
    plot([x2c],[t2c],'b*')
    
    % Chase ---------------------------------------------------------
    t3Ls = t3*ones(size(gs));
    x3Ls = -(M./(gs.^2)) + (1*(gs<=g_s) + (1/(1-gamma))*(gs> g_s)).*(-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gs.^2)).*(t3Ls-1);
    
    for i=1:length(gs)
        plot([x2cLs(i),x3Ls(i)],[t2cLs(i),t3Ls(i)],'r-','linewidth',1.5)
    end
    
    plot([x2c,x3],[t2c,t3],'b-','linewidth',1.5)
    plot([x3],[t3],'b*')
    
    % Sweep, part a ---------------------------------------------------
    c1 = (1/(1-gamma))*(M*Nf+(M/(M-1))*Ns);
    c2 = (1/(1-gamma))*((1/(M-1))*Ns);
    c3 = -M;
    c4 = (1-gamma);
    gm0 = g3;
    tm0 = t3;
    
    t4aLs = t4a*ones(size(gs));
    t4aLs(g_0<gs) = hyp_shock(gamma,c1,c2,c3,c4,gm0,tm0,gs(g_0<gs));
    x4aLs = -(M./(gs.^2)) + (1*(gs<=g_s) + (1/(1-gamma))*(gs> g_s)).*(-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gs.^2)).*(t4aLs-1);
    
    for i=1:1:length(gs)
        plot([x3Ls(i),x4aLs(i)],[t3Ls(i),t4aLs(i)],'r-','linewidth',1.5);
        if g_0<gs(i)
            plot(x4aLs(i),t4aLs(i),'bo')
        end
    end
    
    gmF = g4a;
    gms = linspace(gm0,gmF,1000);
    tms = hyp_shock(gamma,c1,c2,c3,c4,gm0,tm0,gms);
    xms =  -(M./(gms.^2)) + (1/(1-gamma))*(-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gms.^2)).*(tms-1.);
    
    plot(xms,tms,'b-','linewidth',1.5)
    plot([x4a],[t4a],'b*')
    
    % Sweep, part b ---------------------------------------------------
    c1 = (1/(1-gamma))*(M*Nf+(M/(M-1))*Ns);
    c2 = (1/(1-gamma))*((1/(M-1))*Ns);
    c3 = -M;
    c4 = 1;
    gm0 = g4a;
    tm0 = t4a;
    
    t4bLs = t4b*ones(size(gs));
    t4bLs(g_0<gs) = t4aLs(g_0<gs);
    t4bLs(g_s<gs) = hyp_shock(gamma,c1,c2,c3,c4,gm0,tm0,gs(g_s<gs));
    x4bLs = -(M./(gs.^2)) + (1*(gs<=g_s) + (1/(1-gamma))*(gs> g_s)).*(-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gs.^2)).*(t4bLs-1);
    
    for i=1:1:length(gs)
        plot([x4aLs(i),x4bLs(i)],[t4aLs(i),t4bLs(i)],'r-','linewidth',1.5);
        if g_s<gs(i) && gs(i)<g_0
            plot(x4bLs(i),t4bLs(i),'bo')
        end
    end
    
    gmF = g4b;
    gms = linspace(gm0,gmF,1000);
    tms = hyp_shock(gamma,c1,c2,c3,c4,gm0,tm0,gms);
    xms =  -(M./(gms.^2)) + (1/(1-gamma))*(-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gms.^2)).*(tms-1.);
    
    plot(xms,tms,'b-','linewidth',1.5)
    plot([x4b],[t4b],'b*')
    
    % Sweep, part c ---------------------------------------------------
    c1 = (M*Nf+(M/(M-1))*Ns);
    c2 = ((1/(M-1))*Ns);
    c3 = -M;
    c4 = 1/(1-gamma);
    gm0 = g4b;
    tm0 = t4b;
    
    t4cLs = zeros(size(gs));
    t4cLs(g_s<=gs) = t4bLs(g_s<=gs);
    t4cLs(gs<g_s) = hyp_shock(gamma,c1,c2,c3,c4,gm0,tm0,gs(gs<g_s));
    x4cLs = -(M./(gs.^2)) + (1*(gs<=g_s) + (1/(1-gamma))*(gs> g_s)).*(-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gs.^2)).*(t4cLs-1);
    
    for i=1:1:length(gs)
        plot([x4bLs(i),x4cLs(i)],[t4bLs(i),t4cLs(i)],'r-','linewidth',1.5);
        if gs(i)<g_s
            plot(x4cLs(i),t4cLs(i),'bo')
        end
    end
    
    gmF = g4c;
    gms = linspace(gm0,gmF,1000);
    tms = hyp_shock(gamma,c1,c2,c3,c4,gm0,tm0,gms);
    xms =  -(M./(gms.^2)) + (-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gms.^2)).*(tms-1.);
    
    plot(xms,tms,'b-','linewidth',1.5)
    plot([x4c],[t4c],'b*')
    
    axis([1.1*min([-x1,xCs]),1.1*max(xCs),0.,1.1*t4c])
    
end