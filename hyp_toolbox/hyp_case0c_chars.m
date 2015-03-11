function hyp_case0c_chars(xCs,tCs,gCs,M,gamma,Nf,Ns,gs,gmin)
    % Draw characteristics for case 0c
    %   Note that Ns/Nf --> infty (positive slope only) is contained here
    
    % 0c:  a peak forms, then a shock forms on the left
    % RL_LR, LL_LS, L_p, L_RS, L_RR
    
    % ---------------------------------------------------------
    % Key Values
    % ---------------------------------------------------------
    x1 = xCs(1); x2a = xCs(2); x2b = xCs(3); x2c = xCs(4); x4a = xCs(5); x4b = xCs(6);
    t1 = tCs(1); t2a = tCs(2); t2b = tCs(3); t2c = tCs(4); t4a = tCs(5); t4b = tCs(6);
    g1 = gCs(1); g2a = gCs(2); g2b = gCs(3); g2c = gCs(4); g4a = gCs(5); g4b = gCs(6);
    
    g_s = sqrt(M*(M-1)*Nf/Ns+M);
    xLS = -1/((M-1)*Nf/Ns+1);
    xRS =  1/((M-1)*Nf/Ns+1);
    tLL_LS = 1 + (1-gamma)*(M-1)/((M-1)*Nf+Ns);
    xLL_LS = xLS;
    tRL_LR = 1 - 2*((1-gamma)/gamma)/(Nf-Ns);
    xRL_LR = -(2-gamma)/(M*gamma);
    
    % ---------------------------------------------------------
    % Draw the characteristics
    % ---------------------------------------------------------
    
    % Injection -------------------------------------------------
    hyp_inj_chars(M,gamma,t1,gs);
    x1Ls = -(M./(gs.^2));
    x1Rs =  (M./(gs.^2));
    
    % Retreat, Part a -------------------------------------------------
    x2aLs = -(M./(gs.^2)) + ((1/(1-gamma))*(gs<=g_s)+1*(gs> g_s)).*(-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gs.^2))*(t2a-1);
    x2aRs =  (M./(gs.^2)) + ((1/(1-gamma))*(gs> g_s)+1*(gs<=g_s)).*(-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gs.^2))*(t2a-1);
    
    for i=1:length(gs)
        plot([x1Ls(i),x2aLs(i)],[1,t2a],'r-','linewidth',1.5);
        plot([x1Rs(i),x2aRs(i)],[1,t2a],'g-','linewidth',1.5);
    end
    plot([x2a],[t2a],'b*');
    
    % Retreat, Part b -------------------------------------------------
    gp2b = sqrt( (M*(M-1)*Nf/Ns + M) + 2*M*(M-1)*(1-gamma)/(gamma*Ns*(t2b-1)) );
    
    t2bLs = zeros(size(gs));
    t2bLs(gs<=gp2b) = t2b;
    t2bLs(gp2b<gs) = 1 - (2*M*(1-gamma)./(gamma*gs(gp2b<gs).^2)).*((-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gs(gp2b<gs).^2)).^-1);
    
    t2bRs = t2bLs;
    
    x2bLs = -(M./(gs.^2)) + ((1/(1-gamma))*(gs<=g_s)+1*(gs> g_s)).*(-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gs.^2)).*(t2bLs-1);
    x2bRs =  (M./(gs.^2)) + ((1/(1-gamma))*(gs> g_s)+1*(gs<=g_s)).*(-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gs.^2)).*(t2bRs-1);
    for i=1:length(gs)
        plot([x2aLs(i),x2bLs(i)],[t2a,t2bLs(i)],'r-','linewidth',1.5);
        plot([x2aRs(i),x2bRs(i)],[t2a,t2bRs(i)],'g-','linewidth',1.5);
    end
    
    gps = linspace(g2a,gp2b,1000);
    tps = 1 - (2*M*(1-gamma)./(gamma*gps.^2)).*((-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gps.^2)).^-1);
    xps = -(M./(gps.^2)) + (-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gps.^2)).*(tps-1);
    plot(xps,tps,'b-','linewidth',1.5)
    
    plot(x2b,t2b,'b*')
    
    % Retreat, Part c -------------------------------------------------
    c1 = (M*Nf+(M/(M-1))*Ns);
    c2 = ((1/(M-1))*Ns);
    c3 = -M;
    c4 = 1/(1-gamma);
    gm0 = g2b;
    tm0 = t2b;
    
    t2cLs = t2bLs;
    t2cLs(g2b<gs) = hyp_shock(gamma,c1,c2,c3,c4,gm0,tm0,gs(g2b<gs));
    t2cLs(g2c<gs) = 1 - (2*M*(1-gamma)./(gamma*gs(g2c<gs).^2)).*((-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gs(g2c<gs).^2)).^-1);
    
    t2cRs = t2bRs;
    t2cRs(gs<g2c) = t2c;
    t2cRs(g2c<=gs) = 1 - (2*M*(1-gamma)./(gamma*gs(g2c<gs).^2)).*((-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gs(g2c<gs).^2)).^-1);
    
    x2cLs = -(M./(gs.^2)) + ((1/(1-gamma))*(gs<=g_s)+1*(gs> g_s)).*(-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gs.^2)).*(t2cLs-1);
    x2cRs =  (M./(gs.^2)) + ((1/(1-gamma))*(gs> g_s)+1*(gs<=g_s)).*(-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gs.^2)).*(t2cRs-1);
    for i=1:length(gs)
        plot([x2bLs(i),x2cLs(i)],[t2bLs(i),t2cLs(i)],'r-','linewidth',1.5);
        plot([x2bRs(i),x2cRs(i)],[t2bRs(i),t2cRs(i)],'g-','linewidth',1.5);
        if gs(i)>g2b
            plot(x2cLs(i),t2cLs(i),'bo')
        end
        if gs(i)>g2c
            plot(x2cLs(i),t2cLs(i),'bo')
            plot(x2cRs(i),t2cRs(i),'bo')
        end
    end
    
    plot(xLL_LS,tLL_LS,'b*')
    
    gps = linspace(g2a,g2c,1000);
    tps = 1 - (2*M*(1-gamma)./(gamma*gps.^2)).*((-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gps.^2)).^-1);
    xps = -(M./(gps.^2)) + (-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gps.^2)).*(tps-1);
    plot(xps,tps,'b-','linewidth',1.5)
    
    gms = linspace(g2b,g2c,1000);
    tms = hyp_shock(gamma,c1,c2,c3,c4,gm0,tm0,gms);
    xms = -(M./(gms.^2)) + (-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gms.^2)).*(tms-1);
    plot(xms,tms,'b-','linewidth',1.5)
    
    plot([x2c],[t2c],'b*')
    
    % Sweep, part a ---------------------------------------------------
    c1 = (1/(1-gamma))*(M*Nf+(M/(M-1))*Ns);
    c2 = (1/(1-gamma))*((1/(M-1))*Ns);
    c3 = M;
    c4 = 1;
    gm0 = g2c;
    tm0 = t2c;
    gmF = g_s;
    
    t4aRs = zeros(size(gs));
    t4aRs(gs<=g_s) = t4a;
    t4aRs(g_s<gs) = hyp_shock(gamma,c1,c2,c3,c4,gm0,tm0,gs(g_s<gs));
    t4aRs(g2c<gs) = t2cRs(g2c<gs);
    x4aRs = (M./(gs.^2)) + ((1/(1-gamma))*(gs> g_s)+1*(gs<=g_s)).*(-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gs.^2)).*(t4aRs-1);
    
    for i=1:length(gs)
        plot([x2cRs(i),x4aRs(i)],[t2cRs(i),t4aRs(i)],'g-','linewidth',1.5);
        if g_s<gs(i)
            plot(x4aRs(i),t4aRs(i),'bo')
        end
    end
    
    gms = linspace(gm0,gmF,1000);
    tms = hyp_shock(gamma,c1,c2,c3,c4,gm0,tm0,gms);
    xms =  (M./(gms.^2)) + (1/(1-gamma))*(-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gms.^2)).*(tms-1.);
    plot(xms,tms,'b-','linewidth',1.5)
    
    plot([x4a],[t4a],'b*')
    
    % Sweep, part b ---------------------------------------------------
    c1 = (M*Nf+(M/(M-1))*Ns);
    c2 = ((1/(M-1))*Ns);
    c3 = M;
    c4 = 1/(1-gamma);
    gm0 = g4a;
    tm0 = t4a;
    gmF = g4b;
    
    t4bRs = zeros(size(gs));
    t4bRs(g_s<=gs) = t4aRs(g_s<=gs);
    t4bRs(gs<g_s) = hyp_shock(gamma,c1,c2,c3,c4,gm0,tm0,gs(gs<g_s));
    x4bRs = (M./(gs.^2)) + ((1/(1-gamma))*(gs> g_s)+1*(gs<=g_s)).*(-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gs.^2)).*(t4bRs-1);
    
    for i=1:length(gs)
        plot([x4aRs(i),x4bRs(i)],[t4aRs(i),t4bRs(i)],'g-','linewidth',1.5);
        if gs(i)<g4a
            plot(x4bRs(i),t4bRs(i),'bo')
        end
    end
    
    gms = linspace(gm0,gmF,1000);
    tms = hyp_shock(gamma,c1,c2,c3,c4,gm0,tm0,gms);
    xms =  (M./(gms.^2)) + (-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gms.^2)).*(tms-1.);
    plot(xms,tms,'b-','linewidth',1.5)
    
    plot([x4b],[t4b],'b*')
    
    axis([1.1*min([-x1,xCs]),1.1*max(xCs),0.,1.1*t4b])
    
end