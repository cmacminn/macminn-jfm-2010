function hyp_case2_chars(xCs,tCs,gCs,M,gamma,Nf,Ns,gs,gmin)
    % Draw characteristics for case 2
    
    % ---------------------------------------------------------
    % Key Values
    % ---------------------------------------------------------
    x1 = xCs(1); x2a = xCs(2); x2b = xCs(3); x4 = xCs(4);
    t1 = tCs(1); t2a = tCs(2); t2b = tCs(3); t4 = tCs(4);
    g1 = gCs(1); g2a = gCs(2); g2b = gCs(3); g4 = gCs(4);
    
    % ---------------------------------------------------------
    % Draw the characteristics
    % ---------------------------------------------------------
    
    % Injection -------------------------------------------------
    hyp_inj_chars(M,gamma,t1,gs);
    x1Ls = -(M./(gs.^2))*t1;
    x1Rs =  (M./(gs.^2))*t1;
    
    % Retreat, part a -------------------------------------------------
    
    % Draw the left and right fronts
    x2aLs = x1Ls + (1/(1-gamma))*(-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gs.^2))*(t2a-t1);
    x2aRs = x1Rs +               (-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gs.^2))*(t2a-t1);
    for i=1:1:length(gs)
        plot([x1Ls(i),x2aLs(i)],[t1,t2a],'r-','linewidth',1.5);
        plot([x1Rs(i),x2aRs(i)],[t1,t2a],'g-','linewidth',1.5);
    end
    plot([x2a],[t2a],'b*');
    plot([x2a],[t2a],'bo');
    
    
    % Retreat, part b ---------------------------------------------------
    
    % Draw the left and right fronts
    for i=1:1:length(gs)
        if g2a<=gs(i)
            t2bLs(i) = t2a;
            t2bRs(i) = t2a;
        elseif ( g2b<gs(i) && gs(i)<g2a )
            t2bLs(i) = 1 + 2*M./((gamma/(1-gamma))*(M*Nf+(M/(M-1))*Ns-(1/(M-1))*Ns*gs(i)^2));
            t2bRs(i) = 1 + 2*M./((gamma/(1-gamma))*(M*Nf+(M/(M-1))*Ns-(1/(M-1))*Ns*gs(i)^2));
        else
            t2bLs(i) = t2b;
            t2bRs(i) = t2b;
        end
    end
    x2bLs = x2aLs + (1/(1-gamma))*(-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gs.^2)).*(t2bLs-t2a);
    x2bRs = x2aRs +               (-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gs.^2)).*(t2bRs-t2a);
    for i=1:1:length(gs)
        plot([x2aLs(i),x2bLs(i)],[t2a,t2bLs(i)],'r-','linewidth',1.5)
        plot([x2aRs(i),x2bRs(i)],[t2a,t2bRs(i)],'g-','linewidth',1.5)
    end
    gs_lost = (g2b<gs).*(gs<g2a);
    
    % Draw the peak path
    gp0 = g2a;
    gpF = g2b;
    gps = linspace(gp0,gpF,1000);
    tps = 1 + 2*M./((gamma/(1-gamma))*(M*Nf+(M/(M-1))*Ns-(1/(M-1))*Ns*gps.^2));
    xps = (M./(gps.^2)) + (-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gps.^2)).*(tps-1.);
    plot(xps,tps,'b-','linewidth',1.5)
    
    plot(x2bLs(gs_lost~=0),t2bLs(gs_lost~=0),'bo')
    plot(x2bRs(gs_lost~=0),t2bRs(gs_lost~=0),'bo')
    
    plot([x2b],[t2b],'b*')
    
    
    % Sweep ---------------------------------------------------
    
    % Draw the right front
    c1 = (M*Nf+(M/(M-1))*Ns);
    c2 = ((1/(M-1))*Ns);
    c3 = M;
    c4 = 1/(1-gamma);
    gm0 = g2b;
    tm0 = t2b;
    for i=1:1:length(gs)
        if g2b<=gs(i)
            t4Rs(i) = t2bRs(i);
        else
            t4Rs(i) = hyp_shock(gamma,c1,c2,c3,c4,gm0,tm0,gs(i));
        end
    end
    x4Rs = x2bRs + (-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gs.^2)).*(t4Rs-t2bRs);
    for i=1:1:length(gs)
        plot([x2bRs(i),x4Rs(i)],[t2bRs(i),t4Rs(i)],'g-','linewidth',1.5)
    end
    gs_lost = gs<g2b;
    plot(x4Rs(gs_lost~=0),t4Rs(gs_lost~=0),'bo')
    
    % Draw the shock path
    gm0 = g2b;
    tm0 = t2b;
    gmF = 1;
    gms = linspace(gm0,gmF,1000);
    tms = hyp_shock(gamma,c1,c2,c3,c4,gm0,tm0,gms);
    xms =  (M./(gms.^2)) + (-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gms.^2)).*(tms-1.);
    plot(xms,tms,'b-','linewidth',1.5)
    
    plot([x4],[t4],'b*')
    
    axis([1.1*min(x1Ls),1.1*x4,0.,1.1*t4])
    
end