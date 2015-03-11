function hyp_case0a_chars(xCs,tCs,gCs,M,gamma,Nf,Ns,gs,gmin)
    % Draw characteristics for case 0a
    %   Note that Ns/Nf --> infty (positive slope only) is contained here
    
    % 0a:  a shock forms on the left and hits the bottom before a peak forms (no peak)
    % LL_LS, L_LR, L_RL, L_RS, L_RR
    
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
    % Draw the characteristics
    % ---------------------------------------------------------
    
    % Injection -------------------------------------------------
    hyp_inj_chars(M,gamma,t1,gs);
    x1Ls = -(M./(gs.^2));
    x1Rs =  (M./(gs.^2));
    
    % Retreat, Part a -------------------------------------------------
    x2aLs = x1Ls + ((1/(1-gamma))*(gs<=g_s)+1*(gs> g_s)).*(-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gs.^2))*(t2a-1);
    x2aRs = x1Rs + ((1/(1-gamma))*(gs> g_s)+1*(gs<=g_s)).*(-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gs.^2))*(t2a-1);
    
    for i=1:length(gs)
        plot([x1Ls(i),x2aLs(i)],[t1,t2a],'r-','linewidth',1.5);
    end
    for i=length(gs):-1:1
        plot([x1Rs(i),x2aRs(i)],[t1,t2a],'g-','linewidth',1.5);
    end
    plot([x2a],[t2a],'b*');
    
    
    % Retreat, Part b -------------------------------------------------
    c1 = (M*Nf+(M/(M-1))*Ns);
    c2 = ((1/(M-1))*Ns);
    c3 = -M;
    c4 = 1/(1-gamma);
    gm0 = g_s;
    tm0 = tLL_LS;
    
    for i=1:1:length(gs)
        t2bRs(i) = t2b;
    end
    x2bRs = x2aRs + ((1/(1-gamma))*(gs> g_s)+1*(gs<=g_s)).*(-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gs.^2)).*(t2bRs-t2a);
    for i=length(gs):-1:1
        plot([x2aRs(i),x2bRs(i)],[t2a,t2bRs(i)],'g-','linewidth',1.5);
    end
    
    for i=1:1:length(gs)
        if gs(i)<=g_s
            t2bLs(i) = tLL_LS;
        else
            t2bLs(i) = hyp_shock(gamma,c1,c2,c3,c4,gm0,tm0,gs(i));
        end
    end
    x2bLs = x2aLs + ((1/(1-gamma))*(gs<=g_s)+1*(gs> g_s)).*(-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gs.^2)).*(t2bLs-t2a);
    for i=1:1:length(gs)
        plot([x2aLs(i),x2bLs(i)],[t2a,t2bLs(i)],'r-','linewidth',1.5);
    end
    plot(xLL_LS,tLL_LS,'b*')
    
    gms = linspace(g_s,g2b,1000);
    tms = hyp_shock(gamma,c1,c2,c3,c4,g_s,tLL_LS,gms);
    xms = -(M./(gms.^2))*t1 + (-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gms.^2)).*(tms-t1);
    plot(xms,tms,'b-','linewidth',1.5)
    
    plot([x2b],[t2b],'b*')
    
    
    % Chase ---------------------------------------------------------
    x3Rs = x2bRs + ((1/(1-gamma))*(gs> g_s)+1*(gs<=g_s)).*(-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gs.^2))*(t3-t2b);
    for i=length(gs(gs<=g3)):-1:1
        plot([x2bRs(i),x3Rs(i)],[t2b,t3],'g-','linewidth',1.5);
    end
    x3Ls = x3;
    
    plot([x2b,x3],[t2b,t3],'b-','linewidth',1.5)
    
    plot([x3],[t3],'b*')
    
    
    % Sweep ---------------------------------------------------
    
    c1 = (1/(1-gamma))*(M*Nf+(M/(M-1))*Ns);
    c2 = (1/(1-gamma))*((1/(M-1))*Ns);
    c3 = M;
    c4 = 1;
    gm0 = g3;
    tm0 = t3;
    gmF = g_s;
    
    for i=1:1:length(gs)
        if gs(i)<=g_s
            t4aRs(i) = t4a;
        else
            t4aRs(i) = hyp_shock(gamma,c1,c2,c3,c4,gm0,tm0,gs(i));
        end
    end
    x4aRs = x3Rs + ((1/(1-gamma))*(gs> g_s)+1*(gs<=g_s)).*(-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gs.^2)).*(t4aRs-t3);
    
    for i=1:1:length(gs(gs<g3))
        plot([x3Rs(i),x4aRs(i)],[t3,t4aRs(i)],'g-','linewidth',1.5);
    end
    
    gms = linspace(gm0,gmF,1000);
    tms = hyp_shock(gamma,c1,c2,c3,c4,gm0,tm0,gms);
    xms =  (M./(gms.^2)) + (1/(1-gamma))*(-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gms.^2)).*(tms-1.);
    plot(xms,tms,'b-','linewidth',1.5)
    
    plot([x4a],[t4a],'b*')
    
    
    c1 = (M*Nf+(M/(M-1))*Ns);
    c2 = ((1/(M-1))*Ns);
    c3 = M;
    c4 = 1/(1-gamma);
    gm0 = g4a;
    tm0 = t4a;
    gmF = g4b;
    
    for i=1:1:length(gs)
        if gs(i)>=g_s
            t4bRs(i) = t4aRs(i);
        else
            t4bRs(i) = hyp_shock(gamma,c1,c2,c3,c4,gm0,tm0,gs(i));
        end
    end
    x4bRs = x4aRs + ((1/(1-gamma))*(gs> g_s)+1*(gs<=g_s)).*(-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gs.^2)).*(t4bRs-t4aRs);
    
    for i=1:1:length(gs(gs<g3))
        plot([x4aRs(i),x4bRs(i)],[t4aRs(i),t4bRs(i)],'g-','linewidth',1.5);
    end
    
    gms = linspace(gm0,gmF,1000);
    tms = hyp_shock(gamma,c1,c2,c3,c4,gm0,tm0,gms);
    xms =  (M./(gms.^2)) + (-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gms.^2)).*(tms-1.);
    plot(xms,tms,'b-','linewidth',1.5)
    
    plot(x4bRs(gs<g3),t4bRs(gs<g3),'bo')
    plot([x4b],[t4b],'b*')
    
    axis([1.1*min([-x1,xCs]),1.1*max(xCs),0.,1.1*t4b])
    
end