function hyp_case1_chars(xCs,tCs,gCs,M,gamma,Nf,Ns,gs,gmin)
    % Draw characteristics for case 1
    %   Note that Ns/Nf = 0 (flow only, Ruben's solution) is contained here
    
    % ---------------------------------------------------------
    % Key Values
    % ---------------------------------------------------------
    x1 = xCs(1); x2 = xCs(2); x3 = xCs(3); x4 = xCs(4);
    t1 = tCs(1); t2 = tCs(2); t3 = tCs(3); t4 = tCs(4);
    g1 = gCs(1); g2 = gCs(2); g3 = gCs(3); g4 = gCs(4);
    
    % ---------------------------------------------------------
    % Draw the characteristics
    % ---------------------------------------------------------
    
    % Injection -------------------------------------------------
    hyp_inj_chars(M,gamma,t1,gs);
    x1Ls = -(M./(gs.^2))*t1;
    x1Rs =  (M./(gs.^2))*t1;
    
    % Retreat -------------------------------------------------
    g2s = gs;
    
    x2Ls = x1Ls + (1/(1-gamma))*(-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(g2s.^2))*(t2-t1);
    x2Rs = x1Rs +               (-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(g2s.^2))*(t2-t1);
    
    for i=1:length(x2Ls)
        plot([x1Ls(i),x2Ls(i)],[t1,t2],'r-','linewidth',1.5);
    end
    for i=length(x2Rs):-1:1
        plot([x1Rs(i),x2Rs(i)],[t1,t2],'g-','linewidth',1.5);
    end
    plot([x2],[t2],'b*');
    
    
    % Chase ---------------------------------------------------------
    g3s = gs;
    x3Rs = x2Rs + (-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(g3s.^2))*(t3-t2);
    
    plot([x2,x3],[t2,t3],'b-','linewidth',1.5)
    for i=length(x2Rs):-1:1
        plot([x2Rs(i),x3Rs(i)],[t2,t3],'g-','linewidth',1.5);
    end
    plot([x3],[t3],'b*');
    
    
    % Sweep ---------------------------------------------------
    g4s = gs;
    
    c1 = (M*Nf+(M/(M-1))*Ns);
    c2 = ((1/(M-1))*Ns);
    c3 = M;
    c4 = 1/(1-gamma);
    gm0 = g3;
    tm0 = t3;
    gmF = g4;
    gms = linspace(gm0,gmF,1000);
    
    tms = hyp_shock(gamma,c1,c2,c3,c4,gm0,tm0,gms);
    xms =  (M./(gms.^2)) + (-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gms.^2)).*(tms-1.);
    
    t4s_lost = hyp_shock(gamma,c1,c2,c3,c4,gm0,tm0,g4s);
    x4s_lost = x3Rs + (-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(g4s.^2)).*(t4s_lost-t3);
    
    for i=1:length(x4s_lost)
        plot([x3Rs(i),x4s_lost(i)],[t3,t4s_lost(i)],'g-','linewidth',1.5)
    end
    
    plot(xms,tms,'b-','linewidth',1.5)
    plot(x4s_lost,t4s_lost,'bo')
    
    plot([x4],[t4],'b*')
    
end