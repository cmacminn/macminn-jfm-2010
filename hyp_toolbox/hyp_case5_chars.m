function hyp_case5_chars(xCs,tCs,gCs,M,gamma,Nf,Ns,gs,gmin)
    % Draw characteristics for case 5
    %   Note that this case only exists if M<sqrt(2/gamma-1)
    
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
    x2Ls = x1Ls + (1/(1-gamma))*(-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gs.^2))*(t2-t1);
    x2Rs = x1Rs +               (-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gs.^2))*(t2-t1);
    
    for i=1:1:length(gs)
        plot([x1Ls(i),x2Ls(i)],[t1,t2],'r-','linewidth',1.5);
    end
    for i=length(gs):-1:1
        plot([x1Rs(i),x2Rs(i)],[t1,t2],'g-','linewidth',1.5);
    end
    plot([x2],[t2],'b*');
    
    
    % Chase ---------------------------------------------------
    x3Ls = x2Ls + (1/(1-gamma))*(-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gs.^2))*(t3-t2);
    for i=1:length(gs)
        plot([x2Ls(i),x3Ls(i)],[t2,t3],'r-','linewidth',1.5)
    end
    
    % Draw the shock path
    plot([x2,x3],[t2,t3],'b-','linewidth',1.5)
    
    plot([x3],[t3],'bo')
    plot([x3],[t3],'b*')
    
    % Sweep ---------------------------------------------------
    c1 = (1/(1-gamma))*(M*Nf+(M/(M-1))*Ns);
    c2 = (1/(1-gamma))*((1/(M-1))*Ns);
    c3 = -M;
    c4 = (1-gamma);
    gm0 = g3;
    tm0 = t3;
    
    for i=1:1:length(gs)
        if gs(i)==g3
            t4Ls(i) = t3;
        elseif gs(i)<g4
            t4Ls(i) = t4;
        else
            t4Ls(i) = hyp_shock(gamma,c1,c2,c3,c4,gm0,tm0,gs(i));
        end
    end
    x4Ls = x3Ls + (1/(1-gamma))*(-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gs.^2)).*(t4Ls-t3);
    
    for i=1:1:length(gs)
        plot([x3Ls(i),x4Ls(i)],[t3,t4Ls(i)],'r-','linewidth',1.5)
    end
    
    % Draw the shock path
    gms = linspace(g3,g4,1000);
    tms = hyp_shock(gamma,c1,c2,c3,c4,gm0,tm0,gms);
    xms = -(M./(gms.^2)) + (1/(1-gamma))*(-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gms.^2)).*(tms-1);
    plot(xms,tms,'b-','linewidth',1.5)
    
    plot(x4Ls,t4Ls,'bo')
    plot([x4],[t4],'b*')
    
    axis([-1.1*x1,1.1*x4,0.,1.1*t4])
    
end