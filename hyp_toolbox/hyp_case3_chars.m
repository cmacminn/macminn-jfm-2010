function hyp_case3_chars(xCs,tCs,gCs,M,gamma,Nf,Ns,gs,gmin)
    % Draw characteristics for case 3
    %   Note that Ns/Nf = -(M-1) is contained here
    
    % ---------------------------------------------------------
    % Key Values
    % ---------------------------------------------------------
    x1 = xCs(1); x2 = xCs(2); x4 = xCs(3);
    t1 = tCs(1); t2 = tCs(2); t4 = tCs(3);
    g1 = gCs(1); g2 = gCs(2); g4 = gCs(3);
    
    % ---------------------------------------------------------
    % Draw the characteristics
    % ---------------------------------------------------------
    
    % Injection -------------------------------------------------
    hyp_inj_chars(M,gamma,t1,gs);
    x1Ls = -(M./(gs.^2))*t1;
    x1Rs =  (M./(gs.^2))*t1;
    
    % Retreat / Sweep -------------------------------------------------
    
    g2s = gs;
    
    x2Ls = x1Ls + (1/(1-gamma))*(-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(g2s.^2)).*(t2-t1);
    x2Rs = x1Rs +               (-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(g2s.^2)).*(t2-t1);
    
    for i=1:1:length(x2Ls)
        plot([x1Ls(i),x2Ls(i)],[t1,t2],'r-','linewidth',1.5);
    end
    for i=length(x2Rs):-1:1
        plot([x1Rs(i),x2Rs(i)],[t1,t2],'g-','linewidth',1.5);
    end
    plot([x2],[t2],'b*');
    
    
    % Peak path ---------------------------------------------------
    gps = gs(1:end-1);
    
    tps = 1 + 2*M./((gamma/(1-gamma))*(M*Nf+(M/(M-1))*Ns-(1/(M-1))*Ns*gps.^2));
    xps = (M./(gps.^2)) + (-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gps.^2)).*(tps-1);
    
    for i=1:length(xps)
        plot([x2Ls(i),xps(i)],[t2,tps(i)],'r-','linewidth',1.5)
    end
    for i=length(xps):-1:1
            plot([x2Rs(i),xps(i)],[t2,tps(i)],'g-','linewidth',1.5)
    end
    
    plot([x2,fliplr(xps)],[t2,fliplr(tps)],'b-','linewidth',1.5)
    plot(xps,tps,'bo')
    
    plot([x4],[t4],'b*')
    
    axis([1.1*min(x1Ls),1.1*x4,0.,1.1*t4])
end