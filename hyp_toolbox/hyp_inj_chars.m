function hyp_inj_chars(M,gamma,t,gs)
    % Draw characteristics for the injection period, taking Nfi>>Nsi
    
    % ---------------------------------------------------------
    % Key Values
    % ---------------------------------------------------------
    % Injection
    t1 = 1;
    x1 = M;
    
    % ---------------------------------------------------------
    % Draw the characteristics
    % ---------------------------------------------------------
    
    figure(1)
    clf;
    xlabel('\xi');
    ylabel('\tau');
    hold on;
    
    x1Ls = -(M./(gs.^2))*t;
    x1Rs =  (M./(gs.^2))*t;
    xs = [x1Ls,fliplr(x1Rs)];
    
    for i=1:1:length(x1Ls)
        plot([0,x1Ls(i)],[0,t1],'r-','linewidth',1.5);
    end
    for i=length(x1Rs):-1:1
        plot([0,x1Rs(i)],[0,t1],'g-','linewidth',1.5);
    end
    
    axis([-1.1*x1,1.1*x1,0,1.1*t1])
    
end