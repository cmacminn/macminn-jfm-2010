function [xs,gs,xms,gms] = hyp_case3_plume(xCs,tCs,gCs,M,gamma,Nf,Ns,t,gs,gmin)
    % Calculate plume shape for case 3
    %   Note that Ns/Nf = -(M-1) is contained here
    
    % Be quiet?  (Suppress nonessential output?)
    switches_quiet = 1;
    
    % ---------------------------------------------------------
    % Key Values
    % ---------------------------------------------------------
    x1 = xCs(1); x2 = xCs(2); x4 = xCs(3);
    t1 = tCs(1); t2 = tCs(2); t4 = tCs(3);
    g1 = gCs(1); g2 = gCs(2); g4 = gCs(3);
    
    % ---------------------------------------------------------
    % Calculate the plume shape at time t
    % ---------------------------------------------------------
    
    gLs = gs;
    gRs = gs;
    gLms = gs;
    gRms = gs;
    
    if t<0
        
        % ?
        error('You asked for the plume shape for some t<0... ?')
        xs = zeros(1,2*length(gs));
        gs = [gLs,fliplr(gRs)];
        xms = xs;
        gms = gs;
        
    elseif t==0
        
        xs = zeros(1,2*length(gs));
        gs = [gLs,fliplr(gRs)];
        xms = xs;
        gms = gs;
        
    elseif 0<t && t<=t1
        
        % Injection
        if switches_quiet==0
            disp(['Injection:  ' num2str(t) ' < t1 = ' num2str(t1)])
        end
        xLs = -(M./(gs.^2))*t;
        xRs =  (M./(gs.^2))*t;
        xs = [xLs,fliplr(xRs)];
        gs = [gLs,fliplr(gRs)];
        xms = xs;
        gms = gs;
        
    elseif t1<t && t<=t2
        
        % Retreat, part a
        if switches_quiet==0
            disp(['Retreat, part a:  t1 = ' num2str(t1) ' < t = ' num2str(t) ' < t2 = ' num2str(t2)])
        end
        xLs = -(M./(gs.^2)) + (1/(1-gamma))*(-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gs.^2))*(t-1);
        xRs =  (M./(gs.^2)) +               (-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gs.^2))*(t-1);
        xs = [xLs,fliplr(xRs)];
        gs = [gLs,fliplr(gRs)];
        
        xLms = -(M./(gLms.^2));
        xRms = xRs;
        xms = [xLms,fliplr(xRms)];
        gms = [gLms,fliplr(gRms)];
        
    elseif t2<t && t<=t4
        
        % Retreat, part b
        if switches_quiet==0
            disp(['Retreat, part b:  t2 = ' num2str(t2) ' < t = ' num2str(t) ' < t4 = ' num2str(t4)])
        end
        
        gp_t = sqrt(((M-1)/Ns)*((M*Nf+(M/(M-1))*Ns)-2*M*(1-gamma)/(gamma*(t-1))));
        
        gLs(gLs>gp_t) = gp_t;
        gRs(gRs>gp_t) = gp_t;
        xLs = -(M./(gLs.^2)) + (1/(1-gamma))*(-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gLs.^2))*(t-1);
        xRs =  (M./(gRs.^2)) +               (-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gRs.^2))*(t-1);
        
        xs = [xLs,fliplr(xRs)];
        gs = [gLs,fliplr(gRs)];
        
        gp_ts = gLms(gLms>gp_t);
        tps = 1 + 2*M./((gamma/(1-gamma))*(M*Nf+(M/(M-1))*Ns-(1/(M-1))*Ns*gp_ts.^2));
        xLms = -(M./(gLms.^2));
        xRms = xRs;
        xRms(gRms>gp_t) =  (M./(gp_ts.^2)) +               (-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gp_ts.^2)).*(tps-1);
        
        xms = [xLms,fliplr(xRms)];
        gms = [gLms,fliplr(gRms)];
        
    elseif t4<=t
        
        xs = zeros(1,2*length(gs));
        gs = [gLs,fliplr(gRs)];
        
        gp_t = 1;
        gp_ts = gLms(gLms>gp_t);
        tps = 1 + 2*M./((gamma/(1-gamma))*(M*Nf+(M/(M-1))*Ns-(1/(M-1))*Ns*gp_ts.^2));
        xLms = -(M./(gLms.^2));
        xRms =  (M./(gRs.^2)) +               (-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gRms.^2))*(t-1);
        xRms(gRms>gp_t) =  (M./(gp_ts.^2)) +               (-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)./(gp_ts.^2)).*(tps-1);
        
        xms = [xLms,fliplr(xRms)];
        gms = [gLms,fliplr(gRms)];
        
    end
    
end