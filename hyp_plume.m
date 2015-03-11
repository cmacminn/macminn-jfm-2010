function [xss,hss,xmss,hmss,ts] = hyp_plume(M,gamma,Nf,Ns,ts,hs,hmin)
    % Get the plume shape at times ts.
    %   - Returns the positions of thicknesses hs, but hs will be modified according to 
    %       which thicknesses still exist at time t.  hs>max(h(t)) will be set to max(h(t)).
    %   - If hs includes 0 and 1, plume extent and max height will be exact.
    %   - Will also return the profile of residual CO2.
    %   - Plume-shape calculation doesn't directly respect the cut-off thickness hmin,
    %       but t4 comes from hyp_crits which does (esp. for Ns/Nf near -M)
    
    % Add the helper files to the path
    addpath('./hyp_toolbox')
    
    % Tolerance for comparison of floats
    tol = 1E-14;
    
    % Set a default value for hmin
    if exist('hmin')==0
        hmin = 0;
    end
    
    hs_orig = hs;
    
    % Convert hs to gs
    gmin = (M-1)*hmin + 1;
    gs_orig = (M-1)*hs_orig + 1;
    
    [xCs,tCs,hCs,caseName] = hyp_crits(M,gamma,Nf,Ns,hmin);
    gCs = (M-1)*hCs + 1;
    
    if strcmp(caseName,'6b3')
        error(['case ' caseName ' not implemented -- figure out what parameters got you here, and tell Chris.'])
    end
    
    xss = [];
    hss = [];
    xmss = [];
    hmss = [];
    
    for ti=1:length(ts)
        t = ts(ti);
        if t==-1
            t = tCs(end);
            ts(ti) = t;
        end
        
        [xs,gs,xms,gms] = feval(str2func(['hyp_case' caseName '_plume']),xCs,tCs,gCs,M,gamma,Nf,Ns,t,gs_orig,gmin);
        hs = (gs-1)/(M-1);
        hms = (gms-1)/(M-1);
        
        xss = [xss;xs];
        hss = [hss;hs];
        hmss = [hmss;hms];
        xmss = [xmss;xms];
        
    end
    
    % So that we don't have to worry about whether or not these were calculated when plotting
    if isempty(xmss)
        xmss = zeros(size(xss));
        hmss = zeros(size(hss));
    end
    
end