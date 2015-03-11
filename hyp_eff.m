function eff = hyp_eff(M,gamma,Nf,Ns,hmin)
    % Calculate the efficiency factor
    
    % Add the helper files to the path
    addpath('./hyp_toolbox')
    
    % Tolerance for comparison of floats
    tol = 1E-14;
    
    % Set a default value for hmin
    if exist('hmin')==0
        hmin = 0;
    end
    
    [xCs,tCs,hCs,caseName] = hyp_crits(M,gamma,Nf,Ns,hmin);
    
    % Left and right edges of the footprint
    xLm = min([-xCs(1),xCs]);
    xRm = max([-xCs(1),xCs]);
    
    % Length of the footprint
    xfp = abs(xRm-xLm);
    
    % Efficiency
    eff = 2/xfp;
    
end