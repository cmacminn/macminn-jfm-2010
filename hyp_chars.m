function [xCs,tCs,eff] = hyp_chars(M,gamma,Nf,Ns,hmin)
    % Draw characteristics
    %   Characteristics will stop abruptly at the time when max(hs(t))<hmin
    
    % Add the helper files to the path
    addpath('./hyp_toolbox')
    
    % Set a default value for hmin
    if exist('hmin')==0
        hmin = 0;
    end
    
    [xCs,tCs,hCs,caseName] = hyp_crits(M,gamma,Nf,Ns,hmin);
    gCs = (M-1)*hCs+1;
    
    eff = hyp_eff(M,gamma,Nf,Ns,hmin);
    
    % Choose a good range of waves to follow in space-time (i.e., for which to draw the characteristics)
    hs = [0.,logspace(-2.,0.,10)];
    gs = (M-1)*hs + 1;
    gmin = (M-1)*hmin + 1;
    
    figure(gcf)
    clf
    
    feval(str2func(['hyp_case' caseName '_chars']),xCs,tCs,gCs,M,gamma,Nf,Ns,gs,gmin);
    
    xLmost = min([-xCs(1),xCs]);
    xRmost = max([-xCs(1),xCs]);
    axis([1.1*xLmost,1.1*xRmost,0,1.1*max(tCs)])
    
    xlabel('\xi')
    ylabel('\tau')
    drawnow()
    
end