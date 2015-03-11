function [xCs,tCs,hCs,caseName] = hyp_crits(M,gamma,Nf,Ns,hmin)
    % Decide which interval and case this is, and calculate key values
    
    % Add the helper files to the path
    addpath('./hyp_toolbox')
    
    % Tolerance for comparison of floats
    tol = 1E-14;
    
    % Be quiet?  (Suppress nonessential output?)
    switches_quiet = 0;
    
    % Set a default value for hmin
    if exist('hmin')==0
        hmin = 0;
    end
    
    % Enforce an hmin for Ns/Nf near -M if necessary
    if abs(Nf>tol)
        if abs(Ns/Nf-(-M))<1E-4
            hmin = max(1E-5,hmin);
        end
    end
    
    % Convert h to g
    gmin = (M-1)*hmin+1;
    
    Ns01 = 1;
    NsRJ = 0;
    Ns12 = -(M-1)*(2-gamma)/((2-gamma)+M*gamma);
    Ns23 = -(M-1)*(2-gamma)/((2-gamma)+gamma/M);
    NsSC = -(M-1);
    Ns34 = -(M-1)*(2-gamma)/((2-gamma)-gamma/M);
    Ns45 = -(M-1)*(2-gamma)/((2-gamma)-M*gamma);
    Ns56 = -M;
    
    % if Nf==0 and Ns==0
    if abs(Nf)<tol && abs(Ns)<tol
        disp('Both Ns and Nf are 0, so plume does not migrate.')
        xCs = []; tCs = []; hCs = []; caseName = [];
        return
    end
    
    
    if ( Ns56<Ns/Nf && Ns/Nf<=Ns01 )
        % Cases 1-5
        if ( Ns12<=Ns/Nf && Ns/Nf<=Ns01 )
            caseName = '1';
        elseif ( Ns23<Ns/Nf && Ns/Nf<Ns12 )
            caseName = '2';
        elseif ( Ns34<=Ns/Nf && Ns/Nf<=Ns23 )
            if switches_quiet==0
                if Ns/Nf==-(M-1)
                    disp('Special case:  degenerate')
                end
            end
            caseName = '3';
        elseif ( (Ns56-tol)<Ns/Nf && Ns/Nf<Ns34 )
            if abs(Ns/Nf-Ns56)<abs(Ns56*tol)
                hmin = max(hmin,1E-8);
                gmin = (M-1)*hmin+1;
            end
            if sqrt(2/gamma-1)<M
                caseName = '4';
            else
                if Ns45<Ns/Nf
                    caseName = '4';
                else
                    caseName = '5';
                end
            end
        end
    else
        % Cases 0 and 6, including + and - slope only
        
        if Ns/Nf<=Ns56
            if abs(Nf)<tol && switches_quiet==0
                disp('Special case:  negative slope only')
            end
            
            g_s = sqrt(M*(M-1)*Nf/Ns+M); % stationary point
            g_0 = M*(M-1)*Nf/Ns + M; % zero of the flux function
            tRR_RS = 1 - (1-gamma)*(M-1)/((M-1)*Nf+Ns);
            tLR_RL = 1 + 2*((1-gamma)/gamma)/(Nf-Ns);
            
            % Time at which the shock would reach xR0
            c1 = (M*Nf+(M/(M-1))*Ns);
            c2 = ((1/(M-1))*Ns);
            c3 = M;
            c4 = 1/(1-gamma);
            gm0 = g_s;
            tm0 = tRR_RS;
            gmF = g_0;
            tR0 = hyp_shock(gamma,c1,c2,c3,c4,gm0,tm0,gmF);
            
            % Time at which the peak would reach xR0
            tp0 = 1 + (2*(1-gamma)*M/(gamma*g_0^2))*((-(1/(M-1))*Ns+(M*Nf+(M/(M-1))*Ns)/g_0^2)^-1);
            
            % Time at which the shock would reach xRL
            c1 = (M*Nf+(M/(M-1))*Ns);
            c2 = ((1/(M-1))*Ns);
            c3 = M;
            c4 = 1;
            gm0 = g_0;
            tm0 = tR0;
            gmF = M;
            tR_RL = hyp_shock(gamma,c1,c2,c3,c4,gm0,tm0,gmF);
            
            if tRR_RS<tLR_RL
                
                if tR0<tp0
                    
                    if tR_RL<tLR_RL
                        % xRR_RS, xR_R0, xR_RL
                        caseName = '6a';
                    elseif tR0<tLR_RL
                        % xRR_RS, xR_R0, xLR_RL, xR_p
                        caseName = '6b1';
                    else
                        % xRR_RS, xLR_RL, xR_R0, xR_p
                        caseName = '6b2';
                    end
                    
                else
                    % xRR_RS, xLR_RL, xp_R0, xR_p
                    caseName = '6b3';
                end
                
            else
                
                if tR0<tp0
                    % xLR_RL, xRR_RS, xR_R0, xR_p
                    caseName = '6c1';
                else
                    
                    % xLR_RL, xRR_RS, xp_R0, xR_p or xLR_RL, xp_R0, xRR_RS, xR_p
                    % (doesn't matter which)
                    caseName = '6c2';
                end
                
            end
        elseif Ns01<Ns/Nf
            if abs(Nf)<tol && switches_quiet==0
                disp('Special case:  positive slope only')
            end
            
            g_s = sqrt(M*(M-1)*Nf/Ns+M); % stationary point
            tLL_LS = 1 + (1-gamma)*(M-1)/((M-1)*Nf+Ns);
            tRL_LR = 1 - 2*((1-gamma)/gamma)/(Nf-Ns);
            
            % Time at which the shock would reach xLR
            c1 = (M*Nf+(M/(M-1))*Ns);
            c2 = ((1/(M-1))*Ns);
            c3 = -M;
            c4 = 1/(1-gamma);
            gm0 = g_s;
            tm0 = tLL_LS;
            gmF = M;
            tL_LR = hyp_shock(gamma,c1,c2,c3,c4,gm0,tm0,gmF);
            
            if tLL_LS<tRL_LR
                if tL_LR<tRL_LR
                    caseName = '0a';
                else
                    caseName = '0b';
                end
            else
                caseName = '0c';
            end
        end
    end
    
    if switches_quiet==0
        disp(['M = ' num2str(M) ' , Gamma = ' num2str(gamma) ' , Ns/Nf = ' num2str(Ns/Nf) ' -- case ' caseName])
    end
    [xCs,tCs,gCs] = feval(str2func(['hyp_case' caseName '_crits']),M,gamma,Nf,Ns,gmin);
    hCs = (gCs-1)/(M-1);
    
end