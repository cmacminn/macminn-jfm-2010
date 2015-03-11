function [Fs,dFs,d2Fs] = hyp_fluxFn(hs,M,Nf,Ns)
    
    gs = (M-1)*hs+1;
    
    Gs = (M*Nf+((M+1)/(M-1))*Ns) - ((1/(M-1))*Ns)*gs - (M*Nf+(M/(M-1))*Ns)./gs;
    dGs = -((1/(M-1))*Ns) + (M*Nf+(M/(M-1))*Ns)./(gs.^2);
    d2Gs = -2*(M*Nf+(M/(M-1))*Ns)./(gs.^3);
    
    Fs = Gs/(M-1);
    dFs = dGs/((M-1)^2);
    d2Fs = d2Gs/((M-1)^3);
    
end