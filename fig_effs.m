% ---------------------------------------------------------
% Efficiency Factor vs. Ns/Nf
% WARNING: This takes a few minutes.
% ---------------------------------------------------------

clear all; close all; clc;

if exist('./figures/fig_effs.mat')==2
    %% Found previous results -- use them?
    use = input(['Found results file -- use it? (y or n) '],'s');
else
    use = 0;
end
if isempty(use) || strcmp(use,'y')==1
    %% Load existing file
    load('./figures/fig_effs.mat');
else
    
    % M = 10;
    M = 15;
    % M = 30;
    Nf = 1;
    hmin = 0;
    
    gammas = [0.2 0.3 0.4];
    
    for igamma = 1:length(gammas)
        gamma = gammas(igamma)
        
        Ns01 = 1;
        NsRJ = 0;
        Ns12 = -(2-gamma)*(M-1)/((2-gamma)+M*gamma);
        Ns23 = -M*(M-1)*(2-gamma)/(M*(2-gamma)+gamma);
        NsSC = -(M-1);
        Ns34 = -M*(M-1)*(2-gamma)/(M*(2-gamma)-gamma);
        Ns45 = -(2-gamma)*(M-1)/((2-gamma)-M*gamma);
        Ns56 = -M;
        
        N = 200;
        % N = 20;
        if M<sqrt(2/gamma-1)
            Nss = [ -logspace(log10(5*M),log10(-Ns56),N) linspace(Ns56,Ns45,N) linspace(Ns45,Ns34,N) linspace(Ns34,Ns23,N) ...
                    linspace(Ns23,Ns12,N) linspace(Ns12,Ns01,N) logspace(log10(Ns01),log10(5*M),N) ];
        else
            % Nss = [ -logspace(log10(10*M),log10(-Ns56),N) linspace(Ns56,Ns34,N) linspace(Ns34,Ns23,N) ...
            %         linspace(Ns23,Ns12,N) linspace(Ns12,Ns01,N) logspace(log10(Ns01),log10(2*M),N) ];
            Nss = [ -logspace(log10(2*M),log10(-Ns56),2*N) linspace(Ns56,Ns34,2*N) linspace(Ns34,Ns23,2*N) ...
                    linspace(Ns23,Ns12,N) linspace(Ns12,Ns01,N) logspace(log10(Ns01),log10(10),N) ];
        end
        
        % Round to six decimal places to avoid occasional singular values
        Nss = (1E-6)*round((1E6)*Nss);
        
        eval(['effs_hyp_' num2str(igamma) ' = zeros(size(Nss));'])
        for iNs = 1:length(Nss)
            Ns = Nss(iNs);
            eff_hyp = hyp_eff(M,gamma,Nf,Ns,hmin);
            eval(['effs_hyp_' num2str(igamma) '(iNs) = eff_hyp;'])
        end
        eval(['Nss_' num2str(igamma) ' = Nss;'])
        
    end
    
    save('./figures/fig_effs.mat')
    
end


fig1 = figure();
clf;
hold on
xlabel('N_s/N_f')
ylabel('\epsilon')

for igamma = 1:length(gammas)
    gamma = gammas(igamma);
    
    eval(['Nss = Nss_' num2str(igamma) ';'])
    eval(['effs_hyp = effs_hyp_' num2str(igamma) ';'])
    
    figure(fig1)
    plot(Nss/Nf,effs_hyp,'k-','linewidth',1.5)
    axis([-2*M M 0 round(100*1.1*max(gammas)/M)/100])
    axis([min(Nss)/Nf max(Nss)/Nf 0 1.1*max(effs_hyp)])
    
    if M<sqrt(2/gamma-1)
        Nss = Nf*[Ns56 Ns45 Ns34 Ns23 Ns12 Ns01];
    else
        Nss = Nf*[Ns56 Ns34 Ns23 Ns12 Ns01];
    end
    
    effs = [];
    for Ns_i = 1:length(Nss)
        Ns = Nss(Ns_i);
        eff = hyp_eff(M,gamma,Nf,Ns,hmin);
        effs = [effs eff];
    end
    
    figure(fig1)
    plot(Nss/Nf,effs,'.','MarkerSize',14,'Color',[102,102,102]/255)
    drawnow()
end

figure(fig1)
axis([-30,5,0,0.03])
set(gca,'fontsize',16)
set(gcf,'Paperunits','centimeters')
set(gcf,'Paperposition',[1 1 25 15])
print -depsc './figures/fig_effs.eps'
set(gcf,'Paperunits','inches')