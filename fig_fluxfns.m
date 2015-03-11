% ---------------------------------------------------------
%   Plot multiple flux functions
% ---------------------------------------------------------

clear all; close all; clc;
addpath('./hyp_toolbox') % hyp_fluxFns is in the toolbox, not normally for direct use
tol = 1E-14; % Tolerance for comparison of floats

hs = [0.,logspace(-2.,0.,1000)];

Ms = [2 5 10 20 50]
Nf = 1
Ns = 0

fig1 = figure();
xlabel('\eta')
ylabel('F')
axis([0 1 0 1])
set(gca,'xtick',[0,0.2,0.4,0.6,0.8,1.0])
set(gca,'xticklabel',{'0','0.2','0.4','0.6','0.8','1'})
set(gca,'ytick',[0,0.25,0.5,0.75,1.0])
set(gca,'yticklabel',{'-0.0','','','','1'})
hold on
for i=1:length(Ms)
    M = Ms(i);
    Fs = hyp_fluxFn(hs,M,Nf,Ns);
    plot(hs,Fs,'k-','linewidth',1.5)
end
set(gca,'fontsize',16)
set(gcf,'Paperunits','centimeters')
set(gcf,'Paperposition',[1 1 12.5 7.5])
print -depsc './figures/fig_fluxfns_a'
set(gcf,'Paperunits','inches')


Ms = [2 5 10 20 50]
Nf = 0
Ns = 1

fig2 = figure();
xlabel('\eta')
ylabel('F')
axis([0 1 0 0.2])
set(gca,'xtick',[0,0.2,0.4,0.6,0.8,1.0])
set(gca,'xticklabel',{'0','0.2','0.4','0.6','0.8','1'})
set(gca,'ytick',[0,0.05,0.1,0.15,0.2])
set(gca,'yticklabel',{'-0.0','','','','0.2'})
hold on
for i=1:length(Ms)
    M = Ms(i);
    Fs = hyp_fluxFn(hs,M,Nf,Ns);
    plot(hs,Fs,'k-','linewidth',1.5)
end
set(gca,'fontsize',16)
set(gcf,'Paperunits','centimeters')
set(gcf,'Paperposition',[1 1 12.5 7.5])
print -depsc './figures/fig_fluxfns_b'
set(gcf,'Paperunits','inches')


M = 2
Nf = 1
Nss = [-5:1:4];

fig3 = figure();
xlabel('\eta')
ylabel('F')
axis([0 1 -0.5 1.5])
set(gca,'xtick',[0,0.2,0.4,0.6,0.8,1.0])
set(gca,'xticklabel',{'0','0.2','0.4','0.6','0.8','1'})
set(gca,'ytick',[-0.5,0,0.5,1,1.5])
set(gca,'yticklabel',{'-0.5','0','0.5','1','1.5'})
hold on
for i=1:length(Nss)
    Ns = Nss(i);
    Fs = hyp_fluxFn(hs,M,Nf,Ns);
    plot(hs,Fs,'k-','linewidth',1.5)
end
set(gca,'fontsize',16)
set(gcf,'Paperunits','centimeters')
set(gcf,'Paperposition',[1 1 12.5 7.5])
print -depsc './figures/fig_fluxfns_c'
set(gcf,'Paperunits','inches')