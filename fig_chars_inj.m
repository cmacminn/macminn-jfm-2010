% ---------------------------------------------------------
% Injection: generate plots of the profile and characteristics at various times
% ---------------------------------------------------------

clear all; close all; clc;
addpath('./hyp_toolbox') % hyp_inj_chars is in the toolbox, not normally for direct use
tol = 1E-14; % Tolerance for comparison of floats

M = 2;
gamma = 0.5;
Nf = 1;
Ns = 0;

% Characteristics
hs = [0.,logspace(-2.,0.,10)];
gs = (M-1)*hs+1;

xLmost = -1.25*M;
xRmost =  1.25*M;

fig1 = figure(1);
clf;
t = 1;
hyp_inj_chars(M,gamma,t,gs);
axis([xLmost,xRmost,0,2])
xlabel('\xi')
ylabel('\tau')
drawnow()
set(gca,'ytick',[0,0.5,1])
set(gca,'yticklabel',{'0','','1'})
set(gca,'fontsize',16)
set(gcf,'Paperunits','centimeters')
% set(gcf,'Paperposition',[1 1 25 15])
set(gcf,'Paperposition',[1 1 16 10]) % Papers
saveas(fig1,'./figures/fig_inj_chars.fig')
print -depsc './figures/fig_inj_chars.eps'
set(gcf,'Paperunits','inches')


% Plume shape
ts = [0.1,0.5,1];

hmin = 0;
N = 10000;
hs = [0.,logspace(-6.,0.,N)];
[xss,hss,xmss,hmss,ts] = hyp_plume(M,gamma,Nf,Ns,ts,hs,hmin);

for ti = 1:length(ts)
    fig2 = figure(2);
    clf
    set(gca,'fontsize',16)
    xlabel('\xi')
    ylabel('1-\eta')
    t = ts(ti)
    xs = xss(ti,:);
    hs = hss(ti,:);
    fill([xLmost,xLmost,xRmost,xRmost],[1,0,0,1],[255,255,255]/255);
    hold on
    fill([xs(1),xs,xs(end),xs(1)],[1,1-hs,1,1],[121,121,121]/255,'linewidth',1.5)
    plot([0,0],[0,1],'k-')
    axis([xLmost,xRmost,0,1])
    drawnow()
    set(gcf,'Paperunits','centimeters')
    set(gcf,'Paperposition',[1 1 25 3])
    saveas(fig1,'./figures/fig_inj_chars.fig')
    eval(['print -depsc ''./figures/fig_inj_plume_' num2str(ti) ''''])
    eval(['saveas(fig2,''./figures/fig_inj_plume_' num2str(ti) '.fig'')'])
    
end