clear all; close all;

M = 2;
gamma = 0.5;
Nf = 1;
Ns = -0.75;
hmin = 0;

[xCs,tCs,eff] = hyp_chars(M,gamma,Nf,Ns,hmin)
set(gca,'fontsize',16)
set(gcf,'Paperunits','centimeters')
set(gcf,'Paperposition',[1 1 16 10])
print -depsc './figures/fig_chars_migr'
set(gcf,'Paperunits','inches')

N = 10000;
hs = [0.,logspace(-6.,0.,N)];
ts = [1,1.571,2.571,3.672,4.344];
[xss,hss,xmss,hmss] = hyp_plume(M,gamma,Nf,Ns,ts,hs,hmin);

xLmost = min([-xCs(1),xCs]);
xRmost = max(xCs);
for ti = 1:length(ts)
    
    figure(2);
    clf;
    set(gca,'fontsize',16)
    xlabel('\xi')
    ylabel('1-\eta')
    fill(1.1*[xLmost,xLmost,xRmost,xRmost],[1,0,0,1],[1,1,1],'linewidth',1.5);
    axis([1.1*xLmost,1.1*xRmost,0,1])
    hold on
    
    xs = xss(ti,:);
    hs = hss(ti,:);
    xms = xmss(ti,:);
    hms = hmss(ti,:);
    
    fill([xms(1),xms,xms(end),xms(end),xms(1)],[1-hms(1),1-hms,1-hms(end),1,1],[200,200,200]/255,'linewidth',1.5)
    plot(xms,1-hms,'k-','linewidth',1.5)
    fill([xs(1),xs,xs(end),xs(end),xs(1)],[1-hs(1),1-hs,1-hs(end),1,1],[121,121,121]/255,'linewidth',1.5)
    plot([0,0],[0,1],'k-')
    box on;
    if ti>1
        axis off
    end
    drawnow()
    
    set(gcf,'Paperunits','centimeters')
    set(gcf,'Paperposition',[1 1 16 2])
    eval(['print -depsc ''./figures/fig_migr_plume_' num2str(ti) ''''])
    
end