clear all; close all;

M = 3;
gamma = 0;
Nf = 0;
Ns = 0.5;
Ng = 0.3;
Nd = 8E-3;
gammad = 0.3;
Swc = 0;
xRes = 40;

hmin = 1E-3;

[xCs,tCs] = hyp_crits(M,gamma,Nf,Ns,hmin);

N = 10000;
hs = [0.,logspace(-6.,0.,N)];

% Get the plume shape at a selection of times over the whole lifetime, and make a movie
ts = linspace(1,36,100);
% ts = linspace(tCs(1),tCs(end),100);
[xss,hss,xmss,hmss] = hyp_plume(M,gamma,Nf,Ns,ts,hs);
% [xss,hss,hmss,hdss,ts,Vs,Vxss] = nmr_plume(M,gamma,Nf,Ns,Ng,Nd,gammad,Swc,ts,hmin,xRes);

% xLmost = min([-xCs(1),xCs]);
% xRmost = max(xCs);
xLmost = -6;
xRmost = 12;

fig1 = figure();
for ti = 1:length(ts)
    
    disp(['Plotting the plume at time t = ' num2str(ts(ti)) ' (' num2str(ti) '/' num2str(length(ts)) ')'])
    
    figure(fig1)
    hold off
    set(gca,'fontsize',16)
    xlabel('\xi')
    ylabel('1-\eta')
    fill([xLmost,xLmost,xRmost,xRmost],[1,0,0,1],[1,1,1],'linewidth',1);
    axis([1.01*xLmost,1.01*xRmost,0,1.1])
    hold on
    
    if length(xss(:,1))>1
        xs = xss(ti,:);
        xms = xmss(ti,:);
    else
        xs = xss;
        xms = xmss;
    end
    hs = hss(ti,:);
    hms = hmss(ti,:);
    
    hts = hs+hms;
    xs = xs(hts>1E-5);
    hs = hs(hts>1E-5);
    xms = xms(hts>1E-5);
    hms = hms(hts>1E-5);
    if gamma>0
        fill([xms(1),xms,xms(end),xms(end),xms(1)],[1-hms(1),1-hms,1-hms(end),1,1],[200,200,200]/255)
        plot(xms,1-hms,'k-','linewidth')
    end
    fill([xs(1),xs,xs(end),xs(end),xs(1)],[1-hs(1),1-hs,1-hs(end),1,1],[121,121,121]/255)
    plot([0,0],[0,1],'k-')
    box on;
    axis off
    drawnow
    
    set(gcf,'PaperUnits','centimeters')
    set(gcf,'PaperSize',[16 2])
    set(gcf,'PaperPosition',[0 0 16 2])
    print(gcf,['figures/fig_case2_movie_' num2str(ti)],'-djpeg','-r300')
    
end