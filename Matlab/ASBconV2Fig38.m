% ASBconV2Fig38 costs 
% Figure 3.8

clear
N=100;
NN=N+1;
X=10000;
F=1;
%lam=2;
PR=2;
PR2=0;
chi1=1;

%load acn103;
load Rcon103;
load Gbar103;
load mabl103;

RR=Rcon103;
%ma=acn103;
mgbar2=Gbar103;
a=mabl103;

cost=zeros(NN,1);
cost2=zeros(NN,1);
cost3=zeros(NN,1);
ccost3=zeros(X,NN,1);
ccost=zeros(X,NN,1);
ccost2=zeros(X,NN,1);
llim=zeros(X,NN,1);
ulim=zeros(X,NN,1);
llim2=zeros(X,NN,1);
ulim2=zeros(X,NN,1);
llim3=zeros(X,NN,1);
ulim3=zeros(X,NN,1);
lam=zeros(X,1);
chi=zeros(X,1);
M=zeros(X,1);
I=zeros(X,1);
Rstar=zeros(X,1);


for i=1:X
    lam(i)=0.05*i;

    for r=1:NN
        cost(r)=PR*RR(r)+(lam(i)-1)*mgbar2(r);
        cost2(r)=PR2*RR(r)+(lam(i)-1)*mgbar2(r);
        ccost(i,r)=PR*RR(r)+(lam(i)-1)*mgbar2(r); %+ma(r)*F;
        ccost2(i,r)=PR2*RR(r)+(lam(i)-1)*mgbar2(r);
        llim(i,r)=cost(1);
        ulim(i,r)=cost(NN);
        llim2(i,r)=cost2(1);
        ulim2(i,r)=cost2(NN);
        
        
        
    end    
    [M(i), I(i)]=min(cost);

    Rstar(i)=RR(I(i));
end

figure %figure 3.8
tile=tiledlayout(1,4);
tile.Padding='none';
tile.TileSpacing='none';
nexttile
hold on
box on
plot(RR,ccost(40,:)','Color','k')
plot(RR,llim(40,:)','Color','k','LineStyle',':')
plot(RR,ulim(40,:)','Color','k','LineStyle','--')
title('Panel 1: very low cost of ASB (spitting on the sidewalk)')
ylim([llim(40,1)-10 ulim(40,1)+10])
%xlabel('Deterrence resources (R)')
ylabel('E(cost|R)')
legend('Expected costs of crime {\rho}=2 and {\lambda}=2','E(cost|R=0)','E(cost|R=N)','Location','Best')
legend boxoff
hold off
nexttile ([1 2])
hold on
box on
plot(RR,ccost(100,:)','Color','k')
plot(RR,llim(100,:)','Color','k','LineStyle',':')
plot(RR,ulim(100,:)','Color','k','LineStyle','--')
title('Panel 2: Intermediate cost of ASB (interesting case)')
xlabel('Deterrence resources (R)')
%ylabel('E(cost|R)')
legend('Expected costs of crime {\rho}=2 and {\lambda}=5','E(cost|R=0)','E(cost|R=N)','Location','Best')
legend boxoff
hold off
nexttile
hold on
box on
plot(RR,ccost(8000,:)','Color','k')
plot(RR,llim(8000,:)','Color','k','LineStyle',':')
plot(RR,ulim(8000,:)','Color','k','LineStyle','--')
ylim([5000 25000])
title('Panel 3: Very high cost of ASB (nuclear terrorism)')
%xlabel('Deterrence resources (R)')
%ylabel('E(cost|R)')
legend('Expected costs of crime {\rho}=2 and {\lambda}=400','E(cost|R=0)','E(cost|R=N)','Location','Best')
legend boxoff
hold off

for i=1:X
    lam(i)=5; %0.05*i;
    chi(i)=0.02*i;
    if chi(i)<1
        chi(i)=1;
    end

    for r=1:NN
        cost3(r)=PR*RR(r)+(lam(i)-1)*mgbar2(r)+a(r)*chi(i);
        ccost3(i,r)=PR*RR(r)+(lam(i)-1)*mgbar2(r)+a(r)*chi(i);
        llim3(i,r)=cost3(1);
        ulim3(i,r)=cost3(NN);
    end    
end


figure %figure 5.1 with apprehension costs
tile=tiledlayout(1,2);
tile.Padding='none';
tile.TileSpacing='none';
nexttile
hold on
box on
plot(RR,ccost3(40,:)','Color','k')
plot(RR,llim3(40,:)','Color','k','LineStyle',':')
plot(RR,ulim3(40,:)','Color','k','LineStyle','--')
title('Panel 1: \chi = 1')
%ylim([llim3(40,1)-10 ulim3(40,1)+10])
%xlabel('Deterrence resources (R)')
ylabel('E(cost|R)')
legend('Expected costs of crime {\rho}=2 and {\lambda}=5','E(cost|R=0)','E(cost|R=N)','Location','northeast')
hold off
%nexttile ([1 2])
%hold on
%box on
%plot(RR,ccost3(100,:)','Color','k')
%plot(RR,llim3(100,:)','Color','k','LineStyle',':')
%plot(RR,ulim3(100,:)','Color','k','LineStyle','--')
%title('Panel 2: \chi = 2')
%xlabel('Deterrence resources (R)')
%ylabel('E(cost|R)')
%ylim([llim3(100,1)-10 ulim3(100,1)+10])
%legend('Expected costs of crime {\rho}=2 and {\lambda}=5','E(cost|R=0)','E(cost|R=N)','Location','Best')
%hold off
nexttile
hold on
box on
plot(RR,ccost3(500,:)','Color','k')
plot(RR,llim3(500,:)','Color','k','LineStyle',':')
plot(RR,ulim3(500,:)','Color','k','LineStyle','--')
%ylim([llim3(8000,1)-10 ulim3(8000,1)+10])
title('Panel 3: \chi = 10')
%xlabel('Deterrence resources (R)')
%ylabel('E(cost|R)')
legend('Expected costs of crime {\rho}=2 and {\lambda}=5','E(cost|R=0)','E(cost|R=N)','Location','northeast')
hold off
