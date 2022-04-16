% Dynamic Theory of Deterrence and Compliance Fig 3.6
% Fig36.m used data generated from baseline convergence sim: ASBconblV2fig31_32.m, ABSconblV2N.m, ASBconV2Z.m, ASBconV2sig.m, ASBconV2ExpApr.m, ASBconV2hetq.m, ASBconV2unifg.m

clear

load mvbl;
load pmvbl;
load mvZ;
load pmvN;
load mvApr;
load mvsig;
load mvhetq
load mvunifg
N=100;
NN=N+1;
Rcp=zeros(NN,1);
PRcp=zeros(NN,1);
PPRcp=zeros(51,1);
for rcp=1:NN
    Rcp(rcp)=rcp-1;
    PRcp(rcp)=Rcp(rcp)/N;
end
for rcp=1:51
    PPRcp(rcp)=(rcp-1)/51;
end
figure %figure 3.6
tile=tiledlayout(2,3);
tile.Padding='none';
tile.TileSpacing='none';
nexttile
hold on
plot(PRcp,pmvbl,'Color','k','LineWidth',0.5)
plot(PPRcp,pmvN,'Color','k','LineStyle','--')
box on
title('Panel 1: Size of population ')
xlabel('Deterrence resources per agent (R/N)')
%legend('Benchmark N = 100','Comparison N = 50')
%legend boxoff
hold off
nexttile
hold on
plot(Rcp,mvbl,'Color','k','LineWidth',0.5)
plot(Rcp,mvZ,'Color','k','LineStyle','--')
box on
xlabel('Deterrence resources (R)')
title('Panel 2: Length of z-history')
%legend('Benchmark z = 2','Comparison z = 4')
%legend boxoff
ylim([0 N])
hold off
nexttile
hold on
plot(Rcp,mvbl,'Color','k')
plot(Rcp,mvsig,'Color','k','LineStyle','--')
box on
xlabel('Deterrence resources (R)')
title('Panel 3: Tighter distribution of g')
%legend('Benchmark {\sigma}=0.2','Comparison {\sigma}=0.15')
%legend boxoff
ylim([0 N])
hold off
nexttile
hold on
plot(Rcp,mvbl,'Color','k')
plot(Rcp,mvApr,'Color','k','LineStyle','--')
box on
xlabel('Deterrence resources (R)')
title('Panel 4: Different apprehnsion technology')
%legend('Benchmark min apprehension','Comparison exp apprehension')
%legend boxoff
ylim([0 N])
hold off
nexttile
hold on
plot(Rcp,mvbl,'Color','k')
plot(Rcp,mvhetq,'Color','k','LineStyle','--')
box on
xlabel('Deterrence resources (R)')
title('Panel 5: Hetergeneous q')
%legend('Benchmark common q','Comparison het q_{i} = q +/- uniform(-{\sigma},+{\sigma})')
%legend boxoff
ylim([0 N])
hold off
nexttile
hold on
plot(Rcp,mvbl,'Color','k')
plot(Rcp,mvunifg,'Color','k','LineStyle','--')
box on
xlabel('Deterrence resources (R)')
title('Panel 6: Uniform distribution of g')
%legend('Benchmark normal({\mu},{\sigma})','Comparison uniform[0,1.2]')
%legend boxoff
ylim([0 N])
hold off
