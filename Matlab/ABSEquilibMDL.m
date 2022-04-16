clear

N=100;
NN=N+1;
R=NN;
gam=0.8;
mu=0.6;
sig=0.2;
aa=1;
bb=0.25; %0.1765;
F=1;
lam=5;
rho=2;

A=zeros(NN,1);
AA=zeros(NN,1);
RR=zeros(R,1);
ProbEq=zeros(R,NN,1);
%NProbEq=zeros(R,NN,1);
MGB=zeros(NN,1);SDGB=zeros(NN,1);MBB=zeros(NN,1);SDBB=zeros(NN,1);
sdGB=zeros(R,NN,1);sdBB=zeros(R,NN,1);
GProbEq=NaN(R,NN,1);
BProbEq=NaN(R,NN,1);
TProbEq=NaN(R,NN,1);
SProbEq=zeros(R,NN,1);
NEq=zeros(R,1);
GNEq=zeros(R,1);
BNEq=zeros(R,1);
TNEq=zeros(R,1);
v=zeros(NN,1);


EEQv=zeros(R,NN,1);
EvP=zeros(R,NN,1);EvPP=zeros(NN,1);
EEvP=zeros(R,NN,1);
BBEQv=zeros(R,1);
GBEQv=zeros(R,1);
EQmv=zeros(R,1);DEQv=zeros(R,1);mvvbl=zeros(R,1);
pBB=zeros(R,1);pGB=zeros(R,1);gBB=zeros(R,1);gGB=zeros(R,1);
costBB=zeros(R,1);costGB=zeros(R,1);
tBB=zeros(R,1);tGB=zeros(R,1);
BinE=NaN(R,1);TF=zeros(R,1);

load mvbl;
load movbl;
load TGB;
load TBB;

PD=makedist('Normal','mu',mu,'sigma',sig);

for r=1:R
    RR(r)=r-1;
    %EQv=zeros(NN,1);
    for n=1:NN
        v(n)=n-1;
        A(n)=gam*min(1,(RR(r)/v(n)));
        AA(n)=gam*min(1,(RR(r)/(v(n)+1)));
        ProbEq(r,n)=(factorial(N)/(factorial(v(n))*factorial(N-v(n))))*(1-normcdf(A(n)*F,mu,sig))^v(n)*(normcdf(AA(n)*F,mu,sig))^(N-v(n));
        if v(n)<=31
            GProbEq(r,n)=ProbEq(r,n);
        end
        if v(n)>=70
            BProbEq(r,n)=ProbEq(r,n);
        end
        if 31<v(n)<70
            TProbEq(r,n)=ProbEq(r,n);
        end
        EvP(r,n)=abs(v(n)-N*(1-normcdf(F*gam*min(1,RR(r)/v(n)),mu,sig)));
        EvPP(n)=v(n)-N*(1-normcdf(F*gam*min(1,RR(r)/v(n)),mu,sig));
        EEvP(r,n)=N*(1-normcdf(F*gam*min(1,RR(r)/v(n)),mu,sig));
        %EP(n)=N*(1-normcdf(F*gam*min(1,RR(r)/v(n)),mu,sig));
        %if EvP(r,n)<0.5
        %    EQv(n)=N*(1-normcdf(F*gam*min(1,RR(r)/v(n)),mu,sig)); %v(n);
        %end
    end
    EQv=find(EvPP(1:NN-1)<0 & EvPP(2:NN)>0);
    BBEQv(r)=EEvP(r,max(EQv)-1);
    GBEQv(r)=EEvP(r,min(EQv)-1);
    TF(r)=GBEQv(r)-BBEQv(r);
    
    
    %TFl(r)=isempty(find(EvPP(1:NN-1)>0 & EvPP(2:NN)<0));
    %TFh(r)
    if TF(r)<0
        if GBEQv(r)<GBEQv(r-1)
            low=r-1;
        end
        hi=r-1;
        BinE(r)=find(EvPP(1:NN-1)>0 & EvPP(2:NN)<0);
    end
end
for r=1:R

    if RR(r)<low
        GBEQv(r)=0;
    end
    if RR(r)>hi
        BBEQv(r)=0;
    end
    
    

    
    EQmv(r)=TGB(r)*GBEQv(r)+TBB(r)*BBEQv(r);
    mvvbl(r)=mvbl(r);
    DEQv(r)=abs(mvvbl(r)-EQmv(r));
    NEq(r)=sum(ProbEq(r,:));
    GNEq(r)=sum(ProbEq(r,1:31));
    BNEq(r)=sum(ProbEq(r,71:101));
    TNEq(r)=sum(ProbEq(r,32:70));
    SProbEq(r)=sum(ProbEq(r,:));
    pBB(r)=gam*min(1,RR(r)/BBEQv(r));
    pGB(r)=gam*min(1,RR(r)/GBEQv(r));
    tBB(r)=mean(truncate(PD,F*pBB(r),inf));
    tGB(r)=mean(truncate(PD,F*pGB(r),inf));
    gBB(r)=mean(truncate(PD,F*pBB(r),inf))*(1-normcdf(F*pBB(r),mu,sig))*BBEQv(r);
    gGB(r)=mean(truncate(PD,F*pGB(r),inf))*(1-normcdf(F*pGB(r),mu,sig))*GBEQv(r);
    if RR(r)<low
        gGB(r)=NaN;
    end
    if RR(r)>hi
        gBB(r)=NaN;
    end
    costBB(r)=rho*RR(r)+(lam-1)*gBB(r);
    costGB(r)=rho*RR(r)+(lam-1)*gGB(r);
end
MDEQv=mean(DEQv);
SDDEQv=std(DEQv);
[MaxDEQv,I]=max(DEQv);
NProbEq=zeros(R,NN,1);
GBNP=zeros(R,NN,1);
BBNP=zeros(R,NN,1);
vv=zeros(R,NN,1);
Ev=zeros(R,1);sdEv=zeros(R,1);mEv=zeros(R,1);coefv=zeros(R,1);
test=zeros(R,NN,1);test2=zeros(R,NN,1);
zer=zeros(R,1);
for r=1:R
    for n=1:NN
        NProbEq(r,n)=ProbEq(r,n)/SProbEq(r);
        vv(r,n)=v(n)*NProbEq(r,n);
        Ev(r)=Ev(r)+vv(r,n);
        test(r,n)=ProbEq(r,n)/sum(ProbEq(r,1:50));
        test2(r,n)=ProbEq(r,n)/sum(ProbEq(r,51:101));
        GBNP(r,n)=ProbEq(r,n)/sum(ProbEq(r,1:50))*v(n);
        BBNP(r,n)=ProbEq(r,n)/sum(ProbEq(r,51:101))*v(n);
    end
    sdEv(r)=std(vv(r));
    mEv(r)=mean(vv(r));
    coefv(r)=sdEv(r)/mEv(r);
    MGB(r)=sum(GBNP(r,1:50)); %mean(EEvP(r,1:50)); %sum(NProbEq(r,1:50));
    %SDGB(r)=std(GBNP(r,1:50)); %std(EEvP(r,1:50)); %std(NProbEq(r,1:50));
    MBB(r)=sum(BBNP(r,51:101)); %mean(EEvP(r,51:101)); %sum(NProbEq(r,51:101));
    %SDBB(r)=std(BBNP(r,51:101)); %std(EEvP(r,51:101)); %std(NProbEq(r,51:101));
    for n=1:NN
        sdGB(r,n)=((v(n)-MGB(r))^2*test(r,n));
        sdBB(r,n)=((v(n)-MBB(r))^2*test2(r,n));
    end
    SDGB(r)=sqrt((sum(sdGB(r,1:50))));
    SDBB(r)=sqrt((sum(sdBB(r,51:101))));
end
clifeqdata=table(RR(38:45),MGB(38:45),SDGB(38:45),MBB(38:45),SDBB(38:45))


figure % Figure3.1
tile=tiledlayout(4,1);
tile.Padding='none';
tile.TileSpacing='none';
nexttile ([2 1])
hold on
box on
plot(v,ProbEq(39,:),'Color','k');
area(v,TProbEq(39,:),'FaceColor',[0.6 0.6 0.6]);
area(v,GProbEq(39,:),'FaceColor',[0.9 0.9 0.9]);
area(v,BProbEq(39,:),'FaceColor',[0.1 0.1 0.1]);
title('Panel 1: Probability of an equilibrium with v violations, R = 39')
xlabel('Number of violations (v)')
txt={'{\leftarrow} Probability of good eq (v \leq 30)'};
text(18,0.10,txt);
txt1={'{\downarrow} Probability of other eq (30 < v < 70)'};
text(52,0.015,txt1);
txt2={'Probability of bad eq (v \geq 70) {\rightarrow}'} ;
text(75,0.07,txt2);
ylim([0 0.15])
hold off
nexttile
hold on
box on
plot(RR,NEq,'Color','k','LineWidth',1.0)
title('Panel 2: Expected number of equilibria as a function of R')
plot(RR,GNEq,'Color',[0.2 0.2 0.2],'LineStyle','--','LineWidth',1.0) 
plot(RR,BNEq,'Color', [0.6 0.6 0.6],'LineStyle','-.','LineWidth',1.5)
plot(RR,TNEq,'Color', [0.4 0.4 0.4],'LineStyle',':','LineWidth',2.0)
%xlabel('Deterrence resources (R)')
legend('Expected total # of equlibria','Expected # of good equilibria','Expected # of bad equilibria','Expected # of other equilibria')
ylim([0 2.5])
xlim([0 55])
hold off
nexttile
hold on
box on
plot(RR,Ev,'Color','k')
legend('Expected number of total equilibria')
title('Panel 3: Expected violations for a randomly chosen equilibrium')
xlabel('Deterrence resources (R)')
xlim([0 55])
ylim([0 110])
hold off

figure % figure 3.2
hold on
box on
plot(v,v,'Color','k')
plot(v,EEvP(11,:),'Color','k','LineStyle','--','LineWidth',1.0)
plot(v,EEvP(39,:),'Color','k','LineStyle','-.','LineWidth',1.5)
plot(v,EEvP(51,:),'Color','k','LineStyle',':','LineWidth',2.0)
xlabel('Number of violations(v_{t})') 
ylabel('E(v_{t+1}|v_{t})')
legend('E(v_{t+1}|v_{t},R = 10','E(v_{t+1}|v_{t},R = 35','E(v_{t+1}|v_{t},R = 50','location','best')
hold off

figure % figure 3.3
hold on
box on
plot(RR,costBB,'Color','k')
plot(RR,costGB,'Color','k','LineStyle','--','LineWidth',1.0)
legend('Expected cost bad equilibrium','Expected cost good equilibrium') 
xlabel('Deterrence resources (R)')
hold off

figure
hold on
box on
plot(RR,mvvbl,'Color','k')
plot(RR,EQmv,'Color','k','LineStyle','--','LineWidth',1.0)
legend('Mean v given R, convergence sim','Mean v given P(R,v), equilibrium model') 
xlabel('Deterrence resources (R)')
ylabel('Expected violations')
hold off

figure % Figure4.5 the cliff
tile=tiledlayout(3,1);
tile.Padding='none';
tile.TileSpacing='none';
nexttile ([2 1])
hold on 
box on
plot(RR, mvbl,'Color','k')
plot(RR,Ev,'Color','k','LineStyle',':','LineWidth',2.0)
title('Panel 1: Expected violations, dynamic and equilibrium modls')
legend('E(v|R), simulated from dynamic model','E(v|R), from PDF of equilibrium model')
ylim([10 102])
hold off
nexttile ([1 1])
hold on
box on
plot(RR, mvbl-movbl,'Color','k')
plot(RR,zer,'Color','k','LineStyle',':')
title('Panel 2: Difference between mean and modal violations')
xlabel('Deterrence resources (R)')
legend('\Delta = mean-mode','Location','Best')
ylim([-30 22])
hold off

%save C:\Users\kcarlaw\Documents\MATLAB\ASBFINALcodebase\PrEq2 ProbEq -ASCII -DOUBLE;
%save C:\Users\kcarlaw\Documents\MATLAB\ASBFINALcodebase\TPrEq2 TProbEq -ASCII -DOUBLE;
%save C:\Users\kcarlaw\Documents\MATLAB\ASBFINALcodebase\GPrEq2 GProbEq -ASCII -DOUBLE;
%save C:\Users\kcarlaw\Documents\MATLAB\ASBFINALcodebase\BPrEq2 BProbEq -ASCII -DOUBLE;
%save C:\Users\kcarlaw\Documents\MATLAB\ASBFINALcodebase\vEq2 v -ASCII -DOUBLE;