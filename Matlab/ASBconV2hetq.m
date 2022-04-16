%ASBbl3convperV2fig5_1_2_6.m
%Aug, 2021 Baseline, convergence, two bin, perisitence, single R model  

clear

%parameters 

N=100;      %population of agents
NN=N+1;
%M=20;
MM=200;

Block=50000;

%T=MM*Block;     %Time/iteration index for the sim

Rstar=0;      %resources devoted to apprehension
F=1;        %individual cost of apprehnsion given ASA
gam=0.8;
aa=1;
bb=0.25; %0.1765;
mu=0.6;     %mean value of g, individual benefit from ASA
sig=0.2;    %varance of g 
lam=5;     %socail cost conversion of individual ASA for social damage function
rho=2;
Z=2;        %z-history length

BinE=50;

eps=4;

critconv1=0.01;
RR=zeros(NN,1);
mv=zeros(NN,1);
mov=zeros(NN,1);
sdv=zeros(NN,1);
mq=zeros(NN,1);
mA=zeros(NN,1);
coefv=zeros(NN,1);
coefvw=zeros(NN,1);
ma=zeros(NN,1);
sda=zeros(NN,1);
gpiv=zeros(NN,1);
vcon=zeros(NN,Block);OBP=zeros(NN,Block);SBP=zeros(NN,Block);
cost=zeros(NN,1);
cost2=zeros(NN,1);
mgb=zeros(NN,1);
ARATE=zeros(NN,1);
ARATER=NaN(NN,1);
IRATE=zeros(NN,1);
IRATER=NaN(NN,1);
EXCAP=NaN(NN,1);
BH=zeros(NN,1);

edges=zeros(NN+1,1);
for j=1:NN+1
    edges(j)=j-1.5;
end
%vi=normrnd(0,sig,Block,N);

for r=1:NN
    RR(r)=r-1;
   
    c=1;
    b=0;
    cc=0;
    TESTCONV=zeros(MM,1);
    TESTCONV(1)=1;
    while ((cc<1) && (b<MM))
        g=normrnd(mu,sig,Block,N); % unifrnd(0,1,Block,N); %
        vi=unifrnd(-sig,sig,Block,N);
        swi=zeros(Block,1);
        swic=zeros(Block,1);
        vz=zeros(Block,1);
        az=zeros(Block,1);
        q=zeros(Block,1);
        R=zeros(Block,1);
        gbar=zeros(Block,N);
        gbar3=zeros(Block,1);
        A=zeros(Block,1);
        NCw=zeros(Block,1);NCw2=zeros(Block,1);
        vw=zeros(Block,1);aw=zeros(Block,1);
        %vi=zero(N,1);
        ai=zeros(N,1);
        b=b+1;
        BH(r)=b;
        if b<2
            in=Z+1;
        else
            in=1;
        end
        for t=1:Block
            R(t)=RR(r);
            if (b<2)
                if (t<Z+1)
                    for z=1:Z
                        vz(z)=Z*round(unifrnd(0,N));
                        az(z)=Z*round(unifrnd(0,vz(z)));
                    end
                else
                    vz(t)=sum(vw(t-Z:t-1));
                    az(t)=sum(aw(t-Z:t-1));
                end
            else
                if t<1+Z
                    if t<2
                        vz(t)=sum(v(t+(b-1)*Block-Z:t+(b-1)*Block-1));
                        az(t)=sum(a(t+(b-1)*Block-Z:t+(b-1)*Block-1));
                    else
                        vz(t)=sum(v(t+(b-1)*Block-Z-1:(b-1)*Block))+sum(vw(t-(Z-(Z-(t-1))):t-1));
                        az(t)=sum(a(t+(b-1)*Block-Z-1:(b-1)*Block))+sum(aw(t-(Z-(Z-(t-1))):t-1));  
                    end
                else
                    vz(t)=sum(vw(t-Z:t-1));
                    az(t)=sum(aw(t-Z:t-1));
                end
            end
            q(t)=(aa+az(t))/(aa+bb+vz(t));            
            for n=1:N                
                %vi(n)=unifrnd(-sig,sig);
                if (q(t)+vi(t,n))>1
                    vi(t,n)=1-q(t);
                end
                if q(t)+vi(t,n)<0
                    vi(t,n)=0-q(t);
                end
                if (q(t)+vi(t,n))*F<=g(t,n) 
                    vw(t)=vw(t)+1;
                    gbar(t,n)=g(t,n);
                end
            end
            gbar3(t)=sum(gbar(t,:));
            A(t)=gam*min(1,R(t)/vw(t));
            %A(t)=gam*(1-1/(eps^(R(t)/vw(t))));
            aw(t)=binornd(vw(t),A(t));
            NCw(t)=rho*R(t)+(lam-1)*gbar3(t);
        end
        if b<2
            v=vw;
            a=aw;
            gb=gbar3;
            NC=NCw;
            Ap=A;
            qp=q;
        else
            vhold=cat(2,v',vw');
            ahold=cat(2,a',aw');
            gbhold=cat(2,gb',gbar3');
            NChold=cat(2,NC',NCw');
            Aphold=cat(2,Ap',A');
            qphold=cat(2,qp',q');
            v=vhold';
            a=ahold';
            gb=gbhold';
            NC=NChold';
            Ap=Aphold';
            qp=qphold';
        end

        vconv=v(1:b*Block);
        vconvlag=v(1:(b-1)*Block);
        freqvconv=histcounts(vconv(:),edges)/((b)*Block);
        freqvconvlag=histcounts(vconvlag(:),edges)/((b-1)*Block);
        TESTCONV1=zeros(1,NN);
        if b>1
            for j=1:NN
                TESTCONV1(j)=abs(freqvconv(j)-freqvconvlag(j));
            end
            TESTCONV(b)=sum(TESTCONV1);                
        end
        if b > 4
            if (TESTCONV(b)<=critconv1) && (TESTCONV(b-1)<=critconv1) ...
                && (TESTCONV(b-2)<=critconv1) && (TESTCONV(b-3)<=critconv1)...
                && (TESTCONV(b-4)<=critconv1)
                c=b;
                cc=1;
            end
        end
    end

    mv(r)=mean(v);
    mov(r)=mode(v);
    sdv(r)=std(v);
    ma(r)=mean(a);
    sda(r)=std(a);
    mq(r)=mean(qp);
    mA(r)=mean(Ap);
    coefv(r)=sdv(r)/mv(r);
    mgb(r)=mean(gb);
    cost(r)=rho*RR(r)+(lam-1)*mgb(r); 
    cost2(r)=mean(NC);
%    OBP(r,:)=Ap((c-2)*Block+1:(c-1)*Block);
%    SBP(r,:)=qp((c-2)*Block+1:(c-1)*Block);
    
    ARATE(r)=ma(r)/mv(r);
   
    IRATE(r)=min(RR(r),mv(r))/mv(r);
    
    if RR(r)>0
        EXCAP(r)=max(0,(RR(r)-mv(r))/RR(r));
        ARATER(r)=ARATE(r)/RR(r);
        IRATER(r)=IRATE(r)/RR(r);
    end
    
    RR(r)

end

[M, I]=min(cost);
Rst=RR(I);


%figure % Figure 5.2
%tile=tiledlayout(2,1);
%tile.Padding='none';
%tile.TileSpacing='none';
%nexttile
%hold on
%plot(RR,mv,'Color','k','LineStyle','-','LineWidth',0.5)
%plot(RR,mov,'Color','k','LineStyle','--','LineWidth',0.5) 
%box on
%ylabel('Number of violations')
%title('Panel 1: Expected violations in the stationary distribution, E(v|R)')
%legend('Mean violatoins','Modal violations')
%hold off
%nexttile
%hold on
%plot(RR, coefv,'Color','k','LineStyle','-','LineWidth',0.5)
%box on
%title('Panel 2: Coefficient of variation of violations in the statoinary distribution')
%ylabel('\sigma / \mu')
%hold off
%nexttile
%hold on
%plot(RR,gpiv,'Color','k','LineStyle','-','LineWidth',0.5)
%box on
%xlabel('Deterrence resources (R)')
%ylabel('{\phi}(g(R))')
%title('Panel 3: Probability density of {\phi}(g) as a function of R') 
%hold off

%figure %figure 5.6
%tile=tiledlayout(2,1);
%tile.Padding='none';
%tile.TileSpacing='none';
%nexttile
%hold on
%plot(RR,EXCAP,'Color','k','LineStyle','-','LineWidth',0.5)
%plot(RR,IRATE,'Color','k','LineStyle','--','LineWidth',0.5)
%plot(ARATE,'Color','k','LineStyle',':','LineWidth',0.5)
%ylim([0 1.1])
%xlim([0 100])
%box on
%title('Panel 1: Reserve capacity, investigations and apprehensions per unit R')
%legend('Reserve capacity rate','Investigation rate','Apprehension rate','location','northwest')
%hold off
%nexttile
%hold on
%plot(RR,ma,'Color','k','LineStyle','-','LineWidth',0.5)
%plot(RR,(N-mv),'Color','k','LineStyle','--','LineWidth',0.5)
%xlim([0 100])
%box on
%title('Panel 2: Apprehensions and compliance')
%xlabel('Deterrence resources (R)')
%legend('Apprehensions','Compliance','location','northwest')
%hold off

%figure %figure 5.1
%tile=tiledlayout(6,4);
%tile.Padding='none';
%tile.TileSpacing='none';
%nexttile([1 3])
%plot(vcon(38,Block-10000:Block-5000),'Color','k')
%ylim([0 100])
%xlim([0 5000])
%title('Panel 1: R = 37')
%nexttile
%histogram(vcon(38,Block-10000:Block-5000),'Normalization','probability','FaceColor',[0.17 0.17 0.17],'BinWidth',1)
%xlim([0 100])
%ylim([0 0.15])
%xlabel('R = 37')
%view(90,-90)
%nexttile([1 3])
%plot(vcon(40,Block-10000:Block-5000),'Color','k')
%title('Panel 2: R = 39')
%ylim([0 100])
%xlim([0 5000])
%nexttile
%histogram(vcon(40,Block-10000:Block-5000),'Normalization','probability','FaceColor',[0.17 0.17 0.17],'BinWidth',1)
%xlim([0 100])
%ylim([0 0.15])
%xlabel('R = 39')
%view(90,-90)
%nexttile([1 3])
%plot(vcon(41,Block-10000:Block-5000),'Color','k')
%title('Panel 2: R = 40')
%ylim([0 100])
%xlim([0 5000])
%nexttile
%histogram(vcon(41,Block-10000:Block-5000),'Normalization','probability','FaceColor',[0.17 0.17 0.17],'BinWidth',1)
%xlim([0 100])
%ylim([0 0.15])
%xlabel('R = 40')
%view(90,-90)
%nexttile([1 3])
%plot(vcon(42,Block-10000:Block-5000),'Color','k')
%title('Panel 3: R = 41')
%ylim([0 100])
%xlim([0 5000])
%nexttile
%histogram(vcon(42,Block-10000:Block-5000),'Normalization','probability','FaceColor',[0.17 0.17 0.17],'BinWidth',1)
%xlim([0 100])
%ylim([0 0.15])
%xlabel('R = 41')
%view(90,-90)
%nexttile([1 3])
%plot(vcon(43,Block-10000:Block-5000),'Color','k')
%title('Panel 4: R = 42')
%ylim([0 100])
%xlim([0 5000])
%nexttile
%histogram(vcon(43,Block-10000:Block-5000),'Normalization','probability','FaceColor',[0.17 0.17 0.17],'BinWidth',1)
%xlim([0 100])
%ylim([0 0.15])
%xlabel('R = 42')
%view(90,-90)
%nexttile([1 3])
%plot(vcon(47,Block-10000:Block-5000),'Color','k')
%title('Panel 5: R = 46')
%xlabel('Tic (t=0,...,5000)')
%ylim([0 100])
%xlim([0 5000])
%nexttile
%histogram(vcon(47,Block-10000:Block-5000),'Normalization','probability','FaceColor',[0.17 0.17 0.17],'BinWidth',1)
%xlim([0 100])
%ylim([0 0.15])
%xlabel('R = 46')
%ylabel('frequency')
%view(90,-90)
%nexttile([1 4])
%hold on
%box on
%plot(OBP(42,Block-10000:Block-5000),'Color','k')
%plot(SBP(42,Block-10000:Block-5000),'Color','k','Linestyle',':','LineWidth',2.0)
%legend('Objective probability','Subjective probability')
%xlabel('Tic (t=0,...,5000)')
%hold off

%output for figure 3.6 input to CovCostPar.m
save C:\Users\kcarlaw\Documents\MATLAB\ASB2022final\mvhetq mv -ASCII -DOUBLE;



%pmv=mv/N;
%pRR=RR/N;
%output for figure 5.3 from baseline
%save C:\Users\kcarlaw\Documents\MATLAB\ASB2022final\movbl mov -ASCII -DOUBLE;
%save C:\Users\kcarlaw\Documents\MATLAB\ASB2022final\pmvbl pmv -ASCII -DOUBLE;
%save C:\Users\kcarlaw\Documents\MATLAB\ASB2022final\pRR pRR -ASCII -DOUBLE;

%save C:\Users\kcarlaw\Documents\MATLAB\ASB2022final\Rcon103 RR -ASCII -DOUBLE;
%save C:\Users\kcarlaw\Documents\MATLAB\ASB2022final\Gbar103 mgb -ASCII -DOUBLE;
%save C:\Users\kcarlaw\Documents\MATLAB\ASB2022final\mabl103 ma -ASCII -DOUBLE;

%save C:\Users\kcarlaw\Documents\MATLAB\ASB2022final\TGB Tc -ASCII -DOUBLE;
%save C:\Users\kcarlaw\Documents\MATLAB\ASB2022final\TBB Tu -ASCII -DOUBLE;

%save C:\Users\kcarlaw\Documents\MATLAB\ASB2022final\EDGB3 elsc -ASCII -DOUBLE;
%save C:\Users\kcarlaw\Documents\MATLAB\ASB2022final\EDBB3 elsu -ASCII -DOUBLE;

