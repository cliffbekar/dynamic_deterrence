
clear

%parameters
%T=500000000;    %Time/iteration index for the sim
N=100;      %population of agents
NN=N+1;
%Rstar=41;      %fixed resources devoted to apprehension
F=1;        %individual cost of apprehnsion given ASA
gam=0.8;
aa=1;
bb=0.25; %0.1765;
%critical distribution of ASA (anti-scoial act) RN(gi;mu,sig,0,inf)
mu=.6;     %mean value of gi, individual benefit from ASA
sig=0.2;    %varance of gi 
lam=5;     %socail cost conversion of individual ASA for social damage function
PR=2;
Ri=[39,60]; % vectore of R values
RR=1;
Z=2;
Block=10000000;  %number of sim periods in convergence loop

BinE=53;

Evv=zeros(RR,NN,1);
U=zeros(RR,NN,1);
V=zeros(RR,NN,1);
L=zeros(RR,NN,1);
count=zeros(NN,1);
Evvb=zeros(RR,NN,1);
Ub=zeros(RR,NN,1);
Vb=zeros(RR,NN,1);
Lb=zeros(RR,NN,1);
countb=zeros(NN,1);

vp=zeros(NN,1);
for i=1:NN
    vp(i)=i-1;
end

for r=1:RR
    vv=zeros(NN,NN,1);
    vvb=zeros(NN,NN,1);
    Pvv=zeros(NN,NN,1);
    Pvvb=zeros(NN,NN,1);
    
%    for i=1:NN
%        vp(i)=i-1;
        
        q=zeros(Block,1);A=zeros(Block,1);
        vz=zeros(Block,1);az=zeros(Block,1);
        v=zeros(Block,1);a=zeros(Block,1);
        R=zeros(Block,1);
        qb=zeros(Block,1);Ab=zeros(Block,1);
        vzb=zeros(Block,1);azb=zeros(Block,1);
        vb=zeros(Block,1);ab=zeros(Block,1);        
        Rb=zeros(Block,1);
        for t=1:Block
            if (t<=Z)
                for z=1:Z
                    vz(z)=Z*unifrnd(0,N);
                    az(z)=Z*unifrnd(0,vz(z));
                    vzb(z)=Z*unifrnd(0,N);
                    azb(z)=Z*unifrnd(0,vzb(z));

                end
            else
                vz(t)=sum(v(t-Z:t-1));
                az(t)=sum(a(t-Z:t-1));
                vzb(t)=sum(vb(t-Z:t-1));
                azb(t)=sum(ab(t-Z:t-1));
            end
            Rb(t)=Ri(1);
            if t>1
                if v(t-1)<BinE
                    R(t)=Ri(1);
                else
                    R(t)=Ri(2);
                end
            end
            q(t)=(aa+az(t))/(aa+bb+vz(t));
            qb(t)=(aa+azb(t))/(aa+bb+vzb(t));
            g=normrnd(mu,sig,N);
            for n=1:N      
                if q(t)*F<=g(n)
                    v(t)=v(t)+1;
                end
                if qb(t)*F<=g(n)
                    vb(t)=vb(t)+1;
                end                
            end
            %vv(r,t)=v(t);
            A(t)=gam*min(1,R(t)/v(t));
            Ab(t)=gam*min(1,Rb(t)/vb(t));
            %A(t)=gam*(1-alpha^(-(R/v(t))));
            a(t)=binornd(v(t),A(t));
            ab(t)=binornd(vb(t),Ab(t));
            for i=1:NN
                for j=1:NN
                    if t>1
                        if v(t-1)==i-1
                            if v(t)==j-1
                                vv(i,j)=vv(i,j)+1;
                            end
                        end
                    end
                    if t>1
                        if vb(t-1)==i-1
                            if vb(t)==j-1
                                vvb(i,j)=vvb(i,j)+1;
                            end
                        end
                    end
                    
                end
            end

        end
        %vp(i)
    %end
        %TPvv=(vv.*vp)/(Block-1);
    for i=1:NN
        V(r,i)=0;
        Vb(r,i)=0;
        for j=1:NN
            Pvv(i,j)=vv(i,j)/sum(vv(i,:));
            if isnan(Pvv(i,j))
                Pvv(i,j)=0;
            end
            Pvvb(i,j)=vvb(i,j)/sum(vvb(i,:));
            if isnan(Pvvb(i,j))
                Pvvb(i,j)=0;
            end            
        end
    end
    Evv(r,:)=Pvv*vp;
    U(r,:)=Evv(r,:)-vp';
    Evvb(r,:)=Pvvb*vp;
    Ub(r,:)=Evvb(r,:)-vp';
    for i=1:NN
        if U(r,i)>=0
            L(r,i)=1;
        else
            L(r,i)=-1;
        end
        if i>1            
            if L(r,i)-L(r,i-1)~=0
                count(i)=i-1;
            end
        end
        if Ub(r,i)>=0
            Lb(r,i)=1;
        else
            Lb(r,i)=-1;
        end
        if i>1            
            if Lb(r,i)-Lb(r,i-1)~=0
                countb(i)=i-1;
            end
        end 
    end    
end 
%load PrEq2
%load TPrEq2
%load GPrEq2
%load BPrEq2
%load vEq2
[acf,lag]=autocorr(v,10);
[acfb,lagb]=autocorr(vb,10);
xx=[17,17];
yy=[0,17];
xx1=[54,54];
yy1=[0,N];
xx2=[88,88];
yy2=[0,88];


%fig1=figure(1);
%hold on
%box on
%plot(vp, Evv(1,:),'Color','k','LineStyle','--')
%PP=get(gca,'Position');
%plot(vp, vp,'Color','k','LineStyle','-')
%plot(xx,yy,'Color','k','Linestyle',':')
%plot(xx1,yy1,'Color','k','LineStyle',':');
%plot(xx2,yy2,'Color','k','Linestyle',':')
%quiver(vp(1:2:NN),V(1,(1:2:NN)),U(1,(1:2:NN)),V(1,(1:2:NN)),'Color',[0.2 0.2 0.2],'AutoScaleFactor',1.5,'LineWidth',0.9,'MaxHeadSize',0.1)
%xlabel('Violations in current period (v^{t})')
%ylabel('Expected violation in the next period E(v^{t+1}|v^{t})')
%legend('Expected violations given current violations, E(v^{t+1}|v^{t})','Location','northwest') 
%txt={'v^{t}=17'};
%text(15,-3,txt);
%txt1={'v^{t}=54'};
%text(52,-3,txt1);
%txt2={'v^{t}=88'};
%text(86,-3,txt2);
%ylim([-5 N])
%hold off
%fig2=figure(2);
%hold on
%box on
%stem(lag,acf,'Filled','Color','k')
%title('Autocorrelations in violations')
%xlabel('Lags')
%hold off
%insize=0.35;
%[hm, hs]=inset(fig1,fig2,insize);

figure % Figure6.2
tile=tiledlayout(4,1);
tile.Padding='none';
tile.TileSpacing='tight';
nexttile ([2 1])
hold on
box on
plot(vp, Evv(1,:),'Color','k','LineStyle','--')
%PP=get(gca,'Position');
plot(vp, vp,'Color','k','LineStyle','-')
plot(xx,yy,'Color','k','Linestyle',':')
%plot(xx1,yy1,'Color','k','LineStyle',':');
%plot(xx2,yy2,'Color','k','Linestyle',':')
quiver(vp(1:2:NN),V(1,(1:2:NN)),U(1,(1:2:NN)),V(1,(1:2:NN)),'Color',[0.2 0.2 0.2],'AutoScaleFactor',2.5,'LineWidth',0.9,'MaxHeadSize',0.25)
%handaxes=axes('position',[0.15 0.52 0.3 0.3]);
%stem(lag,acf,'Filled','Color','k')
xlabel('Violations in current period (v^{t})')
ylabel('Expected violation in the next period E(v^{t+1}|v^{t})')
legend('E(v^{t+1}|v^{t})','Location','north') 
title('Expected violations given current violations, R_{GB}=39 , R_{BB}=60')
legend boxoff
txt={'v^{t}=17'};
text(15,-3,txt);
txt1={'v^{t}=53'};
text(52,-3,txt1);
txt2={'v^{t}=88'};
text(86,-3,txt2);
ylim([-5 N])
hold off
nexttile ([2 1])
hold on
box on
plot(vp, Evvb(1,:),'Color','k','LineStyle','--')
%PP=get(gca,'Position');
plot(vp, vp,'Color','k','LineStyle','-')
plot(xx,yy,'Color','k','Linestyle',':')
plot(xx1,yy1,'Color','k','LineStyle',':');
plot(xx2,yy2,'Color','k','Linestyle',':')
quiver(vp(1:2:NN),V(1,(1:2:NN)),Ub(1,(1:2:NN)),Vb(1,(1:2:NN)),'Color',[0.2 0.2 0.2],'AutoScaleFactor',1.5,'LineWidth',0.9,'MaxHeadSize',0.15)
%handaxes=axes('position',[0.15 0.52 0.3 0.3]);
%stem(lag,acf,'Filled','Color','k')
xlabel('Violations in current period (v^{t})')
ylabel('Expected violation in the next period E(v^{t+1}|v^{t})')
legend('E(v^{t+1}|v^{t})','Location','north') 
title('Expected violations given current violations, R = 39')
legend boxoff
txt={'v^{t}=17'};
text(15,-3,txt);
txt1={'v^{t}=53'};
text(52,-3,txt1);
txt2={'v^{t}=88'};
text(86,-3,txt2);
ylim([-5 N])
hold off

handaxes1=axes('position',[0.07 0.72 0.22 0.22]);
stem(lag(2:11),acf(2:11),'Filled','Color','k')
title('Autocorrelation in violations')
xlabel('Lags')
xlim([0 10])

handaxes2=axes('position',[0.07 0.22 0.22 0.22]);
stem(lagb(2:11),acfb(2:11),'Filled','Color','k')
title('Autocorrelation in violations')
xlabel('Lags')
xlim([0 10])

handaxes3=axes('position',[0.72 0.62 0.22 0.22]);
histogram(v(1000000:2000000),'Normalization','probability','FaceColor',[0.17 0.17 0.17],'BinWidth',1)
title('frequncy of violations')
xlabel('v')

handaxes4=axes('position',[0.72 0.12 0.22 0.22]);
histogram(vb,'Normalization','probability','FaceColor',[0.17 0.17 0.17],'BinWidth',1)
title('frequncy of violations')
xlabel('v')



%nexttile ([1 1])
%hold on
%box on
%plot(vEq2,PrEq2(39,:),'Color','k');
%area(vEq2,TPrEq2(39,:),'FaceColor',[0.6 0.6 0.6]);
%area(vEq2,GPrEq2(39,:),'FaceColor',[0.9 0.9 0.9]);
%area(vEq2,BPrEq2(39,:),'FaceColor',[0.1 0.1 0.1]);
%title('Panel 1: Probability of an equilibrium with v violations, R = 39')
%xlabel('Number of violations (v)')
%txt={'{\leftarrow} Probability of good eq (v \leq 30)'};
%text(18,0.10,txt);
%txt1={'{\downarrow} Probability of other eq (30 < v < 70)'};
%text(52,0.015,txt1);
%txt2={'Probability of bad eq (v \geq 70) {\rightarrow}'} ;
%text(75,0.07,txt2);
%ylim([0 0.15])
%hold off
