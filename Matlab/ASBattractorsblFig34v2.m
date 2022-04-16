
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
mu=0.6;     %mean value of gi, individual benefit from ASA
sig=0.2;    %varance of gi 
lam=5;     %socail cost conversion of individual ASA for social damage function
PR=2;
R=[39]; % vectore of R values
RR=1;
Z=2;
Block=10000000;  %number of sim periods in convergence loop

BinE=50;

Evv=zeros(RR,NN,1);
U=zeros(RR,NN,1);
V=zeros(RR,NN,1);
L=zeros(RR,NN,1);
count=zeros(NN,1);
vp=zeros(NN,1);

for i=1:NN
    vp(i)=i-1;
end

for r=1:RR
    vv=zeros(NN,NN,1);
    Pvv=zeros(NN,NN,1);
    
%    for i=1:NN
%        vp(i)=i-1;
        
        q=zeros(Block,1);A=zeros(Block,1);
        vz=zeros(Block,1);az=zeros(Block,1);
        v=zeros(Block,1);a=zeros(Block,1);
        for t=1:Block
           if (t<=Z)
                for z=1:Z
                    vz(z)=Z*unifrnd(0,N);
                    az(z)=Z*unifrnd(0,vz(z));
                end
            else
                vz(t)=sum(v(t-Z:t-1));
                az(t)=sum(a(t-Z:t-1));
           end

            q(t)=(aa+az(t))/(aa+bb+vz(t));

            g=normrnd(mu,sig,N);
            for n=1:N      
                if q(t)*F<=g(n)
                    v(t)=v(t)+1;
                end
            end
            %vv(r,t)=v(t);
            A(t)=gam*min(1,R(r)/v(t));
            %A(t)=gam*(1-alpha^(-(R/v(t))));
            a(t)=binornd(v(t),A(t));
            
            for i=1:NN
                for j=1:NN
                    if t>1
                        if v(t-1)==i-1
                            if v(t)==j-1
                                vv(i,j)=vv(i,j)+1;
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
        for j=1:NN
            Pvv(i,j)=vv(i,j)/sum(vv(i,:));
            if isnan(Pvv(i,j))
                Pvv(i,j)=0;
            end
        end
    end
    Evv(r,:)=Pvv*vp;
    U(r,:)=Evv(r,:)-vp';
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
    end
    
    
end 
load PrEq2      %From ASBEquilibMDL.m
load TPrEq2     %From ASBEquilibMDL.m
load GPrEq2     %From ASBEquilibMDL.m
load BPrEq2     %From ASBEquilibMDL.m
load vEq2       %From ASBEquilibMDL.m
[acf,lag]=autocorr(v,10);
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


figure % Figure3.4
tile=tiledlayout(4,1);
tile.Padding='none';
tile.TileSpacing='none';
nexttile ([3 1])
hold on
box on
plot(vp, Evv(1,:),'Color','k','LineStyle','--')
PP=get(gca,'Position');
plot(vp, vp,'Color','k','LineStyle','-')
plot(xx,yy,'Color','k','Linestyle',':')
plot(xx1,yy1,'Color','k','LineStyle',':');
plot(xx2,yy2,'Color','k','Linestyle',':')
quiver(vp(1:2:NN),V(1,(1:2:NN)),U(1,(1:2:NN)),V(1,(1:2:NN)),'Color',[0.2 0.2 0.2],'AutoScaleFactor',1.5,'LineWidth',0.9,'MaxHeadSize',0.1)
%handaxes=axes('position',[0.15 0.52 0.3 0.3]);
%stem(lag,acf,'Filled','Color','k')
xlabel('Violations in current period (v^{t})')
ylabel('Expected violation in the next period E(v^{t+1}|v^{t})')
legend('E(v^{t+1}|v^{t})','Location','east') 
title('Expected violations given current violations, R = 39')
legend boxoff
txt={'v^{t}=17'};
text(15,-3,txt);
txt1={'v^{t}=54'};
text(52,-3,txt1);
txt2={'v^{t}=88'};
text(86,-3,txt2);
ylim([-5 N])
hold off
nexttile ([1 1])
hold on
box on
plot(vEq2,PrEq2(39,:),'Color','k');
area(vEq2,TPrEq2(39,:),'FaceColor',[0.6 0.6 0.6]);
area(vEq2,GPrEq2(39,:),'FaceColor',[0.9 0.9 0.9]);
area(vEq2,BPrEq2(39,:),'FaceColor',[0.1 0.1 0.1]);
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
handaxes2=axes('position',[0.1 0.65 0.25 0.25]);
stem(lag,acf,'Filled','Color','k')
title('Autocorrelation in violations')
xlabel('Lags')
