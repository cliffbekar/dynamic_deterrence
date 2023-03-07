% ASBcon2Rquickfig61.m
% Updated Aug. 30, 2022 K.I. Carlaw

clear

%parameters
T=10000;    %Time/iteration index for the sim
N=100;      %population of agents
F=1;        %individual cost of apprehnsion given ASA
gam=0.8;
aa=1;
bb=0.25; %0.1765;
%critical distribution of ASA (anti-scoial act) RN(gi;mu,sig,0,inf)
mu=0.6;     %mean value of gi, individual benefit from ASA
sig=0.2;    %varance of gi 
lam=5;     %socail cost conversion of individual ASA for social damage function
PR=2;
BinE=50;
Rc=37;
Ru=59;
Rstar=Rc;
Z=2;

load fig61v
v=fig61v;
load fig61bin
bin=fig61bin;
load fig61ARATE
ARATE=fig61ARATE;
load fig61R
R=fig61R;
load fig61A
A=fig61A;
load fig61q
q=fig61q;

%Dimensions of Varialbles

%for r=1:RR
%g=normrnd(mu,sig,T,N);
%R=zeros(T,1);
%R(1)=Rstar;
%R(2)=Rstar;
%R(3)=Rstar;
%R(4)=Rstar;
%v=zeros(T,1);   %number of violators per period
%vz=zeros(T,1);
%a=zeros(T,1);   %number of apprehended violators per period
%az=zeros(T,1);
%A=zeros(T,1); %Aprehnsion probability
%q=zeros(T,1);   %subjective prob of apprehention
%bin=zeros(T,1);ARATE=zeros(T,1);

% Sim loop
%    for t=1:T
%        if t>1
%            if v(t-1)<BinE
%                R(t)=Rc;
%            else
%                R(t)=Ru;
%            end
%        end
%        if (t<=Z)
%            for z=1:Z
%                vz(z)=Z*unifrnd(0,N);
%                az(z)=Z*unifrnd(0,vz(z));
%            end
%        else
%            vz(t)=sum(v(t-Z:t-1));
%            az(t)=sum(a(t-Z:t-1));
%        end
%        q(t)=(aa+az(t))/(aa+bb+vz(t));
%        for n=1:N
%            if q(t)*F<=g(t,n)
%                v(t)=v(t)+1;
%            end
%        end
%        A(t)=gam*min(1,R(t)/v(t));
%        %A(t)=gam*(1-alpha^(-(R/v(t))));
%        a(t)=binornd(v(t),A(t));
%        bin(t)=BinE;
%        ARATE(t)=a(t)/v(t);
%        if isnan(ARATE(t))
%            ARATE(t)=0;
%        end
%    end
%figure 
%histogram(v)
%figure 
%plot([v a])

figure 
tile=tiledlayout(3,1);
tile.Padding='none';
tile.TileSpacing='tight';
nexttile
hold on
box on
colororder({'k','k'})
yyaxis left
%plot(v,'+','MarkerSize',4,'MarkerIndices',1:5:T,'Color','k','LineStyle',"-",'LineWidth',0.5)
plot(v,'Color','k')
plot(bin,'Color','k','LineStyle','-.');
ylim([0 N])
xlim([1400 1500])
ylabel('Number of violations','Color','k')
yyaxis right
%plot(q,'*','MarkerSize',4,'MarkerIndices',1:5:T,'Color','k','LineStyle',"-.",'LineWidth',0.5)
plot(ARATE,'Color','k','LineStyle',"--")
ylim([0 1])
xlim([1400 1500])
ylabel('Apprehension rate','Color','k')
xline(1458,'Color','k','LineStyle',':');
xline(1460,'Color','k','LineStyle',':');
xticks([1400 1410 1420 1430 1440 1450 1460 1470 1480 1490 1500]);
xticklabels({'0','10','20','30','40','50','60','70','80','90','100'})
%xline(2221,'Color','k','LineStyle',':');
%xline(2227,'Color','k','LineStyle',':');
legend('violations (v)','Bin edge = 39','apprehension rate (a/v)','Location','southeast')
legend boxoff
title('Panel 1: Illustrative path of violation and apprehension rate')
hold off
nexttile
hold on
box on
plot(R,'Color','k','LineWidth',1)
ylim([30 65])
xlim([1400 1500])
%xline(1459,'Color','k','LineStyle',':');
%xline(1461,'Color','k','LineStyle',':');
%xline(2221,'Color','k','LineStyle',':');
%xline(2227,'Color','k','LineStyle',':');
ylabel('Enforcement resources (R)')
xticks([1400 1410 1420 1430 1440 1450 1460 1470 1480 1490 1500]);
xticklabels({'0','10','20','30','40','50','60','70','80','90','100'})
txt1='\leftarrow crackdown resources';
txt2='crackdown implemented\rightarrow';
text(1440,58,txt2,'FontSize',12);
%text(2180,30,txt2);
title('Panel 2: Illustrative path of enforcement resource (R) deployment')
hold off 
nexttile
hold on
box on
plot(A,'Color','k')
plot(q,'Color','k','LineStyle','--')
xlim([1400 1500])
legend('Objective probability of apprehension','Subjective probability of apprehension','Location','southeast')
legend boxoff
xlabel('Tics(t)')
%xline(1458,'Color','k','LineStyle',':');
%xline(1460,'Color','k','LineStyle',':');
xticks([1400 1410 1420 1430 1440 1450 1460 1470 1480 1490 1500]);
xticklabels({'0','10','20','30','40','50','60','70','80','90','100'})
title('Panel 3: Illustrative path of objective and subjective probabilities')
hold off



