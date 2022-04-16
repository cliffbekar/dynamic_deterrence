%ASBbl3convperV2fig5_1_2_6.m
%Aug, 2021 Baseline, convergence, two bin, perisitence, single R model  

clear

%parameters 

N=100;      %population of agents
NN=N+1;
%M=20;
MM=200;

Block=20000;

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

eps=4;

dim=3;
MC=150;
cost=zeros(MC,1);
Mcost=zeros(dim,1);SDcost=zeros(dim,1);

RR=[43 44 45];
for i=1:dim
    for m=1:MC
        critconv1=0.01;
        edges=zeros(NN+1,1);
        for j=1:NN+1
            edges(j)=j-1.5;
        end
        c=1;
        b=0;
        cc=0;
        TESTCONV=zeros(MM,1);
        TESTCONV(1)=1;
        while ((cc<1) && (b<MM))
            g=normrnd(mu,sig,Block,N); % unifrnd(0,1,Block,N); %
            swi=zeros(Block,1);
            swic=zeros(Block,1);
            vz=zeros(Block,1);
            az=zeros(Block,1);
            eta=zeros(Block,1);
            R=zeros(Block,1);
            gbar=zeros(Block,N);
            gbar3=zeros(Block,1);
            A=zeros(Block,1);
            NCw=zeros(Block,1);NCw2=zeros(Block,1);
            vw=zeros(Block,1);aw=zeros(Block,1);
            b=b+1;
            %BH(r)=b;
            if b<2
                in=Z+1;
            else
                in=1;
            end
            for t=1:Block
                R(t)=RR(i);
                if (b<2)
                    if (t<Z+1)
                        for z=1:Z
                            vz(z)=Z*unifrnd(0,N);
                            az(z)=Z*unifrnd(0,vz(z));
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
                eta(t)=(aa+az(t))/(aa+bb+vz(t));
                for n=1:N
                    if eta(t)*F<=g(t,n) %g(t+ind(r),n)
                        vw(t)=vw(t)+1;
                        gbar(t,n)=g(t,n); %g(t,n);
                    end
                end
                gbar3(t)=sum(gbar(t,:));
                A(t)=gam*min(1,R(t)/vw(t));
                %A(t)=gam*(1-1/(eps^(R(t)/v(t))));
                aw(t)=binornd(vw(t),A(t));
                NCw(t)=rho*R(t)+(lam-1)*gbar3(t);
            end
            if b<2
                v=vw;
                a=aw;
                gb=gbar3;
                NC=NCw;
                Ap=A;
                qp=eta;
            else
                vhold=cat(2,v',vw');
                ahold=cat(2,a',aw');
                gbhold=cat(2,gb',gbar3');
                NChold=cat(2,NC',NCw');
                Aphold=cat(2,Ap',A');
                qphold=cat(2,qp',eta');
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
        cost(m)=mean(NC);
    end
    Mcost(i)=mean(cost);
    SDcost(i)=std(cost);
    i
end
OPT=[RR' Mcost SDcost];

%save('C:\Users\kcarlaw\Documents\MATLAB\ASB2021\MC150\OPT.txt','OPT','-ascii')