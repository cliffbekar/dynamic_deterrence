%ASBconV2fig31to9.m
%March, 2022 Baseline, convergence, two bin, perisitence, passive R   

clear

%parameters 

N=60;      %population of agents
NN=N+1;
MM=200;

Block=50000;

    %baseline parameterization
F=1;        %saction for ASB
gam=0.8;    % max objective apprehension prob.
aa=1;       %shape poarameter 1 for Beta distribution (Bayesian) 
bb=0.25;    %shape poarameter 2 for Beta distribution (Bayesian)
%mu=0.7;     %mean value of g, individual benefit from ASA
%sig=0.3;    %varance of g 
lam=5;      %socail cost parameter for ASB 
rho=2;      %unit cost of enforcement resources
Z=2;        %z-history length

BinE=50; % Bin edge passive policies

eps=4;      %parameter for exponenetial, objective apprehension function

critconv1=0.01;
RR=zeros(NN,1);
mv=zeros(NN,1);
%mov=zeros(NN,1);
%sdv=zeros(NN,1);
%coefv=zeros(NN,1);
%coefvw=zeros(NN,1);
%ma=zeros(NN,1);
%sda=zeros(NN,1);
%gpiv=zeros(NN,1);
%vcon=zeros(NN,Block);OBP=zeros(NN,Block);SBP=zeros(NN,Block);
cost=zeros(NN,1);
%cost2=zeros(NN,1);
%mgb=zeros(NN,1);
%ARATE=zeros(NN,1);
%ARATER=NaN(NN,1);
%IRATE=zeros(NN,1);
%IRATER=NaN(NN,1);
%EXCAP=NaN(NN,1);
%EI=zeros(NN,1);
%elsc=zeros(NN,1);elsu=zeros(NN,1);pc=zeros(NN,1);pu=zeros(NN,1);Tc=zeros(NN,1);Tu=zeros(NN,1);
%mvc=zeros(NN,1);mvu=zeros(NN,1);sdvc=zeros(NN,1);sdvu=zeros(NN,1);
%SIGGc=zeros(NN,1);SIGGu=zeros(NN,1);
%OPERRc=zeros(NN,1);OPERRu=zeros(NN,1);PrEvc=zeros(NN,1);PrEvu=zeros(NN,1);
BH=zeros(NN,1);

edges=zeros(NN+1,1);
for j=1:NN+1
    edges(j)=j-1.5;
end

X=21;
Y=14;
mmu=zeros(X,1);
ssig=zeros(Y,1);
GenT=zeros(X,Y,1);CritGT=zeros(X,Y,1);
MDROP2=zeros(X,Y,1);
MDROP1=zeros(X,Y,1);
MDROP=zeros(X,Y,1);
mu=0.0*F;
Hmv=zeros(X,NN,1);
for i=1:X
    if i>1
        mu=mu+0.05*F;
    end
    mmu(i)=mu;
    sig=0.1*F;
    for j=1:Y
        %sig(1)=0.05;
        if j>1
            sig=sig+0.05*F;
        end
        ssig(j)=sig;
        CritGT(i,j)=0.35;
        for r=1:NN
            RR(r)=r-1;
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
                vw=zeros(Block,1);aw=zeros(Block,1);ARatw=zeros(Block,1);
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
                    ARatw(t)=aw(t)/vw(t);
                end
                if b<2
                    v=vw;
                    a=aw;
                    gb=gbar3;
                    NC=NCw;
                    Ap=A;
                    qp=eta;
                    ARat=ARatw;
                else
                    vhold=cat(2,v',vw');
                    ahold=cat(2,a',aw');
                    gbhold=cat(2,gb',gbar3');
                    NChold=cat(2,NC',NCw');
                    Aphold=cat(2,Ap',A');
                    qphold=cat(2,qp',eta');
                    ARathold=cat(2,ARat',ARatw');
                    v=vhold';
                    a=ahold';
                    gb=gbhold';
                    NC=NChold';
                    Ap=Aphold';
                    qp=qphold';
                    ARat=ARathold';
                end

                vconv=v(1:b*Block);
                vconvlag=v(1:(b-1)*Block);
                freqvconv=histcounts(vconv(:),edges)/((b)*Block);
                freqvconvlag=histcounts(vconvlag(:),edges)/((b-1)*Block);
                TESTCONV1=zeros(1,NN);
                if b>1
                    for k=1:NN
                        TESTCONV1(k)=abs(freqvconv(k)-freqvconvlag(k));
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
            cost(r)=mean(NC);
            %RR(r)

        end
        if j<2
            Hmv(i,:)=mv(:);
        end
        GenT(i,j)=(mv(1)-mv(2))/N;
        TD=mv(1)-mv(NN);
        DROPn5=zeros(NN,1);
        DROP=zeros(NN,1);
        DROP5=zeros(NN,1);
        for r=1:NN
            if r>5
                DROPn5(r)=(mv(r-3)-mv(r))/N;
            end
            if r>5
                DROP(r)=(mv(r-2)-mv(r))/TD;
            end
            if r>5
                DROP5(r)=(mv(r-3)-mv(r))/TD;
            end
        end
        MDROP2(i,j)=max(DROPn5);
        MDROP1(i,j)=max(DROP);
        MDROP(i,j)=max(DROP5);

        j
    end
    i
end

figure
mesh(MDROP2,'EdgeColor',[0.5 0.5 0.5],'FaceColor',[0.9 0.9 0.9]);
hold on
mesh(CritGT,'EdgeColor',[0.5 0.5 0.5],'FaceColor',[0.3 0.3 0.3]);
%scatter3(2,17,GenT(17,2)',200,'o','filled','MarkerEdgeColor','k','MarkerFaceColor','k')
title('F*${\gamma}$=0.8, DROP=$\frac{E(v|R-5)-E(v|R)}{N}$','Interpreter','latex')
xlabel('{\sigma}','FontSize',25)
xticks([1 2 3 4 5 6 7])
xticklabels({num2str(ssig(1)),num2str(ssig(2)),num2str(ssig(3)),num2str(ssig(4)),num2str(ssig(5)),num2str(ssig(6)),num2str(ssig(7))})
ylabel('{\mu}','FontSize',25)
yticks([1 3 5 7 9 11 13 15 17 19 21])
yticklabels({num2str(mmu(1)),num2str(mmu(3)),num2str(mmu(5)),num2str(mmu(7)),num2str(mmu(9)),num2str(mmu(11)),num2str(mmu(13)),num2str(mmu(15)),num2str(mmu(17)),num2str(mmu(19)),num2str(mmu(21))})
zlabel('Value of DROP')
xlim([1 Y])
ylim([1 X])
hold off

figure
mesh(MDROP2,'EdgeColor',[0.5 0.5 0.5],'FaceColor',[0.9 0.9 0.9]);
hold on
mesh(CritGT,'EdgeColor',[0.5 0.5 0.5],'FaceColor',[0.3 0.3 0.3]);
%scatter3(2,17,GenT(17,2)',200,'o','filled','MarkerEdgeColor','k','MarkerFaceColor','k')
%title('F*${\gamma}$=0.8, DROP=$\frac{E(v|R-5)-E(v|R)}{N}$','Interpreter','latex')
xlabel('${\frac{\sigma}{F*{\gamma}}}$','Interpreter','latex','FontSize',25)
xticks([1 2 3 4 5 6 7])
xticklabels({num2str(ssig(1)/(F*gam)),num2str(ssig(2)/(F*gam)),num2str(ssig(3)/(F*gam)),num2str(ssig(4)/(F*gam)),num2str(ssig(5)/(F*gam)),num2str(ssig(6)/(F*gam)),num2str(ssig(7)/(F*gam))})
ylabel('${\frac{\mu}{F*{\gamma}}}$','Interpreter','latex','FontSize',25)
yticks([1 3 5 7 9 11 13 15 17 19 21])
yticklabels({num2str(mmu(1)/(F*gam)),num2str(mmu(3)/(F*gam)),num2str(mmu(5)/(F*gam)),num2str(mmu(7)/(F*gam)),num2str(mmu(9)/(F*gam)),num2str(mmu(11)/(F*gam)),num2str(mmu(13)/(F*gam)),num2str(mmu(15)/(F*gam)),num2str(mmu(17)/(F*gam)),num2str(mmu(19)/(F*gam)),num2str(mmu(21)/(F*gam))})
zlabel('Value of DROP')
xlim([1 Y])
ylim([1 X])
hold off

%[M, I]=min(cost);
%Rst=RR(I);

%clifeqdata=table(RR(38:45),mvc(38:45),sdvc(38:45),mvu(38:45),sdvu(38:45));

%save C:\Users\kcarlaw\Documents\MATLAB\ASB2022final\mvbl5 mv -ASCII -DOUBLE;

save('C:\Users\kcarlaw\Documents\MATLAB\ASB2022final\fg8DROPn.txt','MDROP2','-ascii');


