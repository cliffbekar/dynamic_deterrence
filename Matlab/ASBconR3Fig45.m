% 3 Bin contingent R strategies to efficiently manage ASB 
%Nov 29, 2021
%optimal BE BE2 Rgb Rtb Rbb

clear
%parameters 
N=100;      %population of agents
NN=N+1;
MM=200;

Block=20000;

F=1;        %individual cost of apprehnsion given ASA
gam=0.8;
aa=1;
bb=0.25; %0.1765;
mu=0.6;     %mean value of gi, individual benefit from ASA
sig=0.2;    %varance of gi 
lam=5;     %socail cost conversion of individual ASA for social damage function
rho=2;
SWrho=5*rho;
mm=0;
Z=2;        %z-history length
X=9;

MC=50;
BnE1=zeros(MC,1);BnE2=zeros(MC,1);RRRgb=zeros(MC,1);RRRtb=zeros(MC,1);RRRbb=zeros(MC,1);CM=zeros(MC,1);
count=zeros(MC,1);

 
    Rgb=30; %round(unifrnd(0,NN));
    Rtb=47; %round(unifrnd(Rgb,NN));
    Rbb=64; %round(unifrnd(Rtb,NN));
    BE=34; %round(unifrnd(0,NN));
    BE2=44; %round(unifrnd(BE,NN));

    crit=1;
    crit2=138.3939;
    RbbM=138.3939+1;

    critconv1=0.01;

    count3=0;
    RR=zeros(NN,1);
    ylin=zeros(NN,1);
    for r=1:NN
        RR(r)=r-1;
        ylin(r)=130.33;
    end

    edges=zeros(NN+1,1);
    for j=1:NN+1
        edges(j)=j-1.5;
    end

        %Search for Bin edge 1
        BinE=zeros(N,1);
        Bcost=zeros(N,1);
        BEd=BE-X;
        BEu=BE+X;
        for i=BEd:BEu
            BinE(i)=i;
            BcostMC=zeros(MC,1);
            for m=1:MC
                c=1;
                b=0;
                cc=0;
                TESTCONV=zeros(MM,1);
                TESTCONV(1)=1;
                while ((cc<1) && (b<MM))
                    g=normrnd(mu,sig,Block,N);
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
                    b=b+1;
                    if b<2
                        in=Z+1;
                    else
                        in=1;
                    end
                    for t=1:Block
                        if t>1
                            if vw(t-1)<BinE(i)
                                R(t)=Rgb;
                            elseif (BinE(i)<=vw(t-1))&&(vw(t-1)<BE2)
                                R(t)=Rtb;
                            elseif BE2<=vw(t-1)
                                R(t)=Rbb;
                            end
                            if R(t)~=R(t-1)
                                swi(t)=abs(Rbb-Rgb); 
                                swic(t)=1;
                            end                    
                        end                
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
                        q(t)=(aa+az(t))/(aa+bb+vz(t));
                        for n=1:N
                            if q(t)*F<=g(t,n) %g(t+ind(r),n)
                                vw(t)=vw(t)+1;
                                gbar(t,n)=g(t,n); %g(t,n);
                            end
                        end
                        gbar3(t)=sum(gbar(t,:));
                        A(t)=gam*min(1,R(t)/vw(t));
                        %A(t)=gam*(1-1/(eps^(R(t)/v(t))));
                        aw(t)=binornd(vw(t),A(t));
                        NCw(t)=rho*R(t)+(lam-1)*gbar3(t)+SWrho*swi(t);
                        NCw2(t)=rho*R(t)+(lam-1)*gbar3(t);
                    end
                    if b<2
                        v=vw;
                        a=aw;
                        gb=gbar3;
                        NC=NCw;
                        NC2=NCw2;
                    else
                        vhold=cat(2,v',vw');
                        ahold=cat(2,a',aw');
                        gbhold=cat(2,gb',gbar3');
                        NChold=cat(2,NC',NCw');
                        NC2hold=cat(2,NC2',NCw2');
                        v=vhold';
                        a=ahold';
                        gb=gbhold';
                        NC=NChold';
                        NC2=NC2hold';
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
                BcostMC(m)=mean(NC2);
            end
            Bcost(i)=mean(BcostMC);
            BinE(i)
        end

        %search for Bin edge 2
        BinE2=zeros(N,1);
        Bcost2=zeros(N,1);
        BE2d=BE2-X;
        BE2u=BE2+X;
        for i=BE2d:BE2u
            BinE2(i)=i;
            Bcost2MC=zeros(MC,1);
            for m=1:MC
                c=1;
                b=0;
                cc=0;
                TESTCONV=zeros(MM,1);
                TESTCONV(1)=1;
                while ((cc<1) && (b<MM))
                    g=normrnd(mu,sig,Block,N);
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
                    b=b+1;
                    if b<2
                        in=Z+1;
                    else
                        in=1;
                    end
                    for t=1:Block
                        if t>1
                            if vw(t-1)<BE
                                R(t)=Rgb;
                            elseif (BE<=vw(t-1))&&(vw(t-1)<BinE2(i))
                                R(t)=Rtb;
                            elseif BinE2(i)<=vw(t-1)
                                R(t)=Rbb;
                            end
                            if R(t)~=R(t-1)
                                swi(t)=abs(Rbb-Rgb); 
                                swic(t)=1;
                            end                    
                        end                
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
                        q(t)=(aa+az(t))/(aa+bb+vz(t));
                        for n=1:N
                            if q(t)*F<=g(t,n) %g(t+ind(r),n)
                                vw(t)=vw(t)+1;
                                gbar(t,n)=g(t,n); %g(t,n);
                            end
                        end
                        gbar3(t)=sum(gbar(t,:));
                        A(t)=gam*min(1,R(t)/vw(t));
                        %A(t)=gam*(1-1/(eps^(R(t)/v(t))));
                        aw(t)=binornd(vw(t),A(t));
                        NCw(t)=rho*R(t)+(lam-1)*gbar3(t)+SWrho*swi(t);
                        NCw2(t)=rho*R(t)+(lam-1)*gbar3(t);
                    end
                    if b<2
                        v=vw;
                        a=aw;
                        gb=gbar3;
                        NC=NCw;
                        NC2=NCw2;
                    else
                        vhold=cat(2,v',vw');
                        ahold=cat(2,a',aw');
                        gbhold=cat(2,gb',gbar3');
                        NChold=cat(2,NC',NCw');
                        NC2hold=cat(2,NC2',NCw2');
                        v=vhold';
                        a=ahold';
                        gb=gbhold';
                        NC=NChold';
                        NC2=NC2hold';
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
                Bcost2MC(m)=mean(NC2);
            end
            Bcost2(i)=mean(Bcost2MC);
            BinE2(i)
        end

        %Search for Rgb
        RRgb=zeros(NN,1);
        Rgbcost=zeros(NN,1);
        Rgbd=Rgb-X;
        Rgbu=Rgb+X;
        for i=Rgbd:Rgbu
            RRgb(i)=i;
            RgbcostMC=zeros(MC,1);
            for m=1:MC
                c=1;
                b=0;
                cc=0;
                TESTCONV=zeros(MM,1);
                TESTCONV(1)=1;
                while ((cc<1) && (b<MM))
                    g=normrnd(mu,sig,Block,N);
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
                    b=b+1;
                    if b<2
                        in=Z+1;
                    else
                        in=1;
                    end
                    for t=1:Block
                        if t>1
                            if vw(t-1)<BE
                                R(t)=RRgb(i);
                            elseif (BE<=vw(t-1))&&(vw(t-1)<BE2)
                                R(t)=Rtb;
                            elseif BE2<=vw(t-1)
                                R(t)=Rbb;
                            end
                            if R(t)~=R(t-1)
                                swi(t)=abs(Rbb-RRgb(i)); 
                                swic(t)=1;
                            end                    
                        end                
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
                        q(t)=(aa+az(t))/(aa+bb+vz(t));
                        for n=1:N
                            if q(t)*F<=g(t,n) %g(t+ind(r),n)
                                vw(t)=vw(t)+1;
                                gbar(t,n)=g(t,n); %g(t,n);
                            end
                        end
                        gbar3(t)=sum(gbar(t,:));
                        A(t)=gam*min(1,R(t)/vw(t));
                        %A(t)=gam*(1-1/(eps^(R(t)/v(t))));
                        aw(t)=binornd(vw(t),A(t));
                        NCw(t)=rho*R(t)+(lam-1)*gbar3(t)+SWrho*swi(t);
                        NCw2(t)=rho*R(t)+(lam-1)*gbar3(t);
                    end
                    if b<2
                        v=vw;
                        a=aw;
                        gb=gbar3;
                        NC=NCw;
                        NC2=NCw2;
                    else
                        vhold=cat(2,v',vw');
                        ahold=cat(2,a',aw');
                        gbhold=cat(2,gb',gbar3');
                        NChold=cat(2,NC',NCw');
                        NC2hold=cat(2,NC2',NCw2');
                        v=vhold';
                        a=ahold';
                        gb=gbhold';
                        NC=NChold';
                        NC2=NC2hold';
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
                RgbcostMC(m)=mean(NC2);
            end
            Rgbcost(i)=mean(RgbcostMC);
            RRgb(i)
        end

        %Search for Rtb
        RRtb=zeros(NN,1);
        Rtbcost=zeros(NN,1);
        Rtbd=Rtb-X;
        Rtbu=Rtb+X;
        for i=Rtbd:Rtbu
            RRtb(i)=i;
            RtbcostMC=zeros(MC,1);
            for m=1:MC
                c=1;
                b=0;
                cc=0;
                TESTCONV=zeros(MM,1);
                TESTCONV(1)=1;
                while ((cc<1) && (b<MM))
                    g=normrnd(mu,sig,Block,N);
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
                    b=b+1;
                    if b<2
                        in=Z+1;
                    else
                        in=1;
                    end
                    for t=1:Block
                        if t>1
                            if vw(t-1)<BE
                                R(t)=Rgb;
                            elseif (BE<=vw(t-1))&&(vw(t-1)<BE2)
                                R(t)=RRtb(i);
                            elseif BE2<=vw(t-1)
                                R(t)=Rbb;
                            end
                            if R(t)~=R(t-1)
                                swi(t)=abs(Rbb-RRgb(i)); 
                                swic(t)=1;
                            end                    
                        end                
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
                        q(t)=(aa+az(t))/(aa+bb+vz(t));
                        for n=1:N
                            if q(t)*F<=g(t,n) %g(t+ind(r),n)
                                vw(t)=vw(t)+1;
                                gbar(t,n)=g(t,n); %g(t,n);
                            end
                        end
                        gbar3(t)=sum(gbar(t,:));
                        A(t)=gam*min(1,R(t)/vw(t));
                        %A(t)=gam*(1-1/(eps^(R(t)/v(t))));
                        aw(t)=binornd(vw(t),A(t));
                        NCw(t)=rho*R(t)+(lam-1)*gbar3(t)+SWrho*swi(t);
                        NCw2(t)=rho*R(t)+(lam-1)*gbar3(t);
                    end
                    if b<2
                        v=vw;
                        a=aw;
                        gb=gbar3;
                        NC=NCw;
                        NC2=NCw2;
                    else
                        vhold=cat(2,v',vw');
                        ahold=cat(2,a',aw');
                        gbhold=cat(2,gb',gbar3');
                        NChold=cat(2,NC',NCw');
                        NC2hold=cat(2,NC2',NCw2');
                        v=vhold';
                        a=ahold';
                        gb=gbhold';
                        NC=NChold';
                        NC2=NC2hold';
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
                RtbcostMC(m)=mean(NC2);
            end
            Rtbcost(i)=mean(RtbcostMC);
            RRtb(i)
        end
        
    
        %search for Rbb
        RRbb=zeros(NN,1);
        Rbbcost=zeros(NN,1);
        Rbbd=Rbb-X;
        Rbbu=Rbb+X;
        for i=Rbbd:Rbbu
            RRbb(i)=i;
            RbbcostMC=zeros(MC,1);
            for m=1:MC
                c=1;
                b=0;
                cc=0;
                TESTCONV=zeros(MM,1);
                TESTCONV(1)=1;
                count1=0;
                count2=0;
                while ((cc<1) && (b<MM))
                    g=normrnd(mu,sig,Block,N);
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
                    b=b+1;
                    if b<2
                        in=Z+1;
                    else
                        in=1;
                    end
                    for t=1:Block
                        if t>1
                            if vw(t-1)<BE
                                R(t)=Rgb;
                            elseif (BE<=vw(t-1))&&(vw(t-1)<BE2)
                                R(t)=Rtb;
                                count1=count1+1;
                            elseif BE2<=vw(t-1)
                                R(t)=RRbb(i);
                                count2=count2+1;
                            end
                            if R(t)~=R(t-1)
                                swi(t)=abs(Rbb-RRgb(i)); 
                                swic(t)=1;
                            end                    
                        end                
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
                        q(t)=(aa+az(t))/(aa+bb+vz(t));
                        for n=1:N
                            if q(t)*F<=g(t,n) %g(t+ind(r),n)
                                vw(t)=vw(t)+1;
                                gbar(t,n)=g(t,n); %g(t,n);
                            end
                        end
                        gbar3(t)=sum(gbar(t,:));
                        A(t)=gam*min(1,R(t)/vw(t));
                        %A(t)=gam*(1-1/(eps^(R(t)/v(t))));
                        aw(t)=binornd(vw(t),A(t));
                        NCw(t)=rho*R(t)+(lam-1)*gbar3(t)+SWrho*swi(t);
                        NCw2(t)=rho*R(t)+(lam-1)*gbar3(t);
                    end
                    if b<2
                        v=vw;
                        a=aw;
                        gb=gbar3;
                        NC=NCw;
                        NC2=NCw2;
                    else
                        vhold=cat(2,v',vw');
                        ahold=cat(2,a',aw');
                        gbhold=cat(2,gb',gbar3');
                        NChold=cat(2,NC',NCw');
                        NC2hold=cat(2,NC2',NCw2');
                        v=vhold';
                        a=ahold';
                        gb=gbhold';
                        NC=NChold';
                        NC2=NC2hold';
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
                RbbcostMC(m)=mean(NC2);
            end
            Rbbcost(i)=mean(RbbcostMC);
            RRbb(i)
        end
        
[M1, I1]=min(Bcost(BEd:BEu));
[M2, I2]=min(Bcost2(BE2d:BE2u));
[Mg, Ig]=min(Rgbcost(Rgbd:Rgbu));
[Mt, It]=min(Rtbcost(Rtbd:Rtbu));
[Mb, Ib]=min(Rbbcost(Rbbd:Rbbu));

[Mx]=max(Rgbcost(Rgbd:Rgbu));


figure % Figure3.4
tile=tiledlayout(4,6);
tile.Padding='none';
tile.TileSpacing='none';
nexttile([2 3])
hold on
box on
plot(BEd:BEu,Bcost(BEd:BEu),'Color','k','LineStyle','-','LineWidth',2)
scatter((BEd+I1-1),Bcost(BEd+I1-1),50,'filled','MarkerFaceColor','k')
scatter((34),129.79,175,'*','MarkerEdgeColor','k')
%scatter(34,130.33,50,'d','MarkerEdgeColor','k')
title('Panel 1: Cost of refined crackdown varying BE around RCD#1')
legend('cost variation','minimum cost','cost RCD#1')
legend boxoff
xlabel('BE')
ylim([Mg-5 Mx+1])
xlim([BEd BEu])
hold off
nexttile([2 3])
hold on
box on
plot(Rgbd:Rgbu,Rgbcost(Rgbd:Rgbu),'Color','k','LineStyle','-','LineWidth',2)
scatter((Rgbd+Ig-1),Rgbcost(Rgbd+Ig-1),50,'filled','MarkerFaceColor','k')
scatter((30),129.79,175,'*','MarkerEdgeColor','k')
%scatter(30,130.33,50,'d','MarkerEdgeColor','k')
title('Panel 2: Cost of refined crackdown varying R_{GB} around RCD#1')
legend('cost variation','minimum cost','cost RCD#1')
legend boxoff
xlabel('R_{GB}')
ylim([Mg-5 Mx+1])
xlim([21 39])
hold off
nexttile([2 2])
hold on
box on
plot(BE2d:BE2u,Bcost2(BE2d:BE2u),'Color','k','LineStyle','-','LineWidth',2)
scatter((BE2d+I2-1),Bcost2(BE2d+I2-1),50,'filled','MarkerFaceColor','k')
scatter((44),129.79,150,'*','MarkerEdgeColor','k')
title('Panel 3: Cost of refined crackdown varying SBE around RCD#1')
legend('cost variation','minimum cost','cost RCD#1')
legend boxoff
xlabel('SBE')
ylim([Mg-5 Mx+1])
xlim([35 53])
hold off
nexttile([2 2])
hold on
box on
plot(Rtbd:Rtbu,Rtbcost(Rtbd:Rtbu),'Color','k','LineStyle','-','LineWidth',2)
scatter((Rtbd+It-1),Rtbcost(Rtbd+It-1),50,'filled','MarkerFaceColor','k')
scatter((47),129.79,150,'*','MarkerEdgeColor','k')
title('Panel 4: Cost of refined crackdown varying R_{BB1} around RCD#1')
legend('cost variation','minimum cost','cost RCD#1')
legend boxoff
xlabel('R_{BB1}')
ylim([Mg-5 Mx+1])
xlim([38 56])
hold off
nexttile([2 2])
hold on
box on
plot(Rbbd:Rbbu,Rbbcost(Rbbd:Rbbu),'Color','k','LineStyle','-','LineWidth',2)
scatter((Rbbd+Ib-1),Rbbcost(Rbbd+Ib-1),50,'filled','MarkerFaceColor','k')
scatter((64),129.79,150,'*','MarkerEdgeColor','k')
%scatter(56,130.33,50,'d','MarkerEdgeColor','k')
title('Panel 5: Cost of refined crackdown varying R_{BB2} around RCD#1')
legend('cost variation','minimum cost','cost RCD#1')
legend boxoff
xlabel('R_{BB2}')
ylim([Mg-5 Mx+1])
xlim([55 73])
hold off

handaxes1=axes('position',[0.07 0.72 0.22 0.22]);
hold on
box on
plot(BEd:BEu,Bcost(BEd:BEu),'Color','k','LineStyle','-')
scatter((BEd+I1-1),Bcost(BEd+I1-1),50,'filled','MarkerFaceColor','k')
scatter((34),129.79,50,'*','MarkerEdgeColor','k')
%scatter(34,130.33,50,'d','MarkerEdgeColor','k')
plot(BEd:BEu,ylin(BEd:BEu),'Color','k','LineStyle',':');
legend('cost variation','minimum cost','cost RCD#1','cost CD#1','Location','southeast')
legend boxoff
xlabel('BE')
ylim([129 131])
xlim([BEd BEu])
hold off

handaxes2=axes('position',[0.57 0.72 0.22 0.22]);
hold on
box on
plot(Rgbd:Rgbu,Rgbcost(Rgbd:Rgbu),'Color','k','LineStyle','-')
scatter((Rgbd+Ig-1),Rgbcost(Rgbd+Ig-1),50,'filled','MarkerFaceColor','k')
scatter((30),129.79,150,'*','MarkerEdgeColor','k')
%scatter(30,130.33,50,'d','MarkerEdgeColor','k')
plot(Rgbd:Rgbu,ylin(Rgbd:Rgbu),'Color','k','LineStyle',':');
legend('cost variation','minimum cost','cost RCD#1','cost CD#1','Location','southeast')
legend boxoff
xlabel('R_{GB}')
ylim([129 131])
xlim([Rgbd Rgbu])
hold off

handaxes3=axes('position',[0.04 0.24 0.2 0.2]);
hold on
box on
plot(BE2d:BE2u,Bcost2(BE2d:BE2u),'Color','k','LineStyle','-')
scatter((BE2d+I2-1),Bcost2(BE2d+I2-1),50,'filled','MarkerFaceColor','k')
scatter((44),129.79,150,'*','MarkerEdgeColor','k')
plot(BE2d:BE2u,ylin(BE2d:BE2u),'Color','k','LineStyle',':');
legend('cost variation','minimum cost','cost RCD#1','cost CD#1','Location','southeast')
legend boxoff
xlabel('SBE')
ylim([129 131])
xlim([BE2d BE2u])
hold off

handaxes4=axes('position',[0.38 0.24 0.2 0.2]);
hold on
box on
plot(Rtbd:Rtbu,Rtbcost(Rtbd:Rtbu),'Color','k','LineStyle','-')
scatter((Rtbd+It-1),Rtbcost(Rtbd+It-1),50,'filled','MarkerFaceColor','k')
scatter((47),129.79,50,'*','MarkerEdgeColor','k')
plot(Rtbd:Rtbu,ylin(Rtbd:Rtbu),'Color','k','LineStyle',':');
legend('cost variation','minimum cost','cost RCD#1','cost CD#1','Location','southeast')
legend boxoff
xlabel('R_{BB1}')
ylim([129 131])
xlim([Rtbd Rtbu])
hold off

handaxes5=axes('position',[0.71 0.24 0.2 0.2]);
hold on
box on
plot(Rbbd:Rbbu,Rbbcost(Rbbd:Rbbu),'Color','k','LineStyle','-')
scatter((Rbbd+Ib-1),Rbbcost(Rbbd+Ib-1),50,'filled','MarkerFaceColor','k')
scatter((64),129.79,50,'*','MarkerEdgeColor','k')
%scatter(56,130.33,50,'d','MarkerEdgeColor','k')
plot(Rbbd:Rbbu,ylin(Rbbd:Rbbu),'Color','k','LineStyle',':');
legend('cost variation','minimum cost','cost RCD#1','cost CD#1','Location','southeast')
legend boxoff
xlabel('R_{BB2}')
ylim([129 131])
xlim([Rbbd Rbbu])
hold off

