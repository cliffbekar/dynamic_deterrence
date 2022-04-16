% ASBconV2Fig38 costs 
% Figure 3.8

%playing with cost function parameters for converged baseline sime

clear
N=70;
NN=N+1;
X=10000;
F=1;
lam=5;
rho=4;
PR2=0;
chi1=1;

playcost=zeros(NN,NN,1);
diagcost=zeros(NN,1);

%load acn103;
load Egg2;
load ERR2;

for i=1:NN
    for j=1:NN
        playcost(i,j)=rho*ER(i,j)+(lam-1)*Egbar(i,j);
        if playcost(i,j)==0
            playcost(i,j)=NaN;
        end
        if i==j
            diagcost(i)=playcost(i,j);
        end
    end
end
M=min(min(playcost));
[Rb,Rg]=find(playcost==M);
Rstar=min(diagcost);

perdiff=(Rstar-M)/M

