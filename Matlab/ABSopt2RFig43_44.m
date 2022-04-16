clear

load costRgb; %Generated from ASBcon2RV6opt.m
load costRbb; %Generated from ASBcon2RV7opt.m
load costb;   %Generated from ASBcon2RV8opt.m
load Ev;      %Generated from ASBcon2RV8opt.m
load Ea;      %Generated from ASBcon2RV8opt.m

dim=71;
D=11;
b=36;Rgb=31;Rbb=58; %CD (BE,RGB,RBB)=(35,30,57)
bl=32;Rgbl=31;Rbbl=27;

%cRgb=zeros(NN,1);cRbb=zeros(NN,1);cbgb=zeros(NN,1);
ccRgb=zeros(D,D,1);ccRbb=zeros(D,D,1);ccb=zeros(D,D,1);
YY=zeros(dim,1);XX=zeros(dim,1);XXb=zeros(dim,1);YYbb=zeros(dim,1);YYgb=zeros(dim,1);
ccostb=NaN(dim,1);Arate=zeros(dim,dim,1);Rcap=zeros(dim,dim,1);AArate=zeros(dim,1);RRcap=zeros(dim,1);
for i=1:dim
    YY(i)=i;
    XX(i)=i;
    XXb(i)=36;
    YYbb(i)=58;
    YYgb(i)=31;
    for j=1:dim
        if i>1
            Rcap(i,j)=max(0,((i-1)-Ev(i,j))/(i-1));
        end
        Arate(i,j)=Ea(i,j)/Ev(i,j);
        if i==j
            ccostb(i)=costb(i,j);
            RRcap(i)=Rcap(i,j);
            AArate(i)=Arate(i,j);
        end        
    end
end
%for i=b-5:b+5
%    for j=Rbb-5:Rbb+5
%        cbgb(j)=costRgb(b,j);
%        cRgb(i)=costRgb(i,Rbb);
%    end
%end

figure
tile=tiledlayout(6,3);
tile.Padding='none';
tile.TileSpacing='none';
nexttile([2 3])
%mesh(costRgb)
mesh(costRgb,'EdgeColor',[0.5 0.5 0.5],'FaceColor',[0.9 0.9 0.9])
%mesh(Rbb-Rbbl:Rbb+Rbbl,b-bl:b+bl,costRgb(b-bl:b+bl,Rbb-Rbbl:Rbb+Rbbl),'EdgeColor',[0.5 0.5 0.5],'FaceColor',[0.9 0.9 0.9])
hold on
box on
plot3(XX(1:71),XXb(1:71),costRgb(b,1:71)','Color','k','LineStyle','--','LineWidth',2.5)
plot3(YYbb(1:71),YY(1:71),costRgb(1:71,Rbb),'Color','k','LineStyle',':','LineWidth',2.5)
scatter3(58,36,costRgb(36,58),200,'^','filled','MarkerFaceColor','k')
ylabel('Bin edge')
xlabel('R_{BB}')
legend('cost','Bin edge=35','R_{BB}=57','optimum CD','Location','Best')
title('Panel 1: cost surface (bin edge, R_{BB}|R_{GB} = 30)')
%zlim([130 140])
xlim([1 71])
ylim([1 71])
hold off
nexttile([2 3])
mesh(costRbb','EdgeColor',[0.5 0.5 0.5],'FaceColor',[0.9 0.9 0.9])
%mesh(Rgb-Rgbl:Rgb+Rgbl,b-bl:b+bl,costRbb(b-bl:b+bl,Rgb-Rgbl:Rgb+Rgbl),'EdgeColor',[0.5 0.5 0.5],'FaceColor',[0.9 0.9 0.9])
hold on
box on
plot3(XX(1:71),YYgb(1:71),costRbb(1:71,Rgb),'Color','k','LineStyle','--','LineWidth',2.5)
plot3(XXb(1:71),YY(1:71),costRbb(b,1:71)','Color','k','LineStyle',':','LineWidth',2.5)
scatter3(36,31,costRbb(36,31)',200,'^','filled','MarkerFaceColor','k')
xlabel('Bin edge')
ylabel('R_{GB}')
legend('cost','bin edge=35','R_{GB}=30','optimum CD','Location','Best')
title('Panel 2: cost surface (Bin edge, R_{GB}|R_{BB} = 57)')
%zlim([130 140])
xlim([1 71])
ylim([1 71])
hold off
nexttile([2 3])
mesh(costb','EdgeColor',[0.5 0.5 0.5],'FaceColor',[0.9 0.9 0.9])
%mesh(1:Rgb+Rgbl,1:Rbb+Rbbl,costb(1:Rbb+Rbbl,Rgb-Rgbl:Rgb+Rgbl),'EdgeColor',[0.5 0.5 0.5],'FaceColor',[0.9 0.9 0.9])
hold on
box on
plot3(YYbb(1:71),YY(1:71),costb(Rbb,1:71),'Color','k','LineStyle',':','LineWidth',2.5)
plot3(YY(1:71),YYgb(1:71),costb(1:71,Rgb)','Color','k','LineStyle','--','LineWidth',2.5)
plot3(XX(1:71),XX(1:71),ccostb(XX(1:71)),'Color','k','LineStyle','-','LineWidth',2.5)
scatter3(58,31,costb(58,31)',200,'^','filled','MarkerEdgeColor','k','MarkerFaceColor','k')
scatter3(44,44,costb(44,44),200,'s','filled','MarkerEdgeColor','k','MarkerFaceColor','k')
xlabel('R_{BB}')
ylabel('R_{GB}')
legend('cost','R_{BB}=57','R_{GB}=30','R_{GB}=R_{BB}=R','optimum CD','optimum passive R*=44','Location','Best')
title('Panel 3: cost surface(R_{GB}, R_{BB}|Bin edge = 35 )')
%zlim([130 140])
xlim([1 71])
ylim([1 71])
hold off

figure
tile=tiledlayout(2,3);
tile.Padding='none';
tile.TileSpacing='none';
nexttile
hold on
box on
plot(XX(1:71),ccostb(XX(1:71)),'Color','k','LineStyle','-','LineWidth',1)
scatter(44,ccostb(44),100,'s','MarkerFaceColor','k','MarkerEdgeColor','k')
xlabel('Deterrence resources (R=R_{GB}=R_{BB})')
xlim([1 71])
ylim([120 320])
hold off
nexttile
hold on
box on
plot(XX(1:Rbb),costb(Rbb,1:Rbb)','Color','k','LineStyle','-','LineWidth',1)
scatter(Rgb-1,costb(Rbb,Rgb-1),100,'^','filled','MarkerFaceColor','k')
title('Panel 1, Cost cross sections, (Bin edge = 35)')
xlabel('Deterrence resources (R_{GB}|R_{BB}=57)')
xlim([1 Rbb])
ylim([120 320])
hold off
nexttile
hold on
box on
plot(YY(Rgb:71),costb(Rgb:71,Rgb),'Color','k','LineStyle','-','LineWidth',1)
scatter(Rbb-1,costb(Rbb-1,Rgb),100,'^','filled','MarkerFaceColor','k')
xlabel('Deterrence resources (R_{BB}|R_{GB}=30)')
xlim([Rgb 71])
ylim([120 320])
hold off
nexttile
hold on
box on
plot(XX(1:71),RRcap(XX(1:71)),'Color','k','LineStyle','-','LineWidth',1)
plot(XX(1:71),AArate(XX(1:71)),'Color','k','LineStyle','--','LineWidth',1)
xlabel('Deterrence resources (R=R_{GB}=R_{BB})')
legend('Reserve capacity','Apprehension rate','Location','northwest')
xlim([1 71])
hold off
nexttile
hold on
box on
plot(XX(1:Rbb),Rcap(Rbb,1:Rbb)','Color','k','LineStyle','-','LineWidth',1)
plot(XX(1:Rbb),Arate(Rbb,1:Rbb)','Color','k','LineStyle','--','LineWidth',1)
title('Panel 2, Efficiecy indexes, (Bin edge = 35)')
xlabel('Deterrence resources (R_{GB}|R_{BB}=57)')
legend('Reserve capacity','Apprehension rate','Location','southeast')
xlim([1 Rbb])
hold off
nexttile
hold on 
box on
plot(YY(Rgb:71),Rcap(Rgb:71,Rgb),'Color','k','LineStyle','-','LineWidth',1)
plot(YY(Rgb:71),Arate(Rgb:71,Rgb),'Color','k','LineStyle','--','LineWidth',1)
xlabel('Deterrence resources (R_{BB}|R_{GB}=30)')
legend('Reserve capacity','Apprehension rate','Location','southeast')
xlim([Rgb 71])
hold off

