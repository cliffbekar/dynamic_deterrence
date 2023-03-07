% ASBconDROPfig4_2.m figure 4.2, Updated: Aug 29, 2022, K.I. Carlaw 

clear

load fg8DROPn.txt
%load oldDROP5.txt
%load mm.txt
%load ss.txt
%load ms.txt

X=21;Y=14;
F=1;
gam=0.8;
CritGT=zeros(X,Y);
mmu=zeros(X,1);
ssig=zeros(Y,1);
mnu(1)=0;
MS=zeros(35,1); ms=zeros(35,1);
for i=1:X
    if i>1
        mmu(i)=mmu(i-1)+0.05*F;
    end
    ssig(1)=0.1;
    for j=1:Y
        if j>1
            ssig(j)=ssig(j-1)+0.05*F;
        end
        CritGT(i,j)=0.35;
    end
end
for i=1:35
    MS(i)=0.35; 
end

mm=[0.0000 0.0500 0.1000 0.1500 0.2000 0.2500 .3000 0.3500 0.4000 0.4500 0.5000 0.5500 0.6000 0.6260 0.6520 0.6700 0.6850 0.7100 0.7200 0.730 0.7350 0.7475 0.76 0.772 0.783 0.792 0.7985 0.8100 0.8150 0.8230 0.8250 0.8320 0.8330 0.8335 0.8350];
sss=[0.4510 0.4667 0.4838 0.4897 0.4940 0.4965 0.4900 0.4810 0.4720 0.4640 0.4445 0.4320 0.4040 0.3885 0.3760 0.3610 0.3550 0.3385 0.3320 0.3270 0.3145 0.3020 0.2875 0.2770 0.263 0.2520 0.2395 0.2270 0.2170 0.2020 0.1910 0.1810 0.1690 0.1590 0.1520];
figure
mesh(ssig,mmu,fg8DROPn,'EdgeColor',[0 0 0],'FaceColor','none');
hold on
mesh(ssig,mmu,CritGT,'EdgeColor',[0.85 0.85 0.85],'FaceColor',[0.95 0.95 0.95]);
plot3(sss,mm,ms,'Color','k','Linewidth',2,'LineStyle','--')
plot3(sss,mm,MS,'Color','k','Linewidth',2)
%title('Panel 1: TDROP','FontSize',16)
%xticks([2 3 4 5 6 7])
xr1={num2str(ssig(2)),num2str(ssig(3)),num2str(ssig(4)),num2str(ssig(5)),num2str(ssig(6)),num2str(ssig(7)),num2str(ssig(8)),num2str(ssig(9)),num2str(ssig(10)),num2str(ssig(11)),num2str(ssig(12)),num2str(ssig(13)),num2str(ssig(14))};
xr2={num2str(ssig(2)/(F*gam)), num2str(ssig(3)/(F*gam)), num2str(ssig(4)/(F*gam)) num2str(ssig(5)/(F*gam)), num2str(ssig(6)/(F*gam)), num2str(ssig(7)/(F*gam)), num2str(ssig(8)/(F*gam)), num2str(ssig(9)/(F*gam)), num2str(ssig(10)/(F*gam)), num2str(ssig(11)/(F*gam)), num2str(ssig(12)/(F*gam)), num2str(ssig(13)/(F*gam)), num2str(ssig(14)/(F*gam))};
xlab=[xr1; xr2];
%xlabs = strtrim(sprintf('%s\\newline%s\\newline%s\n', xlab{:}));
xlabs = strtrim(sprintf('%s\\newline%s\n', xlab{:}));
ax = gca(); 
%ax.XTick = [2 3 4 5 6 7]; 
%ax.XLim = [2 Y];
ax.XTickLabel = xlabs;
ax.TickLabelInterpreter = 'tex';
ax.FontSize = 14;
%yticks([1 3 5 7 9 11 13 15 17 19 21])
yr1={num2str(mmu(1)),num2str(mmu(3)),num2str(mmu(5)),num2str(mmu(7)),num2str(mmu(9)),num2str(mmu(11)),num2str(mmu(13)),num2str(mmu(15)),num2str(mmu(17)),num2str(mmu(19)),num2str(mmu(21))};
yr2={num2str(mmu(1)/(F*gam)),num2str(mmu(3)/(F*gam)),num2str(mmu(5)/(F*gam)),num2str(mmu(7)/(F*gam)),num2str(mmu(9)/(F*gam)),num2str(mmu(11)/(F*gam)),num2str(mmu(13)/(F*gam)),num2str(mmu(15)/(F*gam)),num2str(mmu(17)/(F*gam)),num2str(mmu(19)/(F*gam)),num2str(mmu(21)/(F*gam))};
ylab=[yr1; yr2];
%ylabs = strtrim(sprintf('%s\\newline%s\\newline%s\n', ylab{:}));
ylabs = strtrim(sprintf('%s\\newline%s\n', ylab{:}));
ax = gca(); 
%ax.YTick = [1 3 5 7 9 11 13 15 17 19 21]; 
%ax.YLim = [1 X];
ax.YTickLabel = ylabs;
ax.TickLabelInterpreter = 'tex';
ax.FontSize = 14;
xlabel('${\sigma}$ $\frac{\sigma}{F \cdot {\gamma}}$','Interpreter','latex','FontSize',20)
ylabel('${\mu}$ $\frac{\mu}{F \cdot {\gamma}}$','Interpreter','latex','FontSize',20)
zlabel('DROP')
xlim([0.15 0.55])
ylim([0 1])
zlim([0 1])
box on
hold off


%figure %Figure 3.9
%tile=tiledlayout(1,2);
%tile.Padding='none';
%tile.TileSpacing='none';
%nexttile
%mesh(fg8DROPn,'EdgeColor',[0 0 0],'FaceColor','none');
%hold on
%mesh(CritGT,'EdgeColor',[0.85 0.85 0.85],'FaceColor',[0.95 0.95 0.95]);
%%scatter3(2,17,GenT(17,2)',200,'o','filled','MarkerEdgeColor','k','MarkerFaceColor','k')
%title('Panel 1: TDROP','FontSize',16)
%%xticks([2 3 4 5 6 7])
%xr1={num2str(ssig(2)),num2str(ssig(3)),num2str(ssig(4)),num2str(ssig(5)),num2str(ssig(6)),num2str(ssig(7))};
%xr2={num2str(ssig(2)/(F*gam)), num2str(ssig(3)/(F*gam)), num2str(ssig(4)/(F*gam)) num2str(ssig(5)/(F*gam)), num2str(ssig(6)/(F*gam)), num2str(ssig(7)/(F*gam))};
%xlab=[xr1; xr2];
%%xlabs = strtrim(sprintf('%s\\newline%s\\newline%s\n', xlab{:}));
%xlabs = strtrim(sprintf('%s\\newline%s\n', xlab{:}));
%ax = gca(); 
%ax.XTick = [2 3 4 5 6 7]; 
%ax.XLim = [2 Y];
%ax.XTickLabel = xlabs;
%ax.TickLabelInterpreter = 'tex';
%yticks([1 3 5 7 9 11 13 15 17 19 21])
%yr1={num2str(mmu(1)),num2str(mmu(3)),num2str(mmu(5)),num2str(mmu(7)),num2str(mmu(9)),num2str(mmu(11)),num2str(mmu(13)),num2str(mmu(15)),num2str(mmu(17)),num2str(mmu(19)),num2str(mmu(21))};
%yr2={num2str(mmu(1)/(F*gam)),num2str(mmu(3)/(F*gam)),num2str(mmu(5)/(F*gam)),num2str(mmu(7)/(F*gam)),num2str(mmu(9)/(F*gam)),num2str(mmu(11)/(F*gam)),num2str(mmu(13)/(F*gam)),num2str(mmu(15)/(F*gam)),num2str(mmu(17)/(F*gam)),num2str(mmu(19)/(F*gam)),num2str(mmu(21)/(F*gam))};
%ylab=[yr1; yr2];
%%ylabs = strtrim(sprintf('%s\\newline%s\\newline%s\n', ylab{:}));
%ylabs = strtrim(sprintf('%s\\newline%s\n', ylab{:}));
%ax = gca(); 
%ax.YTick = [1 3 5 7 9 11 13 15 17 19 21]; 
%ax.YLim = [1 X];
%ax.YTickLabel = ylabs;
%ax.TickLabelInterpreter = 'tex';
%xlabel('${\sigma}, \frac{\sigma}{F \cdot {\gamma}}$','Interpreter','latex','FontSize',16)
%ylabel('${\mu}, \frac{\mu}{F \cdot {\gamma}}$','Interpreter','latex','FontSize',16)
%zlabel('Value of DROP')
%%xlim([2 Y])
%%ylim([1 X])
%box on
%hold off
%nexttile
%mesh(oldDROP5,'EdgeColor',[0 0 0],'FaceColor','none');
%hold on
%mesh(CritGT,'EdgeColor',[0.85 0.85 0.85],'FaceColor',[0.95 0.95 0.95]);
%%scatter3(2,17,GenT(17,2)',200,'o','filled','MarkerEdgeColor','k','MarkerFaceColor','k')
%title('Panel 2: RDROP','FontSize',16)
%xr1={num2str(ssig(2)),num2str(ssig(3)),num2str(ssig(4)),num2str(ssig(5)),num2str(ssig(6)),num2str(ssig(7))};
%xr2={num2str(ssig(2)/(F*gam)), num2str(ssig(3)/(F*gam)), num2str(ssig(4)/(F*gam)) num2str(ssig(5)/(F*gam)), num2str(ssig(6)/(F*gam)), num2str(ssig(7)/(F*gam))};
%xlab=[xr1; xr2];
%%xlabs = strtrim(sprintf('%s\\newline%s\\newline%s\n', xlab{:}));
%xlabs = strtrim(sprintf('%s\\newline%s\n', xlab{:}));
%ax = gca(); 
%ax.XTick = [2 3 4 5 6 7]; 
%ax.XLim = [1 Y];
%ax.XTickLabel = xlabs;
%ax.TickLabelInterpreter = 'tex';
%yticks([1 3 5 7 9 11 13 15 17 19 21])
%yr1={num2str(mmu(1)),num2str(mmu(3)),num2str(mmu(5)),num2str(mmu(7)),num2str(mmu(9)),num2str(mmu(11)),num2str(mmu(13)),num2str(mmu(15)),num2str(mmu(17)),num2str(mmu(19)),num2str(mmu(21))};
%yr2={num2str(mmu(1)/(F*gam)),num2str(mmu(3)/(F*gam)),num2str(mmu(5)/(F*gam)),num2str(mmu(7)/(F*gam)),num2str(mmu(9)/(F*gam)),num2str(mmu(11)/(F*gam)),num2str(mmu(13)/(F*gam)),num2str(mmu(15)/(F*gam)),num2str(mmu(17)/(F*gam)),num2str(mmu(19)/(F*gam)),num2str(mmu(21)/(F*gam))};
%ylab=[yr1; yr2];
%%ylabs = strtrim(sprintf('%s\\newline%s\\newline%s\n', ylab{:}));
%ylabs = strtrim(sprintf('%s\\newline%s\n', ylab{:}));
%ax = gca(); 
%ax.YTick = [1 3 5 7 9 11 13 15 17 19 21]; 
%ax.YLim = [1 X];
%ax.YTickLabel = ylabs;
%ax.TickLabelInterpreter = 'tex';
%xlabel('${\sigma}, \frac{\sigma}{F \cdot {\gamma}}$','Interpreter','latex','FontSize',16)
%ylabel('${\mu}, \frac{\mu}{F \cdot {\gamma}}$','Interpreter','latex','FontSize',16)
%zlabel('Value of DROP')
%xlim([2 Y])
%ylim([1 X])
%zlim([0 1])
%box on
%hold off
