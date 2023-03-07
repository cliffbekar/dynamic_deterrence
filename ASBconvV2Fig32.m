
clear

%% Setup the Import Options and import the data
opts = spreadsheetImportOptions("NumVariables", 5);

% Specify sheet and range
opts.Sheet = "Sheet1";
opts.DataRange = "A2:E502";

% Specify column names and types
opts.VariableNames = ["tic1", "v", "ARATE", "p", "q"];
opts.VariableTypes = ["double", "double", "double", "double", "double"];

% Import the data
tbl = readtable("C:\Users\kcarlaw\Documents\MATLAB\ASB2022final\fig32vect.xls", opts, "UseExcel", false);

%% Convert to output type
tic1 = tbl.tic1;
v = tbl.v;
ARATE = tbl.ARATE;
p = tbl.p;
q = tbl.q;

%% Clear temporary variables
clear opts tbl
N=100;
T=length(tic1);
bin=zeros(T,1);
for i=1:T
    bin(i)=50;
end

%ll=4550;
%ul=5050;
figure 
tile=tiledlayout(2,1);
tile.Padding='none';
tile.TileSpacing='none';
nexttile
hold on
box on
colororder({'k','k'})
yyaxis left
%area(tT(4657:4922),NT(4657:4922),'FaceColor',[0.9 0.9 0.9])
plot(tic1,v,'Color','k','LineStyle','-')
plot(tic1,bin,'Color','k','LineStyle','-.')
ylim([0 N])
%xlim([ T])
ylabel('Number of violations','Color','k')
yyaxis right
plot(tic1,ARATE,'Color','k','LineStyle',"--")
ylim([0 1])
%xlim([ll ul])
ylabel('Apprehension rate','Color','k')
xline(4657,'Color','k','LineStyle',':');
xline(4922,'Color','k','LineStyle',':');
legend('Violations','Bin edge = 50','Apprehension rate','Location','south')
legend boxoff
title('Panel 1: Illustrative path of violations and apprehension rate')
txt1={'{\leftarrow} Bad Bin {\rightarrow}'};
text(4750,.65,txt1);
hold off
nexttile
hold on
box on
%area(tT(4657:4922),pT(4657:4922),'FaceColor',[0.9 0.9 0.9])
plot(tic1,p,'Color','k','LineStyle','-')
plot(tic1,q,'Color','k','LineStyle','--')
%xlim([1 T])
ylim([0 1])
xline(4657,'Color','k','LineStyle',':');
xline(4922,'Color','k','LineStyle',':');
legend('Objective probability (p)','Subjective probability (q)','Location','southeast')
legend boxoff
xlabel('Tics(t)')
title('Panel 2: objective (p) and subjective (q) probabilities of apprehension')
%txt2={'{\leftarrow} Bad Bin {\rightarrow}'};
%text(4750,.65,txt2);
hold off
%handaxes1=axes('position',[0.28 0.55 0.16 0.19]);
%hold on
%box on
%colororder({'k','k'})
%yyaxis left
%plot(4648:4667,vf(4648:4667),'Color','k','LineStyle','-')
%plot(bin,'Color','k','LineStyle','-.')
%ylim([0 N])
%yyaxis right
%plot(4648:4667,ARATEf(4648:4667),'Color','k','LineStyle',"--")
%ylim([0 1])
%legend('Violations','Bin edge','Apprehension rate','Location','southeast')
%legend boxoff
%xlim([4648 4667])
%title('Transition into BB')
%hold off
%handaxes2=axes('position',[0.52 0.55 0.16 0.19]);
%hold on
%box on
%colororder({'k','k'})
%yyaxis left
%plot(4913:4932,vf(4913:4932),'Color','k','LineStyle','-')
%plot(bin,'Color','k','LineStyle','-.')
%ylim([0 N])
%yyaxis right
%plot(4913:4932,ARATEf(4913:4932),'Color','k','LineStyle',"--")
%ylim([0 1])
%legend('Violations','Bin edge','Apprehension rate','Location','southwest')
%legend boxoff
%xlim([4913 4932])
%title('Transtion out of BB')
%hold off
%handaxes3=axes('position',[0.28 0.28 0.16 0.19]);
%hold on
%box on
%plot(4648:4667,pf(4648:4667),'Color','k','LineStyle','-')
%plot(4648:4667,qf(4648:4667),'Color','k','LineStyle',"--")
%legend('Objective (p)','Subjective (q)','Location','southeast')
%legend boxoff
%ylim([0 1])
%xlim([4648 4667])
%title('Transition into BB')
%hold off
%handaxes4=axes('position',[0.52 0.28 0.16 0.19]);
%hold on
%box on
%plot(4913:4932,pf(4913:4932),'Color','k','LineStyle','-')
%plot(4913:4932,qf(4913:4932),'Color','k','LineStyle',"--")
%legend('Objective (p)','Subjective (q)','Location','southeast')
%legend boxoff
%ylim([0 1])
%xlim([4913 4932])
%title('Transition into BB')
%title('Transition out of BB')
%hold off
