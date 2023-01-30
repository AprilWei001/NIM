
clear all
close all
clc

cl = [0.789245865397805,0.928357880410425,0.388893363119036;...
        0.0882555694264008,0.704565048234542,0.0175760009204033;...
        0.323204152415515,0.576445456464425,0.638363818800220;...
        0.107400479072226,0.803947292385852,0.912598957116263;...
        0.553509100942747,0.650131247986751,0.924159759082582;...
        0.974333310640643,0.366923881014511,0.959439465792923...
        ];
Data = importdata('InfoAncMafLD.txt');

indNIM = find(Data(:,1) == 1);
indMH  = find(Data(:,1) == 0);

MAF = Data(:,3);
LD = Data(:,4);

%------------cdf----------------
subplot('Position', [0.05 0.16 0.19 0.79]);
h1NIM = cdfplot(MAF(indNIM));
h1NIM.LineWidth = 2;
h1NIM.Color = cl(2,:);
hold on
h1MH = cdfplot(MAF(indMH));

title ''
h1MH.LineWidth = 3;
h1MH.LineStyle = '--';
h1MH.Color = [0.974333310640643,0.366923881014511,0.959439465792923]
xlim([0.001 0.5])
set(gca,'FontSize',25)
xlabel('MAF');
ylabel('Cumulative distribution');
set(gca,'YTick',[0.:0.2:1]);
legend('NIM','MH','Location','southeast')
  set(gca,'linewidth',1)

text(-0.1, 1, 'a', 'FontSize', 35);
subplot('Position', [0.295 0.16 0.19 0.79]);
h2NIM = cdfplot(LD(indNIM));
h2NIM.Color =  cl(2,:)
h2NIM.LineWidth = 2;

hold on
h2MH = cdfplot(LD(indMH));
h2MH.Color = [0.974333310640643,0.366923881014511,0.959439465792923]
h2MH.LineWidth = 3;
h2MH.LineStyle = '--';
set(gca,'YTick',[0.:0.2:1]);
set(gca,'XTick',[200:200:1000]);
title ''
set(gca,'FontSize',25)
xlabel('LD score (UKBB)');
ylabel('Cumulative distribution');
xlim([0 1000])
  set(gca,'linewidth',1)

text(-200, 1, 'b', 'FontSize', 35);

%------------------box------------------
mafBin = prctile( MAF , [0:20:100]); 
ldBin = prctile( LD , [0:20:100]); 


subplot('Position', [0.54 0.16 0.19 0.79]);
AncestryMAF = [];
MAFNew = [];
for i = 1:5;
    indBin = intersect(find(LD <= ldBin(i+1)), find(LD > ldBin(i)));
    indNIMBin = intersect(indBin, indNIM);
    indMHBin = intersect(indBin, indMH);
    AncestryMAF = [AncestryMAF;ones(length(indNIMBin),1)*(i-0.5)];
    AncestryMAF = [AncestryMAF;ones(length(indMHBin),1)*i];
    MAFNew = [MAFNew; MAF(indNIMBin); MAF(indMHBin)];
end

boxplot(MAFNew, AncestryMAF,'positions', [1.7 2.3 3.7 4.3 5.7 6.3 7.7 8.3 9.7 10.3], 'symbol', '','Widths', 0.5);
h = findobj(gca,'Tag','Box');
set(gca,'FontSize',25)
    set(findobj(gca,'type','line'),'linew',1)
  set(gca,'linewidth',1)

for j = 2:2:length(h)
    patch(get(h(j),'XData'), get(h(j),'YData'),  cl(2,:),'FaceAlpha',.5);
end
for j = 1:2:length(h)
    patch(get(h(j),'XData'), get(h(j),'YData'),[0.974333310640643,0.366923881014511,0.959439465792923],'FaceAlpha',.5);
end
set(gca, 'XTick',[2:2:10]);
set(gca,'YTick',[0:0.1:0.5]);
set(gca,'XTickLabel', [{'20'}, {'40'}, {'60'}, {'80'}, {'100'}]);
xlabel('LD score (UKBB) percentile');
ylabel('MAF');
ylim([-0.02 0.5]);

xlim([1 11])
text(-1.23, 0.5, 'c', 'FontSize', 35);
grid on
subplot('Position', [0.785 0.16 0.19 0.79]);
AncestryLD = [];
LDNew = [];
for i = 1:5;
    indBin = intersect(find(MAF <= mafBin(i+1)), find(MAF > mafBin(i)));
    indNIMBin = intersect(indBin, indNIM);
    indMHBin = intersect(indBin, indMH);
    AncestryLD = [AncestryLD;ones(length(indNIMBin),1)*(i-0.5)];
    AncestryLD = [AncestryLD;ones(length(indMHBin),1)*i];
    LDNew = [LDNew; LD(indNIMBin); LD(indMHBin)];
end

boxplot(LDNew, AncestryLD,'positions', [1.7 2.3 3.7 4.3 5.7 6.3 7.7 8.3 9.7 10.3], 'symbol', '','Widths', 0.5);
h = findobj(gca,'Tag','Box');
set(gca,'FontSize',25)
    set(findobj(gca,'type','line'),'linew',1)
  set(gca,'linewidth',1)

for j = 2:2:length(h)
    patch(get(h(j),'XData'), get(h(j),'YData'), cl(2,:),'FaceAlpha',.5);
end
for j = 1:2:length(h)
    patch(get(h(j),'XData'), get(h(j),'YData'), [0.974333310640643,0.366923881014511,0.959439465792923],'FaceAlpha',.5);
end
grid on
set(gca, 'XTick',[2:2:10]);
set(gca,'XTickLabel', [{'20'}, {'40'}, {'60'}, {'80'}, {'100'}]);
xlabel('MAF percentile');
set(gca,'YTick',[0:100:500]);
ylabel('LD score (UKBB)');
ylim([-20 500]);
xlim([1 11])

text(-1.45, 500, 'd', 'FontSize',35);

set(gcf,'PaperPosition',[0 0 26 6])
saveas(1,'ancMafLd.png');
