clear all
close all
clc


cl =[0.0882555694264008,0.704565048234542,0.0175760009204033;...
    0.974333310640643,0.366923881014511,0.959439465792923];

rng(1);

SNPs = zeros(200,1);
indNIM =[10,25,30,35:45,50,52:54,60:65,70:78,90,100,105, 108,110,115:118,125]+20;

SNPs(indNIM)=1;
indMH = find(SNPs==0);
subplot('Position', [0.01 0.1 0.98 0.88])
hold on
scatter(zeros(length(indMH),1)-0.002,indMH,20,'MarkerEdgeColor',cl(2,:), ...
    'MarkerFaceColor',cl(2,:), 'MarkerFaceAlpha', 0.5)
scatter(zeros(length(indNIM),1)+0.002,indNIM,20,'MarkerEdgeColor',cl(1,:), ...
    'MarkerFaceColor',cl(1,:), 'MarkerFaceAlpha', 0.5)
plot([0 0],[0 200],'-k','HandleVisibility','off');
axis off
xlim([-0.01 1])
ylim([0 200])
legend('MH','NIM','Location','westoutside')
set(gca,'FontSize',20)

annotation('textarrow',[0.14 0.2],[0.6 0.6],'Linewidth',2,'Color',[.8 .8 .8])
annotation('textarrow',[0.24 0.3],[0.6 0.6],'Linewidth',2,'Color',[.8 .8 .8])
annotation('textarrow',[0.48 0.54],[0.6 0.6],'Linewidth',2,'Color',[.8 .8 .8])
annotation('textarrow',[0.595 0.655],[0.6 0.6],'Linewidth',2,'Color',[.8 .8 .8])
annotation('textarrow',[0.71 0.77],[0.6 0.6],'Linewidth',2,'Color',[.8 .8 .8])
annotation('textarrow',[0.8 0.86],[0.6 0.6],'Linewidth',2,'Color',[.8 .8 .8])

text(-.05, -10,'QC-ed SNPs','FontSize',18);
%{
scatter(zeros(length(indNIM),1)+0.1,indNIM,20,'MarkerEdgeColor',cl(1,:), ...
    'MarkerFaceColor',cl(1,:), 'MarkerFaceAlpha', 0.5,'HandleVisibility','off')
plot([0.1 0.1],[20 150],'-k','HandleVisibility','off');


text(0.03, 10,'NIM association','FontSize',18);
grid on
%}
indSigNIM =[38:42,50,52, 60:63,75:78,90,100, 108,125]+20;
scatter(zeros(length(indSigNIM),1)+0.12,indSigNIM,20,'MarkerEdgeColor',cl(1,:), ...
    'MarkerFaceColor',cl(1,:), 'MarkerFaceAlpha', 0.5,'HandleVisibility','off')
plot([0.12 0.12],[50 150],'-k','HandleVisibility','off');


text(0.08, 40,'Significant ','FontSize',18);
text(0.08, 30,'NIMs','FontSize',18);


indPrunedNIM =[63,100, 125]+20;
scatter(zeros(length(indPrunedNIM),1)+0.24,indPrunedNIM,20,'MarkerEdgeColor',cl(1,:), ...
    'MarkerFaceColor',cl(1,:), 'MarkerFaceAlpha', 0.5,'HandleVisibility','off')

for i = 1:length(indPrunedNIM)
    text(0.22, indPrunedNIM(i),num2str(i),'FontSize',18);
end
plot([0.24 0.24],[80 150],'-k','HandleVisibility','off');

text(0.2, 60,'Pruned-in','FontSize',18);
text(0.2, 50,'NIMs','FontSize',18);

credibleSNPs = [30,38,40, 45,54,65,78,90,105, 108,110,115,125,70,93,106,107, 120:122, 128]+20;
mergedNIMs=[];
mergedMHs = [];
causalRegion = [];
for i = 1:3
    index = [indPrunedNIM(i) - 30:indPrunedNIM(i)+ 30];
    indMHi = intersect(index,indMH);
    indNIMi = intersect(index, indNIM);
    scatter(zeros(length(indMHi),1)-0.002+0.34+0.015*i,indMHi,20,'MarkerEdgeColor',cl(2,:), ...
    'MarkerFaceColor',cl(2,:), 'MarkerFaceAlpha', 0.5,'HandleVisibility','off')
    scatter(zeros(length(indNIMi),1)+0.002+0.34+0.015*i,indNIMi,20,'MarkerEdgeColor',cl(1,:), ...
    'MarkerFaceColor',cl(1,:), 'MarkerFaceAlpha', 0.5,'HandleVisibility','off')
    plot([0.34+0.015*i 0.34+0.015*i],[min(index) max(index)],'-k','HandleVisibility','off');
    text(0.34+0.015*i, min(index)-10, num2str(i),'FontSize',18,'HandleVisibility','off');
    
    
    
    plot([0.24 ,0.34+0.015*i],[indPrunedNIM(i),indPrunedNIM(i)+30],'--','Color',[0.8 0.8 0.8],...
        'HandleVisibility','off','LineWidth',1.2);
    plot([0.24 ,0.34+0.015*i],[indPrunedNIM(i),indPrunedNIM(i)-30],'--','Color',[0.8 0.8 0.8],...
        'HandleVisibility','off','LineWidth',1.2);

    
    
    indCredibleMHi = intersect(indMHi, credibleSNPs);
    indCredibleNIMi = intersect(indNIMi, credibleSNPs);
    scatter(zeros(length(indCredibleMHi),1)-0.002+0.48+0.015*i,indCredibleMHi,20,'MarkerEdgeColor',cl(2,:), ...
    'MarkerFaceColor',cl(2,:), 'MarkerFaceAlpha', 0.5,'HandleVisibility','off')
    scatter(zeros(length(indCredibleNIMi),1)+0.002+0.48+0.015*i,indCredibleNIMi,20,'MarkerEdgeColor',cl(1,:), ...
    'MarkerFaceColor',cl(1,:), 'MarkerFaceAlpha', 0.5,'HandleVisibility','off')
    plot([0.48+0.015*i 0.48+0.015*i],[min(index) max(index)],'-k','HandleVisibility','off');
    text(0.48+0.015*i, min(index)-10, num2str(i),'FontSize',18,'HandleVisibility','off');
    
    
    
    if length(indCredibleNIMi) > length(indCredibleMHi)
        scatter(zeros(length(indCredibleMHi),1)-0.002+0.616+0.015*i,indCredibleMHi,20,'MarkerEdgeColor',cl(2,:), ...
        'MarkerFaceColor',cl(2,:), 'MarkerFaceAlpha', 0.5,'HandleVisibility','off')
        scatter(zeros(length(indCredibleNIMi),1)+0.002+0.616+0.015*i,indCredibleNIMi,20,'MarkerEdgeColor',cl(1,:), ...
        'MarkerFaceColor',cl(1,:), 'MarkerFaceAlpha', 0.5,'HandleVisibility','off')
        plot([0.616+0.015*i 0.615+0.015*i],[min(index) max(index)],'-k','HandleVisibility','off');
        text(0.616+0.015*i, min(index)-10, num2str(i),'FontSize',18,'HandleVisibility','off');
        mergedNIMs = [mergedNIMs; indCredibleNIMi];
        mergedMHs = [mergedMHs; indCredibleMHi];
        causalRegion = [causalRegion,min(index) max(index) ];

    end
end


text(0.33, 30,'Input data','FontSize',18);

text(0.33, 20,'to SuSiE','FontSize',18);


text(0.44, 30,'Credible sets','FontSize',18);

text(0.6, 30,'Sets with','FontSize',18);
text(0.6, 20,'majority NIM','FontSize',18);

mergedNIMs = unique(mergedNIMs);
mergedMHs = unique(mergedMHs);

scatter(zeros(length(mergedMHs),1)-0.002+0.76,mergedMHs,20,'MarkerEdgeColor',cl(2,:), ...
        'MarkerFaceColor',cl(2,:), 'MarkerFaceAlpha', 0.5,'HandleVisibility','off')
scatter(zeros(length(mergedNIMs),1)+0.002+0.76,mergedNIMs,20,'MarkerEdgeColor',cl(1,:), ...
        'MarkerFaceColor',cl(1,:), 'MarkerFaceAlpha', 0.5,'HandleVisibility','off')
plot([0.76  .76 ],[min(causalRegion) max(causalRegion)],'-k','HandleVisibility','off');

text(0.73, 30,'Merged','FontSize',18);
text(0.73, 20,'region','FontSize',18);
scatter(zeros(length(mergedNIMs),1)+.86,mergedNIMs,20,'MarkerEdgeColor',cl(1,:), ...
        'MarkerFaceColor',cl(1,:), 'MarkerFaceAlpha', 0.5,'HandleVisibility','off')

text(0.825, 30,'Credible','FontSize',18);
text(0.825, 20,'NIM set','FontSize',18);

plot([0.93  .93 ],[min(mergedNIMs) max(mergedNIMs)],'-k','HandleVisibility','off');
text(0.915, 30,'Credible','FontSize',18);
text(0.915, 20,'NIM region','FontSize',18);

set(gcf,'PaperPosition',[0 0 16 6])

text(-0.12, 190,'a','FontSize' ,25);
saveas(1,'finemappingDemo.png');