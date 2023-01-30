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
TYPE = [{'BASELINE'},{'COMMON'},{'RARE'},{'HIGH'},{'LOW'},{'ALL'}];
GWAS = importdata('gwasStat.10.txt');
Susie = importdata('susieStat.10.txt');
GWASType = [];
SuSie = [];


FDPGWAS = GWAS.data(:,3);
FDPSuSie = Susie.data(:,end);
FDPAll = [];
typeAll = [];
subplot('position',  [0.05 , 0.61, 0.35, 0.34])

for i = 1:6
    if i ~=6
        index = find(contains(GWAS.textdata(:,1), TYPE{i})==1);
    else
        index = [1:length(GWAS.data)];
    end
hold on
    boxplot([FDPGWAS(index); FDPSuSie(index)], [ones(length(index),1)*(2*i-0.3); ones(length(index),1)*(2*i+0.3)],'positions', [2*i-0.3 2*i+0.3], 'symbol', '','Widths', 0.5);
    FDPAll = [FDPAll;mean(FDPGWAS(index)), mean(FDPSuSie(index) )];
end
set(findobj(gca,'type','line'),'linew',1)
set(gca,'linewidth',1)


%boxplot(FDPAll, typeAll,'positions', [1.7 2.3 3.7 4.3 5.7 6.3 7.7 8.3 9.7 10.3 11.7 12.3], 'symbol', '','Widths', 0.5);
h = findobj(gca,'Tag','Box');
set(gca,'FontSize',25)

grid on

for j = 1:length(h)
    if round(j/2) == (j/2)
        patch(get(h(j),'XData'), get(h(j),'YData'), cl(2,:),'FaceAlpha',.5);
    else
        patch(get(h(j),'XData'), get(h(j),'YData'),  cl(end,:),'FaceAlpha',.5);
    end
end

legend('Fine mappping','Association testing')

set(gca,'XTick',[2:2:12]);
set(gca,'XTickLabel', TYPE);

ylim([-0.05 0.8])
ylabel('FDP')
text(-1.6, .8,'b','FontSize', 35)
set(gca,'FontSize',20);
xtickangle(80)

finemapping = importdata('credibleNIMRegion.txt');


lenFinemapping = finemapping.data(:,3);
numCredibleNIMs = finemapping.data(:,1);
numNIMs = finemapping.data(:,2);

subplot('Position',  [0.47 , 0.63, 0.22, 0.32] )
histogram(lenFinemapping/1000, 14,'FaceColor',cl(end-1,:),'FaceAlpha', 0.5, 'EdgeColor', cl(end-1,:),'EdgeAlpha',0.5);
set(gca,'FontSize',20);
text(-130, 35,'c','FontSize', 35)
set(gca,'linewidth',1)

xlabel('Credible NIM regions (kb)')
ylabel('# of credible NIM regions')

grid on
subplot('Position',  [0.75 , 0.63, 0.22, 0.32] )

histogram(numCredibleNIMs./numNIMs, [0.1:0.05:1],'FaceColor',cl(end-2,:),'FaceAlpha', 0.5, 'EdgeColor', cl(end-2,:),'EdgeAlpha',0.5);
xlabel('# of credible NIMs/# of NIMs')
ylabel('# of credible NIM regions');
set(gca,'FontSize',20);
grid on
text(-.1, 10,'d','FontSize', 35)
set(gca,'linewidth',1)

rhemc = importdata('summaryRhemc.anc.maf.ld.txt');
pheno = rhemc.textdata(:,1);
totalStat =  cellfun(@(s) strsplit(s,' ('),rhemc.textdata(:,4),'UniformOutput',false);
nimStat = cellfun(@(s) strsplit(s,' ('),rhemc.textdata(:,6),'UniformOutput',false);
H_nim = cellfun(@(s) str2num(s), cellfun(@(s) s(1),nimStat));
se_nim = cellfun(@(s) str2num(s(1:end-1))^2, cellfun(@(s) s(2),nimStat));

H2 = cellfun(@(s) str2num(s), cellfun(@(s) s(1),totalStat));
seH2 = cellfun(@(s) str2num(s(1:end-1))^2, cellfun(@(s) s(2),totalStat));



finemappingPheno = finemapping.textdata(:,1);

L = zeros(length(pheno),1);

for i = 1:length(pheno);
    ind = find(strcmp(finemappingPheno,pheno{i})==1);
    L(i) = length(ind);
end

subplot('Position', [0.06 0.1 0.26 0.34])
histogram(L, 10,'FaceColor',cl(end,:),'FaceAlpha', 0.5, 'EdgeColor', cl(end,:),'EdgeAlpha',0.5);
grid on 
box on
text(-1.8, 50,'e','FontSize', 35)
set(gca,'linewidth',1)

xlabel('# of credible NIM regions')
ylabel('# of phenotypes')
set(gca,'FontSize',20);
subplot('Position', [0.385 0.1 0.26 0.34])
scatter(H2, L,100,'MarkerEdgeColor',cl(end-1,:), 'MarkerFaceColor',cl(end-1,:), 'MarkerFaceAlpha', 0.5);
[rho p] = corr(H2, L, 'type','Spearman')
ylabel('# of credible NIM regions')
xlabel('h^2')
set(gca,'FontSize',20);
grid on 
box on
text(-.12, 10,'f','FontSize', 35)
set(gca,'linewidth',1)

subplot('Position', [0.71 0.1 0.26 0.34])
scatter(H_nim, L,100,'MarkerEdgeColor',cl(end-2,:), 'MarkerFaceColor',cl(end-2,:), 'MarkerFaceAlpha', 0.5);
[rho p] = corr(H_nim, L, 'type','Spearman')
grid on 
box on
xlabel('h^2_{NIM}')
ylabel('# of credible NIM regions')
set(gca,'FontSize',20);
text(-6*10^(-3), 10,'g','FontSize', 35)
set(gca,'linewidth',1)

set(gcf, 'PaperPosition',[0 0 21 13])
saveas(1,'benchmarkSusie_covar.png');
mdl=fitlm([H_nim,H2],L)
function plotFDP(gwas, susie, index)

    data = [gwas(index); susie(index);combined(index)];
    
    cl = [0.789245865397805,0.928357880410425,0.388893363119036;...
        0.0882555694264008,0.704565048234542,0.0175760009204033;...
        0.323204152415515,0.576445456464425,0.638363818800220;...
        0.107400479072226,0.803947292385852,0.912598957116263;...
        0.553509100942747,0.650131247986751,0.924159759082582;...
        0.974333310640643,0.366923881014511,0.959439465792923...
        ];
    n = length(index);
    labels = [repmat({'GWAS'}, n,1); repmat({'Fine mapping'}, n,1)];

    boxplot(data, labels,'positions',[1 2 3],'Widths',.5);
    h = findobj(gca,'Tag','Box');
    for j=1:length(h)
        patch(get(h(j),'XData'),get(h(j),'YData'),cl(end+1-j,:),'FaceAlpha',.5);
    end
    
    xlim([0.5 3.5])
    ylim([-0.05 1])
    %set(h,'Linewidth',2);
    hold on 
    grid on
  xtickangle(80)
end