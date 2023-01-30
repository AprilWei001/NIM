clear all
close all
clc

inSample = importdata('sldscSummary-ukbb.txt');
zDeltaInSample = inSample.data;
inSample = inSample.textdata(2:end,:);
zH2InSample = cellfun(@(s) str2num(s), inSample(:,3));
zNimInSample = cellfun(@(s) str2num(s), inSample(:,5));


rhemc = importdata('summaryRhemc.txt');
zDeltaRhemc = rhemc.data;
rhemc = rhemc.textdata(2:end,:);
zH2Rhemc = cellfun(@(s) str2num(s), rhemc(:,5));
zNimRhemc = cellfun(@(s) str2num(s), rhemc(:,7));

indAnnot = [find(strcmp(rhemc(:,1),'expanded.anc')==1),...
    find(strcmp(rhemc(:,1),'expanded.anc.maf')==1),...
    find(strcmp(rhemc(:,1),'expanded.anc.ld')==1),...
    find(strcmp(rhemc(:,1),'expanded.anc.maf.ld')==1)];

Type = [{'POLY'}, {'COMMON'},{'RARE'},{'HIGH'},{'LOW'},{'ALL'}];
ylimRange = [repmat([-7 5],6,1);...
   repmat([-15 10],6,1);];
panels = 'a':'z';
for i = 6:-1:1
    if i == 6;
        index = [1:60];
    else
        index = find(contains(inSample(:,1),Type{i})==1);
    end
    if i ==1
            subplot('position',  [0.062 + (i-1)*0.156, 0.22, 0.133, 0.36])
            text(-1, 10,'b','FontSize',35)
    else
        subplot('position',  [0.065 + (i-1)*0.156, 0.22, 0.132, 0.36])
    end
    set(gca,'linewidth',1)

    plotLDSC(zDeltaRhemc, index, indAnnot);
    if i==1;
        ylabel('Z-score ($\bf{\hat{\Delta}_{h^2}=0}$)','interpreter','latex');
    else
        set(gca,'YTickLabel',[]);
    end
    set(gca,'FontSize',20);
    xtickangle(90)
    ylim(ylimRange(6+i,:));
end

for i = 6:-1:1    
    if i == 6;
        index = [1:60];
    else
        index = find(contains(inSample(:,1),Type{i})==1);
    end
    subplot('position',  [0.065 + (i-1)*0.156, 0.75, 0.132, 0.2])
       
    plotLDSC(zNimRhemc, index, indAnnot);
    set(gca,'xticklabel',{[]});
    set(gca,'FontSize',20);
    set(gca,'linewidth',1)

    if i==1;
        text(-1, 5,'a','FontSize',35)
        ylabel('Z-score ($\bf{\hat{h^2}_{NIM} = h^2_{NIM}}$)','interpreter','latex');
    end

    if i==1
        title("BASELINE")
    else
        title(Type{i})
        set(gca,'YTickLabel',[]);
    end
    ylim(ylimRange(i,:));
    set(gca,'xticklabel',{[]});

    xtickangle(90)
    
end

set(gcf, 'PaperPosition',[0 0 22 10.5])
saveas(1,'benchmarkRhemc_v3.png');
function plotLDSC(rhemc, index, indAnnot)
    hold on
    for i = 0.55:0.01:4.45
        plot([i i],[-2 2],'-k','LineWidth',10,'Color',[0.93 0.91 0.93]);
    end
    xlim([0.5 4.5])
    data = [rhemc(indAnnot(index,1));...
        rhemc(indAnnot(index,2)); rhemc(indAnnot(index,3));...
        rhemc(indAnnot(index,4))];
    
    cl = [0.789245865397805,0.928357880410425,0.388893363119036;...
        0.0882555694264008,0.704565048234542,0.0175760009204033;...
        0.553509100942747,0.650131247986751,0.924159759082582;...
        0.974333310640643,0.366923881014511,0.959439465792923...
        ];
    n = length(index);
    labels = [ repmat({'Ancestry'}, n, 1); repmat({'Ancestry+MAF'}, n, 1);...
       repmat({'Ancestry+LD'}, n, 1); repmat({'Ancestry+MAF+LD'}, n, 1)];


    boxplot(data, labels,'positions',[1 2 3 4 ], 'Widths',.5);
    set(findobj(gca,'type','line'),'linew',1)
    h = findobj(gca,'Tag','Box');
    for j=1:length(h)
        patch(get(h(j),'XData'),get(h(j),'YData'),cl(end+1-j,:),'FaceAlpha',.6);
    end
    %set(h,'Linewidth',2);
    hold on 
    grid on
end
