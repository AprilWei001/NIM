clear all
close all
clc

outSample = importdata('sldscSummary-1KG.txt');
zDeltaOutSample = outSample.data;
outSample = outSample.textdata(2:end,:);
zH2OutSample = cellfun(@(s) str2num(s), outSample(:,3));
zNimOutSample = cellfun(@(s) str2num(s), outSample(:,5));

inSample = importdata('sldscSummary-ukbb.txt');
zDeltaInSample = inSample.data;
inSample = inSample.textdata(2:end,:);
zH2InSample = cellfun(@(s) str2num(s), inSample(:,3));
zNimInSample = cellfun(@(s) str2num(s), inSample(:,5));

Type = [{'POLY'}, {'COMMON'},{'RARE'},{'HIGH'},{'LOW'},{'ALL'}];

subplot('Position', [0.06 0.2 0.26 0.75])
 plotsldsr(inSample, zDeltaOutSample, zH2InSample, Type)
        ylabel('Z-score ($\hat{\Delta}_{h^2}=0$)','interpreter','latex');
text(-1.5, 30, 'a','FontSize', 35)
legend('S-LDSR (1KG)','S-LDSR (UKBB)','Location','southwest')
subplot('Position', [0.38 0.2 0.26 0.75])
 plotsldsr(inSample, zNimOutSample, zNimInSample, Type)
        ylabel('Z-score ($\hat{h^2}_{NIM} > 0$)','interpreter','latex');
text(-1.5, 4.1, 'b','FontSize', 35)

subplot('Position', [0.71 0.2 0.26 0.75])
  plotsldsr(inSample, zH2OutSample, zNimInSample, Type)
        ylabel('Z-score ($\hat{h^2}=h^2$)','interpreter','latex');
text(-1.5, 17.5, 'c','FontSize', 35)

box on

set(gcf, 'PaperPosition',[0 0 22 7.5])
saveas(1,'benchmarkSLDSR.png');


function plotsldsr(inSample, zDeltaOutSample, zDeltaInSample, Type)
cl = [0.789245865397805,0.928357880410425,0.388893363119036;...
        0.0882555694264008,0.704565048234542,0.0175760009204033;...
        0.323204152415515,0.576445456464425,0.638363818800220;...
        0.107400479072226,0.803947292385852,0.912598957116263;...
        0.553509100942747,0.650131247986751,0.924159759082582;...
        0.974333310640643,0.366923881014511,0.959439465792923...
        ];
        
    hold on
    for i = 0.6:0.1:12.4
        plot([i i],[-2 2],'-k','LineWidth',10,'Color',[0.93 0.91 0.93],'HandleVisibility','off');
    end
    for i = 1:6
    if i == 6;
        index = [1:60];
    else
        index = find(contains(inSample(:,1),Type{i})==1);
    end

    boxplot([zDeltaOutSample(index); zDeltaInSample(index)], [ones(length(index),1)*(2*i-0.3); ones(length(index),1)*(2*i+0.3)],'positions', [2*i-0.3 2*i+0.3], 'symbol', '','Widths', 0.5);
    end

    h = findobj(gca,'Tag','Box');
    set(gca,'FontSize',25)

    grid on

    for j = length(h):-1:1
        if round(j/2) == (j/2)
            patch(get(h(j),'XData'), get(h(j),'YData'), cl(2,:),'FaceAlpha',.5);
        else
            patch(get(h(j),'XData'), get(h(j),'YData'),  cl(end,:),'FaceAlpha',.5);
        end
    end

set(gca,'XTick',[2:2:12]);
xlim([0.5 12.5])
Type{1} = 'BASELINE';
set(gca,'XTickLabel', Type);
set(gca,'FontSize',20);
xtickangle(80)
end



