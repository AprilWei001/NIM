
clear all
close all
clc

rhemc = importdata('summaryRhemc.ALLEUR_NMATCH_sprime_results.anc.maf.ld_order.txt');
group = rhemc.textdata(:,3);
Zscore = [cellfun(@(s) str2num(s), rhemc.textdata(:,5)),...
    cellfun(@(s) str2num(s), rhemc.textdata(:,7)),...
    -rhemc.data]; 

h2Stat = cellfun(@(s) strsplit(s,' ('),rhemc.textdata(:,4),'UniformOutput',false);
h2 = cellfun(@(s) str2num(s), cellfun(@(s) s(1),h2Stat));
V_h2 = cellfun(@(s) str2num(s(1:end-1))^2, cellfun(@(s) s(2),h2Stat));
se_h2 = sqrt(V_h2);
nimStat = cellfun(@(s) strsplit(s,' ('),rhemc.textdata(:,6),'UniformOutput',false);
deltaStat = cellfun(@(s) strsplit(s,' ('),rhemc.textdata(:,10),'UniformOutput',false);
H_nim = cellfun(@(s) str2num(s), cellfun(@(s) s(1),nimStat));
V_nim = cellfun(@(s) str2num(s(1:end-1))^2, cellfun(@(s) s(2),nimStat));
se_nim = sqrt(V_nim);
H_delta = cellfun(@(s) -str2num(s), cellfun(@(s) s(1), deltaStat));
V_delta = cellfun(@(s) str2num(s(1:end-1))^2, cellfun(@(s) s(2), deltaStat));
se_delta = sqrt(V_delta);
uniqueGroup = unique(group);

 cl = [0.789245865397805,0.928357880410425,0.388893363119036;...
        0.0882555694264008,0.704565048234542,0.0175760009204033;...
        0.323204152415515,0.576445456464425,0.638363818800220;...
        0.107400479072226,0.803947292385852,0.912598957116263;...
        0.553509100942747,0.650131247986751,0.924159759082582;...
        0.974333310640643,0.366923881014511,0.959439465792923;...
     0.420935968034598,0.0934387387412213,0.140428451643645;...
        0.993559850546711,0.123618315194761,0.0430377697526367;...
     0.543690308588324,0.00696233562022619,0.458777110820447;...
     0.764791187248605,0.169158920261980,0.936566194437307;...
     0.739498573250049,0.181336000363487,0.316715248729260;...
     0.240066392352203,0.220052169280255,0.870165389849008;...
     0.774815652954463,0.740184497397902,0.921061068011622;...
     0.875016043356001,0.471282901077927,0.705802776525446;...
     0.0671544921529179,0.978191778470600,0.967201985721912;...
     0.695538516796386,0.0480043593654713,0.0620603770294961;...
     0.738637492237587,0.624330338822513,0.396901337597561;...
     0.855044305486398,0.266827179013910,0.233182937809061;...
     0.566311708206001,0.239931052418018,0.870148264940529];

pos = [0.095, 0.68, 0.88, 0.2 ; 0.095, 0.465, 0.88, 0.2; 
    0.095, 0.255, 0.88, 0.2; 0.095, 0.035, 0.88, 0.2 ];

subplot('Position',[0.095, 0.89, 0.88, 0.1])

for i = 1:length(uniqueGroup) 
  h = text(i-0.5, 0,uniqueGroup{i},'FontSize',24,'Color',cl(i,:));
  set(h,'Rotation',30)
end

xlim([0 i+1])
ylim([0 2])
axis off ; 
set(gca,'linewidth',1)


[metaT metaSE] = randomEffectMeta(H_nim,V_nim)
metaT/metaSE
[metaT metaSE] = randomEffectMeta(H_delta,V_delta)
metaT/metaSE
metaResults = zeros(14,2);
for j = 2:3
    subplot('Position', pos(j+1,:))
        set(gca,'linewidth',1)

    hold on
    numGroup=length(uniqueGroup);
    for k = 0.55:0.05:numGroup+0.45
        plot([k k],[-2 2],'-k','LineWidth',5,'Color',[0.93 0.91 0.93])
    end
    for i = 1:numGroup
        ind = find(strcmp(group, uniqueGroup(i))==1);
        groupSize = length(ind);
        if groupSize== 1;
            scatter(i, Zscore(ind,j),100,'MarkerEdgeColor',cl(i,:), 'MarkerFaceColor',cl(i,:) );
        else 
            xPos = [i-0.15: 0.3/(groupSize-1): i+0.15];
            scatter(xPos, Zscore(ind,j),100,'MarkerEdgeColor',cl(i,:), 'MarkerFaceColor',cl(i,:));
        end
        if length(ind) >=4
            if j ==2
                [metaT metaSE] = randomEffectMeta(H_nim(ind),V_nim(ind));
            else
                [metaT metaSE] = randomEffectMeta(H_delta(ind),V_delta(ind));
            end
            metaZ = metaT/metaSE;
            plot([i-0.2 i+0.2], [metaZ metaZ],'-','LineWidth',3,'Color',cl(i,:));
            metaResults(i,j-1) = metaZ;
        end
    end

    set(gca,'XTick',[1:numGroup]);
    set(gca,'XTickLabel', []);
    set(gca,'FontSize',25);
    if j ==2;
        ylabel('Z-score ($\bf{\hat{h^2}_{NIM} = 0}$)','interpreter','latex');
        text(-.7, 5, 'c','FontSize',40);
        ylim([-5 5])
    else
        ylabel('Z-score ($\bf{\hat{\Delta}_{h^2}=0}$)','interpreter','latex');
        text(-.7, 10, 'd','FontSize',40);
        ylim([-10 10])
    end
    
    xlim([0.5, numGroup+0.5])
    box on
    grid on
end
xlabel('Phenotypes')

subplot('Position', pos(1,:))
hold on
xPos = [0.75:(numGroup-.5)/length(h2):numGroup+0.25];

for i = 1:numGroup
    ind = find(strcmp(group, uniqueGroup(i))==1);
    groupSize = length(ind);
    errorbar(xPos(ind), H_nim(ind), se_nim(ind),'o','MarkerSize',12,'MarkerEdgeColor',...
    cl(i,:),'MarkerFaceColor',cl(i,:),  'LineWidth',1,'Color',cl(i,:)) 
end
xlim([0.5, numGroup+0.5])
set(gca,'XTickLabel', []);
set(gca,'FontSize',25);
box on
grid on
    set(gca,'XTick',[1:numGroup]);
    set(gca,'XTickLabel', []);
ylim([-0.01 0.01])
ylabel('$\bf{\hat{h^2}_{NIM}}$','interpreter','latex');
    set(gca,'linewidth',1)

text(-.7, .01, 'a','FontSize',40);
subplot('Position', pos(2,:))
hold on
    set(gca,'linewidth',1)

for i = 1:numGroup
    ind = find(strcmp(group, uniqueGroup(i))==1);
    groupSize = length(ind);
    errorbar(xPos(ind), H_delta(ind), se_delta(ind),'o','MarkerSize',12,'MarkerEdgeColor',...
    cl(i,:),'MarkerFaceColor',cl(i,:), 'LineWidth',1,'Color',cl(i,:)) 
end
set(gca,'FontSize',25);
xlim([0.5, numGroup+0.5])
ylim([-0.01 0.01])
ylabel('$\bf{\hat{h^2}}$)','interpreter','latex');
box on
grid on
set(gca,'XTick',[1:numGroup]);
set(gca,'XTickLabel', []);
text(-.7, .01, 'b','FontSize',40);
set(gcf,'PaperPosition',[0 0 20 20])
  ylabel('$\bf{\hat{\Delta}_{h^2}}$','interpreter','latex');
saveas(1,'figUKBBRhemc_covar_SprimeEUR_phenolabel.png');

function [metaT metaSE] = randomEffectMeta(T,V)
    W = 1./V;
    meanT = W'*T/sum(W);
    Q = W'*(T - meanT).^2;
    dF = length(V) - 1;
    C = sum(W) - W'*W/sum(W);
    if Q > dF
        tao2 = (Q-dF)/C;
    else
        tao2 = 0;
    end
    Vstar = V+tao2;
    Wstar = 1./Vstar;
    metaT = Wstar'*T/sum(Wstar);
    metaV = 1/sum(Wstar);
    metaSE = sqrt(metaV);
end