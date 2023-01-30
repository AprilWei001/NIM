clear all
close all
clc

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
 
rhemc = importdata('summaryRhemc.noNIMPC.anc.maf.ld_order.txt');
group = rhemc.textdata(:,3);
uniqueGroup = unique(group);


h2Stat = cellfun(@(s) strsplit(s,' ('),rhemc.textdata(:,4),'UniformOutput',false);
nimStat = cellfun(@(s) strsplit(s,' ('),rhemc.textdata(:,6),'UniformOutput',false);
deltaStat = cellfun(@(s) strsplit(s,' ('),rhemc.textdata(:,10),'UniformOutput',false);
h2 = cellfun(@(s) str2num(s), cellfun(@(s) s(1),h2Stat));
se_h2 =  cellfun(@(s) str2num(s(1:end-1)), cellfun(@(s) s(2),h2Stat));
H_nim = cellfun(@(s) str2num(s), cellfun(@(s) s(1),nimStat));
se_nim = cellfun(@(s) str2num(s(1:end-1)), cellfun(@(s) s(2),nimStat));
H_delta = cellfun(@(s) -str2num(s), cellfun(@(s) s(1), deltaStat));
se_delta = cellfun(@(s) str2num(s(1:end-1)), cellfun(@(s) s(2), deltaStat));


rhemcAnc = importdata('summaryRhemc.noNIMPC.anc_order.txt');
h2Stat = cellfun(@(s) strsplit(s,' ('),rhemcAnc.textdata(:,4),'UniformOutput',false);
nimStat = cellfun(@(s) strsplit(s,' ('),rhemcAnc.textdata(:,6),'UniformOutput',false);
deltaStat = cellfun(@(s) strsplit(s,' ('),rhemcAnc.textdata(:,10),'UniformOutput',false);
h2Anc = cellfun(@(s) str2num(s), cellfun(@(s) s(1), h2Stat));
se_h2Anc =  cellfun(@(s) str2num(s(1:end-1)), cellfun(@(s) s(2),h2Stat));
H_nimAnc = cellfun(@(s) str2num(s), cellfun(@(s) s(1),nimStat));
se_nimAnc = cellfun(@(s) str2num(s(1:end-1)), cellfun(@(s) s(2),nimStat));
H_deltaAnc = cellfun(@(s) -str2num(s), cellfun(@(s) s(1), deltaStat));
se_deltaAnc = cellfun(@(s) str2num(s(1:end-1)), cellfun(@(s) s(2), deltaStat));
numGroup = length(uniqueGroup);
subplot('Position', [0.01 0.165 0.42 0.75])
hold on
for i = 1:numGroup
    ind = find(strcmp(group, uniqueGroup(i))==1);
    groupSize = length(ind);
    errorbar(h2(ind), h2Anc(ind), se_h2Anc(ind), se_h2Anc(ind), se_h2(ind), se_h2(ind), 's','MarkerSize',10,'MarkerEdgeColor',cl(i,:),'MarkerFaceColor',cl(i,:),'LineWidth',.5,'Color',cl(i,:))
end
grid on
xlabel('$\bf{\hat{h}^2} (Ancestry+MAF+LD)$','Interpreter','latex')
ylabel('$\bf{\hat{h}^2} (Ancestry)$','Interpreter','latex')
set(gca,'FontSize',26)
legend(uniqueGroup,'Location','westoutside')
legend box off
plot([0 0.7], [0 0.7],'--k','LineWidth',1,'HandleVisibility',  'off')
box on
xlim([0 0.7])
ylim([0 0.7])
text(-0.15, 0.7,'a','FontSize',35)
set(gca,'linewidth',1)
subplot('Position', [0.476 0.165 0.22 0.75])
hold on
for i = 1:numGroup
    ind = find(strcmp(group, uniqueGroup(i))==1);
    groupSize = length(ind);
    errorbar(H_nim(ind), H_nimAnc(ind), se_nimAnc(ind), se_nimAnc(ind), se_nim(ind), se_nim(ind), 's','MarkerSize',10,'MarkerEdgeColor',cl(i,:),'MarkerFaceColor',cl(i,:),'LineWidth',.5,'Color',cl(i,:))
end
grid on
xlabel('$\bf{\hat{h}^2_{NIM}} (Ancestry+MAF+LD)$','Interpreter','latex','Position',[-0.00115 -0.0094])
ylabel('$\bf{\hat{h}^2_{NIM}} (Ancestry)$','Interpreter','latex')
set(gca,'FontSize',26)
plot([-0.008 0.008], [-0.008 0.008],'--k','LineWidth',1,'HandleVisibility',  'off')
box on
xlim([-0.008 0.008])
ylim([-0.008 0.008])
text(-0.0105, 0.008, 'b', 'FontSize', 35)
set(gca,'linewidth',1)
subplot('Position', [0.77 0.165 0.22 0.75])
hold on
for i = 1:numGroup
    ind = find(strcmp(group, uniqueGroup(i))==1);
    groupSize = length(ind);
    errorbar(H_delta(ind), H_deltaAnc(ind), se_deltaAnc(ind), se_deltaAnc(ind), se_delta(ind), se_delta(ind), 's','MarkerSize',10,'MarkerEdgeColor',cl(i,:),'MarkerFaceColor',cl(i,:),'LineWidth',.5,'Color',cl(i,:))
end
grid on
xlabel('$\bf{\hat{\Delta}_{h^2_{NIM}}} (Ancestry+MAF+LD)$','Interpreter','latex')
ylabel('$\bf{\hat{\Delta}_{h^2_{NIM }}} (Ancestry)$','Interpreter','latex')
set(gca,'FontSize',26)
plot([-0.025 0.004], [-0.025 0.004],'--k','LineWidth',1,'HandleVisibility',  'off')
box on
xlim([-0.025 0.004])
ylim([-0.025 0.004])
set(gca,'linewidth',1)
text(-0.0328, 0.004,'c','FontSize',35)

set(gcf, 'PaperPosition',[0 0 30 8])
saveas(1,'AncVsAncMafLd_noNIMPC.png');
