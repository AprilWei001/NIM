clear all
close all
clc 


eqtl = importdata('eqtl.txt');
ind_gtexv8 = find(strcmp(eqtl.textdata(:,2),'GTEx/v8')==1);
eqtl_snp = eqtl.textdata(ind_gtexv8,1);
eqtl_snp = cellfun(@(s) strsplit(s,':'), eqtl_snp, 'UniformOutput', false);
eqtl_snp = cellfun(@(s) strcat(s{1},'.', s{2}), eqtl_snp,'UniformOutput', false);
eqtl_tissue = eqtl.textdata(ind_gtexv8,3);
eqtl_gene = eqtl.textdata(ind_gtexv8,4);

unique_eqtl = unique(eqtl_snp);
num_unique_eqtl = length(unique(eqtl_snp))

fid = fopen('supplementary tables/Data S7.csv');

tline = fgetl(fid);
NIM = [];
while(feof(fid) ~= 1)
    tline = strsplit(fgetl(fid),',');
    NIM = [NIM; tline(6:end)'];
end

NIM = unique(NIM);
NIM = NIM(2:end);
numUniqueNIM = length(NIM);
count = [];
for i = 1:length(NIM);
    nim = NIM{i};
    ind = find(strcmp(eqtl_snp, nim)==1);
    tissues = eqtl_tissue(ind,:);
    genes = eqtl_gene(ind,:);
    count = [count; length(unique(tissues)), length(unique(genes))];
end

result = [];
for i = 0:50
    result = [result; length(find(count(:,1) ==i)), length(find(count(:,2)==i))];
end

scatter(find(result(:,1)>0)-1,result(find(result(:,1)>0),1)/numUniqueNIM,'MarkerFaceAlpha', 0.5,'MarkerFaceColor','b')
hold on
scatter(find(result(:,2)>0)-1,result(find(result(:,2)>0),2)/numUniqueNIM,'MarkerFaceAlpha', 0.5,'MarkerFaceColor','r')
legend('tissue', 'gene');
xlabel('# of tissues or genes a given NIM is an eQTL');
ylabel('Proportion of credible NIMs');
set(gca,'FontSize',20);
set(gcf,'PaperPosition',[0 0 7 7])
saveas(1,'eqtl_stat.png');

unique_tissue = unique(eqtl_tissue);
nimPerTissue = [];
for i = 1:length(unique_tissue)
    ind = find(strcmp(eqtl_tissue,unique_tissue{i})==1);
    eqtlThisTissue = unique(eqtl_snp(ind));
    nimPerTissue = [nimPerTissue, length(intersect(eqtlThisTissue, NIM))];
end

hold off

[nimPerTissue ind] = sort(nimPerTissue);
unique_tissue = unique_tissue(ind);
scatter([1:length(unique_tissue)],nimPerTissue,'MarkerFaceAlpha', 0.5,'MarkerFaceColor','r')
xticks([1:length(unique_tissue)])
unique_tissue = cellfun(@(s) strrep(s, '_',' '), unique_tissue,'UniformOutput',false);
xticklabels(unique_tissue)
ylabel('# of unique credible NIMs');
set(gca,'FontSize',20);
set(gcf,'PaperPosition',[0 0 20 8])
saveas(1,'eqtl_tissue.png');

fid2=fopen('eqtl_tissue.txt','w');
for i = 1:length(unique_tissue)
    fprintf(fid2,'%s\t', unique_tissue{i});
    fprintf(fid2,'%d\n', nimPerTissue(i));
end
fclose(fid2);

unique_gene = unique(eqtl_gene);
nimPerGene = [];
for i = 1:length(unique_gene)
    ind = find(strcmp(eqtl_gene,eqtl_gene{i})==1);
    eqtlThisGene = unique(eqtl_snp(ind));
    nimPerGene = [nimPerGene, length(intersect(eqtlThisGene, NIM))];
end


% 
% 
% subplot(1,2,2)
% result = [];
% for i = 0:52
%     result = [result; length(find(nimPerGene ==i))];
% end
% ind = find(result > 0);
% scatter(ind-1,result,'MarkerFaceAlpha', 0.5,'MarkerFaceColor','r')
% 
% xlabel('# of credible NIMs');
% ylabel('Number of Genes');
% 
% set(gca,'FontSize',20);



