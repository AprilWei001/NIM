clear all
close all
clc

inputNIMs = [{'confidentQCedNIMs.txt'}, {'expandedNIM.tags'}];
outputName = [{'confident'},{'expanded'}]

numMafBin = 5;
numLDBin = 5;

for type = 1:2
nimID= importdata(inputNIMs{type});
nimID =sprintfc('%d:%d',nimID);
outname = outputName{type};

ldAll = [];
mafAll =[];
AncestryAll =[];
for chr =1:22
	ld = importdata(strcat('/u/project/sgss/UKBB/data/imp/maf_0.001/ldscores/',num2str(chr),'.score.ld'));
	coor = ld.textdata(2:end,1);
	maf = ld.data(:,3);
	ldScore = ld.data(:,end);
	mafAll = [mafAll;maf];
	ldAll = [ldAll;ldScore];
	[id ind1 ind2] = intersect(coor, nimID);
	row = length(maf);
	Ancestry = [zeros(row,1),ones(row,1)];
	Ancestry(ind1,1) = 1; %NIMs
	Ancestry(ind1,2) = 0; %Humans
	AncestryAll = [AncestryAll; Ancestry];
end

ldResultAll = getQuantileTable(ldAll, numLDBin);
mafResultAll = getQuantileTable(mafAll, numMafBin);

weightAnc = [ones(1,2); 1, 0; 0, sum(AncestryAll(:,1))/sum(AncestryAll(:,2))];
outputMatrixFile(AncestryAll,strcat(outname, '.anc.annot'),'w');
outputMatrixFile(weightAnc,strcat(outname, '.anc.weight'),'w');

outAncMaf = outputNOLProduct(AncestryAll, mafResultAll);
weightAncMaf = [ones(1,10); ones(1,5), zeros(1,5); zeros(1,5), sum(outAncMaf(:, 1:5))./sum(outAncMaf(:, 6:10))];
outputMatrixFile(outAncMaf, strcat(outname, '.anc.maf.annot'),'w');
outputMatrixFile(weightAncMaf,strcat(outname, '.anc.maf.weight'),'w');

outAncLd = outputNOLProduct(AncestryAll, ldResultAll);
weightAncLd = [ones(1,10); ones(1,5), zeros(1,5); zeros(1,5), sum(outAncLd(:, 1:5))./sum(outAncLd(:, 6:10))];
outputMatrixFile(outAncLd,  strcat(outname, '.anc.ld.annot'),'w');
outputMatrixFile(weightAncLd,strcat(outname, '.anc.ld.weight'),'w');

if type == 2
	mafLDResultAll = outputNOLProduct(mafResultAll,ldResultAll);
	outAncMafLd = outputNOLProduct(AncestryAll, mafLDResultAll);
	totalEachCol = sum(outAncMafLd);
	indRemoved = find(totalEachCol < 30);
	weightAncMafLd = [ones(1,50); ones(1,25), zeros(1,25); zeros(1,25), sum(outAncMafLd(:, 1:25))./sum(outAncMafLd(:, 26:50))];
	outAncMafLd(:, indRemoved) = [];
	weightAncMafLd(:, indRemoved) = [];
	outputMatrixFile(indRemoved,'indNearEmptyBins.expanded.anc.maf.ld.txt','w');
	outputMatrixFile(outAncMafLd,  strcat(outname, '.anc.maf.ld.annot'),'w');
	outputMatrixFile(weightAncMafLd,strcat(outname, '.anc.maf.ld.weight'),'w');
end

end

function mat = outputNOLProduct(mat1, mat2)
	mat = [];
	col = size(mat1,2);
	for i = 1:col
		mat = [mat, mat1(:,i).*mat2];
	end
end

function outputMatrixFile(Data,filename,type)
    [row col] = size(Data);
    fid = fopen(filename,type);
    for i = 1:row;
        if (col >1)
            for j = 1:col -1;
                fprintf(fid, '%d ',Data(i,j));
            end
            fprintf(fid,'%d\n',Data(i,j+1));
        else
            fprintf(fid, '%d\n',Data(i,1));
        end
    end
    fclose(fid);
end

function outTable = getQuantileTable(data, numBin)
	bin = quantile(data, [1:numBin]/numBin);
	numRow = size(data,1);
	outTable = zeros(numRow, numBin);
	for i = 1:numBin
		ind = find(data <= bin(i));
		outTable(ind,i) = 1;
	end
	for i = numBin:-1:2
		outTable(:,i) = outTable(:,i) - outTable(:,i-1);
	end
end
