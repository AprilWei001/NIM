clear all
close all
clc

ldAll = [];
mafAll =[];
AncestryAll =[];
nimID= importdata('nimLD.99.tags');
nimID =sprintfc('%d:%d',nimID);
numMafBin = 5;
numLDBin = 5;
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
mafLDResultAll = outputNOLProduct( mafResultAll,ldResultAll);
ancMafLDResultAll = outputNOLProduct(AncestryAll, mafLDResultAll);
totalEachCol = sum(ancMafLDResultAll)
indRemoved = find(totalEachCol < 30)
ancMafLDResultAll(:, indRemoved) = [];

%outputMatrixFile([AncestryAll,ldResultAll],'tags.98.anc.ld.annot','w');
%outputMatrixFile([AncestryAll,mafResultAll], 'tags.98.anc.maf.annot','w');
%outputMatrixFile([AncestryAll,mafResultAll,ldResultAll],'tags.98.anc.maf.ld.annot','w');
%outputMatrixFile(AncestryAll, 'tags.99.nol.anc.annot','w');
%outputMatrixFile(outputNOLProduct(AncestryAll, mafResultAll), 'tags.99.nol.anc.maf.annot','w');
%outputMatrixFile(outputNOLProduct(AncestryAll, ldResultAll),  'tags.99.nol.anc.ld.annot','w');
outputMatrixFile(ancMafLDResultAll,  'tags.99.nol.anc.maf.ld.annot','w');
outputMatrixFile(indRemoved,'indNearEmptyBins.tags.nol.maf.ld.txt','w')
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
