clear all
close all
clc



allNIMs = importdata('../expandedNIM.txt');
allNIMs =sprintfc('%d.%d',allNIMs);
nimDic = containers.Map(allNIMs,[1:length(allNIMs)]);

snpID = importdata('snpID.txt');
snpID = sprintfc('%d.%d', snpID);
allSNPDic = containers.Map(snpID, [1:length(snpID)]);

pcutoff='10';
Pheno = importdata('filename96pheno.txt');
for i = 1:length(Pheno);
	name = Pheno{i};
	if isfile(strcat('susieRange/',name,'.pcutoff.',pcutoff,'.signim.nimrange'))	
		fidOut = fopen(strcat('summary/',name,'.pcutoff.',pcutoff,'.susie.txt'),'w');
		testSegs = importdata(strcat('susieRange/',name,'.pcutoff.',pcutoff,'.signim.nimrange'));	
		total = length(testSegs.data);
		for j = 1:total %this is a temporary fix to deal with the index switch bug. Need to fix later
			chrTest = str2num(testSegs.textdata{j,1});
                	coorEnd = testSegs.data(j);
                	coorStart = strsplit(testSegs.textdata{j,2},'-');
               	 	coorStart = str2num(coorStart{1});
			
			nameSusie = strcat('susie/',name,'.pcutoff.', pcutoff, '.', num2str(j),'.susie'); %this is a temporary fix to deal with the 0 vs 1 index bug
			fid = fopen(nameSusie,'r');
			tline = fgetl(fid);
			CS = [];
			while ischar(tline)
                		if contains(tline, '"X') == 1;
                        		snps = strsplit(tline);
					ind = find(contains(snps, 'X')==1);
					snps = snps(ind);
					snps = cellfun(@(s) s{1}(3:end), cellfun(@(s) strsplit(s, '_'), snps,'UniformOutput', false),'UniformOutput', false);
                			CS = [CS, snps];
				elseif contains(tline, '$coverage')==1;
					tline = fgetl(fid);
					numCSet = length(strsplit(tline)) - 1;
				end
				tline = fgetl(fid);
			end
			fclose(fid);
			
			numSNPSusie  = length(CS);
			numNIMSusie = sum(isKey(nimDic, CS));
			numMHSusie = numSNPSusie - numNIMSusie;

			ind1 = values(allSNPDic, {strcat(num2str(chrTest),'.',num2str(coorStart))});
			ind2 = values(allSNPDic, {strcat(num2str(chrTest),'.',num2str(coorEnd))});
			ind1 = ind1{1};
			ind2 = ind2{1};
			numSNP = ind2 - ind1 + 1;
			numNIM = sum(isKey(nimDic,snpID(ind1:ind2)));
 			numMH = numSNP - numNIM;
			
			output = [chrTest, coorStart, coorEnd, numCSet, numSNP, numNIM, numMH,numNIMSusie, numMHSusie, numSNPSusie];
			fprintf(fidOut, '%d\t', output);
			fprintf(fidOut,'\n');	
		end
		fclose(fidOut);
	end
end
