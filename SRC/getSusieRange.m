clear all
close all
clc


snpID = importdata('snpID.txt');
chrCoor = [0];

for i = 1:22
	ind = find(snpID(:, 1)==i);
	chrCoor = [chrCoor; chrCoor(i) + length(ind)];
end

bin = 10^5;

pcutoff='10';

Name = strsplit(ls(strcat('prunedNIM/*.',pcutoff,'.*.in')));
ind = find(contains(Name,'prune.in')==1);
Name = Name(ind);

for j = 1:length(Name)
        filename = Name{j};
	nimassoc = importdata(filename);
	row = size(nimassoc,1);
	simName = filename(11:end-9);
	outfile= strcat('susieRange/',simName,'.nimrange');
	fid = fopen(outfile,'w');
	for j = 1:row
		chr = nimassoc(j, 1);
		coor = nimassoc(j, 2);
		Coor = snpID(chrCoor(chr)+1:chrCoor(chr+1), 2);
		indLow = find(Coor <= coor - bin);
		if length(indLow) == 0
			coorLow = Coor(1);
		else
			coorLow = Coor(indLow(end));
		end
		indHigh = find(Coor >= coor + bin);
		if length(indHigh) == 0
			coorHigh = Coor(end);
		else
			coorHigh = Coor(indHigh(1));
		end
		fprintf(fid,'%s\n',strcat(num2str(chr),':',num2str(coorLow),'-',num2str(chr),':',num2str(coorHigh)));		
	end
	fclose(fid);
end
