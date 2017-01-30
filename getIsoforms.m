function [altIsoforms, numAltIsoforms, maxIsoform] = getIsoforms (spID, isoNames)

numGenes = length(spID);
altIsoforms = cell(numGenes,1);
numAltIsoforms = zeros(numGenes,1);
maxIsoform = zeros(numGenes,1);
SPmatches = 0;
for i = 1:numGenes
    if ~isempty(spID{i})
        SPmatches = SPmatches + 1;
        ind = find(strncmpi(isoNames,[spID{i} '-'],length(spID{i})+1));
        if ~isempty(ind)
            iso = unique(isoNames(ind));
            numAltIsoforms(i) = length(ind);
            isoNum = zeros(1,length(ind));
            for j = 1:numAltIsoforms(i)
                dashPos = find(iso{j}=='-');
                isoNum(j) = str2double(iso{j}(dashPos+1:length(iso{j})));
            end
            isoNum = sort(isoNum);
            maxIsoform(i) = isoNum(numAltIsoforms(i));
            altIsoforms{i} = isoNum;
        end
    end
end
disp([num2str(sum(numAltIsoforms>0)) ' out of the ' num2str(numGenes) ' reference proteins have atleast one alternative isoform']);

figure
bar([sum(numAltIsoforms==0) sum(numAltIsoforms==1) sum(numAltIsoforms==2) sum(numAltIsoforms==3) sum(numAltIsoforms==4) ...
    sum(numAltIsoforms==5) sum(numAltIsoforms==6) sum(numAltIsoforms==7) sum(numAltIsoforms==8) sum(numAltIsoforms==9) ...
    sum(numAltIsoforms>9)]);
set(gca,'XTick',1:11,'XTickLabel',{'0','1','2','3','4','5','6','7','8','9','\geq 10'});
xlabel('Number of alternative isoforms');
ylabel('Number of reference proteins');
set(gca,'tickDir','out');
box off
