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
