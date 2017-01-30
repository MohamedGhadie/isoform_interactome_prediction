function domPrI = create_domain_protein_matrix (domains, spPFmap, spID)

numGenes = length(spID);
numDom = length(domains);
domPrI = uint8(zeros(numDom,numGenes));
for i = 1:numDom
    sp = spPFmap(strcmpi(spPFmap(:,2),domains{i}),1);
    if ~isempty(sp)
        for j = 1:length(sp)
            geneNum = find(strcmpi(spID,sp{j}),1);
            if ~isempty(geneNum)
                domPrI(i,geneNum) = 1;
            end
        end
    end
end