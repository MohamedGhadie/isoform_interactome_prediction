function domI = create_DDI_matrix (domains, interactions)

numDom = length(domains);
domI = uint8(zeros(numDom,numDom));
for i = 1:size(interactions,1)
    ind1 = find(strcmpi(domains,interactions{i,1}),1);
    ind2 = find(strcmpi(domains,interactions{i,2}),1);
    domI(ind1,ind2) = 1;
    domI(ind2,ind1) = 1;
end