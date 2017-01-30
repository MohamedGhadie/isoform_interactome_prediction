function numCommonPartners = get_num_common_partners (I)

numGenes = size(I,1);
numCommonPartners = zeros(size(I));
for i = 1:numGenes
    for j = i:numGenes
        numCommonPartners(i,j) = sum(I(i,:)&I(j,:));
    end
end
numCommonPartners = numCommonPartners + triu(numCommonPartners,1)';