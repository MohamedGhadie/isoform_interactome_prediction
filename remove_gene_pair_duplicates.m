function [numRepeated, pairs] = remove_gene_pair_duplicates (pairs,genes)
rep = ones(size(pairs,1),1);
rmv = zeros(size(pairs,1),1);
gn = [genes(pairs(:,1)) genes(pairs(:,2))];
numRepeated = 0;
for i = 1:size(pairs,1)
    gn1 = genes(pairs(i,1));
    gn2 = genes(pairs(i,2));
    if ~strcmpi(gn1,'') && ~strcmpi(gn2,'')
        rep_ind = find((strcmpi(gn(:,1),gn1) & strcmpi(gn(:,2),gn2)) | (strcmpi(gn(:,1),gn2) & strcmpi(gn(:,2),gn1)));
        rep(i) = length(rep_ind);
        if rep(i) > 1
            rmv(rep_ind(2:rep(i))) = 1;
            if rmv(i) == 0
                numRepeated = numRepeated + 1;
            end
        end
    end
end
pairs = pairs(~rmv,:);
