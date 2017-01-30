function [GOterms, GOaspects, gene_associations] = load_GO_associations(filename1,filename2)

[~,gene_association_a] = xlsread(filename1);
[~,gene_association_b] = xlsread(filename2);

gene_associations = [gene_association_a; gene_association_b];

GOterms = unique(gene_associations(:,5));
numGOterms = length(GOterms);

GOaspects = char(numGOterms,1);
for i = 1:numGOterms
    aspects = gene_associations(strcmpi(gene_associations(:,5),GOterms{i}),9);
    if sum(strcmpi(aspects,aspects{1})) == length(aspects)
        GOaspects(i) = aspects{1};
    else
        GOaspects(i) = 'N';
    end
end
