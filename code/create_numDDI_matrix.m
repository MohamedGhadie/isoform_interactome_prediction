function numDDImap = create_numDDI_matrix (I, domPrI, domI)

numGenes = size(I,1);
numDDImap = uint8(zeros(numGenes,numGenes));
for i = 1:numGenes
    for j = i:numGenes
        if I(i,j)
            domList1 = find(domPrI(:,i));
            domList2 = find(domPrI(:,j));
            if ~isempty(domList1) && ~isempty(domList2)
                numDDImap(i,j) = sum(sum(domI(domList1,domList2)));
                numDDImap(j,i) = numDDImap(i,j);
            end
        end
    end
end

fprintf('\n');
disp([num2str(sum(sum(numDDImap)>0)) ' proteins in the domain-resolved interactome']);
disp([num2str(sum(sum(triu(numDDImap)>0))) ' PPIs with at least 1 mapping DDI']);
disp([num2str(sum(sum(triu(I))) - sum(sum(triu(numDDImap)>0))) ' PPIs with no mapping DDIs']);
disp([num2str(sum(sum(triu(numDDImap)==1))) ' PPIs with 1 mapping DDI']);
disp([num2str(sum(sum(triu(numDDImap)==2))) ' PPIs with 2 mapping DDIs']);
disp([num2str(sum(sum(triu(numDDImap)==3))) ' PPIs with 3 mapping DDIs']);
disp([num2str(sum(sum(triu(numDDImap)==4))) ' PPIs with 4 mapping DDIs']);
disp([num2str(sum(sum(triu(numDDImap)>4))) ' PPIs with >=5 mapping DDIs']);

figure
bar([sum(sum(triu(I)))-sum(sum(triu(numDDImap)>0)) sum(sum(triu(numDDImap)==1)) ...
    sum(sum(triu(numDDImap)==2)) sum(sum(triu(numDDImap)>2))]);
set(gca,'XTick',1:4,'XTickLabel',{'0','1','2','\geq 3'});
xlabel('Number of DDI annotations');
ylabel('Number of PPIs');
set(gca,'tickDir','out');
box off
