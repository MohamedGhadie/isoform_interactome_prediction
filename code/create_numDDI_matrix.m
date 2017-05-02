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
