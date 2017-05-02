function [numLosingIsoforms, fracLosingIsoforms, numIsoPairs, fracDiffIsoProfilePerGene, avgIsoProfileDistPerGene, allDist] = compare_isoform_interactions (I, domI, domPrI, numDDImap, altIsoforms, numAltIsoforms, isoInterDomains)

numGenes = size(I,1);
inIsoInteractome = zeros(numGenes,1);
numRefProteins = 0;
for i = 1:numGenes
    partners = find(I(i,:));
    partners(partners==i) = [];
    numPartners = length(partners);
    p1domains = find(domPrI(:,i));
    if numPartners > 0
        % pick partners that have a DDI mapping with the canonical of
        % gene i
        for j = numPartners:-1:1
            if numDDImap(i,partners(j)) == 0
                partners(j) = [];
            end
        end
        numPartners = length(partners);
        if numPartners > 0
            inIsoInteractome(i) = 1;
            numRefProteins = numRefProteins + 1;
        end
    end
end

lost = 0;
notlost = 0;
numRefProteins = 0;
numIsoProteins = 0;
numLosingIsoforms = nan(numGenes,1);
fracLosingIsoforms = nan(numGenes,1);
fracDiffIsoProfilePerGene = nan(numGenes,1);
avgIsoProfileDistPerGene = nan(numGenes,1);
allDist = [];
for i = 1:numGenes
    % pick genes with at least one isoform (plus the canonical)
    %if ~isempty(spID{i}) && (maxIsoform(i) > 0)
    if numAltIsoforms(i)>0
        numIso = length(altIsoforms{i});
        partners = find(I(i,:));
        partners(partners==i) = [];
        numPartners = length(partners);
        p1domains = find(domPrI(:,i));
        if numPartners > 0
            % pick partners that have a DDI mapping with the canonical of
            % gene i
            for j = numPartners:-1:1
                if numDDImap(i,partners(j)) == 0
                    partners(j) = [];
                end
            end
            numPartners = length(partners);
            if numPartners > 0
                numRefProteins = numRefProteins + 1;
                numIsoProteins = numIsoProteins + numIso;
                isoformLoss = zeros(numIso,1);
                profiles = [ones(1,numPartners); zeros(numIso,numPartners)];
                for j = 1:numPartners
                    p2domains = find(domPrI(:,partners(j)));
                    p2domPartners = p1domains(sum(domI(p2domains,p1domains),1)>0);
                    numP2domPartners = length(p2domPartners);
                    if numP2domPartners == 0
                        continue
                    end
                    for m = 1:numIso
                        for n = 1:numP2domPartners
                            if ~isempty(find(isoInterDomains{i}{altIsoforms{i}(m)+1}==p2domPartners(n),1))
                                notlost = notlost + 1;
                                profiles(m+1,j) = 1;
                                break
                            elseif n == numP2domPartners && isempty(find(isoInterDomains{i}{altIsoforms{i}(m)+1}==p2domPartners(n),1))
                                lost = lost + 1;
                                isoformLoss(m) = isoformLoss(m) + 1;
                            end
                        end
                    end
                end
                numLosingIsoforms(i) = sum(isoformLoss>0);
                fracLosingIsoforms(i) = numLosingIsoforms(i)/(numIso+1);
                numIs = size(profiles,1);
                dist = zeros(((numIs^2)-numIs)/2,1);
                c = 0;
                for k = 1:numIs-1
                    for q = k+1:numIs
                        c = c + 1;
                        dist(c) = sum(profiles(k,:)~=profiles(q,:))/numPartners;
                    end
                end
                allDist = [allDist; dist];
                fracDiffIsoProfilePerGene(i) = sum(dist>0)/length(dist);
                avgIsoProfileDistPerGene(i) = mean(dist);
            end
        end
    end
end

numIsoPairs = (((numAltIsoforms+1).^2)-(numAltIsoforms+1))/2;
