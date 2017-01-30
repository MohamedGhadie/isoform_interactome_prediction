function [pairs,selected,nodiffTargets,subsetTargets,changeoverTargets] = label_iso_partner_pairs(spID,I,domI,domPrI,isoInterDomains,numDDImap,maxIsoform)

numPairs = 500000;      % estimate for memory allocation
pairs= zeros(numPairs,2);
selected = zeros(numPairs,4);
nodiffTargets = cell(numPairs,1);
subsetTargets = cell(numPairs,1);
changeoverTargets = cell(numPairs,1);
for i = 1:size(I,1)
    spid = spID{i};
    % pick genes with at least two isoforms
    if ~isempty(spid) && (maxIsoform(i) > 0)
        numIso = maxIsoform(i) + 1;
        partners = find(I(i,:));
        partners(partners==i) = [];
        numPartners = length(partners);
        p1domains = find(domPrI(:,i));
        if numPartners > 1
            % pick partners that have a DDI mapping with the canonical of
            % gene i
            for j = numPartners:-1:1
                if numDDImap(i,partners(j)) ==0
                    partners(j) = [];
                end
            end
            numPartners = length(partners);
            if numPartners > 1
                for j = 1:numPartners-1
                    if numPartners > 10
                            [i numPartners j]
                    end
                    p2domains = find(domPrI(:,partners(j)));
                    p2domPartners = p1domains(sum(domI(p2domains,p1domains),1)>0);
                    numP2domPartners = length(p2domPartners);
                    if numP2domPartners == 0
                        continue
                    end
                    isoInteractions = zeros(numIso,2);
                    for m = 1:numIso
                        if ~isempty(isoInterDomains{i}{m})
                            for n = 1:numP2domPartners
                                if ~isempty(find(isoInterDomains{i}{m}==p2domPartners(n),1))
                                    isoInteractions(m,1) = 1;
                                    break
                                end
                            end
                        end
                    end
                    for k = j+1:numPartners
                        p2domains = find(domPrI(:,partners(k)));
                        p2domPartners = p1domains(sum(domI(p2domains,p1domains),1)>0);
                        numP2domPartners = length(p2domPartners);
                        if numP2domPartners == 0
                            continue
                        end
                        isoInteractions(:,2) = zeros(numIso,1);
                        for m = 1:numIso
                            if ~isempty(isoInterDomains{i}{m})
                                for n = 1:numP2domPartners
                                    if ~isempty(find(isoInterDomains{i}{m}==p2domPartners(n),1))
                                        isoInteractions(m,2) = 1;
                                        break
                                    end
                                end
                            end
                        end
                        % record the number of isoforms both partners interact
                        % with for this specific gene
                        numP1interactions = sum(isoInteractions(:,1));
                        numP2interactions = sum(isoInteractions(:,2));
                        numCommonInteractions = sum(isoInteractions(:,1) & isoInteractions(:,2));
                        numDiff = sum(xor(isoInteractions(:,1),isoInteractions(:,2)));
                        nodiff = numCommonInteractions * (numDiff == 0);
                        subset = (numCommonInteractions > 0)* ...
                            ((numCommonInteractions == numP1interactions & numCommonInteractions < numP2interactions) ...
                            + 2*(numCommonInteractions < numP1interactions & numCommonInteractions == numP2interactions));
                        changeover = (numP1interactions > numCommonInteractions) & (numP2interactions > numCommonInteractions);
                        
                        pairInd = find(pairs(:,1)==partners(j) & pairs(:,2)==partners(k),1);
                        if isempty(pairInd)
                            pairInd = find(pairs(:,1)==0,1);
                            if isempty(pairInd)
                                pairInd = size(pairs,1) + 1;
                            end
                        end
                        pairs(pairInd,:) = [partners(j) partners(k)];
                        selected(pairInd,1) = selected(pairInd,1) + (nodiff>0);
                        selected(pairInd,2) = selected(pairInd,2) + (subset>0);
                        selected(pairInd,3) = selected(pairInd,3) + (changeover>0);
                        if nodiff > 0
                            nodiffTargets{pairInd} = [nodiffTargets{pairInd} i];
                        end
                        if subset > 0
                            subsetTargets{pairInd} = [subsetTargets{pairInd} i];
                        end
                        if changeover > 0
                            changeoverTargets{pairInd} = [changeoverTargets{pairInd} i];
                        end
                    end
                end
            end
        end
    end
end

count = find(pairs(:,1)==0,1);
if ~isempty(count)
    numPairs = count - 1;
    pairs = pairs(1:numPairs,:);
    selected = selected(1:numPairs,:);
    nodiffTargets = nodiffTargets(1:numPairs);
    subsetTargets = subsetTargets(1:numPairs);
    changeoverTargets = changeoverTargets(1:numPairs);
end
