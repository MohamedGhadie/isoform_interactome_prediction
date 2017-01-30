function isoInterDomains = get_isoform_interacting_domains (spID, maxIsoform, domains, refDomMap, domAltMap)

% compile interacting domains in each isoform of a gene. Each row in the
% cell array corresponds to one gene. First entry in the row contains
% domains of the reference protein. Entry k+1 contains domains of the kth alternative isoform
numGenes = length(spID);
isoInterDomains = cell(numGenes,1);
for i = 1:numGenes
    spid = spID{i};
    if ~isempty(spid)
        isoInterDomains{i} = cell(1,maxIsoform(i)+1);
        dom = refDomMap(strcmpi(refDomMap(:,1),spid),2);
        dom = unique(dom);
        domNum = [];
        if ~isempty(dom)
            for j = 1:length(dom)
                num = find(strcmpi(domains,dom{j}),1);
                if ~isempty(num)
                    domNum = [domNum; uint16(num)];
                end
            end
        end
        isoInterDomains{i}{1} = domNum;
        
        % do the same for the alternative isoforms
        if maxIsoform(i) > 0
            for k = 1:maxIsoform(i)
                dom = domAltMap(strcmpi(domAltMap(:,2),[spid '-' num2str(k)]),1);
                dom = unique(dom);
                domNum = [];
                if ~isempty(dom)
                    for j = 1:length(dom)
                        num = find(strcmpi(domains,dom{j}),1);
                        if ~isempty(num)
                            domNum = [domNum; uint16(num)];
                        end
                    end
                end
                isoInterDomains{i}{k+1} = domNum;
            end
        end
    else
        isoInterDomains{i} = {};
    end
end
