function [inIsoInteractome, ref_ref_interactions, alt_ref_interactions, numLosingIsoforms, fracLosingIsoforms] = predict_isoform_interactome (spID,domains,PPIs,I,domI,domPrI,numDDImap,altIsoforms,isoInterDomains,interactome,filename)

numGenes = size(I,1);
ref_ref_interactions = cell(sum(sum(triu(I))),4);
fid = fopen(filename,'w');
c = 0;
inIsoInteractome = zeros(numGenes,1);
if strcmpi(interactome,'HI-II-14')
    fprintf(fid,[strjoin({'Isoform_UniProt_ID', 'Partner_UniProt_ID', 'Interaction_Prediction', 'Mapping_DDIs'},'\t') '\n']);
elseif strcmpi(interactome,'IntAct')
    fprintf(fid,[strjoin({'Isoform_UniProt_ID', 'Partner_UniProt_ID', 'Interaction_Prediction', 'Reference_1_IntAct_ID', 'Reference_2_IntAct_ID', 'Mapping_DDIs'},'\t') '\n']);
end
for i = 1:numGenes
    partners = find(I(i,:));
    partners(partners==i) = [];
    numPartners = length(partners);
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
            p1domains = find(domPrI(:,i));
            for j = 1:numPartners
                if i<=partners(j)
                    c = c+1;
                    p2domains = find(domPrI(:,partners(j)));
                    mappingDDIs = cell(1,sum(sum(domI(p1domains,p2domains))));
                    d = 0;
                    for m = 1:length(p1domains)
                        for n = 1:length(p2domains)
                            if domI(p1domains(m),p2domains(n)) > 0
                                d = d + 1;
                                mappingDDIs{d} = [domains{p1domains(m)} '-' domains{p2domains(n)}];
                            end
                        end
                    end
                    ref_ref_interactions(c,:) = {spID{i}, spID{partners(j)}, 'known', mappingDDIs};
                    if strcmpi(interactome,'HI-II-14')
                        fprintf(fid,[strjoin({spID{i}, spID{partners(j)}, 'known', strjoin(mappingDDIs,', ')},'\t') '\n']);
                    elseif strcmpi(interactome,'IntAct')
                        ind = find(strcmpi(PPIs(:,1),spID{i}) & strcmpi(PPIs(:,2),spID{partners(j)}));
                        if length(ind)==1
                            id1 = PPIs{ind,3};
                            id2 = PPIs{ind,4};
                        else
                            ind = find(strcmpi(PPIs(:,1),spID{partners(j)}) & strcmpi(PPIs(:,2),spID{i}));
                            id1 = PPIs{ind,4};
                            id2 = PPIs{ind,3};
                        end
                        fprintf(fid,[strjoin({spID{i}, spID{partners(j)}, 'known', id1, id2, strjoin(mappingDDIs,', ')},'\t') '\n']);
                    end
                end
            end
        end
    end
end
ref_ref_interactions = ref_ref_interactions(1:c,:);
disp([num2str(c) ' reference-reference interactions written to file ' filename]);

lost = 0;
notlost = 0;
numAltIsoforms = arrayfun(@(x) length(x), altIsoforms);
numInteractions = sum(I,2);
alt_ref_interactions = cell(sum(numAltIsoforms.*numInteractions),4);
c = 0;
numLosingIsoforms = nan(numGenes,1);
fracLosingIsoforms = nan(numGenes,1);
for i = 1:numGenes
    numAltIso = length(altIsoforms{i});
    if numAltIso>0
        partners = find(I(i,:));
        partners(partners==i) = [];
        numPartners = length(partners);
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
                p1domains = find(domPrI(:,i));
                isoformLoss = zeros(numAltIso,1);
                for j = 1:numPartners
                    p2domains = find(domPrI(:,partners(j)));
                    p2domPartners = p1domains(sum(domI(p2domains,p1domains),1)>0);
                    numP2domPartners = length(p2domPartners);
                    if numP2domPartners == 0
                        continue
                    end
                    refDDIs = cell(sum(sum(domI(p1domains,p2domains))),2);
                    d = 0;
                    for m = 1:length(p1domains)
                        for n = 1:length(p2domains)
                            if domI(p1domains(m),p2domains(n)) > 0
                                d = d + 1;
                                refDDIs(d,:) = {domains{p1domains(m)}, domains{p2domains(n)}};
                            end
                        end
                    end
                    refDDIs = refDDIs(1:d,:);
                    for m = 1:numAltIso
                        c = c + 1;
                        alt_ref_interactions(c,[1 2]) = {[spID{i} '-' num2str(altIsoforms{i}(m))], spID{partners(j)}};
                        fprintf(fid,[spID{i} '-%d\t' spID{partners(j)} '\t'], altIsoforms{i}(m));
                        if strcmpi(interactome,'IntAct')
                            ind = find(strcmpi(PPIs(:,1),spID{i}) & strcmpi(PPIs(:,2),spID{partners(j)}));
                            if length(ind)==1
                                id1 = PPIs{ind,3};
                                id2 = PPIs{ind,4};
                            else
                                ind = find(strcmpi(PPIs(:,1),spID{partners(j)}) & strcmpi(PPIs(:,2),spID{i}));
                                id1 = PPIs{ind,4};
                                id2 = PPIs{ind,3};
                            end
                        end
                        isoInteractingDom = domains(isoInterDomains{i}{altIsoforms{i}(m)+1});
                        isoDDIs = refDDIs(ismember(refDDIs(:,1),isoInteractingDom),:);
                        if isempty(isoDDIs)
                            lost = lost + 1;
                            isoformLoss(m) = isoformLoss(m) + 1;
                            alt_ref_interactions(c,[3 4]) = {'lost','-'};
                            if strcmpi(interactome,'HI-II-14')
                                fprintf(fid,'lost\t-\n');
                            elseif strcmpi(interactome,'IntAct')
                                fprintf(fid,[strjoin({'lost', id1, id2, '-'},'\t') '\n']);
                            end
                        else
                            notlost = notlost + 1;
                            textDDIs = strjoin(cellfun(@(x,y) [x '-' y], isoDDIs(:,1), isoDDIs(:,2),'UniformOutput',false), ', ');
                            alt_ref_interactions(c,[3 4]) = {'retained',textDDIs};
                            if strcmpi(interactome,'HI-II-14')
                                fprintf(fid,['retained\t' textDDIs '\n']);
                            elseif strcmpi(interactome,'IntAct')
                                fprintf(fid,[strjoin({'retained', id1, id2, textDDIs},'\t') '\n']);
                            end
                        end
                    end
                end
                numLosingIsoforms(i) = sum(isoformLoss>0);
                fracLosingIsoforms(i) = numLosingIsoforms(i)/(numAltIso+1);
            end
        end
    end
end
alt_ref_interactions = alt_ref_interactions(1:c,:);
fclose(fid);
disp([num2str(notlost) ' predicted retained isoform interactions written to file ' filename]);
disp([num2str(lost) ' predicted lost isoform interactions written to file ' filename]);
