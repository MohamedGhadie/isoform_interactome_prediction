function [inIsoInteractome, ref_ref_interactions, alt_ref_interactions, numLosingIsoforms, fracLosingIsoforms] = predict_isoform_interactome (spID,I,domI,domPrI,numDDImap,altIsoforms,isoInterDomains,filename)

numGenes = size(I,1);
ref_ref_interactions = {};
fid = fopen(filename,'w');
c = 0;
inIsoInteractome = zeros(numGenes,1);
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
            for j = 1:numPartners
                if i<=partners(j)
                    c = c+1;
                    ref_ref_interactions{c,1} = spID{i};
                    ref_ref_interactions{c,2} = spID{partners(j)};
                    ref_ref_interactions{c,3} = 'known';
                    fprintf(fid,'%s\t',spID{i});
                    fprintf(fid,'%s\t',spID{partners(j)});
                    fprintf(fid,'known\n');
                end
            end
        end
    end
end
disp([num2str(c) ' reference-reference interactions written to file ' filename]);

lost = 0;
notlost = 0;
alt_ref_interactions = {};
c = 0;
numLosingIsoforms = nan(numGenes,1);
fracLosingIsoforms = nan(numGenes,1);
for i = 1:numGenes
    numAltIso = length(altIsoforms{i});
    if numAltIso>0
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
                isoformLoss = zeros(numAltIso,1);
                for j = 1:numPartners
                    p2domains = find(domPrI(:,partners(j)));
                    p2domPartners = p1domains(sum(domI(p2domains,p1domains),1)>0);
                    numP2domPartners = length(p2domPartners);
                    if numP2domPartners == 0
                        continue
                    end
                    for m = 1:numAltIso
                        c = c + 1;
                        alt_ref_interactions{c,1} = [spID{i} '-' num2str(altIsoforms{i}(m))];
                        alt_ref_interactions{c,2} = spID{partners(j)};
                        fprintf(fid,'%s\t',[spID{i} '-' num2str(altIsoforms{i}(m))]);
                        fprintf(fid,'%s\t',spID{partners(j)});
                        for n = 1:numP2domPartners
                            if ~isempty(find(isoInterDomains{i}{altIsoforms{i}(m)+1}==p2domPartners(n),1))
                                notlost = notlost + 1;
                                alt_ref_interactions{c,3} = 'retained';
                                fprintf(fid,'retained\n');
                                break
                            elseif n == numP2domPartners && isempty(find(isoInterDomains{i}{altIsoforms{i}(m)+1}==p2domPartners(n),1))
                                lost = lost + 1;
                                alt_ref_interactions{c,3} = 'lost';
                                fprintf(fid,'lost\n');
                                isoformLoss(m) = isoformLoss(m) + 1;
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
fclose(fid);
disp([num2str(notlost) ' predicted retained isoform interactions written to file ' filename]);
disp([num2str(lost) ' predicted lost isoform interactions written to file ' filename]);
