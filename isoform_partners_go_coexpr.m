
interactome = 'IntAct';
disp('Loading reference interactome');
if strcmpi(interactome,'HI-II-14')
    interactomeFile = 'HI-II-14.xlsx';
    isoformInteractomeFile = 'HI-II-14_isoform_interactome';
    spEntrezMapFile = 'rolland_spEntrezMap.tab';
    numTimesReported = 1;
    removeSelfInteractions = 1;
    [I, spID, genes] = load_interactome(interactomeFile, [], spEntrezMapFile, numTimesReported, removeSelfInteractions, interactome);
elseif strcmpi(interactome,'IntAct')
    interactomeFile = 'intact.txt';
    isoformInteractomeFile = 'IntAct_isoform_interactome';
    spGeneMapFile = 'intact_spGeneMap.tab.txt';
    numTimesReported = 2;
    removeSelfInteractions = 1;
    [I, spID, genes] = load_interactome(interactomeFile, spGeneMapFile, [], numTimesReported, removeSelfInteractions, interactome);
end
numPPI = sum(sum(triu(I)));
numGenes = size(I,1);

disp('Loading domain-domain interactions');
DDIs = loadDDIs('interaction.xlsx');
domains = unique([DDIs(:,1); DDIs(:,2)]);
numDom = length(domains);
disp('Loading structural domains for reference proteins');
[refDomMap, domRefPos] = load_hmmscan_ref('hmmer_canonical_full.domtab2.txt');
disp('Loading structural domains for alternative isoforms');
% Remove all leading and trailing lines other than raw data from hmmer scan file before loading
% Loading will start from first line
[domAltMap, domAltPos] = load_hmmscan_alt('hmmer_isoforms.domtab2.txt');
disp('Loading reference protein sequences');
[headers, prSeq] = load_ref_protein_seq('human_protein_sequences.tab');
disp('Loading alternative isoform sequences');
isoSeq = load_alt_protein_seq('human_protein_can_iso_sequences.fasta');
isoNames = unique(isoSeq(:,1));
disp('Counting number of alternative isoforms per gene');
[altIsoforms, numAltIsoforms, maxIsoform] = getIsoforms(spID, isoNames);
disp('Identifying isoform interacting domains');
isoInterDomains = get_isoform_interacting_domains(spID, maxIsoform, domains, refDomMap, domAltMap);
disp('Creating domain interaction matrix');
domI = create_DDI_matrix(domains, DDIs);
disp('Creating domain-protein mapping matrix');
domPrI = create_domain_protein_matrix(domains, refDomMap, spID);
disp('Counting number of domain-domain interactions mapping to each protein-protein interaction');
numDDImap = create_numDDI_matrix(I, domPrI, domI);
disp('Calculating number of common partners for each pair of reference proteins');
numCommonPartners = get_num_common_partners(I);
disp('Predicting isoform interactome');
[inIsoInteractome, ref_ref_interactions, alt_ref_interactions, ~, ~] = predict_isoform_interactome(spID,I,domI,domPrI,numDDImap,altIsoforms,isoInterDomains,isoformInteractomeFile);
disp('Comparing isoform interaction profiles');
[numLosingIsoforms, fracLosingIsoforms, numIsoPairs, fracDiffIsoProfilePerGene, avgIsoProfileDistPerGene] = compare_isoform_interactions(I,domI,domPrI,numDDImap,altIsoforms,numAltIsoforms,isoInterDomains);

fprintf('\nPredicted isoform interactome statistics:\n');
disp([num2str(length(unique(ref_ref_interactions(:,1:2)))) ' reference proteins and ' num2str(length(unique(alt_ref_interactions(:,1)))) ' alternative isoforms']);
disp([num2str(size(ref_ref_interactions,1)) ' known reference interactions']);
disp([num2str(sum(strcmpi(alt_ref_interactions(:,3),'retained'))) ' predicted retained isoform interactions']);
disp([num2str(sum(strcmpi(alt_ref_interactions(:,3),'lost'))) ' predicted lost isoform interactions']);
disp([num2str(sum(numAltIsoforms(sum(numDDImap)>0)>0)) ' genes have at least two isoforms (including reference)']);
disp([num2str(sum(numLosingIsoforms>0)) ' genes have at least one isoform losing an interaction']);
fprintf('\n');

% load gene ontology associations
disp('Loading gene ontology associations');
[~,gene_association_a] = xlsread('gene_association.goa_ref_human_a.xlsx');
[~,gene_association_b] = xlsread('gene_association.goa_ref_human_b.xlsx');
gene_associations = [gene_association_a; gene_association_b];

disp('Identifying GO aspect for each GO term');
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

% create protein GO association matrix
disp('Creating protein-GO association matrix');
prGO = uint8(zeros(numGenes,numGOterms));
for i = 1:numGOterms
    ind = find(strcmpi(gene_associations(:,5),GOterms{i}));
    if ~isempty(ind)
        for j = 1:length(ind)
            sp = gene_associations{ind(j),2};
            prGO(strcmpi(spID,sp),i) = 1;
        end
    end
end
clear i j go ind sp aspects

numPrGOterms = sum(prGO,2);
numPrGOFterms = sum(prGO(:,GOaspects=='F'),2);
numPrGOPterms = sum(prGO(:,GOaspects=='P'),2);
numPrGOCterms = sum(prGO(:,GOaspects=='C'),2);

% load gene tissue expression data
disp('Loading gene tissue-expression data');
tissueExpressionFile = 'E-MTAB-513.tsv.txt';
tdfread(tissueExpressionFile);
expr = [adipose adrenal brain breast colon heart kidney leukocyte liver lung lymph_node ovary prostate skeletal_muscle testis thyroid];
expr = log2(expr);
clear adipose adrenal brain breast colon heart kidney leukocyte liver lung lymph_node ovary prostate skeletal_muscle testis thyroid

expr_gene_id = cell(size(Gene_ID,1),1);
expr_gene_names = cell(size(Gene_Name,1),1);
for i = 1:size(Gene_ID,1)
    expr_gene_id{i} = strtrim(Gene_ID(i,:));
    expr_gene_names{i} = strtrim(Gene_Name(i,:));
end
clear i Gene_ID Gene_Name

% identify rows with no expression (-inf) in all tissues
tobeRemoved = [];
for i = 1:size(expr,1)
    if isempty(find(expr(i,:)>-inf,1))
        tobeRemoved = [tobeRemoved i];
    end
end
fprintf('\nThe following rows of the expression matrix have no expression in all tissues:\n');
disp(num2str(tobeRemoved));

% remove rows with no expression in all tissues
expr(tobeRemoved,:) = [];
expr_gene_id(tobeRemoved) = [];
expr_gene_names(tobeRemoved) = [];

% create expression matrix for interacting proteins
disp('Creating tissue expression matrix for all genes in the reference interactome');
numTissues = size(expr,2);
Iexpr = NaN(numGenes,numTissues);
c = [];
for i = 1:numGenes
    if ~isempty(genes{i})
        ind = find(strcmpi(expr_gene_names,genes{i}));
        if length(ind) == 1
            Iexpr(i,:) = expr(ind,:);
        elseif length(ind) > 1
            c = [c i];
        end
    end
end
fprintf('\nThe following genes mapping to more than one expression row were excluded:\n');
disp(num2str(c));
clear i ind tobeRemoved expr_gene_id expr_gene_names c

[pairs,selected,nodiffTargets,subsetTargets,changeoverTargets] = label_iso_partner_pairs(spID,I,domI,domPrI,isoInterDomains,numDDImap,maxIsoform);

category = {'no-difference','subset','changeover'};

nodiff_pairs = pairs(selected(:,1) & ~selected(:,2) & ~selected(:,3),:);
diff_pairs = pairs(~selected(:,1) & (selected(:,2)|selected(:,3)),:);

disp([num2str(size(nodiff_pairs,1)) ' protein pairs interacting with the same subset of isoforms of the same gene']);
disp([num2str(size(diff_pairs,1)) ' protein pairs interacting with different subsets of isoforms of the same gene']);

[numNodiffPairsRepeated, nodiff_pairs] = remove_gene_pair_duplicates(nodiff_pairs,genes);
disp([num2str(numNodiffPairsRepeated) ' protein pairs with identical isoform interaction profiles occur more than once with different uniProt IDs: duplicates removed']);

[numDiffPairsRepeated, diff_pairs] = remove_gene_pair_duplicates(diff_pairs,genes);
disp([num2str(numDiffPairsRepeated) ' protein pairs with different isoform interaction profiles occur more than once with different uniProt IDs: duplicates removed']);

% coexpression and GO similarity of pairs of proteins interacting with the
% same subset of isoforms of the same gene
disp('Calculating GO similarity and coexpression of pairs of proteins with identical isoform interaction profiles');
nodiff_gosim = [];
nodiff_gofsim = [];
nodiff_gopsim = [];
nodiff_gocsim = [];
nodiff_coexpr = [];
for i = 1:size(nodiff_pairs,1)
    if numPrGOterms(nodiff_pairs(i,1)) || numPrGOterms(nodiff_pairs(i,2))
        nodiff_gosim = [nodiff_gosim 1-sum(prGO(nodiff_pairs(i,1),:)~=prGO(nodiff_pairs(i,2),:))/sum(prGO(nodiff_pairs(i,1),:)|prGO(nodiff_pairs(i,2),:))];
    end
    if numPrGOFterms(nodiff_pairs(i,1)) || numPrGOFterms(nodiff_pairs(i,2))
        nodiff_gofsim = [nodiff_gofsim 1-sum(prGO(nodiff_pairs(i,1),GOaspects=='F')~=prGO(nodiff_pairs(i,2),GOaspects=='F'))/sum(prGO(nodiff_pairs(i,1),GOaspects=='F')|prGO(nodiff_pairs(i,2),GOaspects=='F'))];
    end
    if numPrGOPterms(nodiff_pairs(i,1)) || numPrGOPterms(nodiff_pairs(i,2))
        nodiff_gopsim = [nodiff_gopsim 1-sum(prGO(nodiff_pairs(i,1),GOaspects=='P')~=prGO(nodiff_pairs(i,2),GOaspects=='P'))/sum(prGO(nodiff_pairs(i,1),GOaspects=='P')|prGO(nodiff_pairs(i,2),GOaspects=='P'))];
    end
    if numPrGOCterms(nodiff_pairs(i,1)) || numPrGOCterms(nodiff_pairs(i,2))
        nodiff_gocsim = [nodiff_gocsim 1-sum(prGO(nodiff_pairs(i,1),GOaspects=='C')~=prGO(nodiff_pairs(i,2),GOaspects=='C'))/sum(prGO(nodiff_pairs(i,1),GOaspects=='C')|prGO(nodiff_pairs(i,2),GOaspects=='C'))];
    end
    notnan = ~isnan(Iexpr(nodiff_pairs(i,1),:)) & ~isnan(Iexpr(nodiff_pairs(i,2),:));
    num = sum(notnan);
    if num > 7
        [r,~,~,~] = corrcoef([Iexpr(nodiff_pairs(i,1),notnan); Iexpr(nodiff_pairs(i,2),notnan)]','rows','pairwise');
        if ~isnan(r(1,2))
            nodiff_coexpr = [nodiff_coexpr r(1,2)];
        end
    end
end

% coexpression and GO similarity of pairs of proteins interacting with
% different subsets of isoforms of the same gene
disp('Calculating GO similarity and coexpression of pairs of proteins with different isoform interaction profiles');
diff_gosim = [];
diff_gofsim = [];
diff_gopsim = [];
diff_gocsim = [];
diff_coexpr = [];
for i = 1:size(diff_pairs,1)
    if numPrGOterms(diff_pairs(i,1)) || numPrGOterms(diff_pairs(i,2))
        diff_gosim = [diff_gosim 1-sum(prGO(diff_pairs(i,1),:)~=prGO(diff_pairs(i,2),:))/sum(prGO(diff_pairs(i,1),:)|prGO(diff_pairs(i,2),:))];
    end
    if numPrGOFterms(diff_pairs(i,1)) || numPrGOFterms(diff_pairs(i,2))
        diff_gofsim = [diff_gofsim 1-sum(prGO(diff_pairs(i,1),GOaspects=='F')~=prGO(diff_pairs(i,2),GOaspects=='F'))/sum(prGO(diff_pairs(i,1),GOaspects=='F')|prGO(diff_pairs(i,2),GOaspects=='F'))];
    end
    if numPrGOPterms(diff_pairs(i,1)) || numPrGOPterms(diff_pairs(i,2))
        diff_gopsim = [diff_gopsim 1-sum(prGO(diff_pairs(i,1),GOaspects=='P')~=prGO(diff_pairs(i,2),GOaspects=='P'))/sum(prGO(diff_pairs(i,1),GOaspects=='P')|prGO(diff_pairs(i,2),GOaspects=='P'))];
    end
    if numPrGOCterms(diff_pairs(i,1)) || numPrGOCterms(diff_pairs(i,2))
        diff_gocsim = [diff_gocsim 1-sum(prGO(diff_pairs(i,1),GOaspects=='C')~=prGO(diff_pairs(i,2),GOaspects=='C'))/sum(prGO(diff_pairs(i,1),GOaspects=='C')|prGO(diff_pairs(i,2),GOaspects=='C'))];
    end
    notnan = ~isnan(Iexpr(diff_pairs(i,1),:)) & ~isnan(Iexpr(diff_pairs(i,2),:));
    num = sum(notnan);
    if num > 7
        [r,~,~,~] = corrcoef([Iexpr(diff_pairs(i,1),notnan); Iexpr(diff_pairs(i,2),notnan)]','rows','pairwise');
        if ~isnan(r(1,2))
            diff_coexpr = [diff_coexpr r(1,2)];
        end
    end
end

if ~isempty(dir('diffGene_gosim.mat'))
    % load GO similarity of protein partners of products of different genes
    load diffGene_gosim.mat
else
    % calculates GO similarity of all protein pairs that have no
    % common partner. May take several hours
    disp('Calculating GO similarity of pairs of proteins interacting with products of different genes');
    disp('Results will be saved in file diffGene_gosim.mat');
    numNeighbors = sum(I,2);
    maxsum = cumsum(1:numGenes-1);
    maxnum = maxsum(numGenes-1);
    diffGene_gosim = zeros(1,maxnum);
    diffGene_gofsim = zeros(1,maxnum);
    diffGene_gopsim = zeros(1,maxnum);
    diffGene_gocsim = zeros(1,maxnum);
    c = 0;
    cf = 0;
    cp = 0;
    cc = 0;
    Pind = GOaspects=='P';
    Cind = GOaspects=='C';
    for i = 1:numGenes-1
        if numNeighbors(i)>0
            for j = i+1:numGenes
                if numNeighbors(j) > 0
                    if numCommonPartners(i,j) == 0
                        if numPrGOterms(i) || numPrGOterms(j)
                            c = c + 1;
                            diffGene_gosim(c) = 1-sum(prGO(i,:)~=prGO(j,:))/sum(prGO(i,:)|prGO(j,:));
                        end
                        if numPrGOFterms(i) || numPrGOFterms(j)
                            cf = cf + 1;
                            diffGene_gofsim(cf) = 1-sum(prGO(i,GOaspects=='F')~=prGO(j,GOaspects=='F'))/sum(prGO(i,GOaspects=='F')|prGO(j,GOaspects=='F'));
                        end
                        if numPrGOPterms(i) || numPrGOPterms(j)
                            cp = cp + 1;
                            diffGene_gopsim(cp) = 1-sum(prGO(i,Pind)~=prGO(j,Pind))/sum(prGO(i,Pind)|prGO(j,Pind));
                        end
                        if numPrGOCterms(i) || numPrGOCterms(j)
                            cc = cc + 1;
                            diffGene_gocsim(cc) = 1-sum(prGO(i,Cind)~=prGO(j,Cind))/sum(prGO(i,Cind)|prGO(j,Cind));
                        end
                    end
                end
            end
        end
    end
    diffGene_gosim = diffGene_gosim(1:c);
    diffGene_gofsim = diffGene_gofsim(1:cf);
    diffGene_gopsim = diffGene_gopsim(1:cp);
    diffGene_gocsim = diffGene_gocsim(1:cc);
    save diffGene_gosim.mat diffGene_gosim diffGene_gofsim diffGene_gopsim diffGene_gocsim
end

if ~isempty(dir('diffGene_coexpr.mat'))
    % load coexpression of protein partners of products of different genes
    load diffGene_coexpr.mat
else
    disp('Calculating coexpression of pairs of proteins interacting with products of different genes');
    disp('Results will be saved in file diffGene_coexpr.mat');
    numNeighbors = sum(I,2);
    maxsum = cumsum(1:numGenes-1);
    maxnum = maxsum(numGenes-1);
    diffGene_coexpr = zeros(1,maxnum);
    c = 0;
    for i = 1:numGenes-1
        if sum(~isnan(Iexpr(i,:)))>7 && numNeighbors(i)>0
            for j = i+1:numGenes
                if numNeighbors(j) > 0
                    if numCommonPartners(i,j) == 0
                        notnan = ~isnan(Iexpr(i,:)) & ~isnan(Iexpr(j,:));
                        num = sum(notnan);
                        if num > 7
                            mn1 = mean(Iexpr(i,notnan));
                            sd1 = std(Iexpr(i,notnan));
                            mn2 = mean(Iexpr(j,notnan));
                            sd2 = std(Iexpr(j,notnan));
                            if sd1>0 && sd2>0
                                c = c + 1;
                                terms1 = (Iexpr(i,notnan)-mn1)/sd1;
                                terms2 = (Iexpr(j,notnan)-mn2)/sd2;
                                diffGene_coexpr(c) = sum(terms1.*terms2)/(num-1);
                            end
                        end
                    end
                end
            end
        end
    end
    diffGene_coexpr = diffGene_coexpr(1:c);
    save diffGene_coexpr.mat diffGene_coexpr
end

fprintf('\n');
disp(['Mean GO similarity of pairs with identical isoform interaction profiles: ' num2str(mean(nodiff_gosim))]);
disp(['Mean GO similarity of pairs with different isoform interaction profiles: ' num2str(mean(diff_gosim))]);
disp(['Mean GO similarity of protein partners of products of different genes: ' num2str(mean(diffGene_gosim))]);
fprintf('\n');
disp(['Mean molecular function similarity of pairs with identical isoform interaction profiles: ' num2str(mean(nodiff_gofsim))]);
disp(['Mean molecular function similarity of pairs with different isoform interaction profiles: ' num2str(mean(diff_gofsim))]);
disp(['Mean molecular function similarity of protein partners of products of different genes: ' num2str(mean(diffGene_gofsim))]);
fprintf('\n');
disp(['Mean biological process similarity of pairs with identical isoform interaction profiles: : ' num2str(mean(nodiff_gopsim))]);
disp(['Mean biological process similarity of pairs with different isoform interaction profiles: ' num2str(mean(diff_gopsim))]);
disp(['Mean biological process similarity of protein partners of products of different genes: ' num2str(mean(diffGene_gopsim))]);
fprintf('\n');
disp(['Mean cellular component similarity of pairs with identical isoform interaction profiles: ' num2str(mean(nodiff_gocsim))]);
disp(['Mean cellular component similarity of pairs with different isoform interaction profiles: ' num2str(mean(diff_gocsim))]);
disp(['Mean cellular component similarity of protein partners of products of different genes: ' num2str(mean(diffGene_gocsim))]);
fprintf('\n');
disp(['Mean coexpression of pairs with identical isoform interaction profiles: ' num2str(mean(nodiff_coexpr))]);
disp(['Mean coexpression of pairs with different isoform interaction profiles: ' num2str(mean(diff_coexpr))]);
disp(['Mean coexpression of protein partners of products of different genes: ' num2str(mean(diffGene_coexpr))]);
fprintf('\n');

plotdata = [mean(nodiff_gosim) mean(diff_gosim) mean(diffGene_gosim); ...
    mean(nodiff_gofsim) mean(diff_gofsim) mean(diffGene_gofsim); ...
    mean(nodiff_gopsim) mean(diff_gopsim) mean(diffGene_gopsim); ...
    mean(nodiff_gocsim) mean(diff_gocsim) mean(diffGene_gocsim)];
sderror = [std(nodiff_gosim)/sqrt(length(nodiff_gosim)) std(diff_gosim)/sqrt(length(diff_gosim)) std(diffGene_gosim)/sqrt(length(diffGene_gosim)); ...
    std(nodiff_gofsim)/sqrt(length(nodiff_gofsim)) std(diff_gofsim)/sqrt(length(diff_gofsim)) std(diffGene_gofsim)/sqrt(length(diffGene_gofsim)); ...
    std(nodiff_gopsim)/sqrt(length(nodiff_gopsim)) std(diff_gopsim)/sqrt(length(diff_gopsim)) std(diffGene_gopsim)/sqrt(length(diffGene_gopsim)); ...
    std(nodiff_gocsim)/sqrt(length(nodiff_gocsim)) std(diff_gocsim)/sqrt(length(diff_gocsim)) std(diffGene_gocsim)/sqrt(length(diffGene_gocsim))];
figure
hold on
hb = bar(plotdata,0.6,'LineWidth',1.5,'EdgeColor','k');
set(hb(1),'facecolor',[51/255 153/255 1]);
set(hb(2),'facecolor',[0 153/255 76/255]);
set(hb(3),'facecolor',[1 0 0]);
ylim([0 0.35]);
set(gca,'XTick',1:4,'XTickLabel',{'   All three\newlineGO aspects','Molecular\newline function','Biological\newline process','  Cellular\newlinecomponent'},'FontSize',28);
ylabel('GO similarity of pairs of proteins','FontSize',28);
legend({'Pairs of proteins interacting with the same subset of isoforms of the same gene',...
    'Pairs of proteins interacting with different subsets of isoforms of the same gene',...
    'Pairs of proteins interacting with protein products of different genes'},'FontSize',28);
set(gca,'YGrid','on');
set(gca,'tickDir','out');
box off
pause(0.1);
for ib = 1:numel(hb)
    xData = hb(ib).XData+hb(ib).XOffset;
    for j = 1:length(xData)
        line([xData(j) xData(j)], [plotdata(j,ib)-sderror(j,ib) plotdata(j,ib)+sderror(j,ib)],'Color','k','LineWidth',1.5);
    end
end
clear hb ib j plotdata sderror xData

pause(1);
plotdata = [mean(nodiff_coexpr); mean(diff_coexpr); mean(diffGene_coexpr)];
sderror = [std(nodiff_coexpr)/sqrt(length(nodiff_coexpr)); std(diff_coexpr)/sqrt(length(diff_coexpr));...
    std(diffGene_coexpr)/sqrt(length(diffGene_coexpr))];
figure
hold on
hb = bar(plotdata,0.6,'FaceColor','c','LineWidth',1.5,'EdgeColor','k');
ylim([0 0.35]);
pause(0.1);
for ib = 1:numel(hb)
    xData = hb(ib).XData+hb(ib).XOffset;
    for j = 1:length(xData)
        line([xData(j) xData(j)], [plotdata(j,ib)-sderror(j,ib) plotdata(j,ib)+sderror(j,ib)],'Color','k','LineWidth',1.5);
    end
end
set(gca,'XTick',1:3,'XTickLabel',{' Interacting with\newlinethe same subset\newline  of isoforms of\newline the same gene',...
    ' Interacting with\newlinedifferent subsets\newline  of isoforms of\newline the same gene',...
    '   Interacting with\newline  protein products\newline of different genes'},'FontSize',17);
ylabel('Co-expression of pairs of proteins','FontSize',22);
set(gca,'YGrid','on');
set(gca,'tickDir','out');
box off
clear hb ib j xData plotdata sderror

%-------------------------------------------------------------------------
% Calculate p-values using bootstrapping

itr = 100000;
disp('Significance in difference between pairs with identical and different isoform interaction profiles:');
bootMeanSim = bootstrapSim(nodiff_gosim,diff_gosim,itr);
p1 = 2*sum((bootMeanSim(:,1)-bootMeanSim(:,2))>=abs(mean(nodiff_gosim)-mean(diff_gosim)))/size(bootMeanSim,1);
disp(['All GO categories bootstrap test: p = ' num2str(p1)]);
bootMeanSim = bootstrapSim(nodiff_gofsim,diff_gofsim,itr);
p2 = 2*sum((bootMeanSim(:,1)-bootMeanSim(:,2))>=abs(mean(nodiff_gofsim)-mean(diff_gofsim)))/size(bootMeanSim,1);
disp(['Molecular function bootstrap test: p = ' num2str(p2)]);
bootMeanSim = bootstrapSim(nodiff_gopsim,diff_gopsim,itr);
p3 = 2*sum((bootMeanSim(:,1)-bootMeanSim(:,2))>=abs(mean(nodiff_gopsim)-mean(diff_gopsim)))/size(bootMeanSim,1);
disp(['Biological process bootstrap test: p = ' num2str(p3)]);
bootMeanSim = bootstrapSim(nodiff_gocsim,diff_gocsim,itr);
p4 = 2*sum((bootMeanSim(:,1)-bootMeanSim(:,2))>=abs(mean(nodiff_gocsim)-mean(diff_gocsim)))/size(bootMeanSim,1);
disp(['Cellular component bootstrap test: p = ' num2str(p4)]);
bootMeanCoexpr = bootstrapSim(nodiff_coexpr,diff_coexpr,itr);
p5 = 2*sum((bootMeanCoexpr(:,1)-bootMeanCoexpr(:,2))>=abs(mean(nodiff_coexpr)-mean(diff_coexpr)))/size(bootMeanCoexpr,1);
disp(['All GO categories bootstrap test (' num2str(itr) ' resamplings): p = ' num2str(p1)]);
disp(['Molecular function bootstrap test (' num2str(itr) ' resamplings): p = ' num2str(p2)]);
disp(['Biological process bootstrap test (' num2str(itr) ' resamplings): p = ' num2str(p3)]);
disp(['Cellular component bootstrap test (' num2str(itr) ' resamplings): p = ' num2str(p4)]);
disp(['Coexpression bootstrap test  (' num2str(itr) ' resamplings): p = ' num2str(p5)]);

itr = 100;
disp('Significance in difference between pairs with different isoform interaction profiles');
disp('and partners of products of different genes:');
bootMeanSim = bootstrapSim(diff_gosim,diffGene_gosim,itr);
p1 = 2*sum((bootMeanSim(:,1)-bootMeanSim(:,2))>=abs(mean(diff_gosim)-mean(diffGene_gosim)))/size(bootMeanSim,1);
disp(['All GO categories bootstrap test: p = ' num2str(p1)]);
bootMeanSim = bootstrapSim(diff_gofsim,diffGene_gofsim,itr);
p2 = 2*sum((bootMeanSim(:,1)-bootMeanSim(:,2))>=abs(mean(diff_gofsim)-mean(diffGene_gofsim)))/size(bootMeanSim,1);
disp(['Molecular function bootstrap test: p = ' num2str(p2)]);
bootMeanSim = bootstrapSim(diff_gopsim,diffGene_gopsim,itr);
p3 = 2*sum((bootMeanSim(:,1)-bootMeanSim(:,2))>=abs(mean(diff_gopsim)-mean(diffGene_gopsim)))/size(bootMeanSim,1);
disp(['Biological process bootstrap test: p = ' num2str(p3)]);
bootMeanSim = bootstrapSim(diff_gocsim,diffGene_gocsim,itr);
p4 = 2*sum((bootMeanSim(:,1)-bootMeanSim(:,2))>=abs(mean(diff_gocsim)-mean(diffGene_gocsim)))/size(bootMeanSim,1);
disp(['Cellular component bootstrap test: p = ' num2str(p4)]);
bootMeanCoexpr = bootstrapSim(diff_coexpr,diffGene_coexpr,itr);
p5 = 2*sum((bootMeanCoexpr(:,1)-bootMeanCoexpr(:,2))>=abs(mean(diff_coexpr)-mean(diffGene_coexpr)))/size(bootMeanCoexpr,1);
disp(['All GO categories bootstrap test (' num2str(itr) ' resamplings): p = ' num2str(p1)]);
disp(['Molecular function bootstrap test (' num2str(itr) ' resamplings): p = ' num2str(p2)]);
disp(['Biological process bootstrap test (' num2str(itr) ' resamplings): p = ' num2str(p3)]);
disp(['Cellular component bootstrap test (' num2str(itr) ' resamplings): p = ' num2str(p4)]);
disp(['Coexpression bootstrap test  (' num2str(itr) ' resamplings): p = ' num2str(p5)]);
