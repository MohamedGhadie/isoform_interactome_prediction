% This script first predicts the isoform interactome from the reference
% interactome, either HI-II-14 or IntAct interactome. Then it calculates GO
% association similarity and co-expression for pairs of proteins interacting 
% with the same subset of isoforms of the same gene, pairs of proteins
% interacting with different subsets of isoforms of the same gene, and
% pairs of proteins interacting with protein products of different genes.
% GO similarity is calculated as the fraction (Jaccard similarity index) of
% GO terms shared by the two proteins. Co-expression is calculated as
% Pearson's correlation.

% select interactome: HI-II-14 or IntAct
interactome = 'IntAct';

% keep only interactions reported this many times or more
numTimesReported = 2;

% remove self interactions,: 1 for yes, 0 otherwise
removeSelfInteractions = 1;     

% laod processed data (if exists) from .mat file:
% 1 for yes, 0 to process interactome data from scratch
load_processed_data = 1;

% save interactome processed data to .mat file: 
% 1 for yes, 0 otherwise
save_processed_data = 1;

% Calculate p-values for results: 1 for yes, 0 otherwise
calculate_pvalues = 1;

% processed data directory where interactome processed data will be saved
processed_data_dir = 'interactome_processed/';

% load an process interactome data
[I, PPIs, spID, genes, domains, DDIs, refDomMap, domRefPos, domAltMap, domAltPos, prSeq, ...
isoSeq, isoNames, altIsoforms, numAltIsoforms, maxIsoform, isoInterDomains, domI, domPrI, ...
numDDImap, numCommonPartners, inIsoInteractome, ref_ref_interactions, alt_ref_interactions, ...
numLosingIsoforms, fracLosingIsoforms, numIsoPairs, fracDiffIsoProfilePerGene, ...
avgIsoProfileDistPerGene, allDist] ...
= process_interactome(interactome, load_processed_data, save_processed_data, ...
                      numTimesReported, removeSelfInteractions, processed_data_dir);

% display results from domain-resolved interactome and predicted isoform interactome 
display_isoform_interactome_results(I,numDDImap,numAltIsoforms,numIsoPairs,inIsoInteractome,...
                                    ref_ref_interactions,alt_ref_interactions,numLosingIsoforms,...
                                    fracLosingIsoforms,fracDiffIsoProfilePerGene,avgIsoProfileDistPerGene, allDist);

numGenes = size(I,1);

% label pairs of reference proteins based on isoform interaction profiles
[pairs,selected,nodiffTargets,subsetTargets,changeoverTargets] = label_iso_partner_pairs(spID,I,domI,domPrI,isoInterDomains,numDDImap,maxIsoform);

% select pairs with identical isoform interaction profiles
nodiff_pairs = pairs(selected(:,1) & ~selected(:,2) & ~selected(:,3),:);
% select pairs with different isoform interaction profiles
diff_pairs = pairs(~selected(:,1) & (selected(:,2)|selected(:,3)),:);

disp([num2str(size(nodiff_pairs,1)) ' protein pairs interacting with the same subset of isoforms of the same gene']);
disp([num2str(size(diff_pairs,1)) ' protein pairs interacting with different subsets of isoforms of the same gene']);

%-------------------------------------------------------------------------
% load gene ontology associations

if exist('gene_association.goa_ref_human.xlsx', 'file') == 2
    disp('Loading gene ontology associations');
    [~,gene_associations] = xlsread('gene_association.goa_ref_human.xlsx');
% GO data file may be split into two files if too large
elseif (exist('gene_association.goa_ref_human_a.xlsx', 'file') == 2) && (exist('gene_association.goa_ref_human_b.xlsx', 'file') == 2)
    disp('Loading gene ontology associations');
    [~,gene_association_a] = xlsread('gene_association.goa_ref_human_a.xlsx');
    [~,gene_association_b] = xlsread('gene_association.goa_ref_human_b.xlsx');
    gene_associations = [gene_association_a; gene_association_b];
else
    disp('Human gene ontology file (ene_association.goa_ref_human.xlsx) not found. Exiting script');
    return
end

% GO aspect for each term is one of the following 
% F: molecular function
% P: biological process
% C: cellular component
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

% number of GO terms per GO aspect per proteins
numPrGOterms = sum(prGO,2);
numPrGOFterms = sum(prGO(:,GOaspects=='F'),2);
numPrGOPterms = sum(prGO(:,GOaspects=='P'),2);
numPrGOCterms = sum(prGO(:,GOaspects=='C'),2);

%-------------------------------------------------------------------------
% load gene tissue expression data

tissueExpressionFile = 'E-MTAB-513.tsv.txt';
if exist(tissueExpressionFile, 'file') == 2
    disp('Loading gene tissue-expression data');
    tdfread(tissueExpressionFile);
else
    disp('Human gene ontology file (ene_association.goa_ref_human.xlsx) not found. Exiting script');
    return
end
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

%-------------------------------------------------------------------------
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

%-------------------------------------------------------------------------
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

%-------------------------------------------------------------------------
% select file to load GO similarity and co-expression for protein pairs interacting
% protein products of different genes

if strcmpi(interactome,'HI-II-14')
    diffGene_gosim_file = [processed_data_dir 'HI-II-14_diffGene_gosim.mat'];
    diffGene_coexpr_file = [processed_data_dir 'HI-II-14_diffGene_coexpr.mat'];
elseif strcmpi(interactome,'IntAct')
    diffGene_gosim_file = [processed_data_dir 'IntAct_diffGene_gosim.mat'];
    diffGene_coexpr_file = [processed_data_dir 'IntAct_diffGene_coexpr.mat'];
end

% load GO similarity for protein pairs interacting with protein products of
% different genes if data file exists
if (load_processed_data == 1) && (exist(diffGene_gosim_file, 'file') == 2)
    % load GO similarity for all protein pairs interacting with protein products 
    % of different genes
    disp(['Loading GO similarity results for partners of protein products of different genes from file ' diffGene_gosim_file]);
    load(diffGene_gosim_file)
else
    % calculates GO similarity of all protein pairs interacting with protein
    % products of different genes (Takes a long time)
    disp('Calculating GO similarity of pairs of proteins interacting with products of different genes');
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
        if ~mod(i,1000)
            disp(['- Calculated GO similarity for ' num2str(i) ' genes out of ' num2str(numGenes) ' genes']);
        end
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
    if save_processed_data == 1
        save(diffGene_gosim_file, 'diffGene_gosim', 'diffGene_gofsim', 'diffGene_gopsim', 'diffGene_gocsim')
        disp(['GO similarity results for partners of protein products of different genes were saved in file ' diffGene_gosim_file]);
    end
end

% load co-expression for protein pairs interacting with protein products of
% different genes if data file exists
if (load_processed_data == 1) && (exist(diffGene_coexpr_file, 'file') == 2)
    % load coexpression of protein partners of products of different genes
    disp(['Loading co-expression results for partners of protein products of different genes from file ' diffGene_coexpr_file]);
    load(diffGene_coexpr_file)
else
    % calculates co-expression for all protein pairs that have no
    % common partner. Takes a long time
    disp('Calculating tissue co-expression of pairs of proteins interacting with products of different genes');
    numNeighbors = sum(I,2);
    maxsum = cumsum(1:numGenes-1);
    maxnum = maxsum(numGenes-1);
    diffGene_coexpr = zeros(1,maxnum);
    c = 0;
    for i = 1:numGenes-1
        if ~mod(i,1000)
            disp(['- Calculated tissue co-expression for ' num2str(i) ' genes out of ' num2str(numGenes) ' genes']);
        end
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
    if save_processed_data == 1
        save(diffGene_coexpr_file, 'diffGene_coexpr');
        disp(['Coexpression results for partners of protein products of different genes were saved in file ' diffGene_coexpr_file]);
    end
end

%-------------------------------------------------------------------------
% display GO similarity and co-expression results

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

%-------------------------------------------------------------------------
% Plots

% plot GO similarity results
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
hb = bar(plotdata,0.8,'LineWidth',1.5,'EdgeColor','none');
set(hb(1),'facecolor',[51/255 153/255 1]);
set(hb(2),'facecolor',[0 204/255 204/255]);
set(hb(3),'facecolor',[1 102/255 102/255]);
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

% plot co-expression results
pause(1);
plotdata = [mean(nodiff_coexpr); mean(diff_coexpr); mean(diffGene_coexpr)];
sderror = [std(nodiff_coexpr)/sqrt(length(nodiff_coexpr)); std(diff_coexpr)/sqrt(length(diff_coexpr));...
    std(diffGene_coexpr)/sqrt(length(diffGene_coexpr))];
figure
hold on
hb = bar(plotdata,0.6,'EdgeColor','none');
set(hb,'facecolor',[0 204/255 204/255]);
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
% Calculate p-values using bootstrap test

if calculate_pvalues == 1
    itr = 100000;
    disp('Significance in difference between pairs with identical and different isoform interaction profiles:');
    bootMeanSim = bootstrapSim(nodiff_gosim,diff_gosim,itr);
    p1 = 2*sum((bootMeanSim(:,1)-bootMeanSim(:,2))>=abs(mean(nodiff_gosim)-mean(diff_gosim)))/size(bootMeanSim,1);
    disp(['All GO categories bootstrap test (' num2str(itr) ' resamplings): p = ' num2str(p1)]);
    bootMeanSim = bootstrapSim(nodiff_gofsim,diff_gofsim,itr);
    p2 = 2*sum((bootMeanSim(:,1)-bootMeanSim(:,2))>=abs(mean(nodiff_gofsim)-mean(diff_gofsim)))/size(bootMeanSim,1);
    disp(['Molecular function bootstrap test (' num2str(itr) ' resamplings): p = ' num2str(p2)]);
    bootMeanSim = bootstrapSim(nodiff_gopsim,diff_gopsim,itr);
    p3 = 2*sum((bootMeanSim(:,1)-bootMeanSim(:,2))>=abs(mean(nodiff_gopsim)-mean(diff_gopsim)))/size(bootMeanSim,1);
    disp(['Biological process bootstrap test (' num2str(itr) ' resamplings): p = ' num2str(p3)]);
    bootMeanSim = bootstrapSim(nodiff_gocsim,diff_gocsim,itr);
    p4 = 2*sum((bootMeanSim(:,1)-bootMeanSim(:,2))>=abs(mean(nodiff_gocsim)-mean(diff_gocsim)))/size(bootMeanSim,1);
    disp(['Cellular component bootstrap test (' num2str(itr) ' resamplings): p = ' num2str(p4)]);
    bootMeanCoexpr = bootstrapSim(nodiff_coexpr,diff_coexpr,itr);
    p5 = 2*sum((bootMeanCoexpr(:,1)-bootMeanCoexpr(:,2))>=abs(mean(nodiff_coexpr)-mean(diff_coexpr)))/size(bootMeanCoexpr,1);
    disp(['Coexpression bootstrap test  (' num2str(itr) ' resamplings): p = ' num2str(p5)]);
    
    itr = 1000;
    disp('Significance in difference between pairs with different isoform interaction profiles');
    disp('and partners of products of different genes:');
    bootMeanSim = bootstrapSim(diff_gosim,diffGene_gosim,itr);
    p1 = 2*sum((bootMeanSim(:,1)-bootMeanSim(:,2))>=abs(mean(diff_gosim)-mean(diffGene_gosim)))/size(bootMeanSim,1);
    disp(['All GO categories bootstrap test (' num2str(itr) ' resamplings): p = ' num2str(p1)]);
    bootMeanSim = bootstrapSim(diff_gofsim,diffGene_gofsim,itr);
    p2 = 2*sum((bootMeanSim(:,1)-bootMeanSim(:,2))>=abs(mean(diff_gofsim)-mean(diffGene_gofsim)))/size(bootMeanSim,1);
    disp(['Molecular function bootstrap test (' num2str(itr) ' resamplings): p = ' num2str(p2)]);
    bootMeanSim = bootstrapSim(diff_gopsim,diffGene_gopsim,itr);
    p3 = 2*sum((bootMeanSim(:,1)-bootMeanSim(:,2))>=abs(mean(diff_gopsim)-mean(diffGene_gopsim)))/size(bootMeanSim,1);
    disp(['Biological process bootstrap test (' num2str(itr) ' resamplings): p = ' num2str(p3)]);
    bootMeanSim = bootstrapSim(diff_gocsim,diffGene_gocsim,itr);
    p4 = 2*sum((bootMeanSim(:,1)-bootMeanSim(:,2))>=abs(mean(diff_gocsim)-mean(diffGene_gocsim)))/size(bootMeanSim,1);
    disp(['Cellular component bootstrap test (' num2str(itr) ' resamplings): p = ' num2str(p4)]);
    bootMeanCoexpr = bootstrapSim(diff_coexpr,diffGene_coexpr,itr);
    p5 = 2*sum((bootMeanCoexpr(:,1)-bootMeanCoexpr(:,2))>=abs(mean(diff_coexpr)-mean(diffGene_coexpr)))/size(bootMeanCoexpr,1);
    disp(['Coexpression bootstrap test  (' num2str(itr) ' resamplings): p = ' num2str(p5)]);
end
