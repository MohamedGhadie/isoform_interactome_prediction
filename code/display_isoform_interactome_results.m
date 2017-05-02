function [] = display_isoform_interactome_results (I,numDDImap,numAltIsoforms,numIsoPairs,inIsoInteractome,ref_ref_interactions,alt_ref_interactions,numLosingIsoforms,fracLosingIsoforms,fracDiffIsoProfilePerGene,avgIsoProfileDistPerGene, allDist)

%-------------------------------------------------------------------------
% reference interactome

numGenes = size(I,1);
disp([num2str(sum(numAltIsoforms>0)) ' out of the ' num2str(numGenes) ' reference proteins have atleast one alternative isoform']);

figure
bar([sum(numAltIsoforms==0) sum(numAltIsoforms==1) sum(numAltIsoforms==2) sum(numAltIsoforms==3) sum(numAltIsoforms==4) ...
    sum(numAltIsoforms==5) sum(numAltIsoforms==6) sum(numAltIsoforms==7) sum(numAltIsoforms==8) sum(numAltIsoforms==9) ...
    sum(numAltIsoforms>9)]);
set(gca,'XTick',1:11,'XTickLabel',{'0','1','2','3','4','5','6','7','8','9','\geq 10'});
xlabel('Number of alternative isoforms');
ylabel('Number of reference proteins');
set(gca,'tickDir','out');
box off

%-------------------------------------------------------------------------
% domain-resolved interactome

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

%-------------------------------------------------------------------------
% predicted isoform interactome

fprintf('\nPredicted isoform interactome statistics:\n');
disp([num2str(length(unique(ref_ref_interactions(:,1:2)))) ' reference proteins and ' num2str(length(unique(alt_ref_interactions(:,1)))) ' alternative isoforms']);
disp([num2str(size(ref_ref_interactions,1)) ' known reference interactions']);
disp([num2str(sum(strcmpi(alt_ref_interactions(:,3),'retained'))) ' predicted retained isoform interactions']);
disp([num2str(sum(strcmpi(alt_ref_interactions(:,3),'lost'))) ' predicted lost isoform interactions']);
fprintf('\n');

disp([num2str(sum(numAltIsoforms(sum(numDDImap)>0)>0)) ' genes have 2 or more isoforms (including reference)']);
disp([num2str(sum(numLosingIsoforms>0)) ' genes have at least one isoform losing an interaction']);
disp(['Avg. fraction per gene of isoforms losing interactions: ' num2str(mean(fracLosingIsoforms(~isnan(fracLosingIsoforms))))]); 
disp(['Total number of isoform pairs with different interaction profiles: ' num2str(sum(allDist>0))]);
disp(['Total number of isoform pairs with identical interaction profiles: ' num2str(sum(allDist==0))]);
disp(['Avg. fraction per gene of isoform pairs with different profiles: ' num2str(mean(fracDiffIsoProfilePerGene(~isnan(fracDiffIsoProfilePerGene))))]); 
disp(['Avg. fraction per gene of isoform pairs with identical profiles: ' num2str(mean(1-fracDiffIsoProfilePerGene(~isnan(fracDiffIsoProfilePerGene))))]); 
disp(['Avg. pair hamming distance over all genes: ' num2str(mean(allDist))]);

figure
histogram((numAltIsoforms(inIsoInteractome>0)+1));
%xlim([-0.5 12]);
%ylim([0 300]);
xlabel('Number of isoforms');
ylabel('Number of genes in isoform interactome');
box off

figure
histogram(numIsoPairs(inIsoInteractome>0),-0.5:1:700.5);
xlim([-0.5 701.5]);
xlabel('Number of isoform pairs');
ylabel('Number of genes in isoform interactome');
box off

figure
histogram(numLosingIsoforms(~isnan(numLosingIsoforms)));
%xlim([-0.5 12]);
xlabel('Number of isoforms losing interaction');
ylabel('Number of genes with two or more isoforms in isoform interactome');
box off

figure
histogram(fracLosingIsoforms(~isnan(fracLosingIsoforms)),-0.05:0.1:1.05);
xlim([-0.05 1.05]);
xlabel('Fraction of isoforms losing interaction');
ylabel('Number of genes with two or more isoforms in isoform interactome');
box off

figure
histogram(fracDiffIsoProfilePerGene(~isnan(fracDiffIsoProfilePerGene)),-0.05:0.1:1.05);
xlim([-0.05 1.05]);
xlabel('Fraction of isoform pairs with different interaction profiles');
ylabel('Number of genes with two or more isoforms in isoform interactome');
box off

figure
histogram(avgIsoProfileDistPerGene(~isnan(avgIsoProfileDistPerGene)),-0.05:0.1:1.05);
xlim([-0.05 1.05]);
xlabel('Average hamming distance of isoform pair');
ylabel('Number of genes with two or more isoforms in isoform interactome');
box off
