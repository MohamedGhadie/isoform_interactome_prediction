function [I, PPIs, spID, genes, domains, DDIs, refDomMap, domRefPos, domAltMap, domAltPos, prSeq, ...
    isoSeq, isoNames, altIsoforms, numAltIsoforms, maxIsoform, isoInterDomains, domI, domPrI, ...
    numDDImap, numCommonPartners, inIsoInteractome, ref_ref_interactions, alt_ref_interactions, ...
    numLosingIsoforms, fracLosingIsoforms, numIsoPairs, fracDiffIsoProfilePerGene, ...
    avgIsoProfileDistPerGene, allDist] ...
    = process_interactome(interactome, load_processed_data, save_processed_data,...
                          numTimesReported, removeSelfInteractions, processed_data_dir)

processed_data_exists = 0;
if strcmpi(interactome,'HI-II-14')
    if (load_processed_data == 1) && (exist([processed_data_dir 'HI-II-14_data.mat'], 'file') == 2)
        disp(['Loading interactome processed data from file ' processed_data_dir 'HI-II-14_data.mat']);
        processed_data_exists = 1;
        load([processed_data_dir 'HI-II-14_data.mat']);
    else
        disp('Loading and processing interactome data');
        interactomeFile = 'HI-II-14.xlsx';
        isoformInteractomeFile = [processed_data_dir 'HI-II-14_isoform_interactome.txt'];
        spEntrezMapFile = 'HI-II-14_spEntrezMap.tab';
        [I, PPIs, spID, genes] = load_interactome(interactomeFile, [], spEntrezMapFile, numTimesReported, removeSelfInteractions, interactome);
        outputfile = [processed_data_dir 'HI-II-14_data.mat'];
    end
elseif strcmpi(interactome,'IntAct')
    if (load_processed_data == 1) && (exist([processed_data_dir 'IntAct_data.mat'], 'file') == 2)
        disp(['Loading interactome processed data from file ' processed_data_dir 'IntAct_data.mat']);
        processed_data_exists = 1;
        load([processed_data_dir 'IntAct_data.mat']);
    else
        disp('Loading and processing interactome data');
        interactomeFile = 'intact.txt';
        isoformInteractomeFile = [processed_data_dir 'IntAct_isoform_interactome.txt'];
        spGeneMapFile = 'IntAct_spGeneMap.txt';
        [I, PPIs, spID, genes] = load_interactome(interactomeFile, spGeneMapFile, [], numTimesReported, removeSelfInteractions, interactome);
        outputfile = [processed_data_dir 'IntAct_data.mat'];
    end
else
    disp('Invalid interactome name. Exiting script');
    return
end

if processed_data_exists == 0
    numPPI = sum(sum(triu(I)));
    numGenes = size(I,1);
    did3File = '3did_flat.txt';
    domineFile = 'domine_interactions.xlsx';
    processedDDIfile = [processed_data_dir 'DDIs.txt'];
    disp('Loading domain-domain interactions');
    if exist(processedDDIfile, 'file') == 2
        s = tdfread(processedDDIfile);
        DDIs = cell(size(s.dom1,1),2);
        for i = 1:size(s.dom1,1)
            DDIs(i,:) = {strtrim(s.dom1(i,:)), strtrim(s.dom2(i,:))};
        end
        clear dom1 dom2
    elseif (exist(did3File, 'file') == 2) && (exist(domineFile, 'file') == 2)
        DDIs = load_3did_and_domine_DDIs(did3File, domineFile);
        fid = fopen(processedDDIfile,'w');
        fprintf(fid,'dom1\tdom2\n');
        for i = 1:size(DDIs,1)
            fprintf(fid,[DDIs{i,1} '\t' DDIs{i,2} '\n']);
        end
        fclose(fid);
    else
        disp('Domain-domain interaction file not found. Exiting script');
        return
    end
    domains = unique([DDIs(:,1); DDIs(:,2)]);
    numDom = length(domains);
    
    if exist('hmmer_canonical_full.domtab.txt', 'file') == 2
        disp('Loading structural domains for reference proteins');
        [refDomMap, domRefPos] = load_hmmscan_ref('hmmer_canonical_full.domtab.txt');
    else
        disp('Reference-isoform Domain-mapping file (hmmer_isoforms.domtab.txt) not found. Exiting script');
        return
    end
    
    if exist('hmmer_isoforms.domtab.txt', 'file') == 2
        disp('Loading structural domains for alternative isoforms');
        [domAltMap, domAltPos] = load_hmmscan_alt('hmmer_isoforms.domtab.txt');
    else
        disp('Alternative-isoform Domain-mapping file (hmmer_isoforms.domtab.txt) not found. Exiting script');
        return
    end
    
    if exist('human_protein_sequences.tab', 'file') == 2
        disp('Loading reference protein sequences');
        [headers, prSeq] = load_ref_protein_seq('human_protein_sequences.tab');
    else
        disp('Human reference protein sequence file (human_protein_sequences.tab) not found. Exiting script');
        return
    end
    
    if exist('human_protein_can_iso_sequences.fasta', 'file') == 2
        disp('Loading alternative isoform sequences');
        isoSeq = load_alt_protein_seq('human_protein_can_iso_sequences.fasta');
    else
        disp('Human alternative protein sequence file (human_protein_can_iso_sequences.fasta) not found. Exiting script');
        return
    end
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
    numCommonPartners = double(I)*double(I);
    
    disp('Predicting isoform interactome');
    [inIsoInteractome, ref_ref_interactions, alt_ref_interactions, ~, ~] = predict_isoform_interactome(spID,domains,PPIs,I,domI,domPrI,numDDImap,altIsoforms,isoInterDomains,interactome,isoformInteractomeFile);
    
    disp('Comparing isoform interaction profiles');
    [numLosingIsoforms, fracLosingIsoforms, numIsoPairs, fracDiffIsoProfilePerGene, avgIsoProfileDistPerGene, allDist] = compare_isoform_interactions(I,domI,domPrI,numDDImap,altIsoforms,numAltIsoforms,isoInterDomains);
    
    if save_processed_data == 1
        if exist(processed_data_dir, 'dir') ~= 7
            mkdir(processed_data_dir);
        end
        save(outputfile,'interactome','interactomeFile','isoformInteractomeFile', 'numTimesReported', 'removeSelfInteractions',...
            'PPIs', 'I', 'spID', 'genes', 'numPPI', 'numGenes', 'DDIs', 'domains', 'numDom', 'refDomMap',...
            'domRefPos', 'domAltMap', 'domAltPos', 'headers', 'prSeq', 'isoSeq', 'isoNames', 'altIsoforms', 'numAltIsoforms',...
            'maxIsoform', 'isoInterDomains', 'domI', 'domPrI', 'numDDImap', 'numCommonPartners', 'inIsoInteractome',...
            'ref_ref_interactions', 'alt_ref_interactions', 'numLosingIsoforms', 'fracLosingIsoforms', 'numIsoPairs',...
            'fracDiffIsoProfilePerGene', 'avgIsoProfileDistPerGene', 'allDist');
        if strcmpi(interactome,'HI-II-14')
            save(outputfile, 'spEntrezMapFile', '-append');
        elseif strcmpi(interactome,'IntAct')
            save(outputfile, 'spGeneMapFile', '-append');
        end
        disp(['Interactome processed data saved to file ' outputfile]);
    end
end