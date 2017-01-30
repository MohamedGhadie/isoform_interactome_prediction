function isoSeq = load_alt_protein_seq(filename)

isoformData = fastaread(filename);
isoSeq = {};
for i = 1:size(isoformData,1)
    str = strsplit(isoformData(i).Header,'|');
    if ~isempty(strfind(str{2},'-'))
        isoSeq = [isoSeq; {str{2},isoformData(i).Sequence}];
    end
end