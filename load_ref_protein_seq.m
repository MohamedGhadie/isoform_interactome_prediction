function [headers, prSeq] = load_ref_protein_seq(filename)

numLines = 30000;
prSeq = cell(numLines,3);
fid = fopen(filename);
headers = strsplit(fgetl(fid),'\t');
i = 0;
while ~feof(fid)
    i = i + 1;
    str = strsplit(fgetl(fid),'\t');
    prSeq{i,1} = strtrim(str{1});
    prSeq{i,2} = strtrim(str{3});
    prSeq{i,3} = str2double(strtrim(str{2}));
end
fclose(fid);
prSeq = prSeq(1:i,:);