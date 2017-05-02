function DDIs = load_3did_and_domine_DDIs(did3File, domineFile)

DDIs = [];

if exist(domineFile, 'file') == 2
    domineDDIs = load_domine_DDIs(domineFile);
else
    disp(['Domain-domain interaction file ' domineFile ' not found.']);
    return
end

if exist(did3File, 'file') == 2
    DDIs = load_3did_DDIs(did3File);
else
    disp(['Domain-domain interaction file ' did3File ' not found.']);
    return
end

nonoverlapDDIs = zeros(size(domineDDIs,1),1);
for i = 1:size(domineDDIs,1)
    if sum(strcmpi(DDIs(:,1), domineDDIs{i,1}) & strcmpi(DDIs(:,2), domineDDIs{i,2}))==0 ...
            && sum(strcmpi(DDIs(:,1), domineDDIs{i,2}) & strcmpi(DDIs(:,2), domineDDIs{i,1}))==0
        nonoverlapDDIs(i) = 1;
    end
end
DDIs = [DDIs; domineDDIs(nonoverlapDDIs>0,:)];
