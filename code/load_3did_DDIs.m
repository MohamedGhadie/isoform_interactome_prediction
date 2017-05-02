function DDIs = load_3did_DDIs(filename)

numLines = 15000;
numDDI = 0;
DDIs = cell(numLines,2);
fid = fopen(filename);
tline = fgetl(fid);
while ischar(tline)
    if length(tline)>4
        if strcmpi(tline(1:4),'#=ID')
            str = strsplit(tline,'\t');
            if length(str)>1
                dom1 = '';
                dom2 = '';
                for i = 1:length(str)
                    if ~isempty(strfind(str{i},'@Pfam'))
                        PF = strfind(str{i},'PF');
                        dot = find(str{i}=='.');
                        if(length(PF)==1) && (length(dot)==1)
                            if isempty(dom1)
                                dom1 = str{i}(PF:(dot-1));
                            else
                                dom2 = str{i}(PF:(dot-1));
                                break
                            end
                        end
                    end
                end
                if ~isempty(dom1) && ~isempty(dom2)
                    numDDI = numDDI + 1;
                    DDIs{numDDI,1} = dom1;
                    DDIs{numDDI,2} = dom2;
                end
            end
        end
    end
    tline = fgetl(fid);
end
fclose(fid);

DDIs = DDIs(1:numDDI,:);
