function idx = mstrfind(C, S)

if ischar(S)
    S = {S};
end

idx = cellfun(@(s)find(strcmp(C, s)), S, 'UniformOutput', false);
idx(cellfun(@isempty, idx)) = [];
idx = cell2mat(idx);