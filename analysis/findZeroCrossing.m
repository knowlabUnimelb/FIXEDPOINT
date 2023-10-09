function zc = findZeroCrossing(x, fun, howmany)

if nargin == 2
    howmany = 3;
end

switch howmany
    case 3
        
        % Return indices of zero crossings
        zci = @(v) find(v(:).*circshift(v(:), [-1 0]) <= 0);
        
        % Find zero crossings
        zcs{1} = x(zci(fun(1,:)));
        zcs{2} = x(zci(fun(2,:)));
        zcs{3} = x(zci(fun(3,:)));
        
        % Delete start and end
        if numel(zcs{1}) > 2;  zcs{1}([1 end]) = []; end % zcs{1}(zcs{1} < .05) = [];
        if numel(zcs{2}) > 2;  zcs{2}([1 end]) = []; end % zcs{2}(zcs{2} < .05) = [];
        if numel(zcs{3}) > 2;  zcs{3}([1 end]) = []; end % zcs{3}(zcs{3} < .05) = [];
        
        % Find distance between sets of hypothetical points
        hz = allcomb(zcs{1}, zcs{2}, zcs{3}); % Hypothetical set of zero crossings
        for i = 1:size(hz,1)
            d(i,1) = sum(pdist(hz(i,:)')); % Find distances between hypothetical points
        end
        
        % Keep the minimum points
        [m, midx] = min(d);
        zc = hz(midx,:);
        
    case 4

        % Return indices of zero crossings
        zci = @(v) find(v(:).*circshift(v(:), [-1 0]) <= 0);
        
        % Find zero crossings
        zcs{1} = x(zci(fun(1,:)));
        zcs{2} = x(zci(fun(2,:)));
        zcs{3} = x(zci(fun(3,:)));
        zcs{4} = x(zci(fun(4,:)));
        zcs{5} = x(zci(fun(5,:)));
        zcs{6} = x(zci(fun(6,:)));
        
        % Delete start and end
        if numel(zcs{1}) > 2;  zcs{1}([1 end]) = []; end% zcs{1}(zcs{1} < .05) = [];
        if numel(zcs{2}) > 2;  zcs{2}([1 end]) = []; end% zcs{2}(zcs{2} < .05) = [];
        if numel(zcs{3}) > 2;  zcs{3}([1 end]) = []; end% zcs{3}(zcs{3} < .05) = [];
        if numel(zcs{4}) > 2;  zcs{4}([1 end]) = []; end% zcs{3}(zcs{3} < .05) = [];
        if numel(zcs{5}) > 2;  zcs{5}([1 end]) = []; end% zcs{3}(zcs{3} < .05) = [];
        if numel(zcs{6}) > 2;  zcs{6}([1 end]) = []; end% zcs{3}(zcs{3} < .05) = [];
        
        % Find distance between sets of hypothetical points
        hz = allcomb(zcs{1}, zcs{2}, zcs{3}, zcs{4}, zcs{5}, zcs{6}); % Hypothetical set of zero crossings
        for i = 1:size(hz,1)
            d(i,1) = sum(pdist(hz(i,:)')); % Find distances between hypothetical points
        end
        
        % Keep the minimum points
        [m, midx] = min(d);
        zc = hz(midx,:);
end