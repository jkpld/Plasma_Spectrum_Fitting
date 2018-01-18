function [c,ceq,g,geq,initCorrections] = OH_Constraint(p, delta, idx)
% Ensure peaks next to each other are within some range of each other
% g : should be nx3 with the indices and values that can be passed to
% sparse -> sparse(g(:,1),g(:,2),g(:,3))
ceq = [];
geq = [];

% N = numel(A);
p = reshape(p,[],3);
A = p(:,1);
mu = p(:,2);

idx = find(idx);
N = numel(idx);

A = A(idx);
mu = mu(idx);

% Extract the OH peak at 308, which has a tail at longer wavelengths.
idx2 = mu > 307 & mu < 320;
A = A(idx2);
mu = mu(idx2);
idx = idx(idx2);

[~,idx3] = sort(mu);
A = A(idx3);
idx = idx(idx3); 

n = numel(A);
A = A(:);
% constraints
% 1 :  A(i+i) < (1-delta)*A(i)

c = A(2:n) - (1-delta)*A(1:n-1);

% gradient (as sparse matrix)
i = (1:(n-1));
j = 1:n-1;

i = [i, i];
j = idx([j, j+1]);
v = [-(1-delta)*ones(1,n-1), ones(1,n-1)];
% g = sparse(j,i,[v1,v2,v2,v1], N, max(i));

g = [j(:),i(:),v(:)];


if nargout > 4
    initCorrections = zeros(N,1);
    
    An = A;
    for i = 2:n
        An(i) = min(An(i-1)*(1-delta), An(i));
    end

    initCorrections(idx(2:end)) = A(2:end) - An(2:end);
end

end