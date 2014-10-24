function [theNodes, thePolys] = SmolyakGrid(n,node_mu,poly_mu)

if nargin<3, poly_mu = node_mu; end

%%%
% Dimensions
d = numel(n);
ngroups = log2(n-1) + 1;
N = max(ngroups);


%%%
% Node parameter
node_q = max(node_mu);
node_isotropic = (isscalar(node_mu));

if node_isotropic
    node_mu = zeros(1,d);
end

%%%
% Polynomial parameter
poly_q = max(poly_mu);
poly_isotropic = (isscalar(poly_mu));

if poly_isotropic
    poly_mu = zeros(1,d);
end


%%%
% Make grid that identifies node groups, save it in "g"
k = (0:2^(N-1))';
g = zeros(size(k));

p2 = 2;
for it=N:-1:3
    odds = logical(mod(k,p2));
    g(odds) = it;
    k(odds) = 0;
    p2 = 2*p2;
end

g(1) = 2;
g(end) = 2;
g((1+end)/2) = 1;

clear k p2 it N

%%%
% Make disjoint sets
nodeMapping = cell(1,d);
polyMapping = cell(1,d);
nodeIndex = cell(1,d);


for i=1:d
    nodeMapping{1,i} = g(g <= ngroups(i));
    polyMapping{1,i} = sort(nodeMapping{1,i});
    nodeIndex{1,i} = (1:n(i))';
end



%%%
% Set up nodes for first dimension 
nodeSum = nodeMapping{1};
theNodes = nodeIndex{1};
if ~node_isotropic
    isValid = nodeSum <= (node_mu(1) + 1);
    nodeSum = nodeSum(isValid);
    theNodes = theNodes(isValid);
end

%%%
% Set up polynomials for first dimension 
polySum = polyMapping{1};
thePolys = nodeIndex{1};
if ~poly_isotropic
    isValid = polySum <= (poly_mu(1) + 1);
    polySum = polySum(isValid);
    thePolys = thePolys(isValid);
end



%%%
% Compute the grid
for k=2:d
    [theNodes, nodeSum] = ndgrid2(theNodes,nodeSum,nodeMapping{k},k + node_q,node_mu(k));
    [thePolys, polySum] = ndgrid2(thePolys,polySum,polyMapping{k},k + poly_q,poly_mu(k));
end




end



function [newIndices,newSum] = ndgrid2(Indices,groupSum,newGroup,q,qk)

assert(size(Indices,1) == numel(groupSum),...
    'nodeIndices and groupSum must have same number of rows')



[idx,idy] = ndgrid(1:numel(groupSum), 1:numel(newGroup));

idx = idx(:);
idy = idy(:);


newSum = groupSum(idx) + newGroup(idy);

if qk == 0 %isotropic
    isValid = newSum <=q;
else %anisotropic
    isValid = (newSum <=q) & (newGroup(idy) <= qk+1);
end


newSum = newSum(isValid);
newIndices = [Indices(idx(isValid),:)  idy(isValid)];


end





