%% basis class
% Defines a class to represent a multi-dimensional collocation basis
%
% Objects created by this class are of type 'handle', useful for passing the basis by
% reference.
%
% This class is defined to combine _d_ unidimensional bases into a multidimensional one.
% The unidimensional basis is defined as any of these bases:
%
% * basisChebyshev
% * basisSpline
% * basisLinear
%
% Objects of class |basis| have the following properties:
%
% * |d|:    scalar, dimension of basis
% * |n|:    1.d vector of integers, number of nodes per dimension
% * |a|:    1.d vector, lower bounds of variables
% * |b|:    1.d vector, upper bounds of variables
% * |B1|:   1.d array of unidimensional bases
% * |nodes|: matrix with d columns,basis nodes
% * |opts|:  structure with options for chosen type
% * |type|:  basis type, i.e. 'cheb','spli', or 'lin'
%
% Object of class |basis| have the following methods:
%
% * |basis|: class constructor
% * |ExpandBasis|: computes auxilliary fields needed for combining the d 1-basis
% * |Interpolation|: returns the interpolation matrix
% * |nearestNode|: returns the basis node that is closest to input argument
% * |plot|: plots the basis functions for each dimension
% * |disp|: prints a summary of the basis
% * |WarnOutOfBounds|: change behavior of basis if extrapolating
% 
%
% *To date, only basisChebyshev has been implemented*, so calling the function with |type|
% 'spli' or 'lin' returns an error.
%
% Last updated: October 4, 2014.
%
%
% Copyright (C) 2013-2014 Randall Romero-Aguilar
%
% Licensed under the MIT license, see LICENSE.txt







%% 
classdef basis < handle
    
    properties (SetAccess = protected)
        d           % dimension of basis
        n           % number of nodes
        a           % lower bounds
        b           % upper bounds
        B1          % 1.d array of one-dimensional bases
        nodes       % nodes for d-dimensions
        opts        % options for basis of given type
        type        % type of basis
    end
    
    

    methods
        
        %% basis
        function B = basis(...
                n,...           1.d, number of nodes
                a,...           1.d, lower bounds
                b,...           1.d, upper bounds
                opts)%     structure with options parameters
            %%% 
            % Constructor for basis
            %
            %   B = basis(n,a,b,opts)
            % 
            % The following inputs are required:
            %
            % * |n|: scalar or 1.d vector of integers, number of nodes
            % * |a|: 1.d vector, lower bounds
            % * |b|: 1.d vector, upper bounds
            %
            % The optional input |opts| is a structure with optional fields (default
            % values in parenthesis).
            %
            % To choose the type of basis and to set names for its dimensions:
            %
            % * |opts.type|: ('cheb'),'spli', or 'lin'; type of basis
            % * |opts.varnames|: 1.d cell of strings, variable names ({'V1',V2'...})
            % 
            % For Chebyshev basis, the type of nodes is specified by
            %
            % * |opts.nodetype|: 'lobatto', ('gaussian') or 'endpoint'
            %
            % For Splines, the order of the piecewise polynomials is set by
            %
            % * |opts.order|: a positive integer,  (3)
            %
            % For multidimensional basis, the following options control how the
            % unidimensional basis are combined:
            %
            % * |opts.method|:  ('tensor'), 'smolyak', 'complete', 'cluster', or 'zcluster'
            % * |opts.degreeParam|: adjust polynomial degrees for all methods (except
            % 'tensor'), default is 3 for 'smolyak' and max(n-1) otherwise.
            % * |opts.nodeParam|: adjust selection of nodes, for methods 'smolyak',
            % 'cluster', and 'zcluster', default is 3 for 'smolyak', 0 otherwise.
            %
            % As of October 4, 2014 the methods 'cluster' and 'zcluster' are only
            % experimental.
            
            %%%
            % Empty constructor
            if nargin==0
                return
            end
            
            
            %%%
            % Basis dimension
            
            B.d = numel(a);
            
            %%%
            % Use same number of nodes in all dimensions, if n = scalar
            if isscalar(n) 
                n = n*ones(1,B.d);
            end
            
            
            if any(a>=b)
                error('Lower bounds must be less than upper bounds: a < b')
            end
            
            %%% 
            % Default values for optional inputs
            
            if nargin < 4
                opts = struct();
            end
            
            if isfield(opts,'type')
                temptype = lower(opts.type);
                opts.type = temptype(1:4);
            else
                opts.type = 'cheb';
            end
            
            
            if ~isfield(opts,'varnames'), opts.varnames = strcat('V',num2cell(48 + (1:B.d)));end
            
            switch opts.type
                case 'cheb'
                    % Set default type of nodes
                    if ~isfield(opts,'nodetype'), opts.nodetype = 'gaussian'; end
                    
                    % If more than one dimension
                    if B.d >1
                        if ~isfield(opts,'method'), opts.method = 'tensor'; end
                        
                        switch opts.method
                            case 'smolyak'
                                % Adjust number of nodes, as needed
                                n_old = n;
                                n = 2.^ceil(log2(n_old-1))+1;
                                if any(n-n_old)
                                    warning('basis:SmolyakNumberOfNodes','For Smolyak expansion, number of nodes should be n=2^k + 1,for some k = 1,2,...')
                                    fprintf('Adjusting number of nodes\n')
                                    fprintf('\t%6s, %6s\n','Old n','New n')
                                    fprintf('\t%6d, %6d\n',[n_old;n])
                                end
                                
                                % adjust nodetype
                                if  ~strcmp(opts.nodetype,'lobatto')
                                    warning('basis:SmolyakNodetype','nodetype must be ''lobatto'' for Smolyak nodes')
                                    opts.nodetype = 'lobatto';
                                end
                                
                                % set default parameters
                                if ~isfield(opts,'nodeParam'), opts.nodeParam = 3; end
                                if ~isfield(opts,'degreeParam'), opts.degreeParam = 3; end
                                
                            case {'complete','cluster','zcluster'}
                                if ~isfield(opts,'nodeParam'), opts.nodeParam = 0; end
                                if ~isfield(opts,'degreeParam'), opts.degreeParam = max(n)-1; end
                                
                            
                        end
                    end
            end
                    
                    
               
            
            %%%
            % Create the 1-basis
                        
            switch opts.type
                case 'cheb'
                    for i = 1: B.d
                        B1_(i) = basisChebyshev(n(i),a(i),b(i),opts.nodetype,opts.varnames{i});
                    end
                otherwise
                    error('Method not yet implemented')
            end
                        
            %%%
            % Pack values in object
            
            
            B.a = a;
            B.b = b;
            B.n = n;
            B.B1 = B1_;
            B.type = opts.type;
            B.opts.nodetype = B1_(1).nodetype;
            B.opts = opts;
            B.ExpandBasis;
        end %basis
        
        
        %% WarnOutOfBounds
        function WarnOutOfBounds(B,value)
            %%%
            % B.WarnOutOfBounds(bool)
            %
            % If bool = true, B.Interpolation(x) will throw an error if extrapolating
            % (i.e., evaluating outside the interpolation hipercube [B.a, B.b]).
            % Otherwise, the function will continue without warning.
            
            [B.B1.WarnOutOfBounds]  = deal(logical(value));
        end
        
        
        
        %% copy
        function F =copy(B)
            switch class(B)
                case 'basis'
                    F = basis;
                case 'funcApprox'
                    F = funcApprox;
            end
            
            for s = fieldnames(B)'
                s1 = char(s);
                F.(s1) = B.(s1);
            end
        end
        
        %% ExpandBasis
        function ExpandBasis(B)
            %%%
            %   B.ExpandBasis
            %
            % ExpandBasis computes nodes for multidimensional basis and other auxilliary fields required to keep track
            % of the basis (how to combine unidimensional nodes bases). It is called by the constructor method, and as
            % such is not directly needed by the user. |ExpandBasis| updates the following fields in basis |B|:
            %
            % * |B.nodes|: matrix with all nodes, one column by dimension
            % * |B.opts.validPhi|: indices to combine unidimensional bases
            % * |B.opts.validX|: indices to combine unidimensional nodes
            %
            % Combining polynomials depends on value of input |opts.method|:
            %
            % * |'tensor'| takes all possible combinations,
            % * |'smolyak'| computes Smolyak basis, given |opts.degreeParam|,
            % * |'complete'|, |'cluster'|, and |'zcluster'| choose polynomials with degrees not exceeding
            % |opts.degreeParam|
            %
            % Expanding nodes depends on value of field |opts.method|. 
            %
            % * 'tensor' and 'complete' take all possible combinations, 
            % * 'smolyak' computes Smolyak basis, given |opts.nodeParam|
            % * 'cluster' and 'zcluster' compute clusters of the tensor nodes based on |B.opts.nodeParam|
            
            
            %%%
            % Allocate variables
            degs = [B.n]-1; % degree of polynomials
            ldeg = cell(1,B.d);
            
            for i=1:B.d
                ldeg{i} = (0:degs(i))';
            end
            deg_grid = gridmake(ldeg{:});   %degree of polynomials
            idxAll = deg_grid + 1;   % index of elements
            
            
            %%%
            % Expanding bases: 
            
            switch lower(B.opts.method)
                case 'tensor'
                    B.opts.validPhi = deg_grid + 1;
                    
                case {'complete','cluster','zcluster'}
                    deg = sum(deg_grid,2);  % degrees of all possible column products
                    degValid = (deg <= B.opts.degreeParam);
                    idx = deg_grid + 1;  % index of entries
                    B.opts.validPhi = idx(degValid,:);
                    
                case 'smolyak'
                    if B.opts.nodeParam < B.opts.degreeParam
                        warning('Smolyak degree param cannot be bigger than node param; adjusting degree');
                        B.opts.degreeParam = B.opts.nodeParam;
                    end
                    
                    
                    if ~strcmp(B.opts.nodetype,'lobatto')
                        warning('Smolyak expansion requires Lobatto nodes')
                    end
                    
                    updateIDX = false;
                    
                    
                    ngroups = log2(B.n-1)+1;
                    Ngroups = max(ngroups);
                    
                    Nodes = cell(Ngroups,1);
                    Nodes{1,1} = [0,1];
                    
                    for k =Ngroups:-1:2
                        N(k) = 2^(k-1)+1;
                        Nodes{k,1} = [...
                            -cos(pi*((1:N(k))'-1)/(N(k)-1)),... %the nodes (with repetition as k changes)
                            repmat(k,N(k),1)];   % group number
                    end
                    
                    Nodes = cell2mat(Nodes);
                    Nodes(abs(Nodes(:,1))<eps) = 0;
                    
                    
                    nodes1 = cell(1,B.d);
                    groups = cell(1,B.d);
                    sortedGroups = cell(1,B.d);
                    
                    for i=1:B.d
                        validGroup = Nodes(:,2)<= ngroups(i);
                        [nodes1{1,i},idtemp] = unique(Nodes(validGroup,1),'first');
                        groups{1,i} = Nodes(idtemp,2);
                        sortedGroups{1,i} = sort(groups{1,i});
                    end
                    
                    
                    nodesAll = gridmake({B.B1.nodes});
                    xGroupAll = gridmake(groups{:});
                    phiGroupAll = gridmake(sortedGroups{:});
                    
                    validNodes = (sum(xGroupAll,2)<= B.d + B.opts.nodeParam);
                    validPhi = (sum(phiGroupAll,2)<= B.d + B.opts.degreeParam);
                    
                    for i=1:B.d
                        B.nodes(:,i) = nodesAll(validNodes,i);
                    end
                    
                    if updateIDX
                        degs = [B.n]-1; % degree of polynomials
                        ldeg = cell(1,B.d);
                        
                        for i=1:B.d
                            ldeg{i} = (0:degs(i))';
                        end
                        deg_grid = gridmake(ldeg{:});   %degree of polynomials
                        idxAll = deg_grid + 1;   % index of elements
                    end
                    
                    
                    B.opts.validX = idxAll(validNodes,:);  % index of entries
                    B.opts.validPhi = idxAll(validPhi,:);  % index of entries
                    
                    
                otherwise
                    error(...
                        'Expansion type must be ''tensor'', ''complete'', ''smolyak'', ''cluster'', or ''zcluster''.')
            end
            
            
            %%%
            % Expanding nodes
            
            switch lower(B.opts.method)
                case {'tensor', 'complete'}
                    B.nodes = gridmake({B.B1.nodes});
                    B.opts.validX = idxAll;  % index of entries
                    
                case 'smolyak'
                    % done in previous switch
                    
                case {'cluster','zcluster'}
                    H = size(B.opts.validPhi,1) + B.opts.nodeParam;  % number of clusters
                    
                    
                    tempNodes = gridmake({B.B1.nodes});
                    
                    if strcmp(B.opts.method,'cluster')
                        [~,Nodes] = kmeans(tempNodes,H);
                    else
                        IDX = kmeans(zscore(tempNodes),H);
                        Nodes = grpstats(tempNodes,IDX,'mean');
                    end
                    
                    
                    for i=1:B.d
                        B(i).Nodes = Nodes(:,i);
                    end
                    
                    
                otherwise
                    error(...
                        'Expansion type must be ''tensor'', ''complete'', ''smolyak'', ''cluster'', or ''zcluster''.')
            end
        end %ExpandBasis
        
        
        
        
        
        
        
        
        %% disp
        function disp(B)  
            %%%
            % Prints a summary of the basis
            
            switch class(B)
                case 'basis'
                    fprintf('\nBasis of %1.0f dimension(s)\n\n',B.d)
                case 'funcApprox'
                    fprintf('\nFunction approximation for f:A --> B, where\n\t A < R^%d and B < R^%d\n\n',B.d,B.df)
            end
            
            fprintf('\t%-22s %-12s %-16s\n','    Variable','# of nodes','    Interval')
            for i=1:B.d
                fprintf('\t%-22s     %-8.0f [%6.2f, %6.2f]\n',B.opts.varnames{i},B.n(i),...
                    B.a(i),B.b(i))
            end
            
            fprintf('\n\t%-22s %-12s\n','Type of basis:',B.type)
            fprintf('\t%-22s %-12s\n','Type of nodes:',B.opts.nodetype)
            
            if B.d>1
                fprintf('\t%-22s %-12s\n','Expansion method:',B.opts.method)
                fprintf('\t%-22s %-12.0f\n','Total basis nodes:',size(B.opts.validX,1))
                fprintf('\t%-22s %-12.0f\n','Total basis functions:',size(B.opts.validPhi,1))
            end
            fprintf('\n\n')
        end
        
        %% plot
        function plot(B,ii)
            %%%
            %   B.plot(ii)
            %
            % Plots the basis functions. Input ii is optional, to plot only the specified basis functions.
            % Returns a figure with d-subplots (one per dimension).
            
            for i=1:B.d
                if nargin<2
                    ii = 1:B.n(i);
                else
                    ii(ii>B.n(i)) = [];
                end
                subplot(1,B.d,i)
                B.B1(i).plot(ii);
            end
        end
        
        
       
        %% Interpolation
        function Phi = Interpolation(B,...
                x,...      evaluation points
                order,...  order of derivatives
                integrate... integrate if false, derivative if true
                )
            %%%
            %   Phi = B.Interpolation(x,order,integrate)
            %
            % Computes interpolation matrices for basis of d-dimensions. It calls the method |Interpolation| on each
            % unidimensional basis, then combines all of them.
            %
            % Its inputs are 
            %
            % * |x|, k.d matrix or 1.d cell of evaluation points. Defaults to basis nodes. If cell is provided,
            % |Interpolation| is call of values from cell, and then combined (faster than running the method on the
            % combination of nodes directly).
            % * |order|, h.d matrix, for order of derivatives. Defaults to zeros(1,d).
            % * |integrate|, boolean, integrate if true, derivative if false. Defaults to false
            %
            % Output |Phi| returns the interpolation k.g.h array, where g is number of combined basis functions.
            %
            % See also <basis.Evaluate>, <basis.Jacobian> and <basis.Hessian>
            
            
            %%%
            % Check the inputs
            
            if nargin<2; x = [];             end  % evaluate at nodes
            if nargin<3 || isempty(order); order = zeros(1,B.d); end  % no derivative
            if isscalar(order); order = order*ones(1,B.d); end  % use same in all dimension
            if any(order<0), error('order must be nonnegative'),end
            if nargin<4, integrate = false; end
            
            
            %%%
            % Solve for the one-dimensional basis
            
            if B.d==1
                Phi  = B.B1.Interpolation(x,order,integrate);
                return
            end
            
            
            %%%
            % Solve for the multi-dimensional basis
            
            Norders = size(order,1);
            
            if size(order,2)~=B.d
                error('In Interpolation, class ''basis'': order must have d columns')
            end
            
            % Check validity of input x
            if nargin>1  % hasArg(x)
                if (isnumeric(x) && size(x,2)~=B.d && ~isempty(x)),error('In Interpolation, class basis: x must have d columns'); end
                if (iscell(x) && length(x)~=B.d), error('In Interpolation, class basis: x must have d elements'); end
            end
            
            
            %%%
            % HANDLE NODES DIMENSION
            if isempty(x)
                switch B.opts.method
                    case {'cluster' 'zcluster'}
                        x = B.nodes;
                        nrows = size(x,1);
                    otherwise
                        nrows = size(B.opts.validX,1);
                        r = B.opts.validX;    % combine default basis nodes
                end
                
            elseif iscell(x)
                nrows = prod(cellfun(@(v) numel(v),x));
                r = B.opts.validX;    % combine vectors in x
            else % ismatrix(x)
                nrows = size(x,1);
                r = repmat((1:size(x,1))',1,B.d); %use columns of matrix x
            end
            
            
            %%% 
            % HANDLE POLYNOMIALS DIMENSION
            c = B.opts.validPhi;
            ncols = size(c,1);
            
            %%%
            % Preallocate memory
            PHI = zeros(nrows,ncols,Norders,B.d);

            %%%
            % Compute interpolation matrices for each dimension
            for j = 1:B.d
                if isempty(x)
                    Phij = B.B1(j).Interpolation([],order(:,j),integrate);
                    PHI(:,:,:,j) = Phij(r(:,j),c(:,j),:);
                elseif iscell(x)
                    Phij = B.B1(j).Interpolation(x{j},order(:,j),integrate);
                    PHI(:,:,:,j) = Phij(r(:,j),c(:,j),:);
                else
                    Phij = B.B1(j).Interpolation(x(:,j),order(:,j),integrate);
                    PHI(:,:,:,j) = Phij(:,c(:,j),:);
                end
            end
            clear Phij
            
            %%%
            % MULTIPLY individual bases
            Phi = prod(PHI,4);
 
        end % Interpolation
        
        
        
        %% Nearest node
        function [xnode,idx] = nearestNode(B,x)
            %%%
            %   [xnode,idx] = nearestNode(B,x)
            %
            % For a given value x, it finds the nearest node in basis B. This method is useful in simulations, when an
            % initial guess solution at |x| can be given by |xnode|. |idx| is the associated index, e.g. 
            %
            %   xnode = B.nodes(idx,:)'
            
            N = size(B.nodes,1);
            gap = zeros(N,B.d);
            xnode = zeros(B.d,1);
            
            for i = 1:B.d
                gap(:,i) = abs(B.nodes(:,i) - x(i)) / (B.b(i) - B.a(i));
            end
            
            [~, idx] = min(sum(gap,2));
            for i=1:B.d
                xnode(i) = B.nodes(idx,i);
            end
            
        end
        
    end %methods
end %classdef
