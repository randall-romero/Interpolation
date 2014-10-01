%% basis class
% Defines a class to represent a multi-dimensional collocation basis 
%
% Objects created by this class are of type 'handle', useful for passing the basis by reference.
%
% This class is defined to combine _d_ unidimensional bases into a multidimensional one. The unidimensional basis is
% defined as any of these bases:
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
% * |varname|: 	1.d cell of strings, variable names
% * |nodes|: matrix with d columns,basis nodes
% * |opts|: structure with options for chosen type
% * |expansion|: structure with options for combaning the 1-bases
% * |type|:  basis type, i.e. 'cheb','spli', or 'lin'
%
% Object of class |basis| have the following methods:
%
% * |basis|: class constructor
% * |Interpolation|: returns the interpolation matrix
% * |Evaluate|: evaluates the interpolated function
% * |Jacobian|: computes the Jacobian at given values of variables
% * |Hessian|: computes the Hessian at given values of variables
% * |nearestNode|: returns the basis node that is closest to input argument
% * |ExpandBasis|: computes auxilliary fields needed for combining the d 1-basis
% * |plot|: plots the basis functions for each dimension
% * |display|: prints a summary of the basis
% 
%
% *To date, only basisChebyshev has been implemented*, so calling the function with |type| 'spli' or 'lin' returns an
% error.
%
% Implemented by:
%
% * *Randall Romero-Aguilar*
% * |randall.romero@outlook.com|
% * Last updated: September 23, 2014






%% IMPLEMENTATION
classdef basis < handle
    
    properties (SetAccess = protected)
        d           % dimension of basis
        n           % number of nodes
        a           % lower bounds
        b           % upper bounds
        B1          % 1.d array of one-dimensional bases
        varname     % variable names
        nodes       % nodes for d-dimensions
        opts        % options for basis of given type
        expansion   % options to expand basis
        type        % 'cheb','spli', or 'lin'
    end
    
    %% ===== Methods =====

    methods
        
        %% basis
        function B = basis(...
                type,...        string, type of basis
                n,...           1.d, number of nodes
                a,...           1.d, lower bounds
                b,...           1.d, upper bounds
                varname,...     1.d cell of strings, variable names
                opts,...        structure with option for type of basis
                expansion)%     structure with options parameters
            %%% 
            % Constructor for basis
            %
            %   B = basis(type,n,a,b,varname,opts,expansion)
            % 
            % The following inputs are required:
            %
            % * |type|: 'cheb','spli', or 'lin'
            % * |n|: 1.d vector of integers, number of nodes
            % * |a|: 1.d vector, lower bounds
            % * |b|: 1.d vector, upper bounds
            %
            % Optional inputs
            % 
            % * |varname|: 1.d cell of strings, variable names
            % * |opts|: structure, options for specific basis |type|
            % * |expansion|: structure, options for expanding the basis
            %
            % Options passed on structure |opts| depend on basis |type|:
            %
            % * |'cheb'|:  |opts.nodetype| is either 'lobatto','gaussian' (default) or 'endpoint'
            % * |'spli'|:  |opts.order| is order of spline, (default = 3)
            %
            % Options passed on structure |expansion|
            %
            % * |expansion.method|: is either 'tensor' (default), 'smolyak', 'complete', 'cluster', or 'zcluster'
            % * |expansion.degreeParam|: adjust polynomial degrees for all methods (except 'tensor')
            % * |expansion.nodeParam|: adjust selection of nodes, for methods 'smolyak', 'cluster', and 'zcluster'
            %
            % As of September 23, 2014 the methods 'cluster' and 'zcluster' are only experimental.
            
            %%%
            % Basis dimension
            
            B.d = numel(n);
            
            %%%
            % Default values
            
            if nargin < 5 || isempty(varname)
                varname = strcat('V',num2cell(48 + (1:B.d)));
            end
            
            
            switch lower(type)
                case 'cheb'
                    if nargin < 6 || isempty(opts)
                        opts.nodetype = 'gaussian';
                    end
                    
                    if nargin <7
                        expansion.method = 'tensor';
                    end
                    
                    % Adjust number of Smolyak nodes, if selected
                    if strcmp(expansion.method,'smolyak')
                        n_old = n;
                        n = 2.^ceil(log2(n_old-1))+1;
                        if any(n-n_old)
                            warning('For Smolyak expansion, number of nodes adjusted as follows')
                            fprintf('\t%6s, %6s\n','Old n','New n')
                            fprintf('\t%6d, %6d\n',[n_old;n])
                        end
                        
                        if isempty(opts.nodetype) || ~strcmp(opts.nodetype,'lobatto')
                            warning('nodetype must be ''lobatto'' for Smolyak nodes')
                            opts.nodetype = 'lobatto';
                        end
                    end
                    
                otherwise
                    error('method not yet implemented')
            end
            
            %%%
            % Create the 1-basis
                        
            switch lower(type)
                case 'cheb'
                    for i = 1: B.d
                        B1(i) = basisChebyshev(n(i),a(i),b(i),opts.nodetype,varname{i});
                    end
                otherwise
                    error('Method not yet implemented')
            end
                        
            %%%
            % Pack values in object
            B.type = lower(type);
            B.B1 = B1;
            B.a = a;
            B.b = b;
            B.n = n;
            B.varname = varname;
            B.opts.nodetype = B1(1).nodetype;
            B.expansion = expansion;
            B.ExpandBasis;
        end %basis
        
        
        
        
        
        
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
            % * |B.expansion.validPhi|: indices to combine unidimensional bases
            % * |B.expansion.validX|: indices to combine unidimensional nodes
            %
            % Combining polynomials depends on value of input |expansion.method|:
            %
            % * |'tensor'| takes all possible combinations,
            % * |'smolyak'| computes Smolyak basis, given |expansion.degreeParam|,
            % * |'complete'|, |'cluster'|, and |'zcluster'| choose polynomials with degrees not exceeding
            % |expansion.degreeParam|
            %
            % Expanding nodes depends on value of field |expansion.method|. 
            %
            % * 'tensor' and 'complete' take all possible combinations, 
            % * 'smolyak' computes Smolyak basis, given |expansion.nodeParam|
            % * 'cluster' and 'zcluster' compute clusters of the tensor nodes based on |B.expansion.nodeParam|
            
            
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
            
            switch lower(B.expansion.method)
                case 'tensor'
                    B.expansion.validPhi = deg_grid + 1;
                    
                case {'complete','cluster','zcluster'}
                    deg = sum(deg_grid,2);  % degrees of all possible column products
                    degValid = (deg <= B.expansion.degreeParam);
                    idx = deg_grid + 1;  % index of entries
                    B.expansion.validPhi = idx(degValid,:);
                    
                case 'smolyak'
                    if B.expansion.nodeParam < B.expansion.degreeParam
                        warning('Smolyak degree param cannot be bigger than node param; adjusting degree');
                        B.expansion.degreeParam = B.expansion.nodeParam;
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
                    
                    validNodes = (sum(xGroupAll,2)<= B.d + B.expansion.nodeParam);
                    validPhi = (sum(phiGroupAll,2)<= B.d + B.expansion.degreeParam);
                    
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
                    
                    
                    B.expansion.validX = idxAll(validNodes,:);  % index of entries
                    B.expansion.validPhi = idxAll(validPhi,:);  % index of entries
                    
                    
                otherwise
                    error(...
                        'Expansion type must be ''tensor'', ''complete'', ''smolyak'', ''cluster'', or ''zcluster''.')
            end
            
            
            %%%
            % Expanding nodes
            
            switch lower(B.expansion.method)
                case {'tensor', 'complete'}
                    B.nodes = gridmake({B.B1.nodes});
                    B.expansion.validX = idxAll;  % index of entries
                    
                case 'smolyak'
                    % done in previous switch
                    
                case {'cluster','zcluster'}
                    H = size(B(1).expansion.validPhi,1) + B.expansion.nodeParam;  % number of clusters
                    
                    
                    tempNodes = gridmake({B.B1.nodes});
                    
                    if strcmp(B.expansion.method,'cluster')
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
        
        
        
        
        
        
        
        
        %% display
        function display(B)  
            %%%
            % Prints a summary of the basis
                        
            fprintf('\nBasis of %1.0f dimension(s)\n\n',B.d)
            
            fprintf('\t%-22s %-12s %-16s\n','    Variable','# of nodes','    Interval')
            for i=1:B.d
                fprintf('\t%-22s     %-8.0f [%6.2f, %6.2f]\n',B.varname{i},B.n(i),...
                    B.a(i),B.b(i))
            end
            
            fprintf('\n\t%-22s %-12s\n','Type of basis:',B.type)
            fprintf('\t%-22s %-12s\n','Type of nodes:',B.opts.nodetype)
            
            if B.d>1
                fprintf('\t%-22s %-12s\n','Expansion method:',B.expansion.method)
                fprintf('\t%-22s %-12.0f\n','Total basis nodes:',size(B.expansion.validX,1))
                fprintf('\t%-22s %-12.0f\n','Total basis functions:',size(B.expansion.validPhi,1))
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
                switch B.expansion.method
                    case {'cluster' 'zcluster'}
                        x = B.nodes;
                        nrows = size(x,1);
                    otherwise
                        nrows = size(B.expansion.validX,1);
                        r = B.expansion.validX;    % combine default basis nodes
                end
                
            elseif iscell(x)
                nrows = prod(cellfun(@(v) numel(v),x));
                r = B.expansion.validX;    % combine vectors in x
            else % ismatrix(x)
                nrows = size(x,1);
                r = repmat((1:size(x,1))',1,B.d); %use columns of matrix x
            end
            
            
            %%% 
            % HANDLE POLYNOMIALS DIMENSION
            c = B.expansion.validPhi;
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
        
        
        %% Evaluate
        function y = Evaluate(B,...
                coef,...    interpolation coefficiens
                varargin... inputs for interpolation method
                )
            %%%
            %   y = B.Evaluate(coef,x,order,integrate)
            %
            % Evaluates the interpolated functions defined by basis |B| and coefficients |coef|. Inputs |x|, |order| and
            % |integrate| are as required by <Interpolation> method, while |coef| is a g.m matrix (g basis functions and
            % m interpolated functions)
            %
            % Output |y| returns the interpolated functions as a k.m.h array, where k is the number of evaluation
            % points, m the number of functions (number of columns in |coef|), and h is number of order derivatives.
            
            
            if nargin <2, error('Missing ''coef'' input'), end
            Phi = B.Interpolation(varargin{:});
            
            nx = size(Phi,1);  % number of evaluation points
            ny = size(coef,2); % number of evaluated functions
            no = size(Phi,3);  % number of order evaluations
            
            y = zeros(nx,ny,no);
            
            for h = 1:no
                y(:,:,h) = Phi(:,:,h) * coef;
            end
            
            if ny==1  % only one function
                y = squeeze(y);
            end
            
        end %Evaluate
        
        
        
        %% Jacobian
        function [DY,Y] = Jacobian(B,...
                x,...      evaluation points
                coef,...   interpolation coefficients
                index...   optional 1.d boolean indicating if derivative wrt given dimension is required
                )
            %%%
            %   [DY,Y] = B.Jacobian(x,coef,index)
            %
            % Computes the Jacobian of the approximated function f, for f: R^d --> R^m.
            %
            % Inputs:
            %
            % * |x|, k.d matrix of evaluation points.
            % * |coef|, g.m matrix, interpolation coefficients, where g is the number of basis functions and m is the
            % number of interpolated functions.
            % * |index|, 1.d boolean, take partial derivative only wrt variables with index=true. Defaults to true(1,d).
            %
            % Outputs
            %
            % * |DY|, k.m.d1 Jacobian matrix evaluated at |x|, where d1 = sum(index)
            % * |Y|, k.m interpolated function (same as Y = B.Evaluate(coef,x), provided for speed if both Jacobian and
            % funcion value are required).
            %
            % Unlike |Evaluate|, method |Jacobian| requires the |x| argument as an explicit matrix.
            
            
            %%%
            % Check the inputs
            Nfunc = size(coef,2);   %number of interpolated functions
            
            %%%
            % Solve for the one-dimensional basis
            
            if B.d==1
                if nargout == 1
                    Phi  = B.Interpolation(x,1,false);
                    DY = Phi * coef;
                else
                    Phi = B.Interpolation(x,[1;0],false);
                    DY = Phi(:,:,1) * coef;
                    Y = Phi(:,:,2) * coef;
                end
                
                return
            end
            
            %%%
            % Solve for the multi-dimensional basis
            
            % Check validity of input x
            if size(x,2)~=B.d, error('In Jacobian, class basis: x must have d columns'); end
            
            %%%
            % Keep track of required derivatives: Required is logical with true indicating derivative is required
            if nargin<4
                Required = true(1,B.d);
                index = 1:B.d;
            elseif numel(index) < B.d % assume that index have scalars of the desired derivatives
                Required = false(1,B.d);
                Required(index) = true;
            else % assume that index is nonzero for desired derivatives
                Required = logical(index);
                index = find(index);
            end
            
            nRequired = sum(Required);
            
            %%%
            % HANDLE NODES DIMENSION
            Nrows = size(x,1);
            
            %%%
            % HANDLE POLYNOMIALS DIMENSION
            c = B.expansion.validPhi;
            Ncols = size(c,1);
            
            %%%
            % Compute interpolation matrices for each dimension
            
            Phi0 = zeros(Nrows,Ncols,B.d);
            Phi1 = zeros(Nrows,Ncols,B.d);
            
            for k = 1:B.d
                if Required(k)
                    PhiOneDim = B.B1(k).Interpolation(x(:,k),...
                        [0 1],...
                        false);
                    Phi01 = PhiOneDim(:,c(:,k),:);                   
                    Phi0(:,:,k) = Phi01(:,:,1); % function value
                    Phi1(:,:,k) = Phi01(:,:,2); % its derivative
                else
                    PhiOneDim = B.B1(k).Interpolation(x(:,k),...
                        0,...
                        false);
                    Phi0(:,:,k) = PhiOneDim(:,c(:,k));
                end
            end
            
            %%%
            % Compute the Jacobian
            % Preallocate memory
            DY  = zeros(Nrows,Nfunc,nRequired);
            
            % Multiply the 1-dimensional bases
            
            for k=1:nRequired
                Phik = Phi0;
                Phik(:,:,index(k)) = Phi1(:,:,index(k));
                Phi = prod(Phik,3);
                DY(:,:,k) = Phi*coef;
            end
            
            %%%
            % Compute the function if requested
            if nargout > 1
                Y = prod(Phi0,3)*coef;
            end
            
        end % Jacobian
        
        
        
        %% Hessian
        function Hy = Hessian(B,...
                x,...     evaluation points
                coef...   interpolation coefficiens
                )
            %%%
            %   Hy = B.Hessian(x,coef)
            %
            % Computes the Hessian of a function approximated by basis |B| and coefficients |coef|.
            %
            % Its inputs are:
            %
            % * |x|, k.d matrix of evaluation points.
            % * |coef|, g.m matrix, interpolation coefficients.
            %
            % Its output |Hy| returns the k.m.d.d Hessian evaluated at |x|.
            
            
            order = repmat({[0 1 2]'},1,B.d);
            order = gridmake(order{:});
            order = order(sum(order,2)==2,:);
            
            Phi = B.Interpolation(x,order,false);
            
                      
            nx = size(x,1);     % number of evaluated points
            ny = size(coef,2);  % number of evaluated functions
            
            
            %Dy = squeeze(Dy);
            
            Hy = zeros(nx,ny,B.d,B.d);
            
            for k = 1:size(order,1)
                i = find(order(k,:));
                if numel(i)==1
                    Hy(:,:,i,i) = Phi(:,:,k) * coef;%  Dy(:,k);
                else
                    Hy(:,:,i(1),i(2)) = Phi(:,:,k) * coef; %Dy(:,k);
                    Hy(:,:,i(2),i(1)) = Hy(:,:,i(1),i(2)); %Dy(:,k);
                end
            end
        end % Hessian
        
        
        
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

