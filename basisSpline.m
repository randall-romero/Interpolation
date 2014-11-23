%% basisChebyshev class
% Defines a class to represent a univariate Chebyshev basis.  
%
% Objects of class |basisChebyshev| have the following properties:
%
% * |n|:    scalar, number of nodes
% * |a|:    scalar, lower bound
% * |b|:    scalar, upper bound
% * |nodes|: d.1 vector, basis nodes
% * |D|: operators to compute derivative
% * |I|: operators to compute integrals
% * |nodetype|:  node type, i.e. 'gaussian','lobatto', or 'endpoint'
% * |varname|: 	string, variable name. Defaults to 'V0'.
% * |WarnOutOfBounds|: boolean, warns if interpolating outside [a,b] if true. Defaults to 'false'.
%
% Object of class |basisChebyshev| have the following methods:
%
% * |basisChebyshev|: class constructor
% * |SetNodes|: computes the nodes
% * |Interpolation|: returns the interpolation matrix
% * |DiffOperator|: computes operators for derivatives and integrals
% * |Interpolate|: evaluates the interpolated function
% * |nearestNode|: returns the basis node that is closest to input argument
% * |plot|: plots the basis functions
% * |display|: prints a summary of the basis
% * |copy|: makes an independent copy (method inherited from superclass)
% * |CheckValidity|: checks consistency of properties
% * |Reset|: computes nodes and deletes |D| and |I| operators
% * Setter methods for |a|, |b|, |n|, |nodetype|: change respective value and resets the basis.
%
% *NOTE*: Objects created by this class are derived from the superclass
% 'matlab.mixin.Copyable', which is equivalent to a 'handle' but adds a 'copy' method for
% making shallow copies of the basis. For example, in copying basisChebyshev B with 
% 
%   B1 = B
%   B2 = B.copy
% 
% the B1 variable is just a reference to B (then, subsequent changes to B1 also affect B),
% whereas B2 is a value copy of B (thus, changes to B2 will not affect B).
%
%
%
% Last updated: October 6, 2014.
%
% Copyright (C) 2013-2014 Randall Romero-Aguilar
%
% Licensed under the MIT license, see LICENSE.txt

%% 
classdef basisSpline < matlab.mixin.Copyable
    properties (SetAccess = protected)
        nodes       % nodes for own dimension
        D  = cell(0,1)  % operators to compute derivative
        I  = cell(0,1)  % operators to compute integrals
        nodetype = 'cardinal'   % type of nodes        
        k  = 3       % order of spline
    end
    
    properties (SetAccess = public)
        breaks
        varname = 'V0'             % optional name for variable
        WarnOutOfBounds = true    % warns if interpolating outsite [a,b] interval
    end
    
    properties (Dependent)
        n      % number of nodes
        a      % lower limits
        b      % upper limits
    end
    
    
    methods
        
        
        %% basisSpline
        function B = basisSpline(varargin)
                        
            %%% 
            % Constructor for basisChebyshev
            %
            %   B = basisChebyshev(n,a,b,nodetype,varname)
            % 
            % The following inputs are required:
            %
            % * |n|: scalar, number of nodes
            % * |a|: scalar, lower bound
            % * |b|: scalar, upper bound
            %
            % Optional inputs
            % 
            % * |nodetype|: node type, one of 'gaussian' (default), 'endpoint', or 'lobatto'
            % * |varname|: string, variable name. Defaults to 'V0'
            
            %%%
            % If no inputs are provided (needed to preallocate array of basis)
            if nargin==0
                B.breaks = [-1,1];
                return
            end
            
            if ~isscalar(varargin{1})
                if nargin>1 && ~isempty(varargin{2})
                    B.k = varargin{2};
                end
                
                if nargin > 2
                    B.varname = varargin{3};
                end
                
                B.breaks = varargin{1};
            else
                nn = varargin{1};
                aa = varargin{2};
                bb = varargin{3};
                
                
                assert(isscalar(aa),'Lower bound must be a scalar');
                assert(isscalar(bb),'Upper bound must be a scalar');
                assert(aa < bb,'Lower bound must be less than upper bound');
                
                if nargin>3 && ~isempty(varargin{4})
                    B.k = varargin{4};
                end
                
                if nargin > 4
                    B.varname = varargin{5};
                end
                
                
                
                assert(nn + 1 > B.k + 1,'n must be at least order + 1')
                B.breaks = linspace(aa,bb,nn + 1 - B.k);
            end
        end %basisSpline
               
        
        %% Setter and getter methods
        function set.breaks(B,value)
            assert(numel(value)>1,...
                'breakpoint sequence must contain at least two elements');    
            B.breaks = sort(value(:));
            B.SetNodes;
        end
        
        function aa = get.a(B)
            aa = B.breaks(1);
        end
        
        function bb = get.b(B)
            bb = B.breaks(end);
        end
        
        function nn = get.n(B)
            nn = numel(B.breaks) + B.k - 1;
        end
        
        
        
        
        %% SetNodes
        function SetNodes(B)
            %%%
            % B.SetNodes
           
            x=cumsum([B.a*ones(B.k,1);B.breaks;B.b*ones(B.k,1)]);
            x=(x(1+B.k:B.n+B.k) - x(1:B.n))/B.k;
            x(1)   = B.a;
            x(end) = B.b;
            B.nodes = x;
        end %SetNodes
        
        
        
        %% DiffOperator
        function DiffOperator(B,...
                orderD,...            order of derivative
                integrate...          integral if true
                )
            %%%
            % B.DiffOperator(orderD,integrate)
            %
            % Compute operators to differentiate and integrate. This method is needed by |Interpolation| when input
            % |order| is not zero; it should not be called directly by the user.
            %
            % It updates the operators in |B.D| (if integrate is false) or in |B.I| (otherwise). Input |orderD| is the
            % maximum order to compute.
            if nargin < 2, orderD = 1; end
            if nargin < 3, integrate = false; else integrate = logical(integrate); end
            assert(orderD >= 0, 'In DiffOperator: order must be nonnegative')
            
            derivative = ~integrate;            
            if orderD==0, return, end
            
            % Use previously stored values for orders
            if (derivative && length(B.D) >= orderD), return, end
            if (integrate  && length(B.I) >= orderD), return, end
            
            n = B.n;
            k = B.k;
            
            if derivative
                %===============TAKE DERIVATIVE=========
                assert(orderD <= k,'Order of derivative operator cannot be larger than k')
                kk = max(k-1, k-orderD-1);
                augbreaks = [repmat(B.a,kk,1); B.breaks; repmat(B.b,kk,1)];
                temp = k./(augbreaks(k+1:n+k-1) - augbreaks(1:n-1));
                if isempty(B.D)
                    B.D{1} = spdiags([-temp temp],0:1,n-1,n);
                end
                h = length(B.D) + 1;  % index of 1st missing element in B.D
                for i = h:orderD
                    temp=(k+1-i)./(augbreaks(k+1:n+k-i)-augbreaks(i:n-1));
                    B.D{i}=spdiags([-temp temp],0:1,n-i,n+1-i)*B.D{i-1};
                end
            else
                %===============TAKE INTEGRAL=========
                kk=max(k-1,k + orderD-1);
                augbreaks=[repmat(B.a,kk,1); B.breaks; repmat(B.b,kk,1)];
                temp = (augbreaks(kk+2:kk+n+1) - augbreaks(kk-k+1:kk+n-k)) / (k+1);
                if isempty(B.I)
                    B.I{1} = sparse(tril(ones(n+1,1)*temp',-1));
                end
                h = length(B.I)+1;    % index of 1st missing element in B.I
                for i= h:orderD
                    temp=(augbreaks(kk+2:kk+n+i)-augbreaks(kk-k-i+2:kk+n-k))/(k+i);
                    B.I{i}=sparse(tril(ones(n+i,1)*temp',-1))*B.I{i-1};
                end
            end
        end
        
        
        %% Interpolation
        function Phi = Interpolation(B,... ##CompEcon:  chebbas.m
                x,...                      k-vector of the evaluation points
                order,...                  order of derivative
                integrate)%                integrate if true
            %%%
            %   Phi = B.Interpolation(x,order,integrate)
            %
            % Computes interpolation matrices |Phi| for basis |B|.
            %
            % Its inputs are
            %
            % * |x|, k.1 vector of evaluation points. Defaults to basis nodes.
            % * |order|, h.1 vector, for order of derivatives. Defaults to 0.
            % * |integrate|, boolean, integrate if true, derivative if false. Defaults to false
            %
            % Output |Phi| returns the interpolation k.n.h array, where n is number of basis functions.
            %
            % See also <basis.Evaluate>.
            

            if nargin<2 || isempty(x)
                x = B.nodes;  % use basis nodes if x is missing
            else
                x = x(:);   % make sure x is column vector;    
            end
            
            if nargin<3, order = 0; end
            if nargin<4, integrate = false; end
            
          
            if integrate
                sgorder = - order;
            else
                sgorder = order;
            end
            
            assert(max(sgorder) < B.k, 'Derivatives defined for order less than k')
            
            nx = numel(x);
            k = B.k;
            minorder = min(sgorder);
            kaug = k - minorder;
            
            augbreaks = [repmat(B.a,kaug,1); B.breaks; repmat(B.b,kaug,1)];
            ind = lookup(augbreaks,x,3);
            
            % Recursively determine the values of a k-order basis matrix.
            % This is placed in an (m x k+1-order) matrix
            bas = zeros(nx,kaug + 1);
            bas(:,1) = ones(nx,1);
            Phi = ndsparse(numel(order));
            
            for j=1:kaug
                for jj=j:-1:1
                    b0 = augbreaks(ind+jj-j);
                    b1 = augbreaks(ind+jj);
                    temp = bas(:,jj)./(b1-b0);
                    bas(:,jj+1) = (x-b0) .* temp+bas(:,jj+1);
                    bas(:,jj) = (b1-x).*temp;
                end
                % as now contains the order j spline basis
                ii = find((k-j)==sgorder);
                if ~isempty(ii)
                    ii=ii(1);
                    % Put values in appropriate columns of a sparse matrix
                    r = repmat((1:nx)',1,k-sgorder(ii)+1);
                    c = (sgorder(ii)-k : 0) - (sgorder(ii)-minorder);
                    c = repmat(c,nx,1) + repmat(ind,1,k-sgorder(ii)+1);
                    Phi{ii} = sparse(r,c,bas(:,1:k-sgorder(ii)+1),nx,B.n-sgorder(ii));
                    
                    if order(ii)~=0
                        % If needed compute derivative or anti-derivative operator
                        if integrate
                            try
                                Phi{ii} = Phi{ii} * B.I{order(ii)};
                            catch
                                B.DiffOperator(max(order),true);
                                Phi{ii} =  Phi{ii} * B.I{order(ii)};
                            end    
                        else
                            try
                                Phi{ii} = Phi{ii} * B.D{order(ii)};
                            catch
                                B.DiffOperator(max(order));
                                Phi{ii} = Phi{ii} * B.D{order(ii)};
                            end
                        end
                    end
                end
                
            end
        end %Interpolation
        
        
        
        %% Evaluate
        function y = Interpolate(B,...      ##CompEcon:  chebbas.m
                coef,...             coefficient matrix
                varargin...          inputs for Interpolation method
                )
            %%%
            %   y = B.Evaluate(coef,x,order,integrate)
            %
            % Evaluates the interpolated functions defined by basis |B| and coefficients |coef|. Inputs |x|, |order| and
            % |integrate| are as required by <Interpolation> method, while |coef| is a n.m matrix (n basis functions and
            % m interpolated functions)
            %
            % Output |y| returns the interpolated functions as a k.m.h array, where k is the number of evaluation
            % points, m the number of functions (number of columns in |coef|), and h is number of order derivatives.
                        
            
            assert(nargin > 1, 'Missing ''coef'' input')
                        
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
        
        
        
         %% plot
         function plot(B,ii)
             %%%
             %   B.plot(ii)
             %
             % Plots the basis functions. Input ii is optional, to plot only the specified basis functions.
             

            if nargin<2
                ii = 1:B.n;
            else
                ii(ii>B.n) = [];
            end
            x = linspace(B.a,B.b,1+min(120,10*B.n))';
            Phi = B.Interpolation(x);
            Phi = Phi{1};
            hold off
            plot(x,0*x,'k'); hold on
            plot(x,Phi(:,ii),'LineWidth',1.75)
            plot(B.nodes,0*B.nodes,'ro')
            title(B.varname)
            hold off
            
        end
        
        
        
        %% Nearest node
        function [xnode,idx] = nearestNode(B,x)
            %%%
            %   [xnode,idx] = nearestNode(B,x)
            %
            % For a given value x, it finds the nearest node in basis B. This method is useful in simulations, when an
            % initial guess solution at |x| can be given by |xnode|. |idx| is the associated index, e.g. 
            %
            %   xnode = B.nodes(idx,:)'
            
            gap = abs(B.nodes - x);
            [~, idx] = min(gap);
            xnode = B.nodes(idx);
            
        end        
        
        %% display
        function display(B)  
            %%%
            % Prints a summary of the basis
            for i=1:numel(B)
                fprintf('\nSpline Basis of one dimension with order = %d\n\n',B(i).k)
                
                fprintf('\t%-22s %-12s %-16s\n','    Variable','# of nodes','    Interval')
                fprintf('\t%-22s     %-8.0f [%6.2f, %6.2f]\n',B(i).varname,B(i).n,B(i).a,B(i).b)
                fprintf('\n\n')
            end
        end        
        
        
        %% Reset
        function Reset(B)
            %%%
            % Following a change in the basis definition, the nodes must be recomputed and the fields |D| and |I| are no
            % longer valid. This method is called by the setter methods, so the user does not need to call it
            % explicitly.
            %B.CheckValidity;
            B.SetNodes;
            B.I = cell(0,1);
            B.D = cell(0,1);
        end
        
        %% Setter methods
        %
        % The following setter methods allow changing the basis definition (parameters n,a,b,nodetype). Changes to these
%         % parameters require recomputing the nodes and deleting the operators D and I.
%         
%         %%%
%         % Adjusting the lower bound
%         function set.a(B,value)
%             B.a = value;
%             B.Reset;
%         end
%         
%         %%%
%         % Adjusting the upper bound
%         function set.b(B,value)
%             B.b = value;
%             B.Reset;
%         end
%         
%         %%%
%         % Adjusting the number of nodes
%         function set.n(B,value)
%             B.n = value;
%             B.Reset;
%         end
%         
%         %%%
%         % Adjusting the type of nodes
%         function set.nodetype(B,value)
%             B.nodetype = value;
%             B.Reset;
%         end
%         
        
    end %methods
    
    
end %classdef