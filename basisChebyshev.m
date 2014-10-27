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
classdef basisChebyshev < matlab.mixin.Copyable
    properties (SetAccess = protected)
        nodes       % nodes for own dimension
        D           % operators to compute derivative
        I           % operators to compute integrals
    end
    
    properties (SetAccess = public)
        n  =  3     % number of nodes
        a  = -inf     % lower limits
        b  =  inf     % upper limits
        nodetype = 'gaussian'   % type of nodes        
        varname = 'V0'             % optional name for variable
        WarnOutOfBounds = true    % warns if interpolating outsite [a,b] interval
    end
    
    
    
    methods
        
        
        %% basisChebyshev
        function B = basisChebyshev(...
                n,...           scalar, number of nodes
                a,...           scalar, lower bound
                b,...           scalar, upper bound
                nodetype,...    string, 'gaussian','endpoint', or 'lobatto'
                varname)%       string, name for variable
                        
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
            if nargin ==0
                B.n = 3;
                B.a = -1;
                B.b = 1;
                return
            end
            
            
            %%%
            % SET PROPERTIES
            
            B.n = n;
            B.a = a;
            B.b = b;
            
            %%%
            % Set nodetype
            if nargin > 3 && ~isempty(nodetype)
                B.nodetype = lower(nodetype);
            end
            
            %%%
            % Set variable name
            
            if nargin==5 
                B.varname = varname;
            end
            
            %%%
            % Check validity of inputs and set nodes
            B.CheckValidity;
            B.SetNodes;
            
        end %basisChebyshev1
               
        
        
        %% checkValidity
        function CheckValidity(B)
            %%%
            % B.CheckValidity
            %
            % Checks the consistency of some basis B parameters. Called by constructor to ensure the following
            % conditions :
            %
            % * the basis is unidimensional,
            % * a < b, 
            % * n >= 2,
            % * nodetype is either 'gaussian', 'endpoint', or 'lobatto'
        
            assert(isscalar(B.a), 'Parameter ''a'' must be scalar')
            assert(isscalar(B.b), 'Parameter ''b'' must be scalar')
            assert(isscalar(B.n), 'Parameter ''n'' must be scalar')
            
            assert(B.a < B.b, 'Lower bound must be less than upper bound')
            assert(B.n >= 2, 'n must be greater than one')
            
            if ~ismember(lower(B.nodetype),{'gaussian' 'endpoint' 'lobatto'})
                warning(['Unknown nodetype ',B.nodetype,': using gaussian instead'])
                B.nodetype = 'gaussian';
            end
            
        end %CheckValidity
        
        
        
        
        
        %% SetNodes
        function SetNodes(B)
            %%%
            % B.SetNodes
            %
            % Computes the nodes of the basis |B|, depending on |B.nodetype|, which can be 'gaussian', 'endpoint',
            % or 'lobatto'. This method is called by the constructor.
                        
            
            switch B.nodetype
                case 'gaussian'   % Gaussian nodes
                    x = -cos(pi*(1:2*B.n-1)/(2*B.n));
                case 'endpoint'     % Extend nodes to endpoints
                    x = -cos(pi*(1:2*B.n-1)/(2*B.n));
                    x = x/x(end);
                case 'lobatto'     % Lobatto nodes
                    x = - cos(pi*(0:B.n-1)/(B.n-1));
            end
            
            x(abs(x)<eps) = 0;
            B.nodes = B.rescale2_ab(x');
        end %SetNodes
        
        

        function z = rescale2_01(B,x)
            z = (2/(B.b - B.a))* (x-(B.a + B.b)/2);
        end
        
        
        function x = rescale2_ab(B,z)
            x = (B.a + B.b + (B.b-B.a)*z)/2;
        end
        
        
        
        
        
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
            a = B.a;
            b = B.b;
            
            if(n-orderD < 2 && derivative)
                warning('Insufficient nodes: truncating order to = %3d',n-2);
                orderD = n-2;
            end
            
            if derivative
                %===============TAKE DERIVATIVE=========
                if isempty(B.D)
                    [j,i] = meshgrid(1:n,1:n);
                    rc = mod(i+j,2) & (j>i);
                    d = zeros(n,n);
                    d(rc) = (4/(b-a))*(j(rc)-1);
                    d(1,:) = d(1,:)/2;
                    d = sparse(d(1:end-1,:));
                    B.D{1} = d;
                else
                    d = B.D{1};
                end
                
                h = length(B.D) + 1;  % index of 1st missing element in B.D
                for ii = h:orderD
                    B.D{ii} = d(1:n-ii,1:n-ii+1) * B.D{ii-1};
                end
            else
                %===============TAKE INTEGRAL=========
                nn = n + orderD;
                i = (0.25*(b-a))./(1:nn);
                d = sparse(diag(i) + diag(-i(1:nn-2),2));  %good
                d(1,1) = 2*d(1,1);
                
                d0 = ((-1).^(0:nn-1)).*sum(d);
                
                if isempty(B.I)
                    B.I{1} = [d0(1:n);d(1:n,1:n)];
                end
                h = length(B.I)+1;    % index of 1st missing element in B.I
                
                for ii = h:orderD
                    B.I{ii}=[d0(1:n+ii-1);d(1:n+ii-1,1:n+ii-1)] * B.I{ii-1};
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
                hasArg_x = false;
                nx = B.n;
            else
                x = x(:);   % make sure x is column vector;
                hasArg_x = true;
                nx = numel(x);
            end
            
            if nargin<3, order = 0; end
            if nargin<4, integrate = false; end
            
            maxorder = max(order);
            nn = B.n + maxorder * integrate;
            
            % check if out of bounds
            if hasArg_x && B.WarnOutOfBounds
                isBelow =  x < B.a - eps(B.a);
                isAbove =  x > B.b + eps(B.b);
                isOut = isBelow | isAbove;
                
                if any(isOut)
                    xout = x(isOut);
                    error(['OutOfBounds!: Interpolating %s = %f failed because basis defined only in [%4.2f, %4.2f]\n',...
                        'To allow interpolation outside this interval, set basis field ''WarnOutOfBounds'' to false'],...
                        B.varname,xout(1),B.a,B.b)
                end
                
                
            end
            
                 
            % Compute order 0 interpolation matrix
            if hasArg_x || ~strcmp(B.nodetype,'gaussian') % evaluate at arbitrary nodes
                m = length(x);
                z = B.rescale2_01(x);
                bas = zeros(m,nn);
                bas(:,1) = 1;
                bas(:,2) = z;
                z = 2*z;
                for i = 3:nn
                    bas(:,i) = z.*bas(:,i-1) - bas(:,i-2);
                end
            else  % evaluate at usual Gaussian nodes
                temp=((B.n-0.5):-1:0.5)';
                bas  = cos((pi/B.n) * temp *(0:nn-1));
            end
            
            
            % Compute Phi ////////////////////////
            %Phi = cell(length(order),1);
            Phi = zeros(nx,B.n,length(order));
            
            done = nan(length(order),1);
            
            for ii = 1:length(order)
                if any(order(ii)==done)  %done it already
                    %Phi{ii} = Phi{find(order(ii)==done,1)};
                    Phi(:,:,ii) = Phi(:,:,find(order(ii)==done,1));
                else
                    %%%
                    % The operators I and D are saved into the basis when created. The try/catch sequence avoids calling
                    % the DiffOperator method (to speed up execution)
                    
                    if order(ii)==0
                        %Phi{ii} = bas;
                        Phi(:,:,ii) = bas;
                    elseif integrate
                        %Phi{ii} = bas(:,1:B.n+order(ii)) * B.I{order(ii)};
                        try
                            Phi(:,:,ii) = bas(:,1:B.n+order(ii)) * B.I{order(ii)};
                        catch
                            DiffOperator(B,maxorder,integrate);
                            Phi(:,:,ii) = bas(:,1:B.n+order(ii)) * B.I{order(ii)};
                        end
                    else % derivative
                        %Phi{ii} = bas(:,1:B.n-order(ii)) * B.D{order(ii)};
                        try
                            Phi(:,:,ii) = bas(:,1:B.n-order(ii)) * B.D{order(ii)};
                        catch
                            DiffOperator(B,maxorder,integrate);
                            Phi(:,:,ii) = bas(:,1:B.n-order(ii)) * B.D{order(ii)};
                        end
                    end
                    Phi(abs(Phi)<eps(10)) = 0;
                    
                end
                done(ii) = order(ii);
                
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
                fprintf('\nChebyshev Basis of one dimension with %s nodes\n\n',B(i).nodetype)
                
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
            B.CheckValidity;
            B.SetNodes;
            B.I = [];
            B.D = [];
        end
        
        %% Setter methods
        %
        % The following setter methods allow changing the basis definition (parameters n,a,b,nodetype). Changes to these
        % parameters require recomputing the nodes and deleting the operators D and I.
        
        %%%
        % Adjusting the lower bound
        function set.a(B,value)
            B.a = value;
            B.Reset;
        end
        
        %%%
        % Adjusting the upper bound
        function set.b(B,value)
            B.b = value;
            B.Reset;
        end
        
        %%%
        % Adjusting the number of nodes
        function set.n(B,value)
            B.n = value;
            B.Reset;
        end
        
        %%%
        % Adjusting the type of nodes
        function set.nodetype(B,value)
            B.nodetype = value;
            B.Reset;
        end
        
        
    end %methods
    
    
end %classdef