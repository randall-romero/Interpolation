%% interpolator class
% Defines a class to represent an approximated function
%
% Objects created by this class are a subclass of <basis.m basis>, adding fields to
% identify a function and methods to interpolate, compute Jacobian and Hessian.
%
%
% Apart from the properties inherited from <basis.m basis>, objects of class |funcApprox|
% have the following properties:
%
% * |y|: value of interpolated function at basis nodes
% * |c|:  interpolation coefficients
% * |Phi|:   interpolation matrix, evaluated at basis nodes
% * |Phiinv|: inverse of |Phi|
%
% Object of class |funcApprox| have the following methods:
%
% * |funcApprox|: class constructor
% * |updateCoef|: computes interpolation coefficients if |y| is modified
% * |Interpolate|: interpolates the function
% * |Jacobian|: computes the Jacobian of the function
% * |Hessian|: computes the Hessian of the function
%
%
% *To date, only basisChebyshev has been implemented*, so calling the function with |type|
% 'spli' or 'lin' returns an error.
%
% Last updated: November 24, 2014.
%
%
% Copyright (C) 2014 Randall Romero-Aguilar
%
% Licensed under the MIT license, see LICENSE.txt

%%
classdef interpolator < basis
    
    properties (Dependent)
        y    % value of functions at nodes
        c    % interpolation coefficients
        x    % basis nodes
    end
    
    properties %(Access = private)
        fnodes_  % stored values for function at nodes
        coef_    % stored coefficients
        fnodes_is_outdated % if true, calling "y" updates fnodes_ before returning values
        coef_is_outdated   % if true, calling "c" updates coef_ before returning values
    end
    
    properties (SetAccess = protected)
        Phi     % interpolation matrix
        Phiinv  % inverse of interpolation matrix
    end
    
    
    %%
    
    methods
        
        %% funcApprox
        function F = interpolator(B,fnodes)
            %%%
            % Constructor for funcApprox
            %
            %   F = funcApprox(B,fnodes)
            %
            % The inputs are:
            %
            % * |B|: basis object
            % * |fnodes|: matrix, values of the interpolant at basis B nodes
            
            if nargin==0
                return
            end
            
            % If B is not a 'basis', make one
            if ~isa(B,'basis');  B = basis.make(B);  end
            
            % Copy basis properties to new object F
            for s = fieldnames(basis)'
                s1 = char(s);
                F.(s1) = B.(s1);
            end
            
            % Compute interpolation matrix at nodes
            F.Phi = B.Interpolation;
            
            % Compute inverse of interpolation matrix, (not if Spline)
            if strcmp(F.type, 'Chebyshev')
                F.Phiinv = (F.Phi'*F.Phi)\F.Phi';
            end
            
            % Default values
            if nargin < 1
                F.y = zeros(size(F.nodes,1),1);
            else
                F.y = fnodes;
            end
        end
        
        %% SETTERS for "y" and "c"
        %
        %  Since the values of "y" and "c" are mutually dependent, changing one of them requires the other to be
        %  updated by using the interpolation matrix:  y = Phi(x) * c.
        % 
        %  To save time, this updating only takes place when the variable is queried.
        function set.y(F,value)
            assert(size(value,1) == size(F.nodes,1), ...
                'funcApprox:y:size',...
                'New value for y must have %d rows (1 per node)',size(F.nodes,1))
            F.fnodes_ = value;
            F.fnodes_is_outdated = false;
            F.coef_is_outdated = true;
        end
        
        
        function set.c(F,value)
            assert(size(value,1) == size(F.Phi,2), ...
                'funcApprox:c:size',...
                'New value for c must have %d rows (1 per polynomial)',size(F.Phi,2))
            F.coef_ = value;
            F.coef_is_outdated = false;
            F.fnodes_is_outdated = true;
        end
        
        
        %% GETTERS for "y" and "c"
        %
        % If the variable is not outdated, the get methods just returns the stored values. Otherwise, the values are
        % updated using the other variable and then returned.
        function Fval = get.y(F)
            if F.fnodes_is_outdated
                F.fnodes_ =  reshape(F.Phi * reshape(F.coef_,F.n,prod(F.size)),size(F.coef_));
                F.fnodes_is_outdated = false;
            end
            Fval = F.fnodes_;
        end
        
        function Fcoef = get.c(F)
            if F.coef_is_outdated
                switch F.type
                    case 'Chebyshev'
                        F.coef_ = reshape(F.Phiinv * reshape(F.fnodes_,F.n,prod(F.size)),size(F.fnodes_));
                    case 'Spline'
                        F.coef_ = reshape(F.Phi \ reshape(F.fnodes_,F.n,prod(F.size)),size(F.fnodes_));
                end
                F.coef_is_outdated = false;
            end
            Fcoef = F.coef_;
        end
        
        %% Size of the object and getter for "x"
        %
        % The size of the object is defined as the number of functions represented by "y" and "c".
        % "x" is just a shortcut for "nodes",
        
        function sz = size(F)
            if F.fnodes_is_outdated
                sz = size(F.coef_);
            else
                sz = size(F.fnodes_);
            end
            sz(1) = [];
        end
        
        function xx = get.x(F)
            xx = F.nodes;
        end
        
        
        %% Interpolate
        function Y = Interpolate(F,...
                varargin... inputs for interpolation method
                )
            
            %%%
            %   Y = F.Interpolate(x,order,integrate)
            %
            % Interpolates function f: R^d --> R^m at values |x|; optionally computes
            % derivative/integral. It obtains the interpolation matrices by calling
            % |Interpolation| on its basis.
            %
            % Its inputs are
            %
            % * |x|, k.d matrix of evaluation points.
            % * |order|, h.d matrix, for order of derivatives. Defaults to zeros(1,d).
            % * |integrate|, logical, integrate if true, derivative if false. Defaults to
            % false
            %
            % Output |Y| returns the interpolated functions as a k.m.h array, where k is
            % the number of evaluation points, m the number of functions (number of
            % columns in |f.c|), and h is number of order derivatives.
            
            if nargin <2
                Y = F.y;
                return
            end
            
            Phix = F.Interpolation(varargin{:});
            nx = size(Phix,1);  % number of evaluation points
            no = size(Phix,3);  % number of order evaluations
            Y = zeros(nx,F.size,no);
            
            if no==1  %horrible, but necessary because sparse cannot have 3 indices!
                Y = Phix * F.c;
            else
                for h = 1:no
                    Y(:,:,h) = Phix(:,:,h) * F.c;
                end
            end
            
            if F.size==1  % only one function
                Y = squeeze(Y);
            end
            
        end %Evaluate
        
        
        
        %% Jacobian
        function [DY,Y] = Jacobian(F, x, indx,indy, permuted)
            %%%
            %   [DY,Y] = F.Jacobian(x, index)
            %
            % Computes the Jacobian of the approximated function f: R^d --> R^m.
            %
            % Inputs:
            %
            % * |x|, k.d matrix of evaluation points.
            % * |index|, 1.d boolean, take partial derivative only wrt variables with
            % index=true. Defaults to true(1,d).
            %
            % Outputs
            %
            % * |DY|, k.m.d1 Jacobian matrix evaluated at |x|, where d1 = sum(index)
            % * |Y|, k.m interpolated function (same as Y = B.Evaluate(c,x), provided
            % for speed if both Jacobian and
            % funcion value are required).
            
            %%%
            % Restrict function to compute
            if nargin<5 || isempty(indy)
                COEF = F.c;
            elseif isa(indy,'cell')
                COEF = F.c(:,indy{:});
                COEF = reshape(COEF,F.n,[]);
            else
                COEF = F.c(:,indy);
            end
            
            
            
            %%%
            % Solve for the one-dimensional basis
            
            if F.d==1
                if nargout == 1
                    Phix  = F.Interpolation(x,1,false);
                    DY = Phix * COEF;
                else
                    Phix = F.Interpolation(x,[1;0],false);
                    DY = Phix(:,:,1) * COEF;
                    Y = Phix(:,:,2) * COEF;
                end
                
                return
            end
            
            %%%
            % Solve for the multi-dimensional basis
            
            % Check validity of input x
            assert(size(x,2) == F.d, 'In Jacobian, class basis: x must have d columns')
            
            %%%
            % Keep track of required derivatives: Required is logical with true indicating derivative is required
            if nargin<3 || isempty(indx)
                Required = true(1,F.d);
                indx = 1:F.d;
            elseif numel(indx) < F.d % assume that index have scalars of the desired derivatives
                Required = false(1,F.d);
                Required(indx) = true;
            else % assume that index is nonzero for desired derivatives
                Required = logical(indx);
                indx = find(indx);
            end
            
            
            
            
            
            if nargin<5
                permuted = true;
            end
            
            
            
            
            nRequired = sum(Required);
            
            %%%
            % HANDLE NODES DIMENSION
            Nrows = size(x,1);
            
            %%%
            % HANDLE POLYNOMIALS DIMENSION
            C = F.opts.validPhi;
            Ncols = size(C,1);
            
            %%%
            % Compute interpolation matrices for each dimension
            
            Phi0 = zeros(Nrows,Ncols,F.d);
            Phi1 = zeros(Nrows,Ncols,F.d);
            
            for k = 1:F.d
                if Required(k)
                    PhiOneDim = F.B1(k).Interpolation(x(:,k),...
                        [0 1],...
                        false);
                    Phi01 = PhiOneDim(:,C(:,k),:);
                    Phi0(:,:,k) = Phi01(:,:,1); % function value
                    Phi1(:,:,k) = Phi01(:,:,2); % its derivative
                else
                    PhiOneDim = F.B1(k).Interpolation(x(:,k),...
                        0,...
                        false);
                    Phi0(:,:,k) = PhiOneDim(:,C(:,k));
                end
            end
            
            %%%
            % Compute the Jacobian
            % Preallocate memory
            DY  = zeros(Nrows,F.size,nRequired);
            
            % Multiply the 1-dimensional bases
            
            for k=1:nRequired
                Phik = Phi0;
                Phik(:,:,indx(k)) = Phi1(:,:,indx(k));
                Phix = prod(Phik,3);
                DY(:,:,k) = Phix * COEF;
            end
            
            %%%
            % Compute the function if requested
            if nargout > 1
                Y = prod(Phi0,3) * COEF;
            end
            
            
            %%%
            % Permute the Jacobian, if requested
            %
            % Permute the resulting Jacobian so that indices reflect "usual" 2D Jacobian, with 3rd dimension used for
            % different observations.
            %
            % * i = function
            % * j = variable wrt which function is differentiated
            % * k = observation
            
            if permuted
                Y = permute(Y, [2,3,1]);
                DY = permute(DY, [2 3 1]);
            end
            
            
        end % Jacobian
        
        
        
        %% Hessian
        function Hy = Hessian(F,x,indy)
            
            %%%
            %   Hy = F.Hessian(x,c)
            %
            % Computes the Hessian of a function approximated by basis |B| and coefficients |c|.
            %
            % Its input is
            %
            % * |x|, k.d matrix of evaluation points.
            %
            % Its output |Hy| returns the k.m.d.d Hessian evaluated at |x|.
            
            
            if nargin<3 || isempty(indy)
                COEF = F.c;
            else
                COEF = F.c(:,indy);
            end
            
            
            
            
            
            order = repmat({[0 1 2]'},1,F.d);
            order = gridmake(order{:});
            order = order(sum(order,2)==2,:);
            
            Phix = F.Interpolation(x,order,false);
            
            
            nx = size(x,1);     % number of evaluated points
            
            %Dy = squeeze(Dy);
            
            Hy = zeros(nx,F.size,F.d,F.d);
            
            for k = 1:size(order,1)
                i = find(order(k,:));
                if numel(i)==1
                    Hy(:,:,i,i) = Phix(:,:,k) * COEF;%  Dy(:,k);
                else
                    Hy(:,:,i(1),i(2)) = Phix(:,:,k) * COEF; %Dy(:,k);
                    Hy(:,:,i(2),i(1)) = Hy(:,:,i(1),i(2)); %Dy(:,k);
                end
            end
        end % Hessian
        
        %% subsref
        
        function sref = subsref(F,s)
            switch s(1).type
                case '.'
                    sref = builtin('subsref',F,s);
                    return
                case '()'
                    sref = F.Interpolate(s(1).subs{:}); % TODO: Vargout
                    return
                case '{}'
                    if numel(s)==1
                        sref = copy(F);
                        if F.fnodes_is_outdated;
                            sref.c = F.coef_(:,s.subs{:});
                        else
                            sref.y = F.fnodes_(:,s.subs{:});
                        end
                        
                        return
                    else
                        switch s(2).type
                            case '()'
                                error('TBD: should call interpolate with restricted values of coefficients, instead of making copy of F');
                            case '.'
                                switch s(2).subs
                                    case {'y','c'}
                                        if numel(s)==2
                                            s2 = substruct('.',s(2).subs,'()',[{':'} s(1).subs]);
                                            sref = builtin('subsref',F,s2);
                                            return
                                        else
                                            error('Too many index levels');
                                        end
                                    otherwise
                                        sref = builtin('subsref',F,s(2));
                                end
                        end
                    end
            end
        end
        
        %% subsasgn
        function F = subsasgn(F,s,val)
            
            switch s(1).type
                case '.'
                    F = builtin('subsasgn',F,s,val);
                    return
                case '()'
                    if numel(s)==1
                        s2 = substruct('.','y','()',[{':'},s.subs]);
                        F = builtin('subsasgn',F,s2,val);
                        return
                    else
                        switch s(2).type
                            case '.'
                                
                                switch s(2).subs
                                    case {'y','c'}
                                        if numel(s)==2
                                            s2 = substruct('.',s(2).subs,'()',[{':'} s(1).subs]);
                                            F = builtin('subsasgn',F,s2,val);
                                            return
                                        else
                                            error('Too many index levels');
                                        end
                                    otherwise
                                        F = builtin('subsasgn',F,s(2),val);
                                end
                        end
                    end
                case '{}'
                    if numel(s)==1
                        s2 = substruct('.','c','()',[{':'},s.subs]);
                        F = builtin('subsasgn',F,s2,val);
                        return
                    else
                        error('On assigning coefficients with {i,j} indexing, no further indexing allowed')
                    end
                    
                    
                    
                    
            end
        end
        
        
        
        
        
        
        
        
        
    end
    
end