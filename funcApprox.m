%% funcApprox class
% Defines a class to represent an approximated function
%
% Objects created by this class are a subclass of <basis.m basis>, adding fields to
% identify a function and methods to interpolate, compute Jacobian and Hessian.
%
%
% Apart from the properties inherited from <basis.m basis>, objects of class |funcApprox|
% have the following properties:
%
% * |df|:    scalar, number of interpolated functions
% * |fnodes|: value of interpolated function at basis nodes
% * |coef|:  interpolation coefficients
% * |Phi|:   interpolation matrix, evaluated at basis nodes
% * |Phiinv|: inverse of |Phi|
% * |fnames|:  cell of scalars, optional names for interpolated functions
%
% Object of class |funcApprox| have the following methods:
%
% * |funcApprox|: class constructor
% * |updateCoef|: computes interpolation coefficients if |fnodes| is modified
% * |Interpolate|: interpolates the function
% * |Jacobian|: computes the Jacobian of the function
% * |Hessian|: computes the Hessian of the function
% 
%
% *To date, only basisChebyshev has been implemented*, so calling the function with |type|
% 'spli' or 'lin' returns an error.
%
% Last updated: October 4, 2014.
%
%
% Copyright (C) 2014 Randall Romero-Aguilar
%
% Licensed under the MIT license, see LICENSE.txt

%%
classdef funcApprox < basis
 
    properties (Dependent)
        fnodes  % value of functions at nodes
        coef    % interpolation coefficients
        df      % number of functions
    end
    
    properties (Access = private)
       fnodes_  
       coef_
       fnodes_is_outdated
       coef_is_outdated
    end
    
    
    properties     
        fnames  % names for functions
    end
    
    
    properties (SetAccess = protected)
        Phi     % interpolation matrix
        Phiinv  % inverse of interpolation matrix 
    end
    
    
    
    
    methods
        
        %% funcApprox
        function F = funcApprox(B,fnodes,fnames)
            %%% 
            % Constructor for funcApprox
            %
            %   F = funcApprox(B,fnodes,fnames)
            % 
            % The inputs are:
            %
            % * |B|: basis object
            % * |fnodes|: matrix, values of the interpolant at basis B nodes
            % * |fnames|: cell of strings, names of functions (optional)
            
            if nargin==0
                return
            end
            
           if ~isa(B,'basis')
               B = basis.make(B);
           end
            
            
            for s = fieldnames(basis)'
                s1 = char(s);
                F.(s1) = B.(s1);
            end
            
            switch F.type
                case 'Chebyshev'
                    F.Phi = B.Interpolation;
                    F.Phiinv = (F.Phi'*F.Phi)\F.Phi';
                case 'Spline'
                    F.Phi = B.Interpolation;
            end
            
            switch nargin
                case 1
                    F.fnodes = zeros(size(F.nodes,1),1);
                    F.fnames = {'F1'};
                    return;
                case 2
                    F.fnodes = fnodes;
                    F.fnames = strcat('F',num2cell(48 + (1:F.df)));
                    return;
                    
                case 3
                    if ischar(fnames)
                        fnames = {fnames};
                    end
                    
                    if ~isempty(fnodes)
                        F.fnodes = fnodes;
                        F.fnames = fnames;
                    else
                        F.fnames = fnames;
                        F.fnodes = zeros(size(F.nodes,1),1);
                    end
            end
        end
        
        %% SETTERS for fnodes and coef
        function set.fnodes(F,value)
            assert(size(value,1) == size(F.nodes,1), ...
                'funcApprox:fnodes:size',...
                'New value for fnodes must have %d rows (1 per node)',size(F.nodes,1))
            F.fnodes_ = value;
            F.fnodes_is_outdated = false;
            F.coef_is_outdated = true;
        end
        
        
       function set.coef(F,value)
            assert(size(value,1) == size(F.Phi,2), ...
                'funcApprox:coef:size',...
                'New value for coef must have %d rows (1 per polynomial)',size(F.Phi,2))
            F.coef_ = value;
            F.coef_is_outdated = false;
            F.fnodes_is_outdated = true;
        end
        
        
        %% GETTERS for fnodes and coef
        function Fval = get.fnodes(F)
            if F.fnodes_is_outdated
                F.fnodes_ =  F.Phi * F.coef_;
                F.fnodes_is_outdated = false;
            end
            Fval = F.fnodes_;
        end
        
        function Fcoef = get.coef(F)
            if F.coef_is_outdated
                switch F.type
                    case 'Chebyshev'
                        F.coef_ = F.Phiinv * F.fnodes;
                    case 'Spline'
                        F.coef_ = F.Phi \ F.fnodes;
                end
                F.coef_is_outdated = false;
            end
            Fcoef = F.coef_;
        end
        
        function dfval = get.df(F)
            dfval = size(F.fnodes,2);
        end
        
        
        %% Interpolate
        function y = Interpolate(F,...
                varargin... inputs for interpolation method
                )
            
            %%%
            %   y = F.Interpolate(x,order,integrate)
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
            % Output |y| returns the interpolated functions as a k.m.h array, where k is
            % the number of evaluation points, m the number of functions (number of
            % columns in |f.coef|), and h is number of order derivatives.
            
            if nargin <2
                y = F.fnodes;
                return
            end
            
            Phix = F.Interpolation(varargin{:});
            nx = size(Phix,1);  % number of evaluation points
            no = size(Phix,3);  % number of order evaluations
            y = zeros(nx,F.df,no);
            
            if no==1  %horrible, but necessary because sparse cannot have 3 indices!
                y = Phix * F.coef;
            else
                for h = 1:no
                    y(:,:,h) = Phix(:,:,h) * F.coef;
                end
            end
            
            if F.df==1  % only one function
                y = squeeze(y);
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
            % * |Y|, k.m interpolated function (same as Y = B.Evaluate(coef,x), provided
            % for speed if both Jacobian and
            % funcion value are required).
                     
            %%%
            % Restrict function to compute
            if nargin<5 || isempty(indy)
                COEF = F.coef;
            else
                COEF = F.coef(:,indy);
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
            c = F.opts.validPhi;
            Ncols = size(c,1);
            
            %%%
            % Compute interpolation matrices for each dimension
            
            Phi0 = zeros(Nrows,Ncols,F.d);
            Phi1 = zeros(Nrows,Ncols,F.d);
            
            for k = 1:F.d
                if Required(k)
                    PhiOneDim = F.B1(k).Interpolation(x(:,k),...
                        [0 1],...
                        false);
                    Phi01 = PhiOneDim(:,c(:,k),:);
                    Phi0(:,:,k) = Phi01(:,:,1); % function value
                    Phi1(:,:,k) = Phi01(:,:,2); % its derivative
                else
                    PhiOneDim = F.B1(k).Interpolation(x(:,k),...
                        0,...
                        false);
                    Phi0(:,:,k) = PhiOneDim(:,c(:,k));
                end
            end
            
            %%%
            % Compute the Jacobian
            % Preallocate memory
            DY  = zeros(Nrows,F.df,nRequired);
            
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
            %   Hy = F.Hessian(x,coef)
            %
            % Computes the Hessian of a function approximated by basis |B| and coefficients |coef|.
            %
            % Its input is
            %
            % * |x|, k.d matrix of evaluation points.
            %
            % Its output |Hy| returns the k.m.d.d Hessian evaluated at |x|.
            
            
            if nargin<3 || isempty(indy)
                COEF = F.coef;
            else
                COEF = F.coef(:,indy);
            end

            
            
            
                        
            order = repmat({[0 1 2]'},1,F.d);
            order = gridmake(order{:});
            order = order(sum(order,2)==2,:);
            
            Phix = F.Interpolation(x,order,false);
            
            
            nx = size(x,1);     % number of evaluated points
            
            %Dy = squeeze(Dy);
            
            Hy = zeros(nx,F.df,F.d,F.d);
            
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
        
        
    end
    
end