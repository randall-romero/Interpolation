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
 
    properties 
        fnodes  % value of functions at nodes
        fnames  % names for functions
    end
    
    
    properties (SetAccess = protected)
        Phi     % interpolation matrix
        Phiinv  % inverse of interpolation matrix 
        df      % number of functions
        coef    % interpolation coefficients
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
            
            for s = fieldnames(B)'
                s1 = char(s);
                F.(s1) = B.(s1);
            end
            
            F.Phi = B.Interpolation;
            F.Phiinv = (F.Phi'*F.Phi)\F.Phi';
            F.fnodes = fnodes;
            F.df = size(fnodes,2);
            
            if nargin>2
                F.fnames = fnames;
            end
            
        end
        
        %% set.fnodes
        function set.fnodes(F,value)
            %%%
            % Setting new values for fnodes also update the interpolation coefficients
            assert(size(value,1) == size(F.nodes,1), ...
                'funcApprox:fnodes:size',...
                'New value for fnodes must have %d rows (1 per node)',size(F.nodes,1))
            F.fnodes = value;
            F.updateCoef;
        end
        
        
        %% updateCoef
        function updateCoef(F)
            %%% 
            % To update coefficients, the inverse interpolation matrix is premultiplied by
            % the value of the function at the nodes
            F.df = size(F.fnodes,2);
            F.coef = F.Phiinv * F.fnodes;
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
            
            for h = 1:no
                y(:,:,h) = Phix(:,:,h) * F.coef;
            end
            
            if F.df==1  % only one function
                y = squeeze(y);
            end
            
        end %Evaluate
        
        
        
        %% Jacobian
        function [DY,Y] = Jacobian(F, x, index,permuted)
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
            % Solve for the one-dimensional basis
            
            if F.d==1
                if nargout == 1
                    Phix  = F.Interpolation(x,1,false);
                    DY = Phix * F.coef;
                else
                    Phix = F.Interpolation(x,[1;0],false);
                    DY = Phix(:,:,1) * F.coef;
                    Y = Phix(:,:,2) * F.coef;
                end
                
                return
            end
            
            %%%
            % Solve for the multi-dimensional basis
            
            % Check validity of input x
            assert(size(x,2) == F.d, 'In Jacobian, class basis: x must have d columns')
            
            %%%
            % Keep track of required derivatives: Required is logical with true indicating derivative is required
            if nargin<3 || isempty(index)
                Required = true(1,F.d);
                index = 1:F.d;
            elseif numel(index) < F.d % assume that index have scalars of the desired derivatives
                Required = false(1,F.d);
                Required(index) = true;
            else % assume that index is nonzero for desired derivatives
                Required = logical(index);
                index = find(index);
            end
            
            if nargin<4
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
                Phik(:,:,index(k)) = Phi1(:,:,index(k));
                Phix = prod(Phik,3);
                DY(:,:,k) = Phix * F.coef;
            end
            
            %%%
            % Compute the function if requested
            if nargout > 1
                Y = prod(Phi0,3) * F.coef;
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
        function Hy = Hessian(F,x)
            
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
                    Hy(:,:,i,i) = Phix(:,:,k) * F.coef;%  Dy(:,k);
                else
                    Hy(:,:,i(1),i(2)) = Phix(:,:,k) * F.coef; %Dy(:,k);
                    Hy(:,:,i(2),i(1)) = Hy(:,:,i(1),i(2)); %Dy(:,k);
                end
            end
        end % Hessian
        
        
    end
    
end