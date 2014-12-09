classdef ndsparse
    %ndsparse:  A class for multidimensional sparse arrays
    
    properties
        A = cell(0,1) % slice
    end
    
    properties (Dependent)
        nrow
        ncol
    end
    
    
    methods
        function S = ndsparse(nd)
            if nargin<1
                return
            end
            
            if isscalar(nd)
                nd = [nd,1];
            end
            
            S.A = cell(nd);
            
        end
        
        function Sout = subsasgn(S,s,B)
            assert(issparse(B),'Must be a sparse matrix')
            switch s.type
                case '()'
                    id = s.subs;
                    
                    Sout = S;
                    if all(strcmp(id(1:2),':'))
                        othSubs = id(3:end);
                        Sout.A{othSubs{:}} = B;
                    else
                        error('not yet implemented')
                    end
                    
                    
                    
                case '{}'
                    Sout = S;
                    Sout.A{s.subs{:}} = B;
            end    
        end
        
        
        function out = subsref(S,s)
            out = S;
            
            for j=1:numel(s)
                id = s(j).subs;
                
                if isa(out,'ndsparse')
                
                
                switch s(j).type
                    case '()'
                        sss = s(j);
                        sss.subs = id(3:end);
                        out.A = subsref(out.A,sss);
                        
                        id2 = id(1:2);
                        for k=1:numel(out.A)
                            out.A{k} = out.A{k}(id2{:});
                        end
                    case '{}'
                        out = out.A{id{:}};
                    case '.'
                        if strcmp(id,'A')
                            out = out.A;
                        end
                        
                end
                
                else
                    stemp.subs = id;
                    stemp.type = s(j).type;
                    out = subsref(out,stemp);
                end
            
            end
            if isa(out,'ndsparse') && numel(out.A)==1
                out = out.A;
                out = out{1};
            end
            
        end
        
        function sz = size(S,d)
            sz = [size(S.A{1}) size(S.A)];
            
            if nargin>1
                if d>numel(sz)
                    sz = 1;
                else
                    sz = sz(d);
                end
            end
        end
        
        
        function out = prod(S,d)
            assert(d>2,'not yet implemented')
            sz = size(S.A);
                        
            dc = d-2; %cell index
            
            if sz(dc)==1
                out = S;
                return;
            else
                idx = arrayfun(@(n) (1:n)',sz,'UniformOutput',false);
                idx{dc} = 0;
                rhs = num2cell(gridmake(idx{:}));
                rhs(:,dc) = {':'};
                szd = sz(dc);
                if numel(sz) > 2
                    sz(dc) = [];
                else
                    sz(dc) = 1;
                end
                    
                out = ndsparse(sz);
            end
            
            nx = size(rhs,1);
            lhs = rhs;
            lhs(:,dc) = [];
            
            data = S.A;
            for k=1:nx
                temp = data(rhs{k,:});
                temp2 = temp{1};
                for j=2:szd
                    temp2 = temp2 .* temp{j};
                end
                out.A{lhs{k,:}} = temp2;
            end
            
            if numel(out.A)==1
                out = out.A;
                out = out{1};
            end
            
        end
        
        function F = full(S)
            sz = size(S.A);
            
            if prod(sz,2) ==1
                F = full(S.A{1});
            else
                F = zeros([size(S.A{1}),sz]);
                for j=1:sz
                    F{:,:,j} = full(S.A{j});
                end
            end
        end
        
        
        
    end
end

