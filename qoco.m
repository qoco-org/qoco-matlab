% Copyright (c) 2025 Govind M. Chari
% Copyright (c) 2017 Bartolomeo Stellato
%
% This source code is licensed under the BSD 3-Clause License

classdef qoco < handle
    % qoco interface class for QOCO solver
    % This class provides a Matlab interface to the C implementation
    % of the QOCO solver.
    %
    % qoco Properties:
    %   objectHandle - pointer to the C structure of QOCO solver
    %
    % qoco Methods:
    %
    %   setup             - configure solver with problem data
    %   solve             - solve the SOCP
    %
    %   default_settings  - create default settings structure
    %   current_settings  - get the current solver settings structure
    %   update_settings   - update the current solver settings structure
    %
    %   get_dimensions    - get the number of variables and constraints
    
    properties (SetAccess = private, Hidden = true)
        objectHandle % Handle to underlying C instance
    end
    methods(Static)
        %%
        function out = default_settings()
            % DEFAULT_SETTINGS get the default solver settings structure
            out = qoco_mex('default_settings', 'static');
        end
        
        %%
        function out = constant(constant_name)
            % CONSTANT Return solver constant
            %   C = CONSTANT(CONSTANT_NAME) return constant called CONSTANT_NAME
            out = qoco_mex('constant', 'static', constant_name);
        end
        
    end
    methods
        %% Constructor - Create a new solver instance
        function this = qoco(varargin)
            % Construct QOCO solver class
            this.objectHandle = qoco_mex('new', varargin{:});
        end
        
        %% Destructor - destroy the solver instance
        function delete(this)
            % Destroy QOCO solver class
            qoco_mex('delete', this.objectHandle);
        end
        
        %%
        function out = current_settings(this)
            % CURRENT_SETTINGS get the current solver settings structure
            out = qoco_mex('current_settings', this.objectHandle);
        end
        
        %%
        function update_settings(this,varargin)
            % UPDATE_SETTINGS update the current solver settings structure
            
            %second input 'false' means that this is *not* a settings
            %initialization, so some parameter/values will be disallowed
            newSettings = validateSettings(this,false,varargin{:});
            
            %write the solver settings.  C-mex does not check input
            %data or protect against disallowed parameter modifications
            qoco_mex('update_settings', this.objectHandle, newSettings);
            
        end
        
        %%
        function [m,n,p]  = get_dimensions(this)
            % GET_DIMENSIONS get the number of variables and constraints
            
            [m,n,p] = qoco_mex('get_dimensions', this.objectHandle);
            
        end
        
        %%
        function varargout = setup(this, varargin)
            % SETUP configure solver with problem data
            %
            %   setup(n, m, p, P, c, A, b, G, h, l, nsoc, q, options)
            
            nargin = length(varargin);
            
            %dimension checks on user data. Mex function does not
            %perform any checks on inputs, so check everything here
            assert(nargin >= 12, 'incorrect number of inputs');
            [n, m, p, P, c, A, b, G, h, l, nsoc, q] = deal(varargin{1:12});
            
            %
            % Get problem dimensions
            %
            
            if (isempty(P))
                P = sparse(n, n);
            else
                P = sparse(P);
            end
            if (~istriu(P))
                P = triu(P);
            end
            if (isempty(A))
                A = sparse(p, n);
            else
                A = sparse(A);
            end
            if (isempty(G))
                G = sparse(m, n);
            else
                G = sparse(G);
            end
            q = q(:);
            
            
            %
            % Check vector dimensions
            %
            
            assert(length(b) == p, 'Incorrect dimension of b');
            assert(length(h) == m, 'Incorrect dimension of h');
            assert(length(q) == nsoc, 'Incorrect dimension of q');
            
            %make a settings structure from the remainder of the arguments.
            %'true' means that this is a settings initialization, so all
            %parameter/values are allowed.  No extra inputs will result
            %in default settings being passed back
            theSettings = validateSettings(this,true,varargin{13:end});
            
            [varargout{1:nargout}] = qoco_mex('setup', this.objectHandle, n, m, p, P, c, A, b, G, h, l, nsoc, q, theSettings);
        end
        
        %%
        function varargout = solve(this, varargin)
            % SOLVE solve the SOCP
            
            nargoutchk(0,1);  %either return nothing (but still solve), or a single output structure
            [out.x, out.s, out.y, out.z, out.iters, out.setup_time_sec, out.solve_time_sec, out.obj, out.pres, out.dres, out.gap, out.status] = qoco_mex('solve', this.objectHandle);
            if(nargout)
                varargout{1} = out;
            end
            return;
        end
    end
end

function currentSettings = validateSettings(this,isInitialization,varargin)

%get the current settings
if(isInitialization)
    currentSettings = qoco_mex('default_settings', this.objectHandle);
else
    currentSettings = qoco_mex('current_settings', this.objectHandle);
end

%no settings passed -> return defaults
if(isempty(varargin))
    return;
end

%check for structure style input
if(isstruct(varargin{1}))
    newSettings = varargin{1};
    assert(length(varargin) == 1, 'too many input arguments');
else
    newSettings = struct(varargin{:});
end

%get the qoco settings fields
currentFields = fieldnames(currentSettings);

%get the requested fields in the update
newFields = fieldnames(newSettings);

%check for unknown parameters
badFieldsIdx = find(~ismember(newFields,currentFields));
if(~isempty(badFieldsIdx))
    error('Unrecognized solver setting ''%s'' detected',newFields{badFieldsIdx(1)});
end

%everything checks out - merge the newSettings into the current ones
for i = 1:length(newFields)
    currentSettings.(newFields{i}) = double(newSettings.(newFields{i}));
end

end