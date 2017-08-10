classdef MultiframeMotionCoupling< handle
    %MultiframeMotionCoupling
    % Couple all frames successively and solve the super resolution problem
    % EXAMPLE USAGE:
    % Given a video as 4-D matlab array vidLowRes, that we want to upsample
    % by a factor of 4, run
    % MMC = MultiframeMotionCoupling(vidLowRes,'factor',4);
    % MMC.init;
    % MMC.run;
    % vidHighRes = MMC.result1;
    % --- Parameters
    % After constructing the MMC object, you may change the following 
    % parameters directly:
    %
    %%% Algorithm parameters: %%%
    % MMC.factor             % Upsampling factor
    % MMC.alpha              % Regularization weight in super resolution
    % MMC.beta               % Regularization weight in optical flow
    % MMC.kappa              % Regularization weighting   
    % MMC.flowDirection
    % MMC.interpMethod       % Interpolation method
    % MMC.interpKernel       % Custom interpolation kernel
    % MMC.interpAA           % Interpolation AA
    %%% Implementation parameters: %%%
    % MMC.framework          % Choice of framework prost/flexBox/flexBox_vector
    % MMC.comp_mode          % mode preset, i.e. 'fast', 'accurate', 'fastest'
    % MMC.verbose            % verbose > 0 is text output, verbose > 1 is image output
    
    
    %% Class properties:
    properties

        % Problem Structure
        factor@double scalar

        % Output Data
        result1@double            % u (final result) colored
        result2@double            % u-w (temporally stabilized frames) colored
        
        % Solver options
        framework@char            % Choice of framework prost/flexBox
        comp_mode@char            % mode preset, i.e. 'fast' or 'accurate'         
        
        % Parameters
        alpha@double scalar       % weights for regularizer of u
        beta@double scalar        % weights for regularizer of v
        kappa@double scalar       % regularizer adjustment
        flowDirection@char        % 'forward','backward' or 'forward-backward' flow computation and usage
        k@double matrix           % blur kernel
        verbose@double scalar     % verbose mode
        
        % Sampling operator
        interpMethod@char         % Interpolation method
        interpKernel              % Custom interpolation kernel
        interpAA@logical          % Interpolation AA
        
    end
    
    properties (Access = private)
        
        % Input Data
        imageSequenceSmall@double 
        imageSequenceYCbCr@double 
        numFrames@double scalar 
        
        % Main variables
        u@double           % image sequence
        w@double 	       % spatial sequence
        v@double           % flow sequence
        MMCsolver          % main solver object for u subproblem
        
        % Options
        opts               % prost options          
        
        % Parameters
        dimsLarge@double vector          % Output dimensions
        dimsSmall@double vector          % Input dimensions
        h@double scalar                  % time distance heuristic
        tol_pd@double scalar             % primal dual tolerance (set by comp_mode)
        
        
        % Extra Input
        warpOp             % A given warping operator (deactivates flow computation !)
        
        % Fixed primal variables for first frame -> used to test batch video consistency
        u0_frame           % Given primal variable u that should correspond to the first low-res input frame
        w0_frame           % Given primal variable w that should correspond to the first low-res input frame
        
        % Book-keeping
        colorflag@logical scalar          % Is set to 1 for color videos to enable YCbCr treatment
        flowCompute@logical scalar        % no flow needs to be computed if it is given by mechanics or previous runs
        
    end
    
    %%%% Public Methods ---------------------------------------------------------------------------------------------------
    methods 
       %% Initialize object
        function obj = MultiframeMotionCoupling(imageSequenceSmall,varargin)
            vararginParser;
            
            %%% Property standards:
            
            % Input Data
            if ndims(imageSequenceSmall) == 3
                obj.imageSequenceSmall = imageSequenceSmall;
                obj.colorflag = 0;
                obj.imageSequenceYCbCr = 0;
                obj.numFrames = size(imageSequenceSmall,3);
                
            elseif ndims(imageSequenceSmall) == 4
                obj.colorflag = 1;
                obj.numFrames = size(imageSequenceSmall,4);
                
                ycbcrMap = zeros(size(imageSequenceSmall));
                for i = 1:obj.numFrames
                    ycbcrMap(:,:,:,i) = rgb2ycbcr(imageSequenceSmall(:,:,:,i));
                end
                obj.imageSequenceYCbCr = ycbcrMap;
                obj.imageSequenceSmall = squeeze(ycbcrMap(:,:,1,:));
                
            else
                error('Load an actual video sequence');
            end
            
            % Factor options
            if (exist('factor','var'))
                obj.factor = factor;
            else
                obj.factor = 4;
            end
            
            % Problem Structure defaults
            obj.dimsSmall = size(imageSequenceSmall);
            obj.dimsSmall = obj.dimsSmall(1:2);
            obj.dimsLarge = obj.factor*obj.dimsSmall;  % will be overwritten later
            
            
            % Main variables
            %u                  % image sequence
            %v                  % flow sequence
            %k                  % blur kernel
            % will be initialized in init call
            
            % Given flow field:
            if exist('flowField','var')
                if ~isscalar(flowField)
                    obj.v = flowField;
                    obj.flowCompute = 0;
                else
                    obj.flowCompute = 1;
                end
            else
                obj.flowCompute = 1;
            end
            
            % Verbose options
            if (exist('verbose','var'))
                obj.verbose = verbose;
            else
                obj.verbose = 1;
            end
            

            
            % Given warp operator
            if (exist('warpOp','var'))
                if ~isscalar(warpOp)
                    obj.warpOp = warpOp;
                    obj.flowCompute = 0;
                else
                    obj.warpOp = 0;
                end
            else
                obj.warpOp = 0;
            end
            
            % Given fixed primal variables for first frame -> used to test batch video consistency
            if exist('u0_frame','var')
                obj.u0_frame = u0_frame;
            else
                obj.u0_frame = [];
            end
            
            if exist('w0_frame','var')
                obj.w0_frame = w0_frame;
            else
                obj.w0_frame = [];
            end
            
            if (isempty(obj.w0_frame) && ~isempty(obj.u0_frame)) || (isempty(obj.u0_frame) && ~isempty(obj.w0_frame))
                error('Give both fix frames u0 and w0.')
            end
            
            % Solver options
            %MMCsolver                                 % main solver class for u subproblem
            % will be constructed in init_u_... call
            
            % default SR solver:
            obj.framework = 'flexBox';
            obj.comp_mode = 'accurate';
            
            % Will be initialized in init_u calls
            %obj.opts.backend                         % prost backend options
            %obj.opts.opts
            
            % Parameter standards
            obj.alpha  = 0.01;                        % weights for regularizer of u
            obj.beta   = 0.2;                         % weights for regularizer of v
            obj.kappa  = 0.25;
            obj.h      = 1;
            obj.flowDirection = 'backward';
            
            % standard interpolation parameters
            obj.interpMethod  = 'average';           % Interpolation method
            obj.interpKernel  = [];                  % Custom interpolation kernel
            obj.interpAA      = false  ;             % Interpolation AA
            
        end
        
        %% collect inits and validate
        function init(obj)
            
            if obj.verbose > 0
                disp('Initialization...')
            end
            
            %%% Validate all public parameters %%% (ok atleast somewhat)
            MultiframeMotionCoupling.checkposScalar(obj.factor);
            MultiframeMotionCoupling.checkposScalar(obj.alpha);
            MultiframeMotionCoupling.checkposScalar(obj.beta);
            MultiframeMotionCoupling.checkposScalar(obj.kappa);
            
            
            % Possible inputs
            posFlow = {'forward','backward','forward-backward'};
            posFW   = {'prost','flexBox','flexBox-vector'};
            posIM   = {'nearest','bilinear','bicubic','bicubic-0',      ...
                       'lanczos2','lanczos3','pchip',                   ...
                       'custom','average','stride'};
            posComp = {'accurate','fast','fastest'};
            
            validatestring(obj.flowDirection,posFlow);
            validatestring(obj.framework,posFW);
            validatestring(obj.comp_mode,posComp);
            validatestring(obj.interpMethod,posIM);
            
            % Error if flexbox + fastest
            if strcmp(obj.comp_mode,'fastest') && ~strcmp(obj.framework,'prost')
                error(' "Fastest" mode can currently only be used with the "prost" framework');
            end
            
            % Check installed frameworks
            if strcmp(obj.framework,'prost')
                try
                    evalc('prost.list_gpus');
                catch
                    error('prost framework not found on path, or no GPUs found with prost');
                end
            elseif strcmp(obj.framework,'flexBox')
                if ~exist('flexBox','file')==2
                    error(['Cannot find flexBox framework, please install the framework and',...
                          ' add it to your MATLAB path']);
                end
            elseif strcmp(obj.framework,'flexBox-vector')
                if ~exist('flexBox','file')==2
                    error(['Cannot find flexBox framework, please install the framework and',...
                          ' add it to your MATLAB path']);
                end
            end
            
            
            % update sizes to account for factor changes
            obj.dimsLarge = obj.factor*obj.dimsSmall;
            
            % Set tolerances
            if strcmp(obj.comp_mode,'accurate')
                obj.tol_pd = 1e-6;
            else
                obj.tol_pd = 1e-4;
            end
            
            % Set fastest mode
            if strcmp(obj.comp_mode,'fastest')
                if  ~isnan(obj.kappa)
                    warning(['You have chosen the "fastest" computation mode.',newline,...
                             'As such your choice of kappa will be overwritten, ',newline, ...
                             'because this mode uses only additive regularization.'])
                end
                obj.kappa = NaN;
            end
            
            % Call actual flow method
            if obj.flowCompute
                obj.v = zeros([obj.dimsLarge,obj.numFrames-1,2]);
                obj.init_v0;
            end
            
            % Init u with desired framework
            if strcmp(obj.framework,'prost')
                obj.init_u_prost;
            elseif strcmp(obj.framework,'flexBox')
                obj.init_u_flexBox;
            elseif strcmp(obj.framework,'flexBox_vector')
                obj.init_u_flexBox_vector;
            else
                error('Invalid framework string.');
            end
            
            
            if (obj.verbose > 0)
                disp('Initialization finished');
            end
        end 
        
        %% Run the algorithm
        function run(obj)
            if isempty(obj.MMCsolver)
                error('No initialization found, run the "init" method first.');
            end
            % solve u problem
            disp('Solving problem for u');
            obj.solveU;
            
            disp('-------Main Iteration finished!')
            
            
            % Recompute RGB image in color mode
            if obj.colorflag
                obj.recomputeRGB
            else
                obj.result1 = obj.u;
                obj.result2 = obj.w;
            end
            
            disp('Algorithm finished !!');
        end
        
        
    end
    
    %%%% Priavte Methods ---------------------------------------------------------------------------------------------------
    methods (Access = private)
        
        
        %% create prost object for u
        function init_u_prost(obj)
            
            %%%% Create backend and options
            obj.opts.backend = prost.backend.pdhg('stepsize', 'boyd',   ...
                                                   'tau0', 1,           ...
                                                   'sigma0', 1);                       
            if strcmp(obj.comp_mode,'accurate')
                iters = 7500;
            else
                iters = 2500;
            end
            obj.opts.opts = prost.options('max_iters', iters,           ...
                                          'num_cback_calls', 5,         ...
                                          'verbose', true,              ...
                                          'tol_rel_primal', obj.tol_pd, ...
                                          'tol_rel_dual', obj.tol_pd,   ...
                                          'tol_abs_dual', obj.tol_pd,   ...
                                          'tol_abs_primal', obj.tol_pd);               
            
            %%%% Create operators
            
            % Call sampling function to construct matrix representation
            dsOp = samplingOperator(obj.dimsLarge,obj.dimsSmall,        ...
                                    obj.interpMethod,obj.interpKernel,  ...
                                    obj.interpAA);
            
            % Use blur kernel:
            if ~isempty(obj.k) || isscalar(obj.k)
                dsOp = dsOp*RepConvMtx(obj.k,obj.dimsLarge);
            end
            
            if obj.verbose > 0
                disp('Downsampling operator constructed');
            end
            % shorten some notation
            temp = obj.dimsSmall;
            ny=temp(1); nx=temp(2);
            temp = obj.dimsLarge;
            Ny = temp(1);
            Nx = temp(2);
            nc =  obj.numFrames;
            
            % Call warp operator constructor
            if isscalar(obj.warpOp)
                
                % Construct warp in specified direction
                if strcmp(obj.flowDirection,'forward')
                    warpingOp = constructWarpFF(obj.v,'F-I');
                elseif strcmp(obj.flowDirection,'backward')
                    warpingOp = constructWarpFF(obj.v,'I-F');
                elseif strcmp(obj.flowDirection,'forward-backward')
                    warpingOp = constructWarpFB(obj.v);
                end
                % Save to object
                obj.warpOp = warpingOp;
                if obj.verbose > 0
                    disp('Warp operator constructed');
                end
            else
                % Otherwise load given warp operator
                % Make sure that the loaded warp operator comes with the
                % desired flow direction, this is not checked!
                warpingOp = obj.warpOp;
                if obj.verbose > 0
                    disp('Warp operator loaded');
                end
            end
            
            % Compute bicubic estimate
            u_up = zeros(Ny,Nx,nc);
            for i = 1:nc
                u_up(:,:,i) = imresize(obj.imageSequenceSmall(:,:,i),   ...
                                       obj.factor);
            end
            
            
            % Check warp accuracy visually on bicbubic estimate
            if obj.verbose > 1
                u_mv = reshape(warpingOp*u_up(:),Ny,Nx,nc);
                for i = 1:nc-1
                    figure(i), imagesc(u_mv(:,:,i)),
                    title('Visual test of warp accuracy')
                end
            end
            
            %%% Preweight, based on warping accuracy to compute h
            warpPix = sum(abs(warpingOp*u_up(:)))/numel(u_up);
            disp(['Warp energy on bicubic (per pixel): ',num2str(warpPix)]);
            GradOp = spmat_gradient2d(Nx, Ny, nc);
            gradPix = sum(abs(GradOp*u_up(:)))/numel(u_up);
            disp(['Gradient energy on bicubic (per pixel): ',num2str(gradPix)])
            obj.h  = warpPix/gradPix;
            %obj.h = 1;                 % If you want to test without 
                                        % the h, you can manually disable
                                        % at this position.
            
            %%%% Initialize prost variables and operators
            
            % Primal- Core Variables:
            u_vec = prost.variable(Nx*Ny*nc);
            if isscalar(obj.u)
                u_vec.val = u_up(:);
            else
                u_vec.val = obj.u(:);
            end
            w_vec = prost.variable(Nx*Ny*nc);
            if isscalar(obj.w)
                w_vec.val = obj.w(:);
            else
                w_vec.val = u_up(:);
            end
            p = prost.variable(nx*ny*nc);
            g1 = prost.variable(3*Nx*Ny*nc);
            g2 = prost.variable(3*Nx*Ny*nc);
            
            % The following lines construct the primal dual problem 
            %--kappa in (0,1) builds the infimal convolution regularization
            %  (R_spat [IC] R_temp)(u)
            %  as discussed in the paper. This uses a sizable amount of 
            %  memory.
            %--kappa = 1 is mathematically the same algorithm, although the
            %    number of needed dual variables in the implementation is 
            %    halved to decrease memory consumption.
            %--kappa = NaN is the additive regularizer variant discussed in
            %  the paper with 
            %  R(u) = ||\nabla u|| + || W u ||
            %  this is not the best regularizer, but definitely a faster 
            %  and inexpensive option if the accurate regularizer is not
            %  possible.
            %--kappa in (0,1) and u0_frame nonempty is the same algorithm
            %  as option 1 with the frame given in u0_frame (and possibly
            %  also w0_frame used as boundary constraint.
            
            % Connect variables to primal-dual problem
            if obj.kappa ~= 1 && ~isnan(obj.kappa) && isempty(obj.u0_frame)
                
                % Generate min-max problem structure
                obj.MMCsolver = prost.min_max_problem( {u_vec,w_vec}, {p,g1,g2} );
                
                
                % Initialize full downsampling implicitly as kron(id(numFrames,D)
                obj.MMCsolver.add_dual_pair(u_vec, p, prost.block.id_kron_sparse(dsOp,obj.numFrames));
                
                % Add inf-conv duals
                obj.MMCsolver.add_dual_pair(u_vec, g1, prost.block.sparse([warpingOp/obj.h;obj.kappa*spmat_gradient2d(Nx, Ny,nc)]));
                obj.MMCsolver.add_dual_pair(w_vec, g1, prost.block.sparse(-[warpingOp/obj.h;obj.kappa*spmat_gradient2d(Nx, Ny,nc)]));
                obj.MMCsolver.add_dual_pair(w_vec, g2, prost.block.sparse([obj.kappa*warpingOp/obj.h;spmat_gradient2d(Nx, Ny,nc)]));
                
                % Add functions
                obj.MMCsolver.add_function(u_vec,prost.function.sum_1d('ind_box01',1,0,1,0,0));
                obj.MMCsolver.add_function(p, prost.function.sum_1d('ind_box01', 0.5, -0.5, 1, obj.imageSequenceSmall(:), 0)); %l^1
                obj.MMCsolver.add_function(g1, prost.function.sum_norm2(3, false, 'ind_leq0', 1/obj.alpha, 1, 1, 0, 0)); %l^{2,1}
                obj.MMCsolver.add_function(g2, prost.function.sum_norm2(3, false, 'ind_leq0', 1/obj.alpha, 1, 1, 0, 0)); %l^{2,1}
                
            elseif obj.kappa == 1    % Simplify if no infimal convolution is necessary
                
                % Generate min-max problem structure
                obj.MMCsolver = prost.min_max_problem( {u_vec}, {p,g1} );
                
                % Initialize full downsampling implicitly as kron(id(numFrames,D)
                obj.MMCsolver.add_dual_pair(u_vec, p, prost.block.id_kron_sparse(dsOp,obj.numFrames));
                
                % Add regularizer dual
                obj.MMCsolver.add_dual_pair(u_vec, g1, prost.block.sparse([warpingOp/obj.h;spmat_gradient2d(Nx, Ny,nc)]));
                % add functions
                obj.MMCsolver.add_function(u_vec,prost.function.sum_1d('ind_box01',1,0,1,0,0));
                obj.MMCsolver.add_function(p, prost.function.sum_1d('ind_box01', 0.5, -0.5, 1, obj.imageSequenceSmall(:), 0)); %l^1
                obj.MMCsolver.add_function(g1, prost.function.sum_norm2(3, false, 'ind_leq0', 1/obj.alpha/2, 1, 1, 0, 0)); %l^{2,1}
                
            elseif isnan(obj.kappa)    % Do additve TV and Warp
                
                g1 = prost.variable(Nx*Ny*nc);
                g2 = prost.variable(2*Nx*Ny*nc);
                
                % Generate min-max problem structure
                obj.MMCsolver = prost.min_max_problem( {u_vec}, {p,g1,g2} );
                
                % Initialize full downsampling implictely as kron(id(numFrames,D)
                obj.MMCsolver.add_dual_pair(u_vec, p, prost.block.id_kron_sparse(dsOp,obj.numFrames));
                
                % Add regularizer dual
                obj.MMCsolver.add_dual_pair(u_vec, g1, prost.block.sparse(warpingOp/obj.h)); % /obj.h
                obj.MMCsolver.add_dual_pair(u_vec, g2, prost.block.gradient2d(Nx, Ny,nc,false));
                % add functions
                obj.MMCsolver.add_function(u_vec,prost.function.sum_1d('ind_box01',1,0,1,0,0));
                obj.MMCsolver.add_function(p, prost.function.sum_1d('ind_box01', 0.5, -0.5, 1, obj.imageSequenceSmall(:), 0)); %l^1
                obj.MMCsolver.add_function(g1, prost.function.sum_1d('ind_box01', 0.5/obj.alpha, -0.5, 1, 0, 0)); %l^1
                obj.MMCsolver.add_function(g2, prost.function.sum_norm2(2, false, 'ind_leq0', 1/obj.alpha, 1, 1, 0, 0)); %l^{2,1}
                
            elseif obj.kappa ~= 1 && ~isnan(obj.kappa) && ~isempty(obj.u0_frame)
                
                % Add new variables to fix first frame of primal variables
                q1 = prost.variable(Nx*Ny);
                q2 = prost.variable(Nx*Ny);
                
                % Add subvariables to first frame primals
                u1 = prost.sub_variable(u_vec,Nx*Ny);
                w1 = prost.sub_variable(w_vec,Nx*Ny);
                u2 = prost.sub_variable(u_vec,Nx*Ny*(nc-1)); %#ok<NASGU> % these are used internallly!
                w2 = prost.sub_variable(w_vec,Nx*Ny*(nc-1)); %#ok<NASGU>
                
                % Generate min-max problem structure
                obj.MMCsolver = prost.min_max_problem( {u_vec,w_vec}, {p,g1,g2,q1,q2} );
                
                
                % Initialize full downsampling implicitly as kron(id(numFrames,D)
                obj.MMCsolver.add_dual_pair(u_vec, p, prost.block.id_kron_sparse(dsOp,obj.numFrames));
                
                % Add inf-conv duals
                obj.MMCsolver.add_dual_pair(u_vec, g1, prost.block.sparse([warpingOp/obj.h;obj.kappa*spmat_gradient2d(Nx, Ny,nc)]));
                obj.MMCsolver.add_dual_pair(w_vec, g1, prost.block.sparse(-[warpingOp/obj.h;obj.kappa*spmat_gradient2d(Nx, Ny,nc)]));
                obj.MMCsolver.add_dual_pair(w_vec, g2, prost.block.sparse([obj.kappa*warpingOp/obj.h;spmat_gradient2d(Nx, Ny,nc)]));
                
                % Add first frame constraint duals <u,q>
                obj.MMCsolver.add_dual_pair(u1,q1,prost.block.identity);
                obj.MMCsolver.add_dual_pair(w1,q2,prost.block.identity);
                % Add constraint linears <-u0,q>
                obj.MMCsolver.add_function(q1,prost.function.sum_1d('zero', 1, 0, 1, obj.u0_frame(:), 0));
                obj.MMCsolver.add_function(q2,prost.function.sum_1d('zero', 1, 0, 1, obj.w0_frame(:), 0));
                
                
                % Add dual set functions
                obj.MMCsolver.add_function(u_vec,prost.function.sum_1d('ind_box01',1,0,1,0,0));
                obj.MMCsolver.add_function(p, prost.function.sum_1d('ind_box01', 0.5, -0.5, 1, obj.imageSequenceSmall(:), 0)); %l^1
                obj.MMCsolver.add_function(g1, prost.function.sum_norm2(3, false, 'ind_leq0', 1/obj.alpha, 1, 1, 0, 0)); %l^{2,1}
                obj.MMCsolver.add_function(g2, prost.function.sum_norm2(3, false, 'ind_leq0', 1/obj.alpha, 1, 1, 0, 0)); %l^{2,1}
            end
            
        end
        
        
        %% create vector flexBox object for u
        function init_u_flexBox_vector(obj)
            
            %%%% Create operators
            
            % Call sampling function to construct matrix representation
            dsOp = samplingOperator(obj.dimsLarge,obj.dimsSmall,obj.interpMethod,obj.interpKernel,obj.interpAA);
            
            % Use blur kernel:
            if ~isempty(obj.k) || isscalar(obj.k)
                dsOp = dsOp*RepConvMtx(obj.k,obj.dimsLarge);
            end
            
            if obj.verbose > 0
                disp('Downsampling operator constructed');
            end
            % shorten some notation
            temp = obj.dimsLarge;
            Ny = temp(1);
            Nx = temp(2);
            nc =  obj.numFrames;
            
            % Call warp operator constructor
            if isscalar(obj.warpOp)
                
                % Construct warp in specified direction
                if strcmp(obj.flowDirection,'forward')
                    warpingOp = constructWarpFF(obj.v,'F-I');
                elseif strcmp(obj.flowDirection,'backward')
                    warpingOp = constructWarpFF(obj.v,'I-F');
                elseif strcmp(obj.flowDirection,'forward-backward')
                    warpingOp = constructWarpFB(obj.v);
                end
                % Save to object
                obj.warpOp = warpingOp;
                if obj.verbose > 0
                    disp('Warp operator constructed');
                end
            else
                % Otherwise load given warp operator
                % Make sure that the loaded warp operator comes with the
                % desired flow direction
                warpingOp = obj.warpOp;
                if obj.verbose > 0
                    disp('Warp operator loaded');
                end
            end
            
            % Build gradient matrices
            D = spmat_gradient2d(Nx, Ny,nc);
            Dx = D(1:Nx*Ny*nc,:); Dy = D(Nx*Ny*nc+1:end,:);
            
            % Compute bicubic estimate
            u_up = zeros(Ny,Nx,nc);
            for i = 1:nc
                u_up(:,:,i) = imresize(obj.imageSequenceSmall(:,:,i),obj.factor);
            end
            
            
            % Check warp accuracy visually on bicbubic estimate
            if obj.verbose > 1
                u_mv = reshape(warpingOp*u_up(:),Ny,Nx,nc);
                for i = 1:nc-1
                    figure(i), imagesc(u_mv(:,:,i)),title('Visual test of warp accuracy')
                end
            end
            
            %%% Preweight, based on warping accuracy to compute h
            warpPix = sum(abs(warpingOp*u_up(:)))/numel(u_up);
            disp(['Warp energy on bicubic (per pixel): ',num2str(warpPix)]);
            GradOp = spmat_gradient2d(Nx, Ny, nc);
            gradPix = sum(abs(GradOp*u_up(:)))/numel(u_up);
            disp(['Gradient energy on bicubic (per pixel): ',num2str(gradPix)])
            obj.h  = warpPix/gradPix;
            
            %%%% Initialize flexBox variables and operators
            obj.MMCsolver = flexBox;
            obj.MMCsolver.params.tol = obj.tol_pd;
            obj.MMCsolver.params.tryCPP = 1;
            obj.MMCsolver.params.verbose = 1;
            if strcmp(obj.comp_mode,'accurate')
                obj.MMCsolver.params.maxIt = 2500;
            else
                obj.MMCsolver.params.maxIt = 800;
            end
            
            % Add primal variables u and w
            obj.MMCsolver.addPrimalVar([Ny,Nx,nc]); % u
            obj.MMCsolver.addPrimalVar([Ny,Nx,nc]); % w
            
            % Build data term
            obj.MMCsolver.addTerm(L1dataTermOperator(1,kron(speye(nc),dsOp),obj.imageSequenceSmall(:)),1);
            obj.MMCsolver.addTerm(boxConstraint(0,1,[Ny,Nx,nc]),1);
            
            
            % Build infconv regularizer
            % Add first infconv part
            flexA1 = { warpingOp/obj.h, -warpingOp/obj.h, ...
                obj.kappa*Dx,    - obj.kappa*Dx,   ...
                obj.kappa*Dy,    - obj.kappa*Dy        };
            
            obj.MMCsolver.addTerm(L1operatorIso(obj.alpha,2,flexA1),[1,2]);
            % Add second infconv part
            flexA2 = { obj.kappa*warpingOp/obj.h, ...
                Dx,              ...
                Dy                     };
            obj.MMCsolver.addTerm(L1operatorIso(obj.alpha,1,flexA2),2);
            
            % set bicubic starting vector for u and w
            obj.MMCsolver.x{1,1} = u_up(:);
            obj.MMCsolver.x{1,2} = u_up(:);
        end
        %% create flexBox object for u
        function init_u_flexBox(obj)
            
            % shorten some notation
            temp = obj.dimsLarge;
            Ny = temp(1);
            Nx = temp(2);
            N = Nx*Ny;
            nc =  obj.numFrames;
            
            % Failsafe
            if isnan(obj.kappa) || (obj.kappa) == 0
                error('Additive warp and full separation not implemented yet for flexBox.');
            end
            
            %%%% Create operators
            
            % Build downsampling matrix
            if strcmp(obj.interpMethod,'average') && ~obj.interpAA
                dsOp = superpixelOperator(obj.dimsSmall,obj.factor);
            else
                dsOp = samplingOperator(obj.dimsLarge,obj.dimsSmall,    ...
                                        obj.interpMethod,               ...
                                        obj.interpKernel,               ...
                                        obj.interpAA);                  ...
            end
            
            % Add blur matrix
            if ~isempty(obj.k) || isscalar(obj.k)
                %blurOp = convolutionOperator([Ny,Nx],7,sqrt(0.6)); 
                %[blurOP is not supported in flexBoxC++]
                
                %[use this for speed if necessary:]
                %dsOp = dsOp*RepConvMtx(obj.k,obj.dimsLarge);
                % standard:
                dsOp  = concatOperator(dsOp,                            ...
                                       RepConvMtx(obj.k,obj.dimsLarge), ...
                                       'composition');                  
            end
            
            if obj.verbose > 0
                disp('Downsampling operator constructed');
            end
            
            % Construct warp operators and reduced identities
            for i = 1:nc-1
                
                singleField = squeeze(obj.v(:,:,i,:));
                Wi   = warpingOperator(obj.dimsLarge,singleField);
                Id   = ones(N,1);
                
                % Remove out-of range warps
                marker = sum(abs(Wi),2) == 0;
                Wi(marker > 0,:) = 0;
                Id(marker > 0) = 0; 
                
                warps{i} = Wi;  %#ok<AGROW>
                Ids{i}   = Id;  %#ok<AGROW>
            end
            if obj.verbose > 0
                disp('Warp operators constructed');
            end
            
            % Build gradient matrices
            %D = spmat_gradient2d(Nx, Ny,1);
            %Dx = D(1:Nx*Ny,:); Dy = D(Nx*Ny+1:end,:);
            Dx = gradientOperator([Ny,Nx],2);
            Dy = gradientOperator([Ny,Nx],1);
            
            %%%% Compute adaptive parameter h
            
            % Compute bicubic estimate
            u_up = zeros([obj.dimsLarge,nc]);
            for i = 1:nc
                u_up(:,:,i) = imresize(obj.imageSequenceSmall(:,:,i),   ...
                                       obj.factor,'bicubic');
            end
            % Accumulate estimates
            Du = zeros(nc,1);Wu = zeros(nc,1);
            for i = 1:nc-1
                ui = u_up(:,:,i);
                uj = u_up(:,:,i+1);
                Du(i) = sum(abs(Dx*ui(:))+abs(Dy*ui(:)));
                if strcmp(obj.flowDirection,'forward')
                    Wu(i) = sum(abs(-Ids{i}.*uj(:)+warps{i}*ui(:)));
                elseif strcmp(obj.flowDirection,'backward')
                    Wu(i) = sum(abs(Ids{i}.*ui(:)-warps{i}*uj(:)));
                elseif strcmp(obj.flowDirection,'forward-backward')
                    error('todo'); 
                    % todoWu(i) = sum(abs(-Ids{i}*uj(:)+warps{i}*ui(:)));
                else
                    error('');
                end
            end
            ui = u_up(:,:,nc);
            Du(nc) = sum(abs(Dx*ui(:))+abs(Dy*ui(:)));
            gradPix = sum(Du)/N/nc; warpPix = sum(Wu)/N/nc;
            obj.h = warpPix/gradPix;
            if obj.verbose > 0
                disp(['Warp energy on bicubic (per pixel)    : ', num2str(warpPix)]);
                disp(['Gradient energy on bicubic (per pixel): ', num2str(gradPix)]);
            end
            
            %%%% initialize flexBox solver
            
            obj.MMCsolver = flexBox;
            obj.MMCsolver.params.tol = obj.tol_pd;
            obj.MMCsolver.params.tryCPP = 1;
            obj.MMCsolver.params.verbose = 1;
            if strcmp(obj.comp_mode,'accurate')
                obj.MMCsolver.params.maxIt = 2500;
            else
                obj.MMCsolver.params.maxIt = 800;
            end
            
            % Add primal variables u and w
            for i = 1:2*nc
                obj.MMCsolver.addPrimalVar([Ny,Nx]);
            end
            
            % Build data terms
            for i = 1:nc
                
                % Add L1-data
                f_i = obj.imageSequenceSmall(:,:,i);
                obj.MMCsolver.addTerm(L1dataTermOperator(1,dsOp,f_i(:)),i);
                
                % Add box constraint
                obj.MMCsolver.addTerm(boxConstraint(0,1,[Ny,Nx]),i);
            end
            
            % scale Gradients
            kappaDx = concatOperator(diagonalOperator(obj.kappa*ones(N,1)),...
                                     Dx,'composition');
            kappaDy = concatOperator(diagonalOperator(obj.kappa*ones(N,1)),...
                                     Dy,'composition');
            
            % Build nc-1 infconv regularizers
            if strcmp(obj.flowDirection,'forward')
                for i = 1:nc-1
                    % Get warp and identity
                    Wi = warps{i}/obj.h;
                    Id = diagonalOperator(Ids{i}/obj.h);
                    
                    % Add first infconv part
                    flexA1 = {Wi              ,-Id            ,-Wi          ,Id,              ...
                              kappaDx         ,zeroOperator(N),-kappaDx     ,zeroOperator(N), ...
                              kappaDy         ,zeroOperator(N),-kappaDy     ,zeroOperator(N)     };
                    
                    obj.MMCsolver.addTerm(L1operatorIso(obj.alpha,4,flexA1),[i,i+1,nc+i,nc+i+1]);
                    % Add second infconv part
                    Id = diagonalOperator(-obj.kappa*Ids{i}/obj.h);
                    flexA2 = {obj.kappa*Wi  ,Id,...
                              Dx            ,zeroOperator(N), ...
                              Dy            ,zeroOperator(N)      };
                    
                    obj.MMCsolver.addTerm(L1operatorIso(obj.alpha,2,flexA2),[nc+i,nc+i+1]);
                end
                % Build last infconv regularizer
                flexA1 = {kappaDx,-kappaDx, ...
                          kappaDy,-kappaDy      };
                
                obj.MMCsolver.addTerm(L1operatorIso(obj.alpha,2,flexA1),[nc,2*nc]);
                flexA2 = {Dx, ...
                          Dy      };
                
                obj.MMCsolver.addTerm(L1operatorIso(obj.alpha,1,flexA2),2*nc);
                
            elseif strcmp(obj.flowDirection,'backward')
                for i = 1:nc-1
                    % Get warp and identity
                    Wi = warps{i}/obj.h;
                    Id = diagonalOperator(Ids{i}/obj.h);
                    
                    % Add first infconv part
                    flexA1 = {Id              ,-Wi            ,Id          ,-Wi,             ...
                              kappaDx         ,zeroOperator(N),-kappaDx    ,zeroOperator(N), ...
                              kappaDy         ,zeroOperator(N),-kappaDy    ,zeroOperator(N)      };
                    
                    obj.MMCsolver.addTerm(L1operatorIso(obj.alpha,4,flexA1),[i,i+1,nc+i,nc+i+1]);
                    % Add second infconv part
                    Id = diagonalOperator(obj.kappa*Ids{i}/obj.h);
                    
                    flexA2 = {Id            ,-obj.kappa*Wi,...
                              Dx            ,zeroOperator(N), ...
                              Dy            ,zeroOperator(N)      };
                    
                    obj.MMCsolver.addTerm(L1operatorIso(obj.alpha,2,flexA2),[nc+i,nc+i+1]);
                end
                % Build last infconv regularizer
                flexA1 = {kappaDx,-kappaDx, ...
                          kappaDy,-kappaDy      };
                
                obj.MMCsolver.addTerm(L1operatorIso(obj.alpha,2,flexA1),[nc,2*nc]);
                flexA2 = {Dx, ...
                          Dy      };
                
                obj.MMCsolver.addTerm(L1operatorIso(obj.alpha,1,flexA2),2*nc);
            else
                error('todo')
            end
            
            
            %%%% Set bicubic estimation as start vector
            for i = 1:nc
                ui = u_up(:,:,i);
                obj.MMCsolver.x{1,i} = ui(:);
            end
            
            if obj.verbose > 0
                disp('FlexBox object initialized');
            end
            
        end
        
        %% calculate initial velocity fields on low resolution input images and scale them up to target resolution
        function init_v0(obj)
            
            if (obj.verbose > 0)
                disp('Calculating initial velocity fields')
            end
            
            
            for j=1:obj.numFrames-1
                
                % Select correct flow direction
                if strcmp(obj.flowDirection,'forward')
                    uTmpSmall = cat(3,obj.imageSequenceSmall(:,:,j+1),obj.imageSequenceSmall(:,:,j));
                elseif strcmp(obj.flowDirection,'backward')
                    uTmpSmall = cat(3,obj.imageSequenceSmall(:,:,j),obj.imageSequenceSmall(:,:,j+1));
                elseif strcmp(obj.flowDirection,'forward-backward')
                    % alternate forward and backward:
                    if (mod(j,2)==1)
                        uTmpSmall = cat(3,obj.imageSequenceSmall(:,:,j),obj.imageSequenceSmall(:,:,j+1));
                    else
                        uTmpSmall = cat(3,obj.imageSequenceSmall(:,:,j+1),obj.imageSequenceSmall(:,:,j));
                    end
                else
                    error('Invalid flow direction.');
                end
                
                
                motionEstimator = motionEstimatorClass(uTmpSmall,1e-6,obj.beta);
                motionEstimator.verbose = 0;
                
                % Motion estimator energy parameters:
                motionEstimator.dataTerm = 'L1';
                motionEstimator.regularizerTerm ='Huber';
                motionEstimator.doGradientConstancy = 1;
                
                
                % Pyramid scheme parameters:
                motionEstimator.doGaussianSmoothing = 1;
                motionEstimator.medianFiltering = 1;
                motionEstimator.numberOfWarps = 3;
                motionEstimator.steplength = 0.8;
                
                % Run motion estimator
                motionEstimator.init;
                if strcmp(obj.comp_mode,'fast')
                    motionEstimator.steps = motionEstimator.steps(1);
                end
                motionEstimator.runPyramid;
                
                vTmp = motionEstimator.getResult;
                
                % Motion estimator velocity field is x-y instead of m-n.
                % Estimating motion on the low-res frames and interpolating seems to be
                % as efficient as interpolating first and estimating.
                obj.v(:,:,j,1) = imresize(vTmp(:,:,1,2), obj.dimsLarge) * obj.factor;
                obj.v(:,:,j,2) = imresize(vTmp(:,:,1,1), obj.dimsLarge) * obj.factor;
                
                if obj.verbose > 1
                    figure(200+j);imagesc(flowToColorV2(cat(3,obj.v(:,:,j,1),obj.v(:,:,j,2))));axis image;
                end
                drawnow;
                if obj.verbose > 0
                    disp(['Velocity field calculated for frames [',num2str(j),',',num2str(j+1),']'])
                end
                
            end
            drawnow;
        end
        
        
        
        
        %% Call u solver
        function solveU(obj)
            
            if strcmp(obj.framework,'prost')
                % Call prost framework
                prost.solve(obj.MMCsolver, obj.opts.backend, obj.opts.opts);
                
                obj.u = reshape(obj.MMCsolver.primal_vars{1,1}.val,[obj.dimsLarge, obj.numFrames]);
                if obj.kappa ~=1 && ~isnan(obj.kappa)
                    obj.w = reshape(obj.MMCsolver.primal_vars{1,2}.val,[obj.dimsLarge, obj.numFrames]);
                end
            elseif strcmp(obj.framework,'flexBox')
                % Call flexBox framework
                obj.MMCsolver.runAlgorithm;
                
                for i = 1:obj.numFrames
                    ui           = obj.MMCsolver.getPrimal(i);
                    wi           = obj.MMCsolver.getPrimal(obj.numFrames+i);
                    obj.u(:,:,i) = reshape(ui,obj.dimsLarge);
                    obj.w(:,:,i) = reshape(wi,obj.dimsLarge);
                end
            else % Vector valued flexBox
                % Call flexBox
                obj.MMCsolver.runAlgorithm;
                
                % Extract solution
                obj.u = obj.MMCsolver.getPrimal(1);
                obj.w = obj.MMCsolver.getPrimal(2);
            end
            
            
            
            % show solution
            for j=1:obj.numFrames
                if (obj.verbose > 1)
                    figure(100+j);imagesc(obj.u(:,:,j),[0,1]);axis image;colormap(gray);
                    figure(500+j);imagesc(obj.w(:,:,j));axis image;colormap(gray);colorbar;
                    figure(600+j);imagesc(obj.u(:,:,j)-obj.w(:,:,j));axis image;colormap(gray);colorbar;
                end
            end
            drawnow;
            
        end
        
        %% Recompute RGB image
        
        function recomputeRGB(obj)
            % main result
            imageSequenceUp = zeros(obj.dimsLarge(1),obj.dimsLarge(2),3,obj.numFrames);
            for i = 1:obj.numFrames
                imageSequenceUp(:,:,1,i) = obj.u(:,:,i);                                         % actually computed Y
                imageSequenceUp(:,:,2,i) = imresize(obj.imageSequenceYCbCr(:,:,2,i),obj.factor); % bicubic Cb
                imageSequenceUp(:,:,3,i) = imresize(obj.imageSequenceYCbCr(:,:,3,i),obj.factor); % bicubic Cr
                imageSequenceUp(:,:,:,i) = ycbcr2rgb(imageSequenceUp(:,:,:,i));
            end
            obj.result1 = imageSequenceUp;
            
            % u -w
            if obj.kappa ~= 1 && ~isnan(obj.kappa)
                imageSequenceUp = zeros(obj.dimsLarge(1),obj.dimsLarge(2),3,obj.numFrames);
                for i = 1:obj.numFrames
                    imageSequenceUp(:,:,1,i) = obj.u(:,:,i)-obj.w(:,:,i);                            % actually computed Y
                    
                    % correct mean values not guaranteed!
                    imageSequenceUp(:,:,2,i) = imresize(obj.imageSequenceYCbCr(:,:,2,i),obj.factor); % bicubic Cb
                    imageSequenceUp(:,:,3,i) = imresize(obj.imageSequenceYCbCr(:,:,3,i),obj.factor); % bicubic Cr
                    imageSequenceUp(:,:,:,i) = ycbcr2rgb(imageSequenceUp(:,:,:,i));
                end
                obj.result2 = imageSequenceUp;
            end
        end
        
    end
    
    methods (Static)
        function TF = checkposScalar(x)
        % as in the Matlab example :>
            if ~isscalar(x)
                error('Input is not scalar');
            elseif ~isnumeric(x)
                error('Input is not numeric');
            elseif (x <= 0)
                error('Input must be > 0');
            else
                TF = true;
            end
        end
    end
    
end



