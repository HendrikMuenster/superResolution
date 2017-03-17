classdef MultiframeMotionCoupling< handle
    %MultiframeMotionCoupling
    % Couple all frames successively and solve the super resolution problem
    %
    %
    %
    %
    %
    %
    
    %% Class properties:
    properties
        
        % Input Data
        imageSequenceSmall
        imageSequenceYCbCr
        
        % Problem Structure
        numFrames
        factor
        dimsLarge
        dimsSmall
        
        % Main variables
        u                  % image sequence
        w                  % shade sequence
        v                  % flow sequence
        k                  % blur kernel
        
        % Output Data
        result1            % u (final result) colored
        result2            % u-w (temporally stabilized frames) colored
        
        
        % Solver options
        prostU             % main solver class for u subproblem
        
        opts               % prost options
        alpha              % weights for regularizer of u
        beta               % weights for regularizer of v
        kappa              % regularizer adjustment
        h                  % time distance heuristic
        
        verbose
        
        % Sampling operator
        interpMethod       % Interpolation method
        interpKernel       % Custom interpolation kernel
        interpAA           % Interpolation AA
        
        % Extra Input
        warpOp             % A given warping operator (deactivates flow computation !)
        
        % Book-keeping
        colorflag          % Is set to 1 for color videos to enable YCbCr treatment
        flowCompute        % no flow needs to be computed if it is given by mechanics or previous runs
    end
    
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
            
            % Problem Structure defaults
            obj.factor    = 4;
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
            
            % Solver options
            %prostU                                 % main solver class for u subproblem
            % will be constructed in init_u call
            
            
            obj.opts.backend = prost.backend.pdhg('stepsize', 'boyd','tau0', 10, ...
                'sigma0', 0.1);                       % prost backend options
            obj.opts.opts = prost.options('max_iters', 15000, 'num_cback_calls', 5,...
                'verbose', true);%, ...
            %'tol_rel_primal', 1e-10, ...
            %'tol_rel_dual', 1e-10, ...
            %'tol_abs_dual', 1e-10, ...
            %'tol_abs_primal', 1e-10);                      % prost precision options
            
            
            obj.alpha = 0.01;                      % weights for regularizer of u
            obj.beta   = 0.1;                         % weights for regularizer of v
            obj.kappa  = 0.5;
            obj.h      = 1;
            
            % standard interpolation parameters
            obj.interpMethod  = 'bicubic-0';     % Interpolation method
            obj.interpKernel  = [];               % Custom interpolation kernel
            obj.interpAA      = false  ;          % Interpolation AA
            
        end
        
        %% collect inits and validate
        function init(obj)
            
            if obj.verbose > 0
                disp('Initialization...')
            end
            
            % validate and set variables
            
            % update sizes for non-standard factor
            obj.dimsLarge = obj.factor*obj.dimsSmall;
            
            % initialize
            if obj.flowCompute
                obj.v = zeros([obj.dimsLarge,obj.numFrames-1,2]);
            end
            
            % Call actual initialization methods
            if obj.flowCompute
                obj.init_v0;
            end
            
            % Init u
            obj.init_u;
            
            if (obj.verbose > 0)
                disp('Initialization finished');
            end
        end
        
        %% create prost object for u
        function init_u(obj)
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
            % short some notation
            temp = obj.dimsSmall;
            ny=temp(1); nx=temp(2);
            temp = obj.dimsLarge;
            Ny = temp(1);
            Nx = temp(2);
            nc =  obj.numFrames;
            
            % Call warp operator constructor
            if isscalar(obj.warpOp)
                warpingOp = constructWarpFF(obj.v,'F-I');
                obj.warpOp = warpingOp;
                if obj.verbose > 0
                    disp('Warp operator constructed');
                end
            else
                warpingOp = obj.warpOp;
                if obj.verbose > 0
                    disp('Warp operator loaded');
                end
            end
            
            % preweight, based on warping accuracy
            u_up = zeros(Ny,Nx,nc);
            for i = 1:nc
                u_up(:,:,i) = imresize(obj.imageSequenceSmall(:,:,i),obj.factor);
            end
            warpPix = sum(abs(warpingOp*u_up(:)))/numel(u_up);
            disp(['Warp energy on bicubic (per pixel): ',num2str(warpPix)]);
            GradOp = spmat_gradient2d(Nx, Ny, nc);
            gradPix = sum(abs(GradOp*u_up(:)))/numel(u_up);
            disp(['Gradient energy on bicubic (per pixel): ',num2str(gradPix)])
            obj.h  = warpPix/gradPix;
            
            %%% Check warp accuracy visually
            if obj.verbose > 1
                u_mv = reshape(warpingOp*u_up(:),Ny,Nx,nc);
                for i = 1:nc-1
                    figure(i), imagesc(u_mv(:,:,i)),title('Visual test of warp accuracy')
                end
            end
            
            %%%% initialize prost variables and operators
            
            % Primal- Core Variable:
            u_vec = prost.variable(Nx*Ny*nc);
            if isscalar(obj.u)
                u_vec.val = u_up(:);
            else
                u_vec.val = obj.u(:);
            end
            w_vec = prost.variable(Nx*Ny*nc);
            w_vec.val = u_up(:);
            p = prost.variable(nx*ny*nc);
            g1 = prost.variable(3*Nx*Ny*nc);
            g2 = prost.variable(3*Nx*Ny*nc);
            
            % Connect variables to primal-dual problem
            if obj.kappa ~= 1 && ~isnan(obj.kappa)
                
                % Generate min-max problem structure
                obj.prostU = prost.min_max_problem( {u_vec,w_vec}, {p,g1,g2} );
                
                
                % Initialize full downsampling implictely as kron(id(numFrames,D)
                obj.prostU.add_dual_pair(u_vec, p, prost.block.id_kron_sparse(dsOp,obj.numFrames));
                
                % Add inf-conv duals
                obj.prostU.add_dual_pair(u_vec, g1, prost.block.sparse([warpingOp/obj.h;obj.kappa*spmat_gradient2d(Nx, Ny,nc)]));
                obj.prostU.add_dual_pair(w_vec, g1, prost.block.sparse(-[warpingOp/obj.h;obj.kappa*spmat_gradient2d(Nx, Ny,nc)]));
                obj.prostU.add_dual_pair(w_vec, g2, prost.block.sparse([obj.kappa*warpingOp/obj.h;spmat_gradient2d(Nx, Ny,nc)]));
                
                % Add functions
                obj.prostU.add_function(u_vec,prost.function.sum_1d('ind_box01',1,0,1,0,0));
                obj.prostU.add_function(p, prost.function.sum_1d('ind_box01', 0.5, -0.5, 1, obj.imageSequenceSmall(:), 0)); %l^1
                obj.prostU.add_function(g1, prost.function.sum_norm2(3, false, 'ind_leq0', 1/obj.alpha, 1, 1, 0, 0)); %l^{2,1}
                obj.prostU.add_function(g2, prost.function.sum_norm2(3, false, 'ind_leq0', 1/obj.alpha, 1, 1, 0, 0)); %l^{2,1}
                
            elseif obj.kappa == 1    % Simplify if no infimal convolution is necessary
                
                % Generate min-max problem structure
                obj.prostU = prost.min_max_problem( {u_vec}, {p,g1} );
                
                % Initialize full downsampling implictely as kron(id(numFrames,D) if possible
                obj.prostU.add_dual_pair(u_vec, p, prost.block.id_kron_sparse(dsOp,obj.numFrames));
                
                % Add regularizer dual
                obj.prostU.add_dual_pair(u_vec, g1, prost.block.sparse([warpingOp/obj.h;spmat_gradient2d(Nx, Ny,nc)]));
                % add functions
                obj.prostU.add_function(u_vec,prost.function.sum_1d('ind_box01',1,0,1,0,0));
                obj.prostU.add_function(p, prost.function.sum_1d('ind_box01', 0.5, -0.5, 1, obj.imageSequenceSmall(:), 0)); %l^1
                obj.prostU.add_function(g1, prost.function.sum_norm2(3, false, 'ind_leq0', 1/obj.alpha/2, 1, 1, 0, 0)); %l^{2,1}
                
            elseif isnan(obj.kappa)    % Do additve TV and Warp
                
                g1 = prost.variable(Nx*Ny*nc);
                g2 = prost.variable(2*Nx*Ny*nc);
                
                % Generate min-max problem structure
                obj.prostU = prost.min_max_problem( {u_vec}, {p,g1,g2} );
                
                % Initialize full downsampling implictely as kron(id(numFrames,D)
                obj.prostU.add_dual_pair(u_vec, p, prost.block.id_kron_sparse(dsOp,obj.numFrames));
                
                % Add regularizer dual
                obj.prostU.add_dual_pair(u_vec, g1, prost.block.sparse(warpingOp/obj.h)); %*tdist
                obj.prostU.add_dual_pair(u_vec, g2, prost.block.gradient2d(Nx, Ny,nc,false));
                % add functions
                obj.prostU.add_function(u_vec,prost.function.sum_1d('ind_box01',1,0,1,0,0));
                obj.prostU.add_function(p, prost.function.sum_1d('ind_box01', 0.5, -0.5, 1, obj.imageSequenceSmall(:), 0)); %l^1
                obj.prostU.add_function(g1, prost.function.sum_1d('ind_box01', 1/(2*obj.alpha), -0.5, 1, 0, 0)); %l^1
                obj.prostU.add_function(g2, prost.function.sum_norm2(2, false, 'ind_leq0', 1/obj.alpha, 1, 1, 0, 0)); %l^{2,1}
                
                
            end
            
        end
        %% calculate initial velocity fields on low resolution input images and scale them up to target resolution
        function init_v0(obj)
            
            if (obj.verbose > 0)
                disp('Calculating initial velocity fields')
            end
            
            
            for j=1:obj.numFrames-1
                % Select correct flow direction
                uTmpSmall = cat(3,obj.imageSequenceSmall(:,:,j+1),obj.imageSequenceSmall(:,:,j));
                
                motionEstimator = motionEstimatorClass(uTmpSmall,1e-6,obj.beta);
                motionEstimator.verbose = 0;
                
                % motion estimator energy parameters:
                motionEstimator.dataTerm = 'L1';
                motionEstimator.regularizerTerm ='Huber';
                motionEstimator.doGradientConstancy = 1;
                
                
                % pyramid scheme:
                motionEstimator.doGaussianSmoothing = 1;
                motionEstimator.medianFiltering = 1;
                motionEstimator.numberOfWarps = 5;
                motionEstimator.steplength = 0.9;
                
                % run motion estimator
                motionEstimator.init;
                motionEstimator.runPyramid;
                
                vTmp = motionEstimator.getResult;
                
                % motion estimator velocity field is x-y instead of m-n
                % estimating motion on the low-res frames and interpolating seems to be
                % as efficient as interpolating first and estimating
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
        
        
        
        
        %% Call u solver in prost
        function solveU(obj)
            
            %call prost framework
            %tic
            prost.solve(obj.prostU, obj.opts.backend, obj.opts.opts);
            %toc
            obj.u = reshape(obj.prostU.primal_vars{1,1}.val,[obj.dimsLarge, obj.numFrames]);
            if obj.kappa ~=1 && ~isnan(obj.kappa)
                obj.w = reshape(obj.prostU.primal_vars{1,2}.val,[obj.dimsLarge, obj.numFrames]);
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
            imageSequenceUp = zeros(obj.dimsLarge(1),obj.dimsLarge(2),3,obj.numFrames);
            for i = 1:obj.numFrames
                imageSequenceUp(:,:,1,i) = obj.u(:,:,i);                                         % actually computed Y
                imageSequenceUp(:,:,2,i) = imresize(obj.imageSequenceYCbCr(:,:,2,i),obj.factor); % bicubic Cb
                imageSequenceUp(:,:,3,i) = imresize(obj.imageSequenceYCbCr(:,:,3,i),obj.factor); % bicubic Cr
                imageSequenceUp(:,:,:,i) = ycbcr2rgb(imageSequenceUp(:,:,:,i));
            end
            obj.result1 = imageSequenceUp;
            if obj.kappa ~= 1 && ~isnan(obj.kappa)
                imageSequenceUp = zeros(obj.dimsLarge(1),obj.dimsLarge(2),3,obj.numFrames);
                for i = 1:obj.numFrames
                    imageSequenceUp(:,:,1,i) = obj.u(:,:,i)-obj.w(:,:,i);                            % actually computed Y
                    imageSequenceUp(:,:,2,i) = imresize(obj.imageSequenceYCbCr(:,:,2,i),obj.factor); % bicubic Cb
                    imageSequenceUp(:,:,3,i) = imresize(obj.imageSequenceYCbCr(:,:,3,i),obj.factor); % bicubic Cr
                    imageSequenceUp(:,:,:,i) = ycbcr2rgb(imageSequenceUp(:,:,:,i));
                end
                obj.result2 = imageSequenceUp;
            end
        end
        
        
        %% Run the algorithm
        function run(obj)
            
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
    
end



