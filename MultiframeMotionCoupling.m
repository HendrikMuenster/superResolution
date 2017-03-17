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
        framework          % Choice of framework prost/flexBox
        MMCsolver          % main solver object for u subproblem
        
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
            %MMCsolver                                 % main solver class for u subproblem
            % will be constructed in init_u_... call
            
            % default SR solver:
            obj.framework = 'prost';
            
            
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
            
            % Init u with desired framework
            if strcmp(obj.framework,'prost')
                obj.init_u_prost;
            elseif strcmp(obj.framework,'flexBox')
                obj.init_u_flexBox;
            else
                error('Invalid framework string.');
            end
            
            
            if (obj.verbose > 0)
                disp('Initialization finished');
            end
        end
        
        %% create prost object for u
        function init_u_prost(obj)
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
                obj.MMCsolver = prost.min_max_problem( {u_vec,w_vec}, {p,g1,g2} );
                
                
                % Initialize full downsampling implictely as kron(id(numFrames,D)
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
                
                % Initialize full downsampling implictely as kron(id(numFrames,D) if possible
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
                obj.MMCsolver.add_dual_pair(u_vec, g1, prost.block.sparse(warpingOp/obj.h)); %*tdist
                obj.MMCsolver.add_dual_pair(u_vec, g2, prost.block.gradient2d(Nx, Ny,nc,false));
                % add functions
                obj.MMCsolver.add_function(u_vec,prost.function.sum_1d('ind_box01',1,0,1,0,0));
                obj.MMCsolver.add_function(p, prost.function.sum_1d('ind_box01', 0.5, -0.5, 1, obj.imageSequenceSmall(:), 0)); %l^1
                obj.MMCsolver.add_function(g1, prost.function.sum_1d('ind_box01', 1/(2*obj.alpha), -0.5, 1, 0, 0)); %l^1
                obj.MMCsolver.add_function(g2, prost.function.sum_norm2(2, false, 'ind_leq0', 1/obj.alpha, 1, 1, 0, 0)); %l^{2,1}
                
                
            end
            
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
            if isnan(obj.kappa)
                error('Additive warp not implemented yet for flexBox.');
            end
           
            %%%% Create operators
            
            % Build downsampling matrix
            dsOp = samplingOperator(obj.dimsLarge,obj.dimsSmall,obj.interpMethod,obj.interpKernel,obj.interpAA);
            
            % Add blur matrix
            if ~isempty(obj.k) || isscalar(obj.k)
                dsOp = dsOp*RepConvMtx(obj.k,obj.dimsLarge);
            end
            
            if obj.verbose > 0
                disp('Downsampling operator constructed');
            end
            
            % Construct warp operators and reduced identities
            for i = 1:nc-1
                               
                singleField = squeeze(obj.v(:,:,i,:));
                Wi   = warpingOperator(obj.dimsLarge,singleField);
                Id   = speye(N);
                
                % Remove out-of range warps
                marker = sum(abs(Wi),2) == 0;
                Wi(marker > 0,:) = 0;
                Id(marker > 0,:) = 0; %#ok<SPRIX>

                warps{i} = Wi;  %#ok<AGROW>
                Ids{i}   = Id;  %#ok<AGROW>
            end            
           
            % Build gradient matrices
            D = spmat_gradient2d(Nx, Ny,1);
            Dx = D(1:Nx*Ny,:); Dy = D(Nx*Ny+1:end,:);
            
            
            %%%% Compute adaptive parameter h
            
            % Compute bicubic estimate 
            u_up = zeros([obj.dimsLarge,nc]);
            for i = 1:nc
                u_up(:,:,i) = imresize(obj.imageSequenceSmall(:,:,i),obj.factor,'bicubic');
            end
            % Accumulate estimates 
            Du = zeros(nc,1);Wu = zeros(nc,1);
            for i = 1:nc-1
                ui = u_up(:,:,i);
                uj = u_up(:,:,i+1);
                Du(i) = sum(abs(Dx*ui(:))+abs(Dy*ui(:)));
                Wu(i) = sum(abs(-Ids{i}*uj(:)+warps{i}*ui(:)));
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
            obj.MMCsolver.params.tol = 1e-6;
            obj.MMCsolver.params.tryCPP = 1;
            obj.MMCsolver.params.verbose = 1;
            obj.MMCsolver.params.maxIt = 10000;
            
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
            
            % Build nc-1 infconv regularizers
            for i = 1:nc-1
                % Get warp and identity
                Wi = warps{i}/obj.h;
                Id = Ids{i}/obj.h;
                
                % Add first infconv part
                flexA1 = {obj.kappa*Dx    ,zeroOperator(N),-obj.kappa*Dx,   zeroOperator(N), ...
                    obj.kappa*Dy    ,zeroOperator(N),-obj.kappa*Dy,   zeroOperator(N), ...
                    Wi              ,-Id            ,-Wi          ,   Id                   };
                
                obj.MMCsolver.addTerm(L1operatorIso(obj.alpha,4,flexA1),[i,i+1,nc+i,nc+i+1]);
                % Add second infconv part
                flexA2 = {Dx            ,zeroOperator(N), ...
                    Dy            ,zeroOperator(N), ...
                    obj.kappa*Wi  ,-obj.kappa*Id       };
                
                obj.MMCsolver.addTerm(L1operatorIso(obj.alpha,2,flexA2),[nc+i,nc+i+1]); 
            end
            % Build last infconv regularizer
            flexA1 = {obj.kappa*Dx,-obj.kappa*Dx, ...
                      obj.kappa*Dy,-obj.kappa*Dy      };
            
            obj.MMCsolver.addTerm(L1operatorIso(obj.alpha,2,flexA1),[nc,2*nc]);
            flexA2 = {Dx, ...
                      Dy      };
            
            obj.MMCsolver.addTerm(L1operatorIso(obj.alpha,1,flexA2),2*nc);
                    
                  
                    
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
        
        
        
        
        %% Call u solver
        function solveU(obj)
            
            if strcmp(obj.framework,'prost')
                % Call prost framework
                prost.solve(obj.MMCsolver, obj.opts.backend, obj.opts.opts);
                
                obj.u = reshape(obj.MMCsolver.primal_vars{1,1}.val,[obj.dimsLarge, obj.numFrames]);
                if obj.kappa ~=1 && ~isnan(obj.kappa)
                    obj.w = reshape(obj.MMCsolver.primal_vars{1,2}.val,[obj.dimsLarge, obj.numFrames]);
                end
            else
                % Call flexBox framework
                obj.MMCsolver.runAlgorithm;
                
                for i = 1:obj.numFrames
                    ui           = obj.MMCsolver.getPrimal(i);
                    wi           = obj.MMCsolver.getPrimal(obj.numFrames+i);
                    obj.u(:,:,i) = reshape(ui,obj.dimsLarge);
                    obj.w(:,:,i) = reshape(wi,obj.dimsLarge);
                end
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



