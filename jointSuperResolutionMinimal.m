classdef jointSuperResolutionMinimal< handle
    %JOINTSUPERRESOLUTION
    % solve a joint super-resolution and optical flow problem, estimating
    % the optical flow 'v' of an unknown high resolution video 'u' from a low
    % resolution video 'imageSequenceSmall' while jointly recovering the
    % high resolution video and the motion/downsample blur kernels k_i
    %% Class properties:
    properties
        
        % Input Data
        imageSequenceSmall
        imageSequenceYCbCr
        datasetName
        
        % Problem Structure
        numFrames
        factor
        dimsLarge
        dimsSmall
        
        % Main variables
        u                  % image sequence
        w                  % shade sequence
        v                  % flow sequence
        k                  % blur sequence
        
        % Output Data
        result1            % u colored
        result2            % u-w colored
        
        % Ground Truth variables
        gtU
        gtV
        
        % Solver options
        prostU             % main solver class for u subproblem
        mainV              % main solver class for v subproblem
        prostK             % main solver class for k subproblem
        
        opts               % prost options
        alpha1             % weights for regularizer of u
        alpha2             % weights for regularizer of u
        beta               % weights for regularizer of v
        kappa              % regularizer adjustment
        tdist              % time distance heuristic
        
        % kernel estimation (deprecated)
        kOpts
        kernelsize

        verbose
        numMainIt;         % Number of total iterations
        interpMethod       % Interpolation method
        interpKernel       % Custom interpolation kernel
        interpAA           % Interpolation AA
        
        % Extra Input
        warpOp             % A given warping operator (deactivates flow computation !)
        
        % Book-keeping
        currentIt          % notifies solvers of current iteration
        colorflag          % Is set to 1 for color videos to enable YCbCr treatment
        testCase
        flowCompute        % no flow needs to be computed if it is given by mechanics or previous runs
    end
    
    methods
        %% Initialize object
        function obj = jointSuperResolutionMinimal(imageSequenceSmall,varargin)
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
            
            if (exist('datasetName','var'))
                obj.datasetName = datasetName; %#ok<*CPROPLC>
            else
                obj.datasetName = '';
            end
            
            if (exist('testCase','var'))
                obj.testCase = testCase; %#ok<*CPROPLC>
            else
                obj.testCase = 'FB';
            end
            
            % Problem Structure
            
            obj.factor    = 4;
            obj.dimsSmall = size(imageSequenceSmall);
            obj.dimsSmall = obj.dimsSmall(1:2);
            obj.dimsLarge = obj.factor*obj.dimsSmall;  % will be overwritten later if factor is non-standard
            
            
            % Main variables
            %u                  % image sequence
            %v                  % flow sequence
            %k                  % blur sequence
            % will be initialized in init call
            
            % Ground Truth variables
            if (exist('gtU','var'))
                obj.gtU = gtU;
            else
                obj.gtU = 0;
            end
            
            if (exist('gtV','var'))
                obj.gtV = gtV;
            else
                obj.gtV = 0;
            end
            
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
            %mainV                                  % main solver class for v subproblem
            %prostK                                 % main solver class for k subproblem
            % will be constructed in init_u, init_v and init_k calls
            
                                                             
            obj.opts.backend = prost.backend.pdhg('stepsize', 'boyd','tau0', 100, ...
                                                  'sigma0', 0.01);  % prost backend options
            obj.opts.opts = prost.options('max_iters', 12000, 'num_cback_calls', 5, 'verbose', true);  % prost structure options
            obj.numMainIt = 1;                      % Number of total outer iterations
            obj.opts.nsize = 0;                     % radius of local boundary, choose 0 for no local boundaries
            obj.opts.offset = 0.1;                  % offset of local boundary
            
            obj.alpha1 = 0.01;                      % weights for regularizer of u
            obj.alpha2 = 0.01;                      % weights for regularizer of u
            obj.beta   = 0.1;                         % weights for regularizer of v
            obj.kappa  = 0.5;
            obj.tdist  = 1;
            
            
            % kernel estimation
            obj.kOpts.backend = prost.backend.pdhg(...
                'tau0', 100, ...
                'sigma0', 0.01, ...
                'stepsize', 'boyd');                % prost backend options for k
            obj.kOpts.opts = prost.options('max_iters', 2500, ...
                'num_cback_calls', 5, ...
                'verbose', true);                   % prost structure options for k
            obj.kernelsize = 11;                     % standard kernel size
            obj.kOpts.delta = 0.05;                 % l2 penalty on k
            
            %             % starting input kernel
            %             sigmaval = 1/4*sqrt(obj.factor^2-1); %Innerhofer/Pock:
            %             obj.k0 =   exp(-(-3:3).^2 / sigmaval)'*exp(-(-3:3).^2 / sigmaval);
            %             obj.k0 = obj.k0/sum(sum(obj.k0));
            %
            %             obj.lowPass = 'Gaussian';
            %             obj.lowPassParam = sigmaval;
            
            % standard interpolation parameters
            obj.interpMethod  = 'bicubic-0';     % Interpolation method
            obj.interpKernel  = [];               % Custom interpolation kernel
            obj.interpAA      = false  ;          % Interpolation AA
            
            
            % Book-keeping
            obj.currentIt = 1;         % notifies solvers of current iteration
            
            
        end
        
        %% create prost object for u
        function init_u(obj)
            %%%% Create operators
            
            % Call sampling function to construct matrix representation
            dsOp = samplingOperator(obj.dimsLarge,obj.dimsSmall,obj.interpMethod,obj.interpKernel,obj.interpAA);
            %dsOp = superpixelOperator(obj.dimsSmall,obj.factor).matrix;
            
            % Use blur kernel:
            if ~isempty(obj.k) || isscalar(obj.k)
                dsOp = dsOp*RepConvMtx(obj.k,obj.dimsLarge);
            end
            
            % initialize block format
            if obj.numMainIt   ~= 1
                downsamplingOp = kron(speye(obj.numFrames),dsOp);
            else
                downsamplingOp_kron = dsOp;
            end
            if obj.verbose > 0
                disp('downsampling operator constructed');
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
                if strcmp(obj.testCase,'FB')
                    warpingOp = constructWarpFB(obj.v);
                elseif strcmp(obj.testCase,'FMB')
                    warpingOp = constructWarpFMB(obj.v);
                elseif strcmp(obj.testCase,'F-I')
                    warpingOp = constructWarpFF(obj.v,'F-I');
                elseif strcmp(obj.testCase,'I-B')
                    warpingOp = constructWarpFF(obj.v,'I-F');
                end
                obj.warpOp = warpingOp;
            else
                warpingOp = obj.warpOp;
                
            end
            
            if obj.verbose > 0
                disp('warp operator constructed');
            end
            
            % preweight, based on warping accuracy
            u_up = zeros(Ny,Nx,nc);
            for i = 1:nc
                u_up(:,:,i) = imresize(obj.imageSequenceSmall(:,:,i),obj.factor);
            end
            warpPix = sum(abs(warpingOp*u_up(:)))/numel(u_up);
            disp(['warp energy on bicubic (per pixel): ',num2str(warpPix)]);
            GradOp = spmat_gradient2d(Nx, Ny, nc);
            gradPix = sum(abs(GradOp*u_up(:)))/numel(u_up);
            disp(['gradient energy on bicubic (per pixel): ',num2str(gradPix)])
            obj.tdist  = gradPix/warpPix;
            
            %%%% initialize prost variables and operators
            
            % Primal- Core Variable:
            u_vec = prost.variable(Nx*Ny*nc);
            u_vec.val = u_up(:);
            w_vec = prost.variable(Nx*Ny*nc);
            p = prost.variable(nx*ny*nc);
            g1 = prost.variable(3*Nx*Ny*nc);
            g2 = prost.variable(3*Nx*Ny*nc);
            
            % Connect variables to primal-dual problem
            if obj.kappa ~= 1 && ~isnan(obj.kappa)
                
                % Generate min-max problem structure
                obj.prostU = prost.min_max_problem( {u_vec,w_vec}, {p,g1,g2} );
                
                
                % Initialize full downsampling implictely as kron(id(numFrames,D) if possible
                if obj.numMainIt ~= 1
                    obj.prostU.add_dual_pair(u_vec, p, prost.block.sparse(downsamplingOp));
                else
                    obj.prostU.add_dual_pair(u_vec, p, prost.block.id_kron_sparse(downsamplingOp_kron,obj.numFrames));
                end
                
                % Add inf-conv duals
                obj.prostU.add_dual_pair(u_vec, g1, prost.block.sparse([warpingOp*obj.tdist;obj.kappa*spmat_gradient2d(Nx, Ny,nc)]));
                obj.prostU.add_dual_pair(w_vec, g1, prost.block.sparse(-[warpingOp*obj.tdist;obj.kappa*spmat_gradient2d(Nx, Ny,nc)]));
                obj.prostU.add_dual_pair(w_vec, g2, prost.block.sparse([obj.kappa*warpingOp*obj.tdist;spmat_gradient2d(Nx, Ny,nc)]));
                % Add functions

                obj.prostU.add_function(u_vec,prost.function.sum_1d('ind_box01',1,0,1,0,0));
                obj.prostU.add_function(p, prost.function.sum_1d('ind_box01', 0.5, -0.5, 1, obj.imageSequenceSmall(:), 0)); %l^1
                obj.prostU.add_function(g1, prost.function.sum_norm2(3, false, 'ind_leq0', 1/obj.alpha1, 1, 1, 0, 0)); %l^{2,1}
                obj.prostU.add_function(g2, prost.function.sum_norm2(3, false, 'ind_leq0', 1/obj.alpha2, 1, 1, 0, 0)); %l^{2,1}
                
            elseif obj.kappa == 1    % Simplify if no infimal convolution is necessary
                
                % Generate min-max problem structure
                obj.prostU = prost.min_max_problem( {u_vec}, {p,g1} );
                
                % Initialize full downsampling implictely as kron(id(numFrames,D) if possible
                if obj.numMainIt ~= 1
                    obj.prostU.add_dual_pair(u_vec, p, prost.block.sparse(downsamplingOp));
                else
                    obj.prostU.add_dual_pair(u_vec, p, prost.block.id_kron_sparse(downsamplingOp_kron,obj.numFrames));
                end
                
                % Add regularizer dual
                obj.prostU.add_dual_pair(u_vec, g1, prost.block.sparse([warpingOp*obj.tdist;spmat_gradient2d(Nx, Ny,nc)]));
                % add functions
                if obj.opts.nsize ~= 0
                    obj.prostU.add_function(u_vec,prost.function.sum_1d('ind_box01',1./(lmax-lmin),lmin./(lmax-lmin),1,0,0));
                else
                    obj.prostU.add_function(u_vec,prost.function.sum_1d('ind_box01',1,0,1,0,0));
                end
                obj.prostU.add_function(p, prost.function.sum_1d('ind_box01', 0.5, -0.5, 1, obj.imageSequenceSmall(:), 0)); %l^1
                obj.prostU.add_function(g1, prost.function.sum_norm2(3, false, 'ind_leq0', 1/obj.alpha1/2, 1, 1, 0, 0)); %l^{2,1}
            
            elseif isnan(obj.kappa)    % Do standard TV
                
                g1 = prost.variable(Nx*Ny*nc);
                g2 = prost.variable(2*Nx*Ny*nc);
                
                % Generate min-max problem structure
                obj.prostU = prost.min_max_problem( {u_vec}, {p,g1,g2} );
                
                % Initialize full downsampling implictely as kron(id(numFrames,D) if possible
                if obj.numMainIt ~= 1
                    obj.prostU.add_dual_pair(u_vec, p, prost.block.sparse(downsamplingOp));
                else
                    obj.prostU.add_dual_pair(u_vec, p, prost.block.id_kron_sparse(downsamplingOp_kron,obj.numFrames));
                end
                
                % Add regularizer dual
                obj.prostU.add_dual_pair(u_vec, g1, prost.block.sparse(warpingOp*obj.tdist)); %*tdist
                obj.prostU.add_dual_pair(u_vec, g2, prost.block.gradient2d(Nx, Ny,nc,false));
                % add functions
                obj.prostU.add_function(u_vec,prost.function.sum_1d('ind_box01',1,0,1,0,0));
                obj.prostU.add_function(p, prost.function.sum_1d('ind_box01', 0.5, -0.5, 1, obj.imageSequenceSmall(:), 0)); %l^1
                obj.prostU.add_function(g1, prost.function.sum_1d('ind_box01', 1/(2*obj.alpha1), -0.5, 1, 0, 0)); %l^1
                obj.prostU.add_function(g2, prost.function.sum_norm2(2, false, 'ind_leq0', 1/obj.alpha2, 1, 1, 0, 0)); %l^{2,1}
                
                
            end
            
        end
        %% calculate initial velocity fields on low resolution input images and scale them up to target resolution
        function init_v0(obj)
            
            if (obj.verbose > 0)
                disp('Calculating initial velocity fields')
            end
            
            
            for j=1:obj.numFrames-1
                
                
                % Select correct flow direction
                if strcmp(obj.testCase,'FB')
                    if (mod(j,2)==1)%calculate backward flow: v s.t. u_2(x+v)=u_1(x)
                        uTmpSmall = cat(3,obj.imageSequenceSmall(:,:,j),obj.imageSequenceSmall(:,:,j+1));
                    else            %calculate backward flow: v s.t. u_1(x+v)=u_2(x)
                        uTmpSmall = cat(3,obj.imageSequenceSmall(:,:,j+1),obj.imageSequenceSmall(:,:,j));
                    end
                elseif strcmp(obj.testCase,'FMB')
                    uTmpSmall = cat(3,obj.imageSequenceSmall(:,:,j),obj.imageSequenceSmall(:,:,j+1));
                elseif strcmp(obj.testCase,'F-I')
                    uTmpSmall = cat(3,obj.imageSequenceSmall(:,:,j+1),obj.imageSequenceSmall(:,:,j));
                elseif strcmp(obj.testCase,'I-B')
                    uTmpSmall = cat(3,obj.imageSequenceSmall(:,:,j),obj.imageSequenceSmall(:,:,j+1));
                end
                
                
                motionEstimatorLow = motionEstimatorClass(uTmpSmall,1e-6,obj.beta,'doGradientConstancy',1);
                motionEstimatorLow.verbose = 0;
                motionEstimatorLow.init;
                motionEstimatorLow.runPyramid;
                
                vTmp = motionEstimatorLow.getResult;
                
                obj.v(:,:,j,1) = imresize(vTmp(:,:,1,1), obj.dimsLarge) * obj.factor;
                obj.v(:,:,j,2) = imresize(vTmp(:,:,1,2), obj.dimsLarge) * obj.factor;
                
                if obj.verbose > 1
                    figure(200+j);imagesc(flowToColorV2(cat(3,obj.v(:,:,j,1),obj.v(:,:,j,2))));axis image;
                end
                drawnow;
                if obj.verbose > 0
                    disp(['Initial velocity field calculated for frames [',num2str(j),',',num2str(j+1),']'])
                end
                
            end
            drawnow;
        end
        
        %% create blur estimator object in prost
        function init_k(obj)
            %  only to initiate prost block and validate sizes
            
            % create downsampling operator in block format
            downsamplingOp = superpixelOperator(obj.dimsSmall,obj.factor).matrix;
            downsamplingOp = kron(speye(obj.numFrames),downsamplingOp);
            % create DU operator, s.t. DUk -> f <- DKu
            downsampleImageOp = downsamplingOp*imageStackOp(obj.u,obj.kernelsize); % this is kind of a waste of memory ;)
            % short some notation
            temp = obj.dimsSmall;
            ny=temp(1); nx=temp(2);
            nc =  obj.numFrames;
            kx = obj.kernelsize; ky = obj.kernelsize;
            
            %%% Initialize prost variables
            % Primal variable - vector of vectorized kernels of each frame
            k_vec = prost.variable(kx*ky*nc);
            % Dual variables
            p     = prost.variable(nx*ny*nc);   % dual variable for l1 data
            g1    = prost.variable(kx*ky*nc*2); % dual variable for Tikh reg
            
            % Connect variables
            obj.prostK = prost.min_max_problem( {k_vec}, {p,g1} );
            obj.prostK.add_dual_pair(k_vec,p,prost.block.sparse(downsampleImageOp));
            obj.prostK.add_dual_pair(k_vec,g1,prost.block.gradient2d(kx,ky,nc,false));
            
            % Add functions
            obj.prostK.add_function(k_vec,prost.function.sum_ind_simplex(kx*ky,true)); % this is a simplex for each k_i
            obj.prostK.add_function(p,prost.function.sum_1d('ind_box01',0.5,-0.5,1,obj.imageSequenceSmall(:),0)); % l1 data
            obj.prostK.add_function(g1,prost.function.sum_1d('square', 1, 0, 1/obj.kOpts.delta, 0, 0)); % Tikh penalty
        end
        
        %% collect inits and validate
        function init(obj)
            
            if obj.verbose > 0
                disp('Initialization...')
            end
            
            % validate and set variables
            
            % update sizes for non-standard factor
            obj.dimsLarge = obj.factor*obj.dimsSmall;
            
            % initialize u,w,v,k
            obj.u = zeros([obj.dimsLarge,obj.numFrames]);
            obj.w = zeros([obj.dimsLarge,obj.numFrames]);
            if obj.flowCompute
                obj.v = zeros([obj.dimsLarge,obj.numFrames-1,2]);
            end
            %obj.k = zeros(obj.kernelsize^2*obj.numFrames,1);
            
            
            
            % Call actual initialization methods
            if obj.flowCompute
                obj.init_v0;
            end
            
            
           
            
            % Init u
            obj.init_u;
            
            
            % Init a solver for k only if necessary
            if obj.numMainIt > 1
                obj.init_k;
            end
            
            if (obj.verbose > 0)
                disp('Initialization of u, v and k solvers finished');
            end
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
        
        %% Construct u solver update
        function updateProstU(obj)
            % Build new warping operators
            % Adjust parameters
            
            % shorting some notation
            temp = obj.dimsLarge;
            Ny = temp(1);
            Nx = temp(2);
            nc =  obj.numFrames;
            
            % Call warp operator constructor
            %warpingOp = constructWarpFB(obj.v);
            warpingOp = constructWarpFMB(obj.v);
            
            disp(['warp energy on iteration: ',num2str(sum(abs(warpingOp*obj.u(:)))/numel(obj.u))])
            
            
            if obj.verbose > 0
                disp('warping operator updated');
            end
            
            
            % update warping operator in prost
            obj.prostU.data.linop{1,2}{1,4}{1} = [warpingOp*obj.tdist;obj.kappa*spmat_gradient2d(Nx, Ny,nc)];
            obj.prostU.data.linop{1,3}{1,4}{1} = -[warpingOp*obj.tdist;obj.kappa*spmat_gradient2d(Nx, Ny,nc)];
            obj.prostU.data.linop{1,4}{1,4}{1} = [obj.kappa*warpingOp*obj.tdist;spmat_gradient2d(Nx, Ny,nc)];
            
        end
        
        %% Call v solver from motionEstimationGUI (which wraps flexBox)
        function solveV(obj)
            
            % Build input into solveV
            uw = obj.u-obj.w; % this is a slight simplification, as +|Grad(u-w)| is not represented
            
            for j=1:obj.numFrames-1
                if (mod(j,2)==1)    %calculate backward flow: v s.t. u_2(x+v)=u_1(x)
                    uTmp = cat(3,uw(:,:,j),uw(:,:,j+1));
                else                %calculate backward flow: v s.t. u_1(x+v)=u_2(x)
                    uTmp = cat(3,uw(:,:,j+1),uw(:,:,j));
                end
                
                motionEstimator = motionEstimatorClass(uTmp,1e-6,obj.beta,'doGradientConstancy',1);
                motionEstimator.regularizerTerm = 'huber';
                motionEstimator.verbose = 0;
                motionEstimator.init;
                motionEstimator.runPyramid;
                
                vTmp =  motionEstimator.getResult;
                obj.v(:,:,j,1) = vTmp(:,:,1,1);
                obj.v(:,:,j,2) = vTmp(:,:,1,2);
                
                if obj.verbose > 0
                    disp(['Updated high-res velocity field calculated for frames [',num2str(j),',',num2str(j+1),']']);
                end
                
            end
            
            if obj.verbose > 1
                for j = 1:obj.numFrames-1
                    figure(200+j);imagesc(flowToColorV2(cat(3,obj.v(:,:,j,1),obj.v(:,:,j,2))));axis image;
                end
            end
            
            
        end
        
        %% Call appropriate solver for blur problem
        function solveBlur(obj)
            
            % Update DU Operator
            downsamplingOp = kron(speye(obj.numFrames),superpixelOperator(obj.dimsSmall,obj.factor).matrix);
            downsampleImageOp = downsamplingOp*imageStackOp(obj.u,obj.kernelsize); % this is kind of a waste of memory ;)
            
            % update DU operator in prost
            obj.prostK.data.linop{1,1}{1,4}{1} = downsampleImageOp;
            
            % call prost framework
            tic
            prost.solve(obj.prostK,obj.kOpts.backend,obj.kOpts.opts);
            toc
            obj.k  = reshape(obj.prostK.primal_vars{1,1}.val,obj.kernelsize,obj.kernelsize,obj.numFrames);
            if obj.verbose > 1
                % draw subplots of blur kernels results
                msize = ceil(sqrt(obj.numFrames+1));
                maxis = max(obj.k(:));
                figure(300),
                for i = 1:obj.numFrames
                    subplot(msize,msize,i);
                    imagesc(obj.k(:,:,i)); caxis([0,maxis]);
                    title(['Blur in frame ',num2str(i)])
                end
                subplot(msize,msize,obj.numFrames+1)
                imagesc(obj.kOpts.initKx'*obj.kOpts.initKy); caxis([0,maxis]);
                title('Initial blur');
                drawnow;
            end
            
            
            %Update operator in prostU object
            DK = [];
            for i = 1:obj.numFrames % this needs to be done smarter at some point
                tmpOp = RepConvMtx(obj.k(:,:,i),obj.dimsLarge); %
                tmpOp = superpixelOperator(obj.dimsSmall,obj.factor).matrix*tmpOp;
                DK    = blkdiag(DK,tmpOp);
            end
            obj.prostU.data.linop{1,1}{1,4}{1} = DK;
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
        
        
        %% Calculate Error margins if ground truths are given
        function [psnrErr,ssimErr, psnrV] = calculateErrors(obj)
            psnrErr = -1;
            psnrV = -1;
            
            %image error
            if (~isscalar(obj.gtU))
                if obj.colorflag
                    % 4D PSNR
                    psnrErr = psnr(obj.result2(20:end-20,20:end-20,:,:),obj.gtU(20:end-20,20:end-20,:,1:obj.numFrames));
                    
                    for j = 1:3 % SSIM computation takes a while ...
                        ssimErr(j) = ssim(squeeze(obj.result2(20:end-20,20:end-20,j,:)),squeeze(obj.gtU(20:end-20,20:end-20,j,1:obj.numFrames))); %#ok<AGROW>
                    end
                    ssimErr = mean(ssimErr); % mean over all color channels
                else
                    % 3D PSNR
                    psnrErr = psnr(obj.result2(20:end-20,20:end-20,:),obj.gtU(20:end-20,20:end-20,1:obj.numFrames));
                    
                    % SSIM computation takes a while ...
                    ssimErr = ssim(obj.result2(20:end-20,20:end-20,:),obj.gtU(20:end-20,20:end-20,1:obj.numFrames));
                end
                
                disp(['PSNR (central patch): ',num2str(psnrErr),' dB']);
                disp(['SSIM (central patch): ',num2str(ssimErr),' ']);
                
            end
            
            %flow error
            if (~isscalar(obj.gtV))
                % 4D PSNR
                psnrV = psnr(obj.v(20:end-20,20:end-20,:,:),obj.gtV(20:end-20,20:end-20,:,:));
                disp(['PSNR of flow field(central patch): ',num2str(psnrV),' dB']);
            end
        end
        
        %% Run the algorithm
        function run(obj)
            
            for i=1:obj.numMainIt
                
                obj.currentIt = i;
                % solve u problem
                disp('Solving problem for u');
                obj.solveU;
                
                % solve problem v
                if obj.currentIt ~= obj.numMainIt  % ;) minimize Eqyptian activity during testing
                    disp('Solving problem for v');
                    obj.solveV;
                end
                
                % update warping operators and parameters for u problem
                if obj.currentIt ~= obj.numMainIt
                    obj.updateProstU;
                end
                
                % update blur kernels
                if obj.currentIt ~= obj.numMainIt
                    disp('Solving blur problem');
                    obj.solveBlur;
                end
                
                disp(['-------Main Iteration ',num2str(i), ' finished !'])
                
                % Produce preliminary output for testing
                if obj.currentIt ~= obj.numMainIt
                    obj.recomputeRGB;
                    centralSlice = ceil(obj.numFrames/2);
                    outImage = obj.result1(20:end-20,20:end-20,:,centralSlice);
                    psnrVal = psnr(outImage,obj.gtU(20:end-20,20:end-20,:,centralSlice)) %#ok<NOPRT>
                    imwrite(outImage,['../outFolder/',obj.datasetName,'/','8_SuperResolutionOurs_P_final_it_',num2str(i),'_PSNR_',num2str(psnrVal,4),'.tif']);
                    outImage = obj.result2(20:end-20,20:end-20,:,centralSlice);
                    psnrVal = psnr(outImage,obj.gtU(20:end-20,20:end-20,:,centralSlice));
                    imwrite(obj.result2(20:end-20,20:end-20,:,centralSlice),['../outFolder/',obj.datasetName,'/SRvariations/','8_SuperResolutionOurs_P_final_2_its_u-w',num2str(psnrVal,4),'.tif']);
                    WriteToVideo(obj.result1,['../outFolder/',obj.datasetName,'/videos/','8_SuperResolutionOurs_P_final_2_its_.avi']);
                end
            end
            
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



