classdef jointSuperResolutionJonas < handle
    %JOINTSUPERRESOLUTION
    % solve a joint super-resolution and optical flow problem, estimating
    % the optical flow 'v' of an unknown high resolution video 'u' from a low
    % resolution video 'imageSequenceSmall' while jointly recovering the
    % high resolution video.
    %% Class properties:
    properties
        
        % Input Data
        imageSequenceSmall
        datasetName
        
        % Problem Structure
        numFrames
        factor
        dimsLarge
        dimsSmall
        
        % Main variables
        u                  % image sequence
        v                  % flow sequence
        k                  % blur sequence
        
        % Ground Truth variables
        gtU
        gtV
        
        % Solver options
        prostU             % main solver class for u subproblem
        mainV              % main solver class for v subproblem
        prostK             % main solver class for k subproblem
        
        opts               % prost options
        regU, regUq, regTGV% various u regularizer options
        regV               % v regularizer options
        alpha              % weights for regularizer of u
        beta               % weights for regularizer of v
        gamma              % coupling intensity (= WarpingOp factor)
        kOpts              % k regularizer and backend options
        kernelsize         % targeted kernel size for each k_i
        verbose
        profiler           % initialize profiler GUI for init and run methods
        numMainIt;         % Number of total iterations
        
        
        
        % Currently unused properties
        % vForward
        % vExists
        % numWarpTerm
        
        
        % Book-keeping
        currentIt          % notifies solvers of current iteration
    end
    
    methods
        %% Initialize object
        function obj = jointSuperResolutionJonas(imageSequenceSmall,varargin)
            vararginParser;
            
            %%% Property standards:
            
            % Input Data
            obj.imageSequenceSmall = imageSequenceSmall;
            if (exist('datasetName','var'))
                obj.datasetName = datasetName;
            else
                obj.datasetName = '';
            end
            
            % Problem Structure
            obj.numFrames = size(imageSequenceSmall,3);
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
            
            % Solver options
            %prostU                                 % main solver class for u subproblem
            %mainV                                  % main solver class for v subproblem
            %prostK                                 % main solver class for k subproblem
            % will be constructed in init_u, init_v and init_k calls
            
            obj.numMainIt = 1;                      % Number of total outer iterations
            obj.opts.backend = prost.backend.pdhg(...
                'tau0', 100, ...
                'sigma0', 0.01, ...
                'stepsize', 'boyd');                % prost backend options
            obj.opts.opts = prost.options('max_iters', 500, ...
                'num_cback_calls', 5, ...
                'verbose', true);                   % prost structure options
            
            obj.regU = 'TV';                        % u regularizer options
            obj.regUq = 1;                          % u regularizer options
            obj.regTGV = 0.3;                       % u regularizer options
            obj.regV = 'Huber';                     % v regularizer options
            
            obj.alpha = 0.01*ones(obj.numMainIt,1); % weights for regularizer of u
            obj.beta = 0.1*ones(obj.numMainIt,1);   % weights for regularizer of v
            obj.gamma = ones(obj.numMainIt,1);      % coupling intensity (= WarpingOp factor)
            
            obj.kOpts.backend = prost.backend.pdhg(...
                'tau0', 100, ...
                'sigma0', 0.01, ...
                'stepsize', 'boyd');                % prost backend options
            obj.kOpts.opts = prost.options('max_iters', 500, ...
                'num_cback_calls', 5, ...
                'verbose', true);                   % prost structure options
            obj.kernelsize = 5;                     % standard kernel size
            obj.kOpts.delta = 0.05*ones(obj.numMainIt,1);   % l2 penalty on k
            obj.kOpts.zeta  = 0.01*ones(obj.numMainIt,1);   % Tikh penalty on k
            obj.kOpts.initKx =  [0.05 0.15 0.6 0.15 0.05];  % initial separable kernel (x)
            obj.kOpts.initKy = [0.05 0.15 0.6 0.15 0.05];   % initial separable kernel (y)
            
            
            if (exist('verbose','var'))
                obj.verbose = verbose;
            else
                obj.verbose = 1;
            end
            
            
            % unused
            %vForward
            %vExists
            %numWarpTerm
            
            
            % Book-keeping
            obj.currentIt = 1;         % notifies solvers of current iteration
            
            
        end
        
        %% create prost object for u
        function init_u(obj)
            %%%% Create operators
            
            % create downsampling operator
            downsamplingOp = superpixelOperator(obj.dimsSmall,obj.factor).matrix;
            % and incorporate blur:
            downsamplingOp =downsamplingOp*createSparseMatrixFrom1dSeperableKernel( ...
                                  obj.kOpts.initKx,obj.kOpts.initKy, obj.dimsLarge(1),obj.dimsLarge(2), 'Neumann');                   
            % initialize block format
            downsamplingOp = kron(speye(obj.numFrames),downsamplingOp);
           
           
            % short some notation
            temp = obj.dimsSmall;
            ny=temp(1); nx=temp(2);
            temp = obj.dimsLarge;
            Ny = temp(1);
            Nx = temp(2);
            nc =  obj.numFrames;
            
            
            % create warping operator from given v
            for i = 1:nc-1
                % extract flow field
                singleField = squeeze(obj.v(:,:,i,:));
                
                % create warping operator forward and backward
                warp = warpingOperator(obj.dimsLarge,singleField);
                idOp = speye(Nx*Ny);
                
                % find out of range warps in each of the operators and set the corresponding line in the other operator also to zero
                marker = sum(abs(warp),2) == 0;
                warp(marker > 0,:) = 0;
                idOp(marker > 0,:) = 0; 
                
                %Build the sparse total warping operator
                if i == 1
                    warpingOp = sparse([],[],[],Nx*Ny*(nc-1),Nx*Ny*nc,2*nnz(warp)*nc); % TODO: check preallocation accuracy
                end
                
                if (mod(i,2)==1)
                    %add warping operator forward
                    warpingOp((Nx*Ny*(i-1)+1):(Nx*Ny*i), (Nx*Ny*(i-1)+1):(Nx*Ny*i)) = -idOp; %#ok<*SPRIX> %;)
                    warpingOp((Nx*Ny*(i-1)+1):(Nx*Ny*i), (Nx*Ny*i+1):(Nx*Ny*(1+i))) = warp;
                else
                    %add warping operator backward
                    warpingOp((Nx*Ny*(i-1)+1):(Nx*Ny*i), (Nx*Ny*(i-1)+1):(Nx*Ny*i)) = warp ; 
                    warpingOp((Nx*Ny*(i-1)+1):(Nx*Ny*i), (Nx*Ny*i+1):(Nx*Ny*(1+i))) = -idOp;
                end
            end
            % use weighting parameter gamma
            warpingOp = obj.gamma(obj.currentIt)*warpingOp;
            
            
            %%%% initialize prost variables and operators
            
            % Primal- Core Variable:
            u_vec = prost.variable(Nx*Ny*nc);
            
            % Primal - Varying regularizer primals:
            if strcmp(obj.regU,'TGV')
                w_vec = prost.variable(Nx*Ny*nc*2); % additional primal variable for TGV reformulation
            end
            
            % Dual - Core Variables
            p = prost.variable(nx*ny*nc);           % dual variable for fitting to data
            q = prost.variable(Nx*Ny*(nc-1));       % dual variable for warping between data
            
            
            % Dual - Varying regularizer duals:
            if strcmp(obj.regU,'TV')
                g = prost.variable(2*Nx*Ny*nc);     % dual variable for TV regularization
            elseif strcmp(obj.regU,'TGV')
                g1 = prost.variable(2*Nx*Ny*nc);     % dual variable for TGV regularization (to u)
                g2 = prost.variable(4*Nx*Ny*nc);     % dual variable for TGV regularization (to w)
            else
                error('invalid regularization type')
            end
            
            %%%% Connect variables to a prost min-max problem and add functions
            if strcmp(obj.regU,'TV')
                obj.prostU = prost.min_max_problem( {u_vec}, {q,p,g} );
                % Core
                obj.prostU.add_dual_pair(u_vec, p, prost.block.sparse(downsamplingOp));
                obj.prostU.add_dual_pair(u_vec, q, prost.block.sparse(warpingOp));
                % Regularizer
                obj.prostU.add_dual_pair(u_vec, g, prost.block.gradient2d(Nx,Ny,nc,false));
                
                % add functions
                % Core
                obj.prostU.add_function(p, prost.function.sum_1d('ind_box01', 0.5, -0.5, 1, obj.imageSequenceSmall(:), 0)); %l^1
                obj.prostU.add_function(q, prost.function.sum_1d('ind_box01', 0.5, -0.5, 1, 0, 0)); %l^1
                % Regularizer
                if obj.regUq == 1
                    obj.prostU.add_function(g, prost.function.sum_norm2(2, false, 'ind_leq0', 1/obj.alpha(obj.currentIt), 1, 1, 0, 0)); %l^{2,1}
                else
                    obj.prostU.add_function(g, ...
                        prost.function.conjugate(...
                        prost.function.sum_norm2(...
                        2, false, 'lq', ...
                        1,0,obj.alpha(obj.currentIt),0,0, obj.regUq))); % l^\regUq
                end
                
            elseif strcmp(obj.regU,'TGV')
                obj.prostU = prost.min_max_problem( {u_vec,w_vec}, {q,p,g1,g2} );
                % Core
                obj.prostU.add_dual_pair(u_vec, p, prost.block.sparse(downsamplingOp));
                obj.prostU.add_dual_pair(u_vec, q, prost.block.sparse(warpingOp));
                % Regularizer
                obj.prostU.add_dual_pair(w_vec, g1, prost.block.identity(-1));
                obj.prostU.add_dual_pair(w_vec, g2, prost.block.sparse(spmat_gradient2d(Nx, Ny, 2*nc)));
                obj.prostU.add_dual_pair(u_vec, g1, prost.block.gradient2d(Nx,Ny,nc,false));
                
                % Add functions
                % Core
                obj.prostU.add_function(p, prost.function.sum_1d('ind_box01', 0.5, -0.5, 1, obj.imageSequenceSmall(:), 0)); %l^1
                obj.prostU.add_function(q, prost.function.sum_1d('ind_box01', 0.5, -0.5, 1, 0, 0)); %l^1
                % Regularizer
                if obj.regUq == 1
                    obj.prostU.add_function(g1, prost.function.sum_norm2(2, false, 'ind_leq0', 1/obj.alpha(obj.currentIt), 1, 1, 0, 0)); %l^{2,1}
                    obj.prostU.add_function(g2, prost.function.sum_norm2(4, false, 'ind_leq0', 1/(obj.alpha(obj.currentIt)*obj.regTGV), 1, 1, 0, 0)); %l^{2,1}
                else
                    obj.prostU.add_function(g1, ...
                        prost.function.conjugate(...
                        prost.function.sum_norm2(...
                        2, false, 'lq', ...
                        1,0,obj.alpha(obj.currentIt),0,0, obj.regUq))); % l^\regUq
                    obj.prostU.add_function(g2, ...
                        prost.function.conjugate(...
                        prost.function.sum_norm2(...
                        4, false, 'lq', ...
                        1,0,obj.alpha(obj.currentIt)*obj.regTGV,0,0, obj.regUq))); % l^\regUq
                end
            end
            
            
            % update solver backend if a nonconvex regularizer is chosen
            if obj.regUq < 1
                    gamma_factor = 1/sqrt(prod(obj.dimsLarge))*10; % experimental :>
                    obj.opts.backend = prost.backend.pdhg('stepsize', 'alg2', ...
                             'residual_iter', -1, ...
                             'alg2_gamma', gamma_factor);
                    obj.opts.opts = prost.options('max_iters', 1500, ...         
                                        'num_cback_calls', 0, ...
                                        'verbose', true); 
                    %( initialization is actually the bigger problem than a
                    %  few thousand extra iterations )
            end
            
        end
        
        %% calculate initial velocity fields on low resolution input images and scale them up to target resolution
        function init_v0(obj)
            if (obj.verbose > 0)
                disp('Calculating initial velocity fields')
            end
            
            
            for j=1:obj.numFrames-1
                if (mod(j,2)==1)%calculate backward flow: v s.t. u_2(x+v)=u_1(x)
                    uTmpSmall = cat(3,obj.imageSequenceSmall(:,:,j),obj.imageSequenceSmall(:,:,j+1));
                else%calculate backward flow: v s.t. u_1(x+v)=u_2(x)
                    uTmpSmall = cat(3,obj.imageSequenceSmall(:,:,j+1),obj.imageSequenceSmall(:,:,j));
                end
                
                motionEstimatorLow = motionEstimatorClass(uTmpSmall,1e-6,obj.beta(1),'doGradientConstancy',1);
                motionEstimatorLow.verbose = 0;%obj.verbose;
                motionEstimatorLow.init;
                motionEstimatorLow.runPyramid;
                
                vTmp = motionEstimatorLow.getResult;
                
                obj.v(:,:,j,1) = imresize(vTmp(:,:,1,1), obj.dimsLarge) * obj.factor;
                obj.v(:,:,j,2) = imresize(vTmp(:,:,1,2), obj.dimsLarge) * obj.factor;
                
                if (obj.verbose > 0)
                    figure(200+j);imagesc(flowToColorV2(cat(3,obj.v(:,:,j,1),obj.v(:,:,j,2))));axis image;
                end
                %todo: add comparison with ground-truth flow
                
                
                if (obj.verbose > 0)
                    disp(['Initial velocity field calculated for frames [',num2str(j),',',num2str(j+1),']'])
                end
                
            end
        end
        
        %% create real motion estimator object in motionEstimatorClass
        %(this seems actually unnecessary due to reinit in solveV ?)
        function init_v(obj)
            
            obj.mainV = motionEstimatorClass(obj.u,1e-6,obj.beta,'doGradientConstancy',1);
            obj.mainV.verbose = 1;
            obj.mainV.init;
            %obj.vExists = 0;
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
            g1    = prost.variable(kx*ky*nc);   % dual variable for l2 norm reg
            g2    = prost.variable(kx*ky*nc*2); % dual variable for Tikh reg
            
            % Connect variables
            obj.prostK = prost.min_max_problem( {k_vec}, {p,g1,g2} ); 
            obj.prostK.add_dual_pair(k_vec,p,prost.block.sparse(downsampleImageOp));
            obj.prostK.add_dual_pair(k_vec,g1,prost.block.identity(1));
            obj.prostK.add_dual_pair(k_vec,g2,prost.block.gradient2d(kx,ky,nc,false));
            
            % Add functions
            obj.prostK.add_function(k_vec,prost.function.sum_ind_simplex(kx*ky,true)); % this should be a simplex for each k_i
            obj.prostK.add_function(p,prost.function.sum_1d('ind_box01',0.5,-0.5,1,obj.imageSequenceSmall(:),0)); % l1 data
            obj.prostK.add_function(g1,prost.function.sum_1d('square', 1, 0, 1/obj.kOpts.delta(1), 0, 0)); % l2 penalty
            obj.prostK.add_function(g2,prost.function.sum_1d('square', 1, 0, 1/obj.kOpts.zeta(1), 0, 0)); % Tikh penalty
        end
        
        %% collect inits and validate
        function init(obj)
            
            % validate and set variables
            
            % update sizes for non-standard factor
            obj.dimsLarge = obj.factor*obj.dimsSmall;
            
            % re-initialize u,v
            obj.u = zeros([obj.dimsLarge,obj.numFrames]);
            obj.v = zeros([obj.dimsLarge,obj.numFrames,2]);
            obj.k = zeros(obj.kernelsize^2*obj.numFrames,1);
            
            if length(obj.alpha) ~= obj.numMainIt
                warning('u regularization weights wrong, iteration number extended or alpha padded')
                if obj.numMainIt < length(obj.alpha)
                    obj.numMainIt = length(obj.alpha);
                else
                    obj.alpha(end+1:obj.numMainIt) = obj.alpha(end);
                end
            end
            if length(obj.beta) ~= obj.numMainIt
                warning('v regularization weights wrong, iteration number extended or beta padded')
                if obj.numMainIt < length(obj.beta)
                    obj.numMainIt = length(obj.beta);
                else
                    obj.alpha(end+1:obj.numMainIt) = obj.beta(end);
                end
            end
            if length(obj.gamma) ~= obj.numMainIt
                warning('warp weights wrong, iteration number extended or gamma padded')
                if obj.numMainIt < length(obj.gamma)
                    obj.numMainIt = length(obj.gamma);
                else
                    obj.gamma(end+1:obj.numMainIt) = obj.gamma(end);
                end
            end
            if length(obj.kOpts.delta) ~= obj.numMainIt
                warning('k l2 regularization weights wrong, iteration number extended or delta padded')
                if obj.numMainIt < length(obj.kOpts.delta)
                    obj.numMainIt = length(obj.kOpts.delta);
                else
                    obj.kOpts.delta(end+1:obj.numMainIt) = obj.kOpts.delta(end);
                end
            end
            if length(obj.kOpts.zeta) ~= obj.numMainIt
                warning('k Tikh regularization weights wrong, iteration number extended or zeta padded')
                if obj.numMainIt < length(obj.kOpts.zeta)
                    obj.numMainIt = length(obj.kOpts.zeta);
                else
                    obj.kOpts.zeta(end+1:obj.numMainIt) = obj.kOpts.zeta(end);
                end
            end
            
            if obj.profiler  % start profiler gui
                profile -memory on
                setpref('profiler', 'showJitLines', true); % I wonder if this still works for 2016a
                profile on
            end
            
            % Call actual initialization methods
            obj.init_v0;
            obj.init_v;
            obj.init_u;
            obj.init_k;
            
            if obj.profiler
                profile off
                profile viewer
                pause();
            end
            
            
            
            if (obj.verbose > 0)
                disp('Initialization of u, v and k solvers finished');
            end
        end
        
        %% Call u solver in prost
        function solveU(obj)
            
            %call prost framework
            tic
            %prost_('set_gpu', 0)
            prost.solve(obj.prostU, obj.opts.backend, obj.opts.opts);
            toc
            obj.u = reshape(obj.prostU.primal_vars{1,1}.val,[obj.dimsLarge(1), obj.dimsLarge(2), obj.numFrames]);
            

            %%% sanity check:
            %obj.u = call_prost_wrapper(obj);



            % show solution
            for j=1:obj.numFrames
                if (obj.verbose > 0)
                    figure(100+j);imagesc(obj.u(:,:,j),[0,1]);axis image;colormap(gray);drawnow
                end
            end
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
            
            
            % create warping operator from given v
            for i = 1:nc-1
                % extract flow field
                singleField = squeeze(obj.v(:,:,i,:));
                
                % create warping operator forward and backward
                warp = warpingOperator(obj.dimsLarge,singleField);
                idOp = speye(Nx*Ny);
                
                % find out of range warps in each of the operators and set the corresponding line in the other operator also to zero
                marker = sum(abs(warp),2) == 0;
                warp(marker > 0,:) = 0;
                idOp(marker > 0,:) = 0;
                
                %Build the sparse total warping operator
                if i == 1
                    warpingOp = sparse([],[],[],Nx*Ny*(nc-1),Nx*Ny*nc,2*nnz(warp)*(nc-1)); % TODO: check preallocation accuracy
                end
                
                if (mod(i,2)==1)
                    %add warping operator forward
                    warpingOp((Nx*Ny*(i-1)+1):(Nx*Ny*i), (Nx*Ny*(i-1)+1):(Nx*Ny*i)) = -idOp;
                    warpingOp((Nx*Ny*(i-1)+1):(Nx*Ny*i), (Nx*Ny*i+1):(Nx*Ny*(1+i))) = warp;
                else
                    %add warping operator backward
                    warpingOp((Nx*Ny*(i-1)+1):(Nx*Ny*i), (Nx*Ny*(i-1)+1):(Nx*Ny*i)) = warp ; % flexbox updates operator and operatorT here, dunno if this is correctly represented ...
                    warpingOp((Nx*Ny*(i-1)+1):(Nx*Ny*i), (Nx*Ny*i+1):(Nx*Ny*(1+i))) = -idOp;
                end
            end
            
            % update warping operator in prost
             obj.prostU.data.linop{1,2}{1,4}{1} = warpingOp;

            
            % update alpha for next main iteration
            if obj.currentIt ~= obj.numMainIt
                if obj.regUq == 1
                    obj.prostU.data.prox_fstar{1,3}{1,5}{1,4}{1} = 1/obj.alpha(obj.currentIt+1);
                else
                    obj.prostU.data.prox_fstar{1,2}{1,5}{1,1}{1,5}{1,4}{3} = obj.alpha(obj.currentIt+1);
                end
            end
            % the question is however where the second TGV parameter can be found ...
            
        end
        
        %% Call v solver from motionEstimationGUI (which wraps flexBox)
        function solveV(obj)
            for j=1:obj.numFrames-1
                if (mod(j,2)==1)    %calculate backward flow: v s.t. u_2(x+v)=u_1(x)
                    uTmp = cat(3,obj.u(:,:,j),obj.u(:,:,j+1));
                else                %calculate backward flow: v s.t. u_1(x+v)=u_2(x)
                    uTmp = cat(3,obj.u(:,:,j+1),obj.u(:,:,j));
                end
                
                motionEstimator = motionEstimatorClass(uTmp,1e-6,obj.beta(obj.currentIt),'doGradientConstancy',1);
                motionEstimator.regularizerTerm = obj.regV;
                motionEstimator.verbose = 0;%obj.verbose; 
                motionEstimator.init;
                motionEstimator.runPyramid;
                
                vTmp =  motionEstimator.getResult;
                obj.v(:,:,j,1) = vTmp(:,:,1,1);
                obj.v(:,:,j,2) = vTmp(:,:,1,2);
                
                if obj.verbose
                    disp(['Updated high-res velocity field calculated for frames [',num2str(j),',',num2str(j+1),']']);
                end
                
            end
            
            if (obj.verbose > 0)
                figure(200+j);imagesc(flowToColorV2(cat(3,obj.v(:,:,j,1),obj.v(:,:,j,2))));axis image;
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
            if obj.verbose 
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

            
            % Update operator in prostU object
            DK = [];
            for i = 1:obj.numFrames % this needs to be done smarter at some point
                %tmpOp = writeKernelToSparseDownsamplingMatrix(obj.k(:,:,i),1,obj.dimsLarge(1),obj.dimsLarge(2)); % there is a bug in here somewhere
                tmpOp = RepConvMtx(obj.k(:,:,i),obj.dimsLarge); %
                tmpOp = superpixelOperator(obj.dimsSmall,obj.factor).matrix*tmpOp;
                DK    = blkdiag(DK,tmpOp);
            end
            obj.prostU.data.linop{1,1}{1,4}{1} = DK;
        end
        
        %% Calculate Error margins if ground truths are given
        function [errorU,errorV] = calculateErrors(obj)
            errorU = -1;
            errorV = -1;
            
            if size(obj.gtU) ~= size(obj.u)
                warning('Wrong ground truth data given, error estimation terminated')
                return
            end
            
            
            %image error
            if (~isscalar(obj.gtU))
                errorU = 0;
                for j=1:obj.numFrames
                    %squared error
                    err = (obj.u(:,:,j) - obj.gtU(:,:,j)).^2;
                    
                    %cut out inner window
                    err = err(20:end-20,20:end-20);
                    
                    %figure(124);imagesc(err);colorbar;
                    %figure(123);imagesc(err);colorbar;
                    %sum up
                    err = sum(sqrt(err(:))) / numel(err);
                    errorU = errorU + err;
                end
                errorU = errorU / obj.numFrames; %average;
                
                psnrErr = psnr(obj.u(20:end-20,20:end-20,:),obj.gtU(20:end-20,20:end-20,1:obj.numFrames));
                ssimErr = ssim(obj.u(20:end-20,20:end-20,:),obj.gtU(20:end-20,20:end-20,1:obj.numFrames));
                
                disp(['PSNR (central patch): ',num2str(psnrErr),' dB']);
                disp(['SSIM (central patch): ',num2str(ssimErr),' ']);
                
                % save psnr value and parameters, just in case
                save(['./parameters/',obj.datasetName,'_',num2str(psnrErr,4),'.mat'],'psnrErr','ssimErr');
                regU = obj.regU; regUq  = obj.regUq; regTGV = obj.regTGV; alpha = obj.alpha; regV = obj.regV; beta= obj.beta; gamma = obj.gamma;  %#ok<NASGU,PROP>
                save(['./parameters/',obj.datasetName,'_',num2str(psnrErr,4),'.mat'],'regU','regUq', ...
                    'regTGV','alpha','regV','beta','gamma','-append');
                
            end
            
            %flow error
            if (~isscalar(obj.gtV))
                errorV = 0;
                for j=1:obj.numFrames-1
                    field = squeeze(obj.v(:,:,j,:));
                    fieldGT = squeeze(obj.gtV(:,:,j,:));
                    
                    figure(121);imagesc(flowToColorV2(field));
                    figure(122);imagesc(flowToColorV2(fieldGT));
                    
                    err = absoluteError(field,fieldGT);
                    errorV = errorV + err;
                end
                errorV = errorV / (obj.numFrames-1);
            end
        end
        
        %% Run the algorithm
        function run(obj)
            
            if obj.profiler
                profile clear
                profile on
            end
            
            for i=1:obj.numMainIt
                
                obj.currentIt = i;
                
                %% solve u problem
                disp('Solving problem for u');
                obj.solveU;
                
                %% solve problem v
                if obj.currentIt ~= obj.numMainIt  % ;) minimize Eqyptian activity during testing
                    disp('Solving problem for v');
                    obj.solveV;
                end
                
                %% update warping operators and parameters for u problem
                obj.updateProstU;
                
                %% update blur kernels       
                disp('Solving blur problem');
                obj.solveBlur;

                
                disp(['-------Main Iteration ',num2str(i), ' finished !'])
            end
            
            
            
            if obj.profiler
                profile off
                profile view
            end
            disp('Algorithm finished !!');
        end
        
    end
    
end

