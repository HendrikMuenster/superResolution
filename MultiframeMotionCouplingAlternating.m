classdef MultiframeMotionCouplingAlternating< MultiframeMotionCoupling
    %MultiframeMotionCoupling
    % Couple all frames successively and solve the super resolution problem
    % alternate over high-res optical flow estimation and super-resolution
    %
    %
    %
    %
    %
    
    %% Class properties:
    properties
        outerIts            % number of outer iterations
        currentIt           % current outer iteration
        psnrVals            % psnrVals for each frame per iteration
        imageSequenceLarge  % for psnr check
       
    end
    
    %% Methods: 
    methods
        
        %% decorate constructor
        function obj = MultiframeMotionCouplingAlternating(imageSequenceSmall,imageSequenceLarge,varargin)
            obj@MultiframeMotionCoupling(imageSequenceSmall);
            obj.imageSequenceLarge = imageSequenceLarge;
            obj.outerIts = 5;
            obj.currentIt = 1;
            obj.psnrVals = 0;
        end
        
        %% decorate init
        function init(obj)
            
            if ~strcmp(obj.framework,'prost')
                error('Alternating only in prost for now');
            end
            obj.psnrVals = zeros(obj.numFrames,obj.outerIts);
            init@MultiframeMotionCoupling(obj);
        end
        
        %% solve high-res optical flow problem
        % is actually init_v, but applied to u instead of imageSequenceSmall
        % todo for rewrite ...
        function solveV(obj) 
            if (obj.verbose > 0)
                disp('Calculating new velocity fields')
            end
            
            
            for j=1:obj.numFrames-1
                
                % Select correct flow direction
                if strcmp(obj.flowDirection,'forward')
                    uTmp = cat(3,obj.u(:,:,j+1),obj.u(:,:,j));
                elseif strcmp(obj.flowDirection,'backward')
                    uTmp = cat(3,obj.u(:,:,j),obj.u(:,:,j+1));
                elseif strcmp(obj.flowDirection,'forward-backward')
                    % alternate forward and backward:
                    if (mod(j,2)==1)
                        uTmp = cat(3,obj.u(:,:,j),obj.u(:,:,j+1));
                    else
                        uTmp = cat(3,obj.u(:,:,j+1),obj.u(:,:,j));
                    end
                else
                    error('Invalid flow direction.');
                end
                
                
                motionEstimator = motionEstimatorClass(uTmp,1e-6,obj.beta);
                motionEstimator.verbose = 0;
                
                % Motion estimator energy parameters:
                motionEstimator.dataTerm = 'L1';
                motionEstimator.regularizerTerm ='Huber';
                motionEstimator.doGradientConstancy = 1;
                
                
                % Pyramid scheme parameters:
                motionEstimator.doGaussianSmoothing = 1;
                motionEstimator.medianFiltering = 1;
                motionEstimator.numberOfWarps = 5;
                motionEstimator.steplength = 0.9;
                
                % Run motion estimator
                motionEstimator.init;
                motionEstimator.runPyramid;
                
                vTmp = motionEstimator.getResult;
                
                % Motion estimator velocity field is x-y instead of m-n.
                obj.v(:,:,j,1) = vTmp(:,:,1,2);
                obj.v(:,:,j,2) = vTmp(:,:,1,1);
                
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
        
        %% update values in optimization problem
        function updateMMCsolver(obj)
            
            % shorten some notation
            temp = obj.dimsLarge;
            Ny = temp(1);
            Nx = temp(2);
            nc =  obj.numFrames;
            
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
                    disp('New warp operator constructed');
                end
            
            % keep h fixed for new iterations !?!
            % h = obj.h
            
            % Create new operator blocks and insert into MMCsolver
            if obj.kappa ~= 1 && ~isnan(obj.kappa) && isempty(obj.u0_frame)
                u_g1_op = [warpingOp/obj.h;obj.kappa*spmat_gradient2d(Nx, Ny,nc)];
                w_g1_op = -[warpingOp/obj.h;obj.kappa*spmat_gradient2d(Nx, Ny,nc)];
                w_g2_op = [obj.kappa*warpingOp/obj.h;spmat_gradient2d(Nx, Ny,nc)];
                
                obj.MMCsolver.data.linop{1, 2}{1, 4}{1, 1}  = u_g1_op;
                obj.MMCsolver.data.linop{1, 3}{1, 4}{1, 1}  = w_g1_op;
                obj.MMCsolver.data.linop{1, 4}{1, 4}{1, 1}  = w_g2_op;
            elseif obj.kappa == 1
                u_g1_op = [warpingOp/obj.h;spmat_gradient2d(Nx, Ny,nc)];
                
                obj.MMCsolver.data.linop{1, 2}{1, 4}{1, 1}  = u_g1_op;
                
            elseif isnan(obj.kappa) 
                obj.MMCsolver.data.linop{1, 2}{1, 4}{1, 1}  = warpingOp/obj.h;
                
            elseif obj.kappa ~= 1 && ~isnan(obj.kappa) && ~isempty(obj.u0_frame)
                u_g1_op = [warpingOp/obj.h;obj.kappa*spmat_gradient2d(Nx, Ny,nc)];
                w_g1_op = -[warpingOp/obj.h;obj.kappa*spmat_gradient2d(Nx, Ny,nc)];
                w_g2_op = [obj.kappa*warpingOp/obj.h;spmat_gradient2d(Nx, Ny,nc)];
                
                obj.MMCsolver.data.linop{1, 2}{1, 4}{1, 1}  = u_g1_op;
                obj.MMCsolver.data.linop{1, 3}{1, 4}{1, 1}  = w_g1_op;
                obj.MMCsolver.data.linop{1, 4}{1, 4}{1, 1}  = w_g2_op;   
            end
        end
        
        
        %% run multiple iterations
        function run(obj)
            
            % outer iteration, order is v0-u-v-u-v ...
            % k is not updated
            for jj=1:obj.outerIts
               disp('Solving problem for u');
               obj.currentIt = jj;
               
               obj.solveU;
               obj.solveV;
               obj.updateMMCsolver;
               
               obj.getPSNR;
               disp(['-------Main Iteration ',num2str(jj),' finished!'])
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
        
        
        %% save psnr values of iteration
        function getPSNR(obj)
            
            obj.recomputeRGB;
            for ii = 1:obj.numFrames
                outImage = obj.result1(20:end-20,20:end-20,:,ii);
                obj.psnrVals(ii,obj.currentIt) = round(psnr(outImage,obj.imageSequenceLarge(20:end-20,20:end-20,:,ii)),2);
                %ssimErr(ii) = round(ssim(outImage,imageSequenceLarge(20:end-20,20:end-20,:,ii)),3);
            end
            
            disp(['PSNR in current iterate is on average:',num2str(mean(obj.psnrVals(:,obj.currentIt)))]);
            
        end
    end
end



