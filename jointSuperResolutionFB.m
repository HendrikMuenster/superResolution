classdef jointSuperResolutionFB < handle
    %JOINTSUPERRESOLUTION Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
         numFrames;
         dimsLarge
         dimsSmall
         factor
         alpha %weight for tv u
         beta %weight for tv v
         betaStart %initial weight for tv v
         numMainIt;
         imageSequenceSmall
         mainU
         mainV
         numWarpTerm
         u %image sequence
         v %flow sequence
         vForward
         
         verbose
         
         gtU
         gtV
         
         vExists
    end
    
    methods
        function obj = jointSuperResolutionFB(imageSequenceSmall,alpha,beta,factor,varargin)
            vararginParser;
            
            obj.imageSequenceSmall = imageSequenceSmall;
            
            obj.alpha = alpha;
            obj.beta = beta;
            obj.factor = factor;
            obj.dimsSmall = size(imageSequenceSmall);
            obj.dimsSmall = obj.dimsSmall(1:2);
            
            obj.dimsLarge = obj.factor*obj.dimsSmall;
            
            obj.numMainIt = 3;
            obj.numFrames = size(imageSequenceSmall,3);
            
            obj.u = zeros([obj.dimsLarge,obj.numFrames]);
            obj.v = zeros([obj.dimsLarge,obj.numFrames,2]);
            
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
            
            if (exist('verbose','var'))
                obj.verbose = verbose;
            else
                obj.verbose = 1;
            end
            
            if (exist('betaStart','var'))
                obj.betaStart = betaStart;
            else
                obj.betaStart = beta;
            end
        end
        
        function init(obj)
            if (obj.verbose > 0)
                disp('Calculating initial velocity fields')
            end
            
            %% calculate initial velocity fields on low resolution input images and scale them up to target resolution
            for j=1:obj.numFrames-1
                if (mod(j,2)==1)%calculate backward flow: v s.t. u_2(x+v)=u_1(x)
                    uTmpSmall = cat(3,obj.imageSequenceSmall(:,:,j),obj.imageSequenceSmall(:,:,j+1));
                else%calculate backward flow: v s.t. u_1(x+v)=u_2(x)
                    uTmpSmall = cat(3,obj.imageSequenceSmall(:,:,j+1),obj.imageSequenceSmall(:,:,j));
                end
                
                motionEstimatorLow = motionEstimatorClass(uTmpSmall,1e-6,obj.betaStart,'doGradientConstancy',1);
                motionEstimatorLow.verbose = obj.verbose;
                motionEstimatorLow.init;
                motionEstimatorLow.runPyramid;
                
                vTmp = motionEstimatorLow.getResult;
                
                obj.v(:,:,j,1) = imresize(vTmp(:,:,1,1), obj.dimsLarge) * obj.factor;
                obj.v(:,:,j,2) = imresize(vTmp(:,:,1,2), obj.dimsLarge) * obj.factor;

                if (obj.verbose > 0)
                    figure(200+j);imagesc(flowToColorV2(cat(3,obj.v(:,:,j,1),obj.v(:,:,j,2))));axis image;
                end
                %todo: add comparison with ground-truth flow
            end
            
            if (obj.verbose > 0)
                disp('Initial velocity fields calculated')
            end
            
            
            
            %% init real motion estimator
            obj.mainV = motionEstimatorClass(obj.u,1e-6,obj.beta,'doGradientConstancy',1);
            obj.mainV.verbose = 1;
            obj.mainV.init;
            obj.vExists = 0;
            %% create flexBox object for u
            obj.mainU = flexBox;

            obj.mainU.params.showPrimals = 100;
            obj.mainU.params.tol = 1e-6;
            obj.mainU.params.verbose = 1;
            obj.mainU.params.tryCPP = 1;

    
            %initialize downsampling operator
            downsamplingOp = superpixelOperator(obj.dimsSmall,obj.factor);

            %add for each frame data term and tv term
            for i=1:obj.numFrames
                %add primal variable for each image
                obj.mainU.addPrimalVar(obj.dimsLarge);

                %add data term
                obj.mainU.addTerm(L1dataTermOperator(1,downsamplingOp.matrix,obj.imageSequenceSmall(:,:,i)),i);

                %add tv term for each primal var
                obj.mainU.addTerm(L1gradientIso(obj.alpha,obj.dimsLarge),i);
            end

            %connect subsequent images by optical flow term
            for i=1:obj.numFrames - 1
                %extract flow field
                singleField = squeeze(obj.v(:,:,i,:));

                %create warping operator forward and backward
                warp = warpingOperator(obj.dimsLarge,singleField);

                %find out of range warps in each of the operators and set the
                %corresponding line in the other operator also to zero
                marker = sum(abs(warp),2) == 0;
                
                idOp = speye(size(warp,1));

                warp(marker > 0,:) = 0;
                idOp(marker > 0,:) = 0;

                if (mod(i,2)==1)
                    %add warping operator to flexbox
                    obj.mainU.addTerm(L1operatorAniso(1,2,{-idOp,warp}),[i,i+1]);
                    obj.numWarpTerm(i) = numel(obj.mainU.duals);
                else
                    %add warping operator to flexbox
                    obj.mainU.addTerm(L1operatorAniso(1,2,{warp,-idOp}),[i,i+1]);
                    obj.numWarpTerm(i) = numel(obj.mainU.duals);              
                end
            end
            
            if (obj.verbose > 0)
                disp('Initialization finished');
            end
        end
        
        function solveU(obj)
            obj.mainU.runAlgorithm;

            % extract solution and store into matrix u
            for j=1:obj.numFrames
                obj.u(:,:,j) = obj.mainU.getPrimal(j);
                if (obj.verbose > 0)
                    figure(100+j);imagesc(obj.u(:,:,j),[0,1]);axis image;colormap(gray);drawnow
                end
            end
        end
        
        
        
        function updateFlexBoxU(obj)
            for j=1:obj.numFrames-1
                %extract flow field
                singleField = squeeze(obj.v(:,:,j,:));

                %create warping operator forward and backward
                warp = warpingOperator(obj.dimsLarge,singleField);
                idOp = speye(size(warp,1));

                %find out of range warps in each of the operators and set the
                %corresponding line in the other operator also to zero
                marker = sum(abs(warp),2) == 0;
                

                warp(marker > 0,:) = 0;
                idOp(marker > 0,:) = 0;
                
                if (mod(j,2)==1)
                    obj.mainU.duals{obj.numWarpTerm(j)}.operator{1} = -idOp;
                    obj.mainU.duals{obj.numWarpTerm(j)}.operatorT{1} = -idOp;
                    
                    obj.mainU.duals{obj.numWarpTerm(j)}.operator{2} = warp;
                    obj.mainU.duals{obj.numWarpTerm(j)}.operatorT{2} = warp';
                else
                    obj.mainU.duals{obj.numWarpTerm(j)}.operator{1} = warp;
                    obj.mainU.duals{obj.numWarpTerm(j)}.operatorT{1} = warp';
                    
                    obj.mainU.duals{obj.numWarpTerm(j)}.operator{2} = -idOp;
                    obj.mainU.duals{obj.numWarpTerm(j)}.operatorT{2} = -idOp;
                end
            end
        end
           
        function solveV(obj)
            for j=1:obj.numFrames-1
                if (mod(j,2)==1)%calculate backward flow: v s.t. u_2(x+v)=u_1(x)
                    uTmp = cat(3,obj.u(:,:,j),obj.u(:,:,j+1));
                else%calculate backward flow: v s.t. u_1(x+v)=u_2(x)
                    uTmp = cat(3,obj.u(:,:,j+1),obj.u(:,:,j));
                end
                
                motionEstimatorLow = motionEstimatorClass(uTmp,1e-6,obj.beta,'doGradientConstancy',1);
                motionEstimatorLow.verbose = obj.verbose;
                motionEstimatorLow.init;
                motionEstimatorLow.runPyramid;
                
                vTmp = motionEstimatorLow.getResult;
                
                obj.v(:,:,j,1) = imresize(vTmp(:,:,1,1), obj.dimsLarge);
                obj.v(:,:,j,2) = imresize(vTmp(:,:,1,2), obj.dimsLarge);
            end
        end
        
        function [errorU,errorV] = calculateErrors(obj)
            errorU = -1;
            errorV = -1;
            
            %image error
            if (~isscalar(obj.gtU))
                % extract solution and store into matrix u
                for j=1:obj.numFrames
                    obj.u(:,:,j) = obj.mainU.getPrimal(j);
                end
                
                errorU = 0;
                for j=1:obj.numFrames
                    %squared error
                    err = (obj.u(:,:,j) - obj.gtU(:,:,j)).^2;

                    %cut out inner window
                    err = err(20:end-20,20:end-20);
                    
                    figure(124);imagesc(err);colorbar;
                    figure(123);imagesc(err);colorbar;
                    %sum up
                    err = sum(sqrt(err(:))) / numel(err);

                    ssimErr = ssim(obj.u(20:end-20,20:end-20,j),obj.gtU(20:end-20,20:end-20,j))
                    
                    errorU = errorU + err;
                end
                errorU = errorU / obj.numFrames; %average;
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
        
        function run(obj)
            useFlexboxForUupdate = false;
            updateBlurKernels = true;
            
            for i=1:obj.numMainIt
                %% solve u problem
                disp('Solving problem for u');
                obj.solveU;

                %% solve problem v
                disp('Solving problem for v');
                obj.solveV;

                %% update warping operators in flexBox
                obj.updateFlexBoxU;
                
                %% update blur kernels
                disp('Lerning optimal blur kernels');
                if updateBlurKernels
                    nr = 1;
                    kernelSize = 4;
                    alph=0.65; %learning rate
                    for k=1:size(obj.mainU.duals,2)
                        if isa(obj.mainU.duals{k}, 'L1dataTermOperator')
                            [A,~] = optimalBlurAndDownsamplingKernel(obj,kernelSize, nr);
                            obj.mainU.duals{k}.operator{1} = alph*A+ (1-alph)*obj.mainU.duals{k}.operator{1};
                            nr = nr+1;
                        end
                    end
                end           
            end
            
            disp('Finished');
        end
        
    end
    
end


