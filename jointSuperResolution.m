classdef jointSuperResolution < handle
    %JOINTSUPERRESOLUTION Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
         numFrames
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
         
         verbose
         
         gtU
         gtV
         
         vExists
    end
    
    methods
        function obj = jointSuperResolution(imageSequenceSmall,alpha,beta,factor,varargin)
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
                obj.betaStart = 0.01;
            end
        end
        
        function init(obj)
            if (obj.verbose > 0)
                disp('Calculating initial velocity fields')
            end
            %% calculate initial velocity fields on low resolution input images and scale them up to target resolution
            motionEstimatorLow = motionEstimatorClass(obj.imageSequenceSmall,1e-6,obj.betaStart);
            motionEstimatorLow.verbose = obj.verbose;
            motionEstimatorLow.init;
            motionEstimatorLow.runPyramid;
            
            if (obj.verbose > 0)
                disp('Initial velocity fields calculated')
            end
            
            vTmp = motionEstimatorLow.getResult;
            
            for j=1:obj.numFrames-1
                obj.v(:,:,j,1) = imresize(vTmp(:,:,j,1), obj.dimsLarge) * obj.factor;
                obj.v(:,:,j,2) = imresize(vTmp(:,:,j,2), obj.dimsLarge) * obj.factor;

                if (obj.verbose > 0)
                    figure(200+j);imagesc(flowToColorV2(cat(3,obj.v(:,:,j,1),obj.v(:,:,j,2))));axis image;
                end
                %todo: add comparison with ground-truth flow
            end
            
            %% init real motion estimator
            obj.mainV = motionEstimatorClass(obj.u,1e-6,obj.beta);
            obj.mainV.verbose = 1;
            obj.mainV.init;
            obj.vExists = 0;
            %% create flexBox object for u
            obj.mainU = flexBox;

            obj.mainU.params.showPrimals = 100;
            obj.mainU.params.tol = 1e-6;
            obj.mainU.params.verbose = 1;
            obj.mainU.params.tryCPP = 1;
    
            %initialize downsampling operator - HACKED BY MICHAEL:
            %ADDITIONAL BLUR BEFORE DOWNSAMPLING
            downsamplingOp = superpixelOperator(obj.dimsSmall,obj.factor);
            temp = obj.dimsSmall*obj.factor;
            downsamplingOp.matrix =downsamplingOp*createSparseMatrixFrom1dSeperableKernel([0.1 0.2 0.4 0.2 0.1], [0.1 0.2 0.4 0.2 0.1], temp(1),temp(2), 'Neumann');
            
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

                %extract flow field between i and i+1
                singleField = squeeze(obj.v(:,:,i,:));

                %create warping operator forward and backward
                warpF = warpingOperator(obj.dimsLarge,0.5*singleField);
                warpB = warpingOperator(obj.dimsLarge,-0.5*singleField);

                %find out of range warps in each of the operators and set the
                %corresponding line in the other operator also to zero
                idF = sum(abs(warpF),2) == 0;
                idB = sum(abs(warpB),2) == 0;

                marker = idF + idB;
                marker = marker > 0;

                warpF(marker,:) = 0;
                warpB(marker,:) = 0;

                %add warping operator to flexbox
                obj.mainU.addTerm(L1operatorAniso(1,2,{warpB,-warpF}),[i,i+1]);
                obj.numWarpTerm(i) = numel(obj.mainU.duals);
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
                singleField = squeeze(obj.v(:,:,j,:));

                warpF = warpingOperator(obj.dimsLarge,0.5*singleField);
                warpB = warpingOperator(obj.dimsLarge,-0.5*singleField);

                idF = sum(abs(warpF),2) == 0;
                idB = sum(abs(warpB),2) == 0;

                marker = idF + idB;
                marker = marker > 0;

                warpF(marker,:) = 0;
                warpB(marker,:) = 0;

                obj.mainU.duals{obj.numWarpTerm(j)}.operator{1} = warpB;
                obj.mainU.duals{obj.numWarpTerm(j)}.operatorT{1} = warpB';

                obj.mainU.duals{obj.numWarpTerm(j)}.operator{2} = -warpF;
                obj.mainU.duals{obj.numWarpTerm(j)}.operatorT{2} = -warpF';
            end
        end
           
        function solveV(obj)
%             for j=1:obj.numFrames
%                 obj.u(:,:,j) = obj.mainU.getPrimal(j);
%             end
%             
            obj.mainV.resetImages(obj.u);
            
            if (~obj.vExists)
                obj.mainV.runPyramid;
                obj.vExists = 1;
            else
                obj.mainV.runLevel(1);
            end
            
            obj.v = obj.mainV.getResult;
            %[obj.v,~] = motionEstimationPyramid(obj.u,size(obj.u),1e-5,obj.beta,'L1TVOpticalFlowNonlinear',7,'adjustStepsize',0,'gradientConstancy',0);
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
            for i=1:obj.numMainIt
                %% solve u problem
                useFlexboxForUupdate = false;
                updateBlurKernels = true;
                
                disp('Solving problem for u');
                if useFlexboxForUupdate
                    obj.solveU;
                else
                    regParam = 0.02 + (obj.numMainIt-i)*0.04;
                    warpingTermParam = 0.1;
                    useConvexModel = true;
                    obj.u = solveSuperresolutionWithProst(obj,regParam,warpingTermParam,useConvexModel);
                end
                

                %% solve problem v
                disp('Solving problem for v');
                obj.solveV;

                %% update warping operators in flexBox
                disp('Updating operators in flexBox');
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


