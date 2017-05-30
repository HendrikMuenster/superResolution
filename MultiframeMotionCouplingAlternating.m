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
       
    end
    
    %% Methods: 
    methods
        
        %% decorate init
        function init(obj)
            
            if ~strcmp(obj.framework,'flexBox')
                error('Alternating only in flexBox for now');
            end
            init@MultiframeMotionCoupling(obj);
        end
        
        %% solve high-res optical flow problem
        % is actually init_v, but applied to u instead of imageSequenceSmall
        % todo for rewrite ...
        function solveV(obj) 
            if (obj.verbose > 0)
                disp('Calculating initial velocity fields')
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
            
            % keep h fixed for new iterations !?!
            % h = obj.h
        end
        
        
        % run multiple iterations
        function run(obj)
            
        end
    end
end



