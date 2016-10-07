function v = solveSuperresolutionWithProst(obj,regParam,warpingTermParam,useConvexModel)

temp = obj.dimsSmall;
ny=temp(1); nx=temp(2);
temp = obj.dimsLarge;
Ny = temp(1);
Nx = temp(2);
nc =  obj.numFrames;

subSaplingoperator = obj.mainU.duals{1}.operator{1};
f = obj.imageSequenceSmall;

temp = obj.mainU.duals{5}.operator{1};
warpingOperator = sparse([],[],[],Nx*Ny*(nc-1),Nx*Ny*nc,2*nnz(temp)*nc); 
for i=1:nc-1
    warpingOperator((Nx*Ny*(i-1)+1):(Nx*Ny*i), (Nx*Ny*(i-1)+1):(Nx*Ny*i)) = obj.mainU.duals{2*nc+i}.operator{1};
    warpingOperator((Nx*Ny*(i-1)+1):(Nx*Ny*i), (Nx*Ny*i+1):(Nx*Ny*(1+i))) = obj.mainU.duals{2*nc+i}.operator{2};
end
warpingOperator = warpingTermParam*warpingOperator;

%primal variable
u = prost.variable(Nx*Ny*nc);

%dual variables 
p = prost.variable(nx*ny*nc); %dual variable for fitting to data
q = prost.variable(Nx*Ny*(nc-1)); %dual variable for warping between data
g = prost.variable(2*Nx*Ny*nc);% dual variable for TV regularization
g2 = prost.variable(2*Nx*Ny*nc); %dual variable for additional quadratic regularization

prob = prost.min_max_problem( {u}, {q,p,g,g2} );

prob.add_dual_pair(u, p, prost.block.id_kron_sparse(subSaplingoperator,nc)); 
prob.add_dual_pair(u, q, prost.block.sparse(warpingOperator)); 
prob.add_dual_pair(u, g, prost.block.gradient2d(Nx,Ny,nc,false)); 
prob.add_dual_pair(u, g2, prost.block.gradient2d(Nx,Ny,nc,false)); 

prob.add_function(p, prost.function.sum_1d('ind_box01', 0.5, -0.5, 1, f(:), 0)); %ell^1
prob.add_function(q, prost.function.sum_1d('ind_box01', 0.5, -0.5, 1, 0, 0)); %ell^1
%prob.add_function(q, prost.function.sum_1d('square', 1, 0, 1, 0,
%0));%ell^2-squared
prob.add_function(g2, prost.function.sum_1d('square', 1, 0, 10, 0, 0));%ell^2-squared
if useConvexModel
    callbacks = 5;
	prob.add_function(g, prost.function.sum_norm2(2, false, 'ind_leq0', 1/regParam, 1, 1, 0, 0)); %ell^{2,1}
     backend = prost.backend.pdhg(...
        'tau0', 100, ...
        'sigma0', 0.01, ...
        'stepsize', 'boyd');
else
    callbacks = 0;
    alpha = 0.5;
    prob.add_function(g, ...
                      prost.function.conjugate(...
                          prost.function.sum_norm2(... 
                              2, false, 'lq', ...
                              1,0,regParam,0,0, alpha))); % l^\alpha

    % specify solver options
    backend = prost.backend.pdhg('stepsize', 'alg2', ...
                             'residual_iter', -1, ...
                             'alg2_gamma', 0.1);
end
                         

opts = prost.options('max_iters', 500, ...
                     'num_cback_calls', callbacks, ...
                     'verbose', true);

tic;
result = prost.solve(prob, backend, opts);
toc;
v = reshape(u.val, [Ny,Nx,nc]);
