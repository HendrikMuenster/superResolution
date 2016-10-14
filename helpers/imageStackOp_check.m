%
% imageStackOpSanity
%
clearvars

u = randi(10,1000,1000,3);
k = [0,1,0;1,4,1;0,1,0];
Ku = imfilter(u,k,'replicate','conv');
U_sparse = imageStackOp(u,3);

 
Kuh = U_sparse*[k(:);k(:);k(:)];

sum(abs((Kuh-Ku(:))))


%%% check also writeKernel...

K = writeKernelToSparseDownsamplingMatrix(k,1,1000,1000);
K = blkdiag(K,K,K);

D  = superpixelOperator([250,250],4).matrix;
D = blkdiag(D,D,D);
Ku2 = D*K*u(:);
DKu = D*Ku(:);
sum(abs((Ku2-DKu)))


%%

DK = [];
for i = 1:3
    tmpOp = writeKernelToSparseDownsamplingMatrix(k,1,1000,1000);
    tmpOp = superpixelOperator([250,250],4).matrix*tmpOp;
    DK    = blkdiag(DK,tmpOp);
end

Ku3 = DK*u(:);
sum(abs((Ku2-Ku3)))

%%
DK2 = writeKernelToSparseDownsamplingMatrix(k,1,200,250);
%%
k = [0,1,0;1,4,1;0,1,0];
u = rand(250,500);
tmpOp = convmtx2(k,25,500);
blerg = conv2(u,k)';
%%%

%A(i-1,i):A(i+1,i) = k(1:3,2)
%A(i-m-1,i):A(i-m+1,i) = k(1:3,1)
%A(i+m-1,i):A(i+m+1,i) = k
