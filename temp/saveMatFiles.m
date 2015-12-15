
version = 'act';
ind = 450;
trajLen = int2str(ind);
FinvVec = ilc.Finv(:);
save(['Finv_',version,'_',trajLen,'.txt'],'FinvVec','-ascii');
M = zeros(size(ilc.Ql));
dimu = 7;
M(end-2*dimu+1:end,end-2*dimu+1:end) = 1;
Mat = pinv(ilc.F' * M * ilc.F, 0.05) * (ilc.F' * M);
FinvmayerVec = Mat(:);
save(['Finv_',version,'_mayer_',trajLen,'.txt'],'FinvmayerVec', '-ascii');
Hvec = ilc.H(:);
save(['H_',version,'_',trajLen,'.txt'],'Hvec','-ascii');