function [genesigs,samplesigs] = nnmf_unmix_quick(Data,factors)

warning('off','all')

stream0 = RandStream('mt19937ar','Seed',0);

%display('nnmf');
opt = statset('MaxIter',10,'Display','off');
[genesigs,samplesigs] = nnmf(Data,factors,'options',opt,'replicates',20,'algorithm','mult');
opt = statset('MaxIter',1000,'Display','off');
[genesigs,samplesigs] = nnmf(Data,factors,'options',opt,'w0',genesigs,'h0',samplesigs);
m = norm(Data);
%display(sprintf(' Data MSE          : %2.2f',m))
e = norm(genesigs*samplesigs-Data);
%display(sprintf(' Reconstruction MSE: %2.2f',e))

for i = size(genesigs,2):-1:1;
	if sum(genesigs(:,i)~=0)==0
		genesigs(:,i) = [];
		samplesigs(i,:) = [];
	end
end


