function [mixMat, mixDval] = D2D3_MIX_Rand(A, q, DesignInfo)
	
	% Avoid the number of changed column over M3
	if q >= DesignInfo.M3 - DesignInfo.M3Fix,
		q = DesignInfo.M3 - DesignInfo.M3Fix;
	end
	if q >= DesignInfo.M2,
		q = DesignInfo.M2;
	end
	
	% randi() needs statistical toolbox, so we don't use it
	RandPool = setdiff(DesignInfo.CandInd, A);
	ColInd_A = randperm(DesignInfo.M3 - DesignInfo.M3Fix, q);
	ColInd_Pool = randperm(DesignInfo.M2, q);
	
	A(ColInd_A) = RandPool(ColInd_Pool);
	
	mixMat = A;
	mixDval = D_obj(A, DesignInfo.X2, DesignInfo.Hadamard);


% Subfunctions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Dval = D_obj(D3ind, X2, Hadmat)
	
	Dval = round(det(X2'*Hadmat(:, D3ind)*Hadmat(:, D3ind)'*X2));
	
return