function [mixCore, mixDval] = D1_MIX_Rand(A, q, DesignInfo, fixInd, mixType)
	
	% Avoid the number of changed column over nTwoFI
	if q > DesignInfo.nTwoFI,
		q = DesignInfo.nTwoFI;
	end
	
	TAR = unique(A(randperm(DesignInfo.nTwoFI, q),:));
	if size(A,1) == 1,
		RandPool = setdiff(1:(DesignInfo.nRuns-1), [unique(A), fixInd]);
	else,
		RandPool = setdiff(1:(DesignInfo.nRuns-1), [unique(A)', fixInd]);
	end
	ColInd_A = randperm(length(RandPool), length(TAR));
	
	mixCore = A;
	for j = 1:length(TAR),
		mixCore(A == TAR(j)) = RandPool(ColInd_A(j));
	end
	%mixCore
	[X2, Qual] = X2Qualifier(mixCore, DesignInfo.Hadamard, DesignInfo.nRuns, DesignInfo.nTwoFI, 2);
	if mixType == 2, mixDval = D_obj(fixInd, X2, DesignInfo.Hadamard); end
	if mixType == 3, mixDval = D_obj2(mixCore, fixInd, X2, DesignInfo.Hadamard, DesignInfo.nFactors); end

% Subfunctions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X2, Qual] = X2Qualifier(Core, Hadmat, n, e, Qadj)
	
	X2 = zeros(n, e);
	for j = 1:e,
		X2(:,j) = Hadmat(:, Core(j,1)).*Hadmat(:, Core(j,2)); 
	end
	C = abs([X2'*Hadmat(:, unique(Core)), X2'*X2 - diag(n*ones(1, e))]);
	Qual = length(unique(C)) - Qadj;
	
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Dval = D_obj(D3ind, X2, Hadmat)
	
	Dval = round(det(X2'*Hadmat(:, D3ind)*Hadmat(:, D3ind)'*X2));
	
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Dval = D_obj2(Core, D2ind, X2, Hadmat, m)
	
	if size(Core, 1) == 1,
	mat = horzcat(ones(size(Hadmat,1),1), Hadmat(:, [unique(Core), D2ind]), X2);
	else
	mat = horzcat(ones(size(Hadmat,1),1), Hadmat(:, [unique(Core)', D2ind]), X2);
	end
	Dval = round(det(mat'*mat)/(size(Hadmat,1)^(m-size(X2,2)+1)));
	
return