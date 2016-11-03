function [mixCore, mixDval] = D1_MIX_AB(A, B, q, DesignInfo, fixInd, mixType)
	
	% Avoid the number of changed column over nTwoFI
	if q > DesignInfo.nTwoFI,
		q = DesignInfo.nTwoFI;
	end
	
	[X2, Qual] = X2Qualifier(A, DesignInfo.Hadamard, DesignInfo.nRuns, DesignInfo.nTwoFI, 2);
	
	% Deletion Process
	Dvals = zeros(1, (DesignInfo.nTwoFI));
	for i = 1:DesignInfo.nTwoFI,
		if mixType == 2, 
			Dvals(i) = D_obj(fixInd, X2(:,setdiff(1:DesignInfo.nTwoFI, i)), DesignInfo.Hadamard); 
		end
		if mixType == 3, 
			Dvals(i) = D_obj2(A(setdiff(1:size(A,1),i),:), fixInd, X2(:,setdiff(1:DesignInfo.nTwoFI, i)), DesignInfo.Hadamard, DesignInfo.nFactors - 1); 
		end
	end
	[val, loc] = sort(Dvals, 'descend'); % 20140921
	TAR = unique(A(loc(1:q), :));
	
	% Addition Process
	B_uni = unique(B);
	B_comb = nchoosek(B_uni(randperm(DesignInfo.M1)), length(TAR));
	parfor i = 1:size(B_comb,1);
		A_temp = A;
		for j = 1:length(TAR),
			A_temp(A == TAR(j)) = B_comb(i,j);
		end
		if length(unique(A_temp)) < DesignInfo.M1,
			mixDvals(i) = 0;
		else,
			[X2, Qual] = X2Qualifier(A_temp, DesignInfo.Hadamard, DesignInfo.nRuns, DesignInfo.nTwoFI, 2);
			if mixType == 2, mixDvals(i) = D_obj(fixInd, X2, DesignInfo.Hadamard); end
			if mixType == 3, mixDvals(i) = D_obj2(A_temp, fixInd, X2, DesignInfo.Hadamard, DesignInfo.nFactors); end
		end
	end
	
	[mixDval, loc] = max(mixDvals);
	mixCore = A;
	for j = 1:length(TAR),
		mixCore(A == TAR(j)) = B_comb(loc,j);
	end
	%mixCore
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