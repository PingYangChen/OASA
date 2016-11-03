function [mixMat, mixDval] = D2D3_MIX_AB(A, B, q, DesignInfo)
	
	% Avoid the number of changed column over M3
	if q >= DesignInfo.M3 - DesignInfo.M3Fix,
		q = DesignInfo.M3 - DesignInfo.M3Fix;
	end
	
	% Deletion Process
	for qInd = 1:q,
		Dvals = zeros(1, (DesignInfo.M3 - DesignInfo.M3Fix - qInd + 1));
		% Delete the i-th index of A, evaluate the D value
		A_minus = zeros((DesignInfo.M3 - DesignInfo.M3Fix - qInd + 1), (DesignInfo.M3 - DesignInfo.M3Fix - qInd));
		for i = 1:(DesignInfo.M3 - DesignInfo.M3Fix - qInd + 1),
			DeleteInd = setdiff(1:(DesignInfo.M3 - DesignInfo.M3Fix - qInd + 1), i);
			A_minus(i,:) = A(DeleteInd);
			Dvals(i) = D_obj([A(DeleteInd), DesignInfo.D3FixInd], DesignInfo.X2, DesignInfo.Hadamard);
		end
		% Find the best index of A making the best D for A(:,-i)
		[BestVal, loc] = max(Dvals);
		% Delete the best index then go next deletion process
		A = A_minus(loc,:);
	end
	
	% Addition Process
	for qInd = 1:q,
		Dvals = zeros(1, (DesignInfo.M3 - DesignInfo.M3Fix - qInd + 1));
		% Insert the i-th index of B, evaluate the D value
		for i = 1:(DesignInfo.M3 - DesignInfo.M3Fix - qInd + 1),
			if sum(A == B(i)) > 0,
				Dvals(i) = 0;
			else
				Dvals(i) = D_obj([A, B(i), DesignInfo.D3FixInd], DesignInfo.X2, DesignInfo.Hadamard);
			end
		end
		% Find the best index of B making the best D for [R, B(:,i)]
		[BestVal, loc] = max(Dvals);
		% Add the best index then go next addition process
		A = [A, B(loc)];
		B = B(setdiff(1:(DesignInfo.M3 - DesignInfo.M3Fix - qInd + 1), loc));
	end
	
	mixMat = A;
	mixDval = BestVal;
	

% Subfunctions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Dval = D_obj(D3ind, X2, Hadmat)
	
	Dval = round(det(X2'*Hadmat(:, D3ind)*Hadmat(:, D3ind)'*X2));
	
return