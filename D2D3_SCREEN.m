function D3Res = D2D3_SCREEN(DesignInfo, AlgInfo)

	DesignInfo.CandInd = setdiff(1:(DesignInfo.nRuns - 1), unique(DesignInfo.Core));
	
	% Screening Step
	A = (DesignInfo.X2'*DesignInfo.Hadamard(:, DesignInfo.CandInd)).^2;
	colLen = sum(A,1);
	
	Screen = setdiff(find(colLen == max(colLen)), find(colLen == min(colLen)));
	DesignInfo.D3FixInd = DesignInfo.CandInd(Screen);
	DesignInfo.M3Fix = length(DesignInfo.D3FixInd);
	
	% Three scenarios after Screening Step
	if DesignInfo.M3 - DesignInfo.M3Fix <= 0, % D3 is fully identified
		
		DesignInfo.D3FixInd = DesignInfo.D3FixInd(randperm(DesignInfo.M3Fix, DesignInfo.M3));
		D3Res.GlobalHist = D_obj(DesignInfo.D3FixInd, DesignInfo.X2, DesignInfo.Hadamard);
		D3Res.MixCounter = [];
		D3Res.LocalHist = [];
		D3Res.AllHist = [];
		D3Res.ind.D2 = setdiff(DesignInfo.CandInd, DesignInfo.D3FixInd);
		D3Res.ind.D3 = DesignInfo.D3FixInd;
		
	elseif DesignInfo.M3 - DesignInfo.M3Fix == 1, % One column of D3 have not been identified
		
		DesignInfo.CandInd = setdiff(DesignInfo.CandInd, DesignInfo.D3FixInd);
		
		Dvals = zeros(1, (DesignInfo.M2 + DesignInfo.M3 - DesignInfo.M3Fix));
		for i = 1:length(DesignInfo.CandInd),
			Dvals(i) = D_obj([DesignInfo.D3FixInd, DesignInfo.CandInd(i)], DesignInfo.X2, DesignInfo.Hadamard);	
		end
		[BestD, loc] = max(Dvals);
		D3Augment = DesignInfo.CandInd(loc);
		
		D3Res.GlobalHist = BestD;
		D3Res.LocalHist = [];
		D3Res.AllHist = [];
		D3Res.MixCounter = [];
		
		D3Res.ind.D2 = setdiff(DesignInfo.CandInd, DesignInfo.CandInd(loc));
		D3Res.ind.D3 = [DesignInfo.D3FixInd, DesignInfo.CandInd(loc)];
		
	else % Two more columns of D3 have not been identified
		
		DesignInfo.CandInd = setdiff(DesignInfo.CandInd, DesignInfo.D3FixInd);
		D3Res = D2D3_SIBSSD(DesignInfo, AlgInfo);
		
	end


% Subfunctions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Dval = D_obj(D3ind, X2, Hadmat)
	
	Dval = round(det(X2'*Hadmat(:, D3ind)*Hadmat(:, D3ind)'*X2));
	
return