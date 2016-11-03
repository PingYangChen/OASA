
function D1Res = D1_SIBSSD(DesignInfo, AlgInfo, D3Res, mixType)

	if mixType == 2, 
		fixInd = D3Res.ind.D3; 
		fixNum = DesignInfo.M3;
		AlgInfo.D1 = AlgInfo.D1f3;
	end
	if mixType == 3, 
		fixInd = D3Res.ind.D2; 
		fixNum = DesignInfo.M2;
		AlgInfo.D1 = AlgInfo.D1f2;
	end
	CoreVertex 	= unique(DesignInfo.CoreFig);
	
	D1Res.AllHist = zeros(AlgInfo.D1.Iter + 1, AlgInfo.D1.Swarm);
	D1Res.LocalHist = zeros(AlgInfo.D1.Iter + 1, AlgInfo.D1.Swarm);
	D1Res.GlobalHist = zeros(1, AlgInfo.D1.Iter + 1);
	D1Res.MixCounter = zeros(AlgInfo.D1.Iter, 3);
	% D1 Updating Procedure
	% Initialize the Swarm from the rest columns of DesignInfo.Hadamard matrix
	Core_set  = zeros([size(DesignInfo.CoreFig), AlgInfo.D1.Swarm]);
	DvalRaw = zeros(1, AlgInfo.D1.Swarm);

	CandInd = setdiff(1:(DesignInfo.nRuns - 1), fixInd);
	for i = 1:AlgInfo.D1.Swarm, 
		D1_ColInd = randperm(DesignInfo.nRuns - 1 - fixNum, DesignInfo.M1);
		C_temp = DesignInfo.CoreFig;
		for j = 1:DesignInfo.M1,
			C_temp(DesignInfo.CoreFig == CoreVertex(j)) = CandInd(D1_ColInd(j));  
		end
		Core_set(:,:,i) = C_temp;
		[X2_set(:,:,i), Qual] = X2Qualifier(C_temp, DesignInfo.Hadamard, DesignInfo.nRuns, DesignInfo.nTwoFI, 2);
		if mixType == 2, DvalRaw(i) = D_obj(fixInd, X2_set(:,:,i), DesignInfo.Hadamard); end
		if mixType == 3, DvalRaw(i) = D_obj2(C_temp, fixInd, X2_set(:,:,i), DesignInfo.Hadamard, DesignInfo.nFactors); end
	end
			
	% Save the information for PSO process
	D1Res.AllHist(1,:) = DvalRaw;
	
	D1Res.LB_set = Core_set;
	D1Res.LocalHist(1,:) = DvalRaw;
	[BestVal, loc] = max(DvalRaw);
	D1Res.GB_mat = Core_set(:,:,loc);
	D1Res.GlobalHist(1) = BestVal;
	% Set Space for Moved-core
	New_Core_set = Core_set;
	
	for ttD1 = 1:AlgInfo.D1.Iter,
	
		New_Dvals = zeros(1, AlgInfo.D1.Swarm);
		% Mix Step
		parfor i = 1:AlgInfo.D1.Swarm,	
			% Mix with LB
			[mixLB_temp, D_LB_temp] = D1_MIX_AB(Core_set(:,:,i), D1Res.LB_set(:,:,i), AlgInfo.D1.qLB, DesignInfo, fixInd, mixType);
			DmixLB(i) = D_LB_temp;
			mixLB(:,:,i) = mixLB_temp;
			% Mix with GB
			[mixGB_temp, D_GB_temp] = D1_MIX_AB(Core_set(:,:,i), D1Res.GB_mat, AlgInfo.D1.qGB, DesignInfo, fixInd, mixType);
			DmixGB(i) = D_GB_temp;
			mixGB(:,:,i) = mixGB_temp;	
		end

		% Move Step
		EXwGB = (DmixGB > DmixLB) & (DmixGB > DvalRaw);
		EXwLB = (DmixLB > DvalRaw) & ((DmixGB < DmixLB) | (DmixGB < DvalRaw));
		EXwRC = find(ones(1,AlgInfo.D1.Swarm) - EXwGB - EXwLB);
		
		New_Core_set(:,:,EXwGB) = mixGB(:,:,EXwGB);
		New_Dvals(EXwGB) = DmixGB(EXwGB);
		D1Res.MixCounter(ttD1, 1) = sum(EXwGB);
		
		New_Core_set(:,:,EXwLB) = mixLB(:,:,EXwLB);
		New_Dvals(EXwLB) = DmixLB(EXwLB);
		D1Res.MixCounter(ttD1, 2) = sum(EXwLB);
		
		for i = 1:AlgInfo.D1.Swarm,
			if any(i == EXwRC),
				[mixRC, DmixRC] = D1_MIX_Rand(Core_set(:,:,i), AlgInfo.D1.qRC, DesignInfo, fixInd, mixType);
				New_Core_set(:,:,i) = mixRC;
				New_Dvals(i) = DmixRC;
			end
		end

		% Update LB
		D1Res.LocalHist(ttD1+1,:) = D1Res.LocalHist(ttD1,:);
		UpdateLB = New_Dvals > D1Res.LocalHist(ttD1+1,:);
		D1Res.LocalHist(ttD1+1, UpdateLB) = New_Dvals(UpdateLB);
		D1Res.LB_set(:,:,UpdateLB) = New_Core_set(:,:,UpdateLB);
		% Update GB
		[BestVal, loc] = max(D1Res.LocalHist(ttD1+1,:));
		D1Res.GB_mat = D1Res.LB_set(:,:,loc);
		D1Res.GlobalHist(ttD1+1) = BestVal;
		% Save the history value
		D1Res.AllHist(ttD1+1,:) = New_Dvals;
		D1Res.MixCounter(ttD1, 3) = AlgInfo.D1.Swarm - D1Res.MixCounter(ttD1, 1) - D1Res.MixCounter(ttD1, 2);
		% Update Swarms
		Core_set = New_Core_set;
		DvalRaw = New_Dvals;
		
	end
	
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