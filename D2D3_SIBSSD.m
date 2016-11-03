function D3Res = D2D3_SIBSSD(DesignInfo, AlgInfo)
	
	% Construct historical record of D value, local best D value and global D value for each AlgInfo.D3.Swarm and iteration
	D3Res.AllHist = zeros(AlgInfo.D3.Iter+1, AlgInfo.D3.Swarm);
	D3Res.LocalHist = zeros(AlgInfo.D3.Iter+1, AlgInfo.D3.Swarm);
	D3Res.GlobalHist = zeros(1, AlgInfo.D3.Iter+1);
	% MixCounter is a record, for each iteration, for the number of Swarms move to
	% Column 1: mixGB
	% Column 2: mixLB
	% Column 3: mixRC
	D3Res.MixCounter = zeros(AlgInfo.D3.Iter, 3);
	
	% Initialize the Swarm from the rest columns of DesignInfo.Hadamard matrix
	D3_set = zeros(AlgInfo.D3.Swarm, DesignInfo.M3 - DesignInfo.M3Fix);
	DvalRaw = zeros(1, AlgInfo.D3.Swarm);
	
	for i = 1:AlgInfo.D3.Swarm,
		D3_ColInd = randperm(DesignInfo.M2 + DesignInfo.M3 - DesignInfo.M3Fix, DesignInfo.M3 - DesignInfo.M3Fix);
		D3_set(i,:) = DesignInfo.CandInd(D3_ColInd);
		% Compute D value for each D2 combining with D1 and X2
		DvalRaw(i) = D_obj([DesignInfo.CandInd(D3_ColInd), DesignInfo.D3FixInd], DesignInfo.X2, DesignInfo.Hadamard); 
	end
	
	% Save the information for PSO process
	D3Res.AllHist(1,:) = DvalRaw;
	
	D3Res.LB_set = D3_set;
	D3Res.LocalHist(1,:) = DvalRaw;

	[BestVal, loc] = max(DvalRaw);
	D3Res.GB_mat = D3_set(loc,:);
	D3Res.GlobalHist(1) = BestVal;

	% Set Space for Moved-D3
	New_D3_set = D3_set;
	
	% Start DPSO updating
	for ttD3 = 1:AlgInfo.D3.Iter,
				
		New_Dvals = zeros(1,AlgInfo.D3.Swarm);
		% Mix Step
		parfor i = 1:AlgInfo.D3.Swarm,	
			% Mix with LB
			[mixLB_temp, D_LB_temp] = D2D3_MIX_AB(D3_set(i,:), D3Res.LB_set(i,:), AlgInfo.D3.qLB, DesignInfo);
			DmixLB(i) = D_LB_temp;
			mixLB(i,:) = mixLB_temp;
			% Mix with GB
			[mixGB_temp, D_GB_temp] = D2D3_MIX_AB(D3_set(i,:), D3Res.GB_mat, AlgInfo.D3.qGB, DesignInfo);
			DmixGB(i) = D_GB_temp;
			mixGB(i,:) = mixGB_temp;	
		end
		
		% Move Step
		EXwGB = (DmixGB > DmixLB) & (DmixGB > DvalRaw);
		EXwLB = (DmixLB > DvalRaw) & ((DmixGB < DmixLB) | (DmixGB < DvalRaw));
		EXwRC = find(ones(1,AlgInfo.D3.Swarm) - EXwGB - EXwLB);
		
		New_D3_set(EXwGB,:) = mixGB(EXwGB,:);
		New_Dvals(EXwGB) = DmixGB(EXwGB);
		D3Res.MixCounter(ttD3, 1) = sum(EXwGB);
		
		New_D3_set(EXwLB,:) = mixLB(EXwLB,:);
		New_Dvals(EXwLB) = DmixLB(EXwLB);
		D3Res.MixCounter(ttD3, 2) = sum(EXwLB);
		
		for i = 1:AlgInfo.D3.Swarm,
			if any(i == EXwRC),
				[mixRC, DmixRC] = D2D3_MIX_Rand(D3_set(i,:), AlgInfo.D3.qRC, DesignInfo);
				New_D3_set(i,:) = mixRC;
				New_Dvals(i) = DmixRC;
			end
		end
		
		% Update LB
		D3Res.LocalHist(ttD3+1,:) = D3Res.LocalHist(ttD3,:);
		UpdateLB = New_Dvals > D3Res.LocalHist(ttD3+1,:);
		D3Res.LocalHist(ttD3+1, UpdateLB) = New_Dvals(UpdateLB);
		D3Res.LB_set(UpdateLB,:) = New_D3_set(UpdateLB,:);
		% Update GB
		[BestVal, loc] = max(D3Res.LocalHist(ttD3+1,:));
		D3Res.GB_mat = D3Res.LB_set(loc,:);
		D3Res.GlobalHist(ttD3+1) = BestVal;
		% Save the history value
		D3Res.AllHist(ttD3+1,:) = New_Dvals;
		D3Res.MixCounter(ttD3, 3) = AlgInfo.D3.Swarm - D3Res.MixCounter(ttD3, 1) - D3Res.MixCounter(ttD3, 2);
		% Update Swarms
		D3_set = New_D3_set;
		DvalRaw = New_Dvals;
	end
	
	D3Res.ind.D2 = setdiff(DesignInfo.CandInd, D3Res.GB_mat);
	D3Res.ind.D3 = [D3Res.GB_mat DesignInfo.D3FixInd];

	
% Subfunctions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Dval = D_obj(D3ind, X2, Hadmat)
	
	Dval = round(det(X2'*Hadmat(:, D3ind)*Hadmat(:, D3ind)'*X2));
	
return
	