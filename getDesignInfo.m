function DesignInfo = getDesignInfo()
	
	DesignInfo = struct( ...
	...
		'Hadamard',		[], 		... % The Hadamard matrix
		'nRuns',		[], 			... % The number of run size of design
		'nFactors',		[], 		... % The number of factors of design
		'CoreFig',		[] 			... % The matrix of two factor interactions
	);
	