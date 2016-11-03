function [bound, bdType] = getBound(D3Res, DesignInfo)

	%DesignInfo.nTwoFI		= size(DesignInfo.CoreFig, 1);
	%DesignInfo.M1			  = length(unique(DesignInfo.CoreFig));
	%DesignInfo.M2 			= DesignInfo.nFactors - length(unique(DesignInfo.CoreFig));
	%DesignInfo.M3 			= DesignInfo.nRuns - DesignInfo.nFactors - 1;
	
	%D1 = DesignInfo.Hadamard(:,unique(DesignInfo.Core));
	%D2 = DesignInfo.Hadamard(:,D3Res.ind.D2);
	%D3 = DesignInfo.Hadamard(:,D3Res.ind.D3);
	[X2, Qual] = X2Qualifier(DesignInfo.Core, DesignInfo.Hadamard, DesignInfo.nRuns, DesignInfo.nTwoFI, 2);
	Amat = DesignInfo.Hadamard(:, D3Res.ind.D3)'*X2;
	Alen = sum(Amat .* Amat, 1);
	
	if mod(DesignInfo.nRuns, 8) == 4,
		LVAL = ((DesignInfo.nRuns/4)^2 - (DesignInfo.nRuns - DesignInfo.M3 - 3));
		if ((DesignInfo.nTwoFI < 4) && (all(Alen == 16*LVAL))),
			bdType = 'Sharp';
			if DesignInfo.nTwoFI == 1,
				bound = 16*LVAL;
			end
			if DesignInfo.nTwoFI == 2,
				check = mod(DesignInfo.M3, 2); 
				switch check,
					case 0
						bound = 16*16*LVAL*LVAL;
					case 1
						bound = 16*16*(LVAL*LVAL - 1);
				end
			end
			if DesignInfo.nTwoFI == 3,
				check = mod(DesignInfo.M3, 4); 
				switch check,
					case 0
						bound = 16*16*16*LVAL*LVAL*LVAL; 
					case 1
						bound = 16*16*16*(LVAL*LVAL*LVAL - 3*LVAL + 2); 
					case 2
						bound = 16*16*16*(LVAL*LVAL*LVAL - 4*LVAL);
					case 3
						bound = 16*16*16*(LVAL*LVAL*LVAL - 3*LVAL - 2); 
				end
			end
		else,
			bdType = 'General';
			bound = (16*LVAL)^DesignInfo.nTwoFI;
		end
	elseif mod(DesignInfo.nRuns, 8) == 0,
		bdType = 'General';
		bound = DesignInfo.nRuns^(2*DesignInfo.nTwoFI);
	else
		error('The run size should be a multiple of 4 \n');
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

