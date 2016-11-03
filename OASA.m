function output = OASA(DesignInfo, AlgInfo)
	
	% In This Version:
	%
	
	fprintf('--------------------------------- \n');
	fprintf('Orthogonal Array Search Algorithm \n');
	% Start recording time
	tic;
	D3Res = [];
	% Set design informations 
	DesignInfo.nTwoFI		= size(DesignInfo.CoreFig, 1);
	DesignInfo.M1				= length(unique(DesignInfo.CoreFig));
	DesignInfo.M2 			= DesignInfo.nFactors - length(unique(DesignInfo.CoreFig));
	DesignInfo.M3 			= DesignInfo.nRuns - DesignInfo.nFactors - 1;
	
	fprintf(strcat('Runs: ', num2str(DesignInfo.nRuns), ...
								 '; Factors: ', num2str(DesignInfo.nFactors), '; 2fi Str: \n'));
	DesignInfo.CoreFig
	
	output.DvalHist = zeros(1, AlgInfo.Iter);
	output.GlobalHist = zeros(1, AlgInfo.Iter);
	output.StopAt = AlgInfo.Iter;
	output.AcceptRate = ones(1, AlgInfo.Iter);
	% Compute the Upper Bound for the Target Design
	% V is the number of two connected edges on the graph.
%	v = 0;
%	for i = 1:M1,
%		k = sum(CoreFig(:) == CoreVertex(i));
%		v = v + 0.5*k*(k-1);
%	end
		
	fprintf('Initializing... ');
	% Generate Qualified D1 and X2
	CoreVertex 	= unique(DesignInfo.CoreFig);
	
	DesignInfo.Core = DesignInfo.CoreFig;
	[DesignInfo.X2, Qual] = X2Qualifier(DesignInfo.Core, DesignInfo.Hadamard, DesignInfo.nRuns, DesignInfo.nTwoFI, 2);
	Qual = Qual*AlgInfo.X2Qual; % Decide to tune the first D1
	X2time = 0;
	Qadj = 2;
	while Qual > 0,
		X2time = X2time + 1;
		D1choose = randperm(DesignInfo.nRuns - 1);
		for j = 1:DesignInfo.M1,
			DesignInfo.Core(DesignInfo.CoreFig == CoreVertex(j)) = D1choose(j);  
		end
		[DesignInfo.X2, Qual] = X2Qualifier(DesignInfo.Core, DesignInfo.Hadamard, DesignInfo.nRuns, DesignInfo.nTwoFI, Qadj);
		if X2time >= AlgInfo.X2tune,
			Qadj = Qadj + 1;
			X2time = 0;
		end
	end
	
	% Run Approach
	D3Res = D2D3_SCREEN(DesignInfo, AlgInfo);	
	[output.Bound, output.bdType] = getBound(D3Res, DesignInfo);
	fprintf(['Bound Type: ', output.bdType, '\n']);
	fprintf('Efficiency: %d. \n', D3Res.GlobalHist(end)*100/output.Bound);
	
	saveBest = D3Res;
	saveBest.Core = DesignInfo.Core;
	
	% Save the initial state
	Init.Core = DesignInfo.Core;
	Init.D2   = D3Res.ind.D2;
	Init.Dval = D3Res.GlobalHist(end);
	
	if AlgInfo.Iter == 0,
		output.GlobalHist = saveBest.GlobalHist(end);
	end
	
	StopFlag = 0;
	if D3Res.GlobalHist(end) >= output.Bound,
		StopFlag = 1;
		output.StopAt = 0;
		output.GlobalHist = saveBest.GlobalHist(end);
		fprintf('Initial status obtain 100 %% efficiency.  Shutting down.\n');
	end
	
	tt = 0;

	% Start DPSO updating
	while (StopFlag == 0) && tt < AlgInfo.Iter,
	
		tt = tt + 1;
		
		%unique(DesignInfo.Core)'
		%D3Res.ind.D2
		%D3Res.ind.D3
		
		if DesignInfo.M2 > DesignInfo.M1,
			D1ResMixD2 = D1_SIBSSD(DesignInfo, AlgInfo, D3Res, 2);
		else, 
			D1ResMixD2.GlobalHist = 0;
		end
		
		if DesignInfo.M3 > DesignInfo.M1,
			D1ResMixD3 = D1_SIBSSD(DesignInfo, AlgInfo, D3Res, 3);
		else, 
			D1ResMixD3.GlobalHist = 0;
		end
		
		[val, loc] = max([D1ResMixD2.GlobalHist(end), D1ResMixD3.GlobalHist(end)]/D3Res.GlobalHist(end));
		
		if val > 1,
		
			switch loc,
			
				case 1,
					DesignInfo.Core = D1ResMixD2.GB_mat;
					[DesignInfo.X2, Qual] = X2Qualifier(DesignInfo.Core, DesignInfo.Hadamard, DesignInfo.nRuns, DesignInfo.nTwoFI, 2);
					
				case 2,% Switch with D3 index
					output.AcceptRate(tt) = 2;
					DesignInfo.Core = D1ResMixD3.GB_mat;	
					[DesignInfo.X2, Qual] = X2Qualifier(DesignInfo.Core, DesignInfo.Hadamard, DesignInfo.nRuns, DesignInfo.nTwoFI, 2);					
			end
			
		else, 
			
			% If None of D1ResMixD2 or D1ResMixD3 is better than the original D1, switch with random column 
			output.AcceptRate(tt) = 3; % 3 stand for switch with random columns, which need DPSO to find new best design
			DesignInfo.Core = DesignInfo.CoreFig;
			
			Qual = 1; 
			Qadj = 2;
			X2time = 0;
			if DesignInfo.M2 < DesignInfo.M1,
				nTemp = DesignInfo.M1 - DesignInfo.M2;
			end
			if DesignInfo.M3 < DesignInfo.M1,
				nTemp = DesignInfo.M3;
			end
			
			while Qual > 0,
			
				X2time 	= X2time + 1;
				if DesignInfo.M2 >= DesignInfo.M1 & DesignInfo.M3 >= DesignInfo.M1,
					nTemp 	= randperm(DesignInfo.M1 - 1, 1);
				end
				
				RandInd = [D3Res.ind.D3(randperm(DesignInfo.M3, nTemp)), ...
							D3Res.ind.D2(randperm(DesignInfo.M2, DesignInfo.M1 - nTemp))];
				RandInd = RandInd(randperm(length(RandInd)));
				for j = 1:DesignInfo.M1,
					DesignInfo.Core(DesignInfo.CoreFig == CoreVertex(j)) = RandInd(j);  
				end
				[DesignInfo.X2, Qual] = X2Qualifier(DesignInfo.Core, DesignInfo.Hadamard, DesignInfo.nRuns, DesignInfo.nTwoFI, Qadj);
				Qual = Qual*AlgInfo.X2Qual;
				if X2time >= AlgInfo.X2tune,
					Qadj = Qadj + 1;
					X2time = 0;
				end
			end
			
		end

		% update D3
		D3Res = D2D3_SCREEN(DesignInfo, AlgInfo);	
				
		if D3Res.GlobalHist(end) >= max([Init.Dval, output.DvalHist]),
			saveBest = D3Res;
			saveBest.Core = DesignInfo.Core;	
		end
			
		output.DvalHist(tt) = D3Res.GlobalHist(end);
		output.GlobalHist(tt) = max([Init.Dval, output.DvalHist]);
		
		if output.GlobalHist(tt) >= output.Bound,
			StopFlag = 1;
			output.StopAt = tt;
			fprintf('Obtain 100 %% efficiency. At: %d th iteration.  Shutting down. \n', tt);
		end		
		
	end	
	
	if output.StopAt == AlgInfo.Iter,
		fprintf('Obatain max. iteration. %d %% efficiency.  Shutting down.  \n', saveBest.GlobalHist(end).*100/output.Bound);
	end
	
	% plots
%	if output.StopAt > 0,
%		PlotSummary = figure(1);
%		plot(0, Init.Dval/output.Bound, 'og', 'MarkerFaceColor', 'g')%, 'MarkerEdgeColor', 'k');
%		hold on;
%		plot(1:output.StopAt, output.GlobalHist(1:output.StopAt)./output.Bound, 'o-b', 'MarkerFaceColor', 'b', 'MarkerSize', 3);
%		hold off;
%		ylim([min(cat(2, Init.Dval, output.DvalHist(1:output.StopAt))./output.Bound) 1]);
%		set(gca, 'YGrid', 'on');
%		%title('Overall','interpreter','latex');
%		xlabel('Iteration','interpreter','latex');
%		ylabel('Ratio to the upper bound of $\det\mathbf{A}^T\mathbf{A}$','interpreter','latex');
%		legend('Initial Status', 'DualDPSO', 'Location', 'Best');
%		if length(AlgInfo.save.Overall) > 0,
%			saveas(PlotSummary, [pwd, '/', AlgInfo.save.Overall, '.eps'], 'epsc');
%			saveas(PlotSummary, [pwd, '/', AlgInfo.save.Overall, '.fig'], 'fig');
%			saveas(PlotSummary, [pwd, '/', AlgInfo.save.Overall, '.png'], 'png');
%		end
%		
%	else,
%		PlotSummary = figure(1); plot(0);
%	end
%	if length(saveBest.GlobalHist) == AlgInfo.D3.Iter + 1,
%		PlotD3Res = figure(2);
%		plot(0:AlgInfo.D3.Iter, saveBest.GlobalHist./output.Bound, 'LineWidth', 3, 'Color', 'blue')
%		hold on;
%		plot(0:AlgInfo.D3.Iter, saveBest.LocalHist./output.Bound, 'LineStyle', '--', 'LineWidth', 1, 'Color', 'black')
%		hold off;
%		%title('The last D3 update procedure','interpreter','latex');
%		ylim([min([saveBest.GlobalHist, min(saveBest.LocalHist(:))])/output.Bound 1]);
%		xlabel('Iteration','interpreter','latex');
%		ylabel('Ratio to the upper bound of $\det\mathbf{A}^T\mathbf{A}$','interpreter','latex');
%		legend('GolbalBest','LocalBest', 'Location', 'Best')
%		if length(AlgInfo.save.D3Res) > 0,
%			saveas(PlotD3Res, [pwd, '/', AlgInfo.save.D3Res, '.eps'], 'epsc');
%			saveas(PlotD3Res, [pwd, '/', AlgInfo.save.D3Res, '.fig'], 'fig');
%			saveas(PlotD3Res, [pwd, '/', AlgInfo.save.D3Res, '.png'], 'png');
%		end
%	end

	
	% Save outputs
	if output.StopAt > 0,
		output.BestD 	= output.GlobalHist(output.StopAt);
	else,
		output.BestD 	= output.GlobalHist;
	end
	output.Core 	= saveBest.Core;
	output.BestOA 	= horzcat(DesignInfo.Hadamard(:,unique(output.Core)), DesignInfo.Hadamard(:,saveBest.ind.D2));
	
	[X2, Qual] 		= X2Qualifier(output.Core, DesignInfo.Hadamard, DesignInfo.nRuns, DesignInfo.nTwoFI, 2);
	output.mat 		= struct('D1', DesignInfo.Hadamard(:,unique(output.Core)), 'D2', DesignInfo.Hadamard(:,saveBest.ind.D2), ...
						'D3', DesignInfo.Hadamard(:,saveBest.ind.D3), 'X2', X2);
	output.ind 		= struct('D1', unique(output.Core), 'D2', saveBest.ind.D2, 'D3', saveBest.ind.D3);
	output.init 	= Init;
	output.D3Res 	= saveBest;
	% Stop recording time
	output.time 	= toc;
	fprintf('============== End ============== \n\n');
	% End of the Approach
	

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