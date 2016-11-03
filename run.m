clear all
format long g
% 
DesignInfo = getDesignInfo();
DesignInfo.nRuns = 36;
DesignInfo.nFactors = 18;
CoreSet = TwoFIstruct(3);
DesignInfo.CoreFig = CoreSet(:,:,2); 

Row_1 = [1,1,-1,-1,1,1,1,1,-1,1,-1,1,-1,-1,-1,-1,1,1,-1];
hadmat_20 = zeros(20, 20-1);
for i = 1:20,
	if i == 1,
	hadmat_20(i,:) = Row_1;
	elseif i == 20,
	hadmat_20(i,:) = (-1).*ones(1,20-1);
	else,
	hadmat_20(i,:) = Row_1([i:(20-1), 1:(i-1)]);
	end
end

%MAT = horzcat(ones(20,1), hadmat_20);
%MAT = ConvertHadamard(DesignInfo.nRuns, 'http://neilsloane.com/hadamard/had.24.pal.txt');
%MAT = ConvertHadamard(DesignInfo.nRuns, 'http://neilsloane.com/hadamard/had.28.will.txt');
%MAT = ConvertHadamard(DesignInfo.nRuns, 'http://neilsloane.com/hadamard/had.28.pal2.txt');
%MAT = ConvertHadamard(DesignInfo.nRuns, 'http://neilsloane.com/hadamard/had.32.pal.txt');
%MAT = ConvertHadamard(DesignInfo.nRuns, 'http://neilsloane.com/hadamard/had.36.will.txt');
MAT = ConvertHadamard(DesignInfo.nRuns, 'http://neilsloane.com/hadamard/had.36.pal2.txt');
%MAT = ConvertHadamard(DesignInfo.nRuns, 'http://neilsloane.com/hadamard/had.52.will.txt');
DesignInfo.Hadamard = MAT(:,2:end);
%DesignInfo.Hadamard = csvread('H20-Q.csv');

AlgInfo = getAlgInfo(); % Initialize Algorithm Settings
AlgInfo.Iter = 50; % Total OASA Iterations

% Settings for The Updating [D2,D3] Step
AlgInfo.D3.Swarm = 8; % Swarm Size
AlgInfo.D3.Iter = 20;	% Number of Iterations
AlgInfo.D3.qGB = 1; 	% Number of Exchanged Columns with Global Best
AlgInfo.D3.qLB = 2; 	% Number of Exchanged Columns with Local Best
AlgInfo.D3.qRC = 2; 	% Number of Exchanged Columns with Randomly Chosen Columns in Candidate Set

% Settings for The Updating D1 Step (under fixed D2)
AlgInfo.D1f2.Swarm = 32;	% Swarm Size
AlgInfo.D1f2.Iter = 20;		% Number of Iterations
AlgInfo.D1f2.qGB = 1; 		% Number of Exchanged Columns with Global Best
AlgInfo.D1f2.qLB = 1; 		% Number of Exchanged Columns with Local Best
AlgInfo.D1f2.qRC = 2; 		% Number of Exchanged Columns with Randomly Chosen Columns in Candidate Set

% Settings for The Updating D1 Step (under fixed D3)
AlgInfo.D1f3.Swarm = 32; 	% Swarm Size
AlgInfo.D1f3.Iter = 20;		% Number of Iterations
AlgInfo.D1f3.qGB = 1; 		% Number of Exchanged Columns with Global Best
AlgInfo.D1f3.qLB = 1; 		% Number of Exchanged Columns with Local Best
AlgInfo.D1f3.qRC = 2; 		% Number of Exchanged Columns with Randomly Chosen Columns in Candidate Set


% Run OASA
output = OASA(DesignInfo, AlgInfo)