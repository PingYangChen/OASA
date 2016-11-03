% Initialize Target Design Settings
DesignInfo = getDesignInfo();
% Set Number of Runs
DesignInfo.nRuns = 36;			
% Set Number of Factors
DesignInfo.nFactors = 18;		
% Choose Type of Two Factor Interaction
CoreSet = TwoFIstruct(3);		
DesignInfo.CoreFig = CoreSet(:,:,2); 
% Get Hadamard Matrix Based on Run Size (Please visit http://neilsloane.com/hadamard/)
MAT = ConvertHadamard(DesignInfo.nRuns, 'http://neilsloane.com/hadamard/had.36.will.txt');
%MAT = ConvertHadamard(DesignInfo.nRuns, 'http://neilsloane.com/hadamard/had.36.pal2.txt'); % another Hadamard matrix
DesignInfo.Hadamard = MAT(:,2:end); % the 1st column is for intercept term

% Note that, for 20 runs, we have four Hadamard matrices in our repository, for example,
% DesignInfo.Hadamard = csvread('H20-Q.csv');

% Initialize Algorithm Settings
AlgInfo = getAlgInfo(); 
AlgInfo.Iter = 50; 		% Total OASA Iterations

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

