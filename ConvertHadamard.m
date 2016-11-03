
function Hadamard = ConvertHadamard(nRuns, SymbolHadamardURL)
%% Convert Hadamard matrix from symbol to logical number on http://neilsloane.com/hadamard/
SymbolHadamard = urlread(SymbolHadamardURL);
Hadamard = zeros(nRuns,nRuns);
for i = 1:nRuns,
	temp_str = SymbolHadamard((1+(nRuns+1)*(i-1)):((nRuns+1)*i));
	temp_vec = zeros(1,nRuns);
	for j = 1:nRuns,
		switch temp_str(j)
			case '+', temp_vec(j) = 1;
			case '-', temp_vec(j) = -1;
		end
	end
	% Force the first column to 1's (for intercept)
	if temp_vec(1) < 0,
		temp_vec = (-1).*temp_vec;
	end
	Hadamard(i,:) = temp_vec;
end

