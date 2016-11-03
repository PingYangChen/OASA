function Core = TwoFIstruct(e)

	if e > 4,
		error('We consider the situation that the number of 2fi is no more than 4');
	end

	switch e
		
		case 1,
			
			Core = [1 2];
		
		case 2,
		
			Core = cat(3, [1 2; 3 4], [1 2; 2 3]);
			
		case 3,
		
			Core = cat(3, [1 2; 3 4; 5 6], [1 2; 2 3; 4 5], [1 2; 1 3; 1 4], [1 2; 2 3; 3 4], [1 2; 2 3; 1 3]);
			
		case 4,
		
			Core = cat(3, [1 2; 3 4; 5 6; 7 8], [1 2; 2 3; 4 5; 6 7], [1 2; 1 3; 1 4; 5 6], [1 2; 2 3; 3 4; 5 6], ...
				[1 2; 2 3; 1 3; 4 5], [1 2; 2 3; 4 5; 5 6], [1 2; 1 3; 1 4; 1 5], [1 2; 2 3; 2 4; 4 5], ...
				[1 2; 2 3; 3 4; 4 5], [1 2; 2 3; 1 3; 3 4], [1 2; 2 3; 3 4; 1 4]);
		
	end
