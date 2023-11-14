classdef controller
	%CONTROLLERS Summary of this class goes here
	%   Detailed explanation goes here
	
	properties
		Ts_ctrl (1,1) {mustBePositive, mustBeNumeric, mustBeFinite} ...
			= 5e-4;			% Controller Ts

		Nic;
		
	end

	properties(Dependent)

		Ad_ctrl;
		Bd_ctrl;

	end
	
	methods
		function obj = controller(inputArg1,inputArg2)
			%CONTROLLERS Construct an instance of this class
			%   Detailed explanation goes here
			%obj.Property1 = inputArg1 + inputArg2;
		end
		
		function outputArg = method1(obj,inputArg)
			%METHOD1 Summary of this method goes here
			%   Detailed explanation goes here
			%outputArg = obj.Property1 + inputArg;
		end
	end
end

