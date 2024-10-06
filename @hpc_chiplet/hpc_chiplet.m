classdef hpc_chiplet < thermal_model & power_model & perf_model & handle
    %HPC_CHIPLET Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        vd (1,1) {mustBePositive, mustBeInteger} ...
			= 3;		% Voltage/Power Domains
		VDom {mustBeNonnegative, mustBeNumericOrLogical, mustBeNonempty, mustBeLessThanOrEqual(VDom,1)} ...
			= [1 0 0; 1 0 0; 0 1 0; 0 1 0; 1 0 0; 1 0 0; 0 1 0; 0 1 0; 0 0 1; 0 0 1; 0 0 1; 0 0 1];		% Structure of the Voltage Domains Nc x vd
		

		%measure_noise = 1;
		sensor_noise (1,1) {mustBeNonnegative, mustBeNumericOrLogical} ...
			= 1;
		%T_noise_max = 1.0;
		sensor_noise_amplitude (3,1) {mustBeNonnegative, mustBeNumeric, mustBeNonempty, mustBeFinite} ...
			= [1.0; 1.0; 1.0];	% Process, Voltage, Temperature amplitudes.
    end
    
    methods
        function obj = hpc_chiplet()
            %HPC_CHIPLET Construct an instance of this class
            %   Detailed explanation goes here
            
        end
        
        function [lvd, lVDom] = domain_give_model(obj)  
            lvd = obj.vd;
            lVDom = obj.VDom;
        end
    end

	%% Dependent Variables
	methods
    	function obj = set.VDom(obj, val)
			cmp = [obj.Nc obj.vd];
			if all(size(val) == cmp) && ~isempty(val)
				obj.VDom = val;
			else
				warning("[LAB Error]: VDom wrong input size.");
				disp(strcat("[LAB] VDom required size:", num2str(cmp(2))));
			end
    	end
    end
end

