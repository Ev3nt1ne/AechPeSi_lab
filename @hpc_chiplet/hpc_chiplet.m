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

    methods(Static)
        function g=gaussian_filter(Filter_size, sigma)
			%size=5; %filter size, odd number
			size=Filter_size;
			g=zeros(size,size); %2D filter matrix
			%sigma=2; %standard deviation
			%gaussian filter
			for i=-(size-1)/2:(size-1)/2
    			for j=-(size-1)/2:(size-1)/2
        			x0=(size+1)/2; %center
        			y0=(size+1)/2; %center
        			x=i+x0; %row
        			y=j+y0; %col
        			g(y,x)=exp(-((x-x0)^2+(y-y0)^2)/2/sigma/sigma);
    			end
			end
			%normalize gaussian filter
			sum1=sum(g);
			sum2=sum(sum1);
			g=g/sum2;

			%plot 3D
			%g1=Gaussian_filter(50,2);
			%g2=Gaussian_filter(50,7);
			%g3=Gaussian_filter(50,11);
			%figure(1);
			%subplot(1,3,1);surf(g1);title('filter size = 50, sigma = 2');
			%subplot(1,3,2);surf(g2);title('filter size = 50, sigma = 7');
			%subplot(1,3,3);surf(g3);title('filter size = 50, sigma = 11');
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

