classdef mpc_hpc < controller
	%MPC Summary of this class goes here
	%   Detailed explanation goes here
	
	properties
		%% MPC Controller
		Q;
		R;
		Rs;
        Rt;
		Nhzn = 3;
		Ctx;
		Ctu;
		Cty;
		
		umin;
		umax;
		xmin;
		xMax;
		ymin;
		yMax;


		Ts_obs (1,1) {mustBePositive, mustBeNumeric, mustBeFinite} ...
			= 1e-3;			% Observer Ts

		Obs_poles = [0.8 0.1];		% Poles of the Luemberg Observer
        Qcov;
        Rcov;
        Gw;

		save_solver_stats = 0;

        LK;

	end

	properties(Dependent)
		Ad_obs;
		Bd_obs;
	end

	properties(SetAccess=protected, GetAccess=public)		
		mpc_ctrl;					% Persistent variable to optimize controller (prev: Controller)
		failed;
		xlplot;
		tmpc;
		output_mpc;
        Tobs;

		solver_stats;
	end
	
	methods
		function obj = mpc_hpc()
			%MPC Construct an instance of this class
			%   Detailed explanation goes here

			% check Yalmip
			try 
				yalmip('clear')
			catch
				addpath(genpath([pwd filesep 'YALMIP']));
				savepath
			end
			% check OSQP
			try
				osqp;
			catch
				if obj.osunix
					%{
					LIBCURL_PATH = "/lib/x86_64-linux-gnu/";
					path1 = getenv('LD_LIBRARY_PATH');			  % Store existing path
					path = strcat(LIBCURL_PATH, ':', path1);      % Add compatible libcurl to path
					setenv('LD_LIBRARY_PATH', path);              % Set the path
					cd osqp-matlab
					!git submodule update --init --recursive
					make_osqp
					addpath(pwd)
					savepath
					cd ../
					%}
					websave('install_osqp.m','https://raw.githubusercontent.com/osqp/osqp-matlab/master/package/install_osqp.m');
					install_osqp;
				else
					%TODO Windows
				end
			end
			
		end
	end

	methods(Abstract=true)
		[obj] = setup_mpc(obj)
		[uout] = call_mpc(obj)
	end
end

