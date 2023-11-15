classdef black_wolf < mpc
	%UNTITLED Summary of this class goes here
	%   Detailed explanation goes here
	
	properties

	end
	
	methods
		function obj = black_wolf()
			%UNTITLED Construct an instance of this class
			%   Detailed explanation goes here

		end		
	end

	methods
		function obj = setup_mpc(obj, hc)
			
			% It's good practice to start by clearing YALMIPs internal database 
			% Every time you call sdpvar etc, an internal database grows larger
			yalmip('clear')

			x = sdpvar(repmat(obj.Ns,1,obj.Nhzn+2),repmat(1,1,obj.Nhzn+2));
			u = sdpvar(repmat(obj.Nc,1,obj.Nhzn),repmat(1,1,obj.Nhzn));
			
			ot = sdpvar(obj.Nc,1); %TODO: decide if single or with horizon
			w = sdpvar(obj.Nc,1); %TODO: decide if single or with horizon
			h0v = sdpvar(obj.Nc,1);

			uo = sdpvar(obj.Nc,1);

			%sdpvar ly_xref;
			ly_uref = sdpvar(obj.Ni_c,1);
			%sdpvar ly_yref;
			ly_usum = sdpvar(obj.Nhzn,obj.vd+1);
			%ly_Cty = dpvar(repmat(obj.Nout,1,obj.Nhzn),repmat(1,1,obj.Nhzn));

			Adl_mpc = obj.Ad_mpc;
			Bdl_mpc = obj.Bd_mpc;
			
			Bu = Bdl_mpc(:,1:obj.Nc);
			Bd = Bdl_mpc(:,obj.Nc+1:end);
			
			[h0, h1, h2] = obj.pws_ls_approx();			
			
			C2 = eye(obj.Ns);
			for k=2:hc.full_model_layers
				C2([k:hc.full_model_layers:end]) = 0;
			end
			C2([end-hc.add_states+1:end], :) = 0;
			
			constraints = [];
			constraints = [x{2} == (Adl_mpc + C2*h1)*x{1}-(C2*h1*273.15*ones(hc.Ns,1))+(Bu+(hc.Cc'*h2))*(u0*w)+Bd*ot + hc.Cc*(h0 + h0v)];

			for k = 1 : obj.Nhzn
					constraints = [constraints, x{k+2} == (Adl_mpc + C2*h1)*x{k+1}-(C2*h1*273.15*ones(hc.Ns,1))+(Bu+(hc.Cc'*h2))*(u{k}.*w)+Bd*ot + hc.Cc*(h0 + h0v)];

					%TODO: Ct I did only for y, and only for max
					if (~isinf(obj.umin))
						%constraints = [constraints, u{k}.*w*(h2+1) + h1*(obj.C(1:obj.Nc,:)*(x{k+1}-273.15)) + h0 >= obj.umin];
					end
					if (~isinf(obj.uMax))
						%constraints = [constraints, u{k}.*w*(h2+1) + h1*(obj.C(1:obj.Nc,:)*(x{k+1}-273.15)) + h0 <= obj.uMax - obj.Ctu(k,:)'];
					end
					if (~isinf(obj.xmin))
						%constraints = [constraints, x{k+2} >= obj.xmin];
					end
					if (~isinf(obj.xMax))
						%constraints = [constraints, x{k+2} <= obj.xMax];
					end
					if (~isinf(obj.ymin))
						%constraints = [constraints, obj.C)*x{k+2} >= obj.ymin];
					end
					if (~isinf(obj.yMax))
						constraints = [constraints, hc.C*x{k+2} <= obj.yMax - obj.Cty(k,:)'];
					end
					%if (~isinf(obj.usum(1)))
						constraints = [constraints,sum(u{k}.*w*(h2+1) + h1*(hc.Cc*(x{k+1}-273.15)) + h0v) <= ly_usum(k,1) - sum(obj.Ctu(k,:),2)];
						%TODO should I also add the "resulting" thing. i.e.
						%	the above formula with x{k+2}?
						%	maybe no, because I use horizon for this
					%end
					%if (~isinf(obj.usum(2:end)))
						%constraints = [constraints, obj.VDom'*(u{k}.*w*(h2+1) + h1*(obj.C(1:obj.Nc,:)*(x{k+1}-273.15)) + h0) <= ly_usum(k,2:end)' - (obj.Ctu(k,:)*obj.VDom)'];
					%end
				end
			
			objective = 0;

			for k = 1:obj.Nhzn
				objective = objective + (u{k}-ly_uref)'*obj.R*(u{k}-ly_uref) + u{k}'*obj.R2*(u{k}) + x{k+2}'*obj.Q*x{k+2};
			end
			
			%ops = sdpsettings('verbose',1,'solver','quadprog', 'usex0',1);
			ops = sdpsettings('verbose',1,'solver','osqp', 'usex0',0); %You have specified an initial point, but the selected solver (OSQP) does not support warm-starts through YALMIP
			ops.quadprog.TolCon = 1e-2;
			ops.quadprog.TolX = 1e-5;
			ops.quadprog.TolFun = 1e-3;
			ops.convertconvexquad = 0;
			%ops.quadprog.MaxPCGIter = max(1, ops.quadprog.MaxPCGIter * 3);
			ops.quadprog.MaxIter = 50;
			obj.Controller = optimizer(constraints,objective,ops,{x{1},u0,w,ot,h0v,ly_uref,ly_usum},{u{1}, x{3}})
			
		end %lin_mpc_setup
		function uout = call_mpc(obj, x, u0, w, ot, h0v, uref, usum)
			
			%{
			if nargin < 5 || ...
				error("Needs current x, u0, w, ot");
			end
			if nargin < 7 || ...
				usum = obj.usum;
			end
			if nargin < 8 || ...
				uref = obj.uref;
			end
			
			if size(usum,1) < obj.Nhzn
				usum(size(usum,1)+1:obj.Nhzn,:) = repmat(usum(size(usum,1),:), obj.Nhzn - size(usum,1),1);
			end
			%}
			
			[uout, problem,~,~,optimizer_object] = obj.Controller({x, u0, w, ot, h0v, uref, usum});

			% Analyze error flags
			if problem
				warning(yalmiperror(problem))
			end			
		end %lin_mpc
	end
end

