
%{ 
% Missing:

quantum_us = 50;
customColormap;
observable;

properties(SetAccess=protected, GetAccess=public)		
	Controller;					% Persistent variable to optimize controller
	polyFV_opt;
	
	% wl:
	% How much? coefficient
	dur_c = 5/6;
	min_secondary_space = 0.075;
	
	% model
	 wl_index;
	 qt_storage;
	 V_s;
	 F_s;
	 A_s;
	 B_s;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% In the constructor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%obj.x_init = ones(obj.Ns,1)*obj.temp_amb;
% TODO: remove this after leakage is fixed:
leak_store = obj.exp_leakage;
obj.exp_leakage = 0;
%obj.min_pw_red = obj.power_compute(ones(obj.Nc,1)*obj.F_min,ones(obj.Nc,1)*obj.V_min,(obj.core_crit_temp)*ones(obj.Nc,1),[zeros(obj.Nc,obj.ipl-1) ones(obj.Nc,1)],ones(obj.Nc,1));
obj.min_pw_red = obj.power_compute(ones(obj.Nc,1)*obj.F_min,ones(obj.Nc,1)*obj.V_Max,(obj.core_crit_temp)*ones(obj.Nc,1),[zeros(obj.Nc,obj.ipl-1) ones(obj.Nc,1)],ones(obj.Nc,1));
obj.exp_leakage = leak_store;

%}

%%

a = thermal_model;

%a.VDom = [1, NaN, inf, 'l'];

%a.Ac_nom = [1 1;1 1];

a.C = ones(a.Nout, a.Ns);

%% Things to put in the constructor
Ac_nom
Bc_nom
Ac_true
Bc_true
C
D
param_dev_per

% things to ALSO fix in set
RC_floorplan_length
CPw_floorplan_length
R_floorplan_material
C_floorplan_material
sensor_active

%
add_states
t_comp
wid_comp
len_comp


%%
hpc = hpc_lab;

hpc.t_init = hpc.temp_amb*ones(hpc.Ns,1);
hpc.sim_tm_autonomous()

%%
tnA = hpc.Ac_nom;

hpc2 = hpc_system();

% Cores:
hpc2.Nc = 9; %72;
hpc2.Nh = 3; %9;			% Number of rows
hpc2.Nv = 3; %8;			% Number of cols

% Version of the Thermal Model
hpc2.thermal_model_ver = 0;
% Exponential Leakage
hpc2.exp_leakage = 1;
% Noise and Variation
hpc2.model_variation = 1;
hpc2 = hpc2.create_core_pw_noise();
hpc2.measure_noise = 1;

% External Ambient Temperature
temp_amb = 25.0 + 273.15;

% Create Thermal Model:
hpc2 = hpc2.create_thermal_model();

toA = hpc2.Ac_nom;

dif = tnA - toA


difB = hpc.Bc_nom - hpc2.Bc_nom

%%

aa = [
	

    0.0706    0.2000    0.5882
    0.0803    0.2175    0.5760
    0.0900    0.2351    0.5637
    0.0997    0.2526    0.5515
    0.1094    0.2702    0.5392
    0.1191    0.2877    0.5270
    0.1288    0.3053    0.5147
    0.1385    0.3228    0.5025
    0.1482    0.3404    0.4902
    0.1579    0.3579    0.4779
    0.1676    0.3755    0.4657
    0.1774    0.3930    0.4534
    0.1871    0.4106    0.4412
    0.1968    0.4281    0.4289
    0.2065    0.4457    0.4167
    0.2162    0.4632    0.4044
    0.2259    0.4808    0.3922
    0.2356    0.4983    0.3799
    0.2453    0.5159    0.3676
    0.2550    0.5334    0.3554
    0.2647    0.5510    0.3431
    0.2744    0.5685    0.3309
    0.2841    0.5861    0.3186
    0.2938    0.6036    0.3064
    0.3035    0.6212    0.2941
    0.3132    0.6387    0.2819
    0.3229    0.6563    0.2696
    0.3326    0.6738    0.2574
    0.3424    0.6914    0.2451
    0.3521    0.7089    0.2328
    0.3618    0.7265    0.2206
    0.3715    0.7440    0.2083
    0.3812    0.7616    0.1961
    0.3909    0.7791    0.1838
    0.4006    0.7967    0.1716
    0.4103    0.8142    0.1593
    0.4200    0.8318    0.1471
    0.4297    0.8493    0.1348
    0.4394    0.8669    0.1225
    0.4491    0.8844    0.1103
    0.4588    0.9020    0.0980
    0.4714    0.9042    0.0981
    0.4840    0.9065    0.0981
    0.4966    0.9088    0.0982
    0.5092    0.9111    0.0982
    0.5218    0.9134    0.0983
    0.5343    0.9156    0.0983
    0.5469    0.9179    0.0984
    0.5595    0.9202    0.0984
    0.5721    0.9225    0.0984
    0.5847    0.9248    0.0985
    0.5973    0.9270    0.0985
    0.6098    0.9293    0.0986
    0.6224    0.9316    0.0986
    0.6350    0.9339    0.0987
    0.6476    0.9362    0.0987
    0.6602    0.9384    0.0988
    0.6728    0.9407    0.0988
    0.6854    0.9430    0.0989
    0.6979    0.9453    0.0989
    0.7105    0.9476    0.0990
    0.7231    0.9498    0.0990
    0.7357    0.9521    0.0990
    0.7483    0.9544    0.0991
    0.7609    0.9567    0.0991
    0.7735    0.9590    0.0992
    0.7860    0.9612    0.0992
    0.7986    0.9635    0.0993
    0.8112    0.9658    0.0993
    0.8238    0.9681    0.0994
    0.8364    0.9704    0.0994
    0.8490    0.9726    0.0995
    0.8616    0.9749    0.0995
    0.8741    0.9772    0.0995
    0.8867    0.9795    0.0996
    0.8993    0.9818    0.0996
    0.9119    0.9840    0.0997
    0.9245    0.9863    0.0997
    0.9371    0.9886    0.0998
    0.9497    0.9909    0.0998
    0.9622    0.9932    0.0999
    0.9748    0.9954    0.0999
    0.9874    0.9977    0.1000
    1.0000    1.0000    0.1000
    0.9954    0.9891    0.0989
    0.9908    0.9782    0.0977
    0.9862    0.9672    0.0966
    0.9816    0.9563    0.0954
    0.9770    0.9454    0.0943
    0.9724    0.9345    0.0931
    0.9678    0.9236    0.0920
    0.9632    0.9126    0.0908
    0.9586    0.9017    0.0897
    0.9540    0.8908    0.0885
    0.9494    0.8799    0.0874
    0.9448    0.8690    0.0862
    0.9402    0.8580    0.0851
    0.9356    0.8471    0.0839
    0.9310    0.8362    0.0828
    0.9264    0.8253    0.0816
    0.9218    0.8144    0.0805
    0.9172    0.8034    0.0793
    0.9126    0.7925    0.0782
    0.9080    0.7816    0.0770
    0.9034    0.7707    0.0759
    0.8989    0.7598    0.0747
    0.8943    0.7489    0.0736
    0.8897    0.7379    0.0724
    0.8851    0.7270    0.0713
    0.8805    0.7161    0.0701
    0.8759    0.7052    0.0690
    0.8713    0.6943    0.0678
    0.8667    0.6833    0.0667
    0.8621    0.6724    0.0655
    0.8575    0.6615    0.0644
    0.8529    0.6506    0.0632
    0.8483    0.6397    0.0621
    0.8437    0.6287    0.0609
    0.8391    0.6178    0.0598
    0.8345    0.6069    0.0586
    0.8299    0.5960    0.0575
    0.8253    0.5851    0.0563
    0.8207    0.5741    0.0552
    0.8161    0.5632    0.0540
    0.8115    0.5523    0.0529
    0.8069    0.5414    0.0517
    0.8023    0.5305    0.0506
    0.7977    0.5195    0.0494
    0.7931    0.5086    0.0483
    0.7885    0.4977    0.0471
    0.7839    0.4868    0.0460
    0.7793    0.4759    0.0448
    0.7747    0.4649    0.0437
    0.7701    0.4540    0.0425
    0.7655    0.4431    0.0414
    0.7609    0.4322    0.0402
    0.7563    0.4213    0.0391
    0.7517    0.4103    0.0379
    0.7471    0.3994    0.0368
    0.7425    0.3885    0.0356
    0.7379    0.3776    0.0345
    0.7333    0.3667    0.0333
    0.7287    0.3557    0.0322
    0.7241    0.3448    0.0310
    0.7195    0.3339    0.0299
    0.7149    0.3230    0.0287
    0.7103    0.3121    0.0276
    0.7057    0.3011    0.0264
    0.7011    0.2902    0.0253
    0.6966    0.2793    0.0241
    0.6920    0.2684    0.0230
    0.6874    0.2575    0.0218
    0.6828    0.2466    0.0207
    0.6782    0.2356    0.0195
    0.6736    0.2247    0.0184
    0.6690    0.2138    0.0172
    0.6644    0.2029    0.0161
    0.6598    0.1920    0.0149
    0.6552    0.1810    0.0138
    0.6506    0.1701    0.0126
    0.6460    0.1592    0.0115
    0.6414    0.1483    0.0103
    0.6368    0.1374    0.0092
    0.6322    0.1264    0.0080
    0.6276    0.1155    0.0069
    0.6230    0.1046    0.0057
    0.6184    0.0937    0.0046
    0.6138    0.0828    0.0034
    0.6092    0.0718    0.0023
    0.6046    0.0609    0.0011
    0.6000    0.0500         0
    0.5947    0.0494         0
    0.5894    0.0488         0
    0.5841    0.0482         0
    0.5788    0.0476         0
    0.5735    0.0471         0
    0.5682    0.0465         0
    0.5629    0.0459         0
    0.5576    0.0453         0
    0.5524    0.0447         0
    0.5471    0.0441         0
    0.5418    0.0435         0
    0.5365    0.0429         0
    0.5312    0.0424         0
    0.5259    0.0418         0
    0.5206    0.0412         0
    0.5153    0.0406         0
    0.5100    0.0400         0
    0.5047    0.0394         0
    0.4994    0.0388         0
    0.4941    0.0382         0
    0.4888    0.0376         0
    0.4835    0.0371         0
    0.4782    0.0365         0
    0.4729    0.0359         0
    0.4676    0.0353         0
    0.4624    0.0347         0
    0.4571    0.0341         0
    0.4518    0.0335         0
    0.4465    0.0329         0
    0.4412    0.0324         0
    0.4359    0.0318         0
    0.4306    0.0312         0
    0.4253    0.0306         0
    0.4200    0.0300         0
    0.4147    0.0294         0
    0.4094    0.0288         0
    0.4041    0.0282         0
    0.3988    0.0276         0
    0.3935    0.0271         0
    0.3882    0.0265         0
    0.3829    0.0259         0
    0.3776    0.0253         0
    0.3724    0.0247         0
    0.3671    0.0241         0
    0.3618    0.0235         0
    0.3565    0.0229         0
    0.3512    0.0224         0
    0.3459    0.0218         0
    0.3406    0.0212         0
    0.3353    0.0206         0
    0.3300    0.0200         0
    0.3247    0.0194         0
    0.3194    0.0188         0
    0.3141    0.0182         0
    0.3088    0.0176         0
    0.3035    0.0171         0
    0.2982    0.0165         0
    0.2929    0.0159         0
    0.2876    0.0153         0
    0.2824    0.0147         0
    0.2771    0.0141         0
    0.2718    0.0135         0
    0.2665    0.0129         0
    0.2612    0.0124         0
    0.2559    0.0118         0
    0.2506    0.0112         0
    0.2453    0.0106         0
    0.2400    0.0100         0
    0.2347    0.0094         0
    0.2294    0.0088         0
    0.2241    0.0082         0
    0.2188    0.0076         0
    0.2135    0.0071         0
    0.2082    0.0065         0
    0.2029    0.0059         0
    0.1976    0.0053         0
    0.1924    0.0047         0
    0.1871    0.0041         0
    0.1818    0.0035         0
    0.1765    0.0029         0
    0.1712    0.0024         0
    0.1659    0.0018         0
    0.1606    0.0012         0
    0.1553    0.0006         0
    0.1500         0         0
]

