%% Demonstration script for new NVU model

%% Version control
% This is the latest version of the NVU model: Version 2.0
% K+, NO, Astrocytic Ca2+, TRPV4, ECS, Neuron

%% Construct NVU
% The NVU consists of a number of submodules, implemented as MATLAB
% classes, presently an astrocyte, a lumped SMC/EC model, and a model of
% the wall mechanics. 
%
% The parameters of each of these submodules are
% specified when the modules are constructed, here in the call to NVU, see
% for example the |SMCEC| part of the |NVU| call below:
%
% Options for the ODE solver (currently |ode15s|) are provided by
% specifying the |odeopts| parameter. The code works fine with default
% tolerances.
clear all



odeopts = odeset('RelTol', 1e-04, 'AbsTol', 1e-04, 'MaxStep', 0.5, 'Vectorized', 1);

XLIM1 = 0; XLIM2 = 250;
FIG_NUM = 1;


NEURONAL_START      = 20;      % Start of neuronal stimulation
NEURONAL_END        = 40;      % End of neuronal stimulation 
CURRENT_STRENGTH    = 0.014;  % Strength of current input in mA/cm2  (0.009 for subthreshold, 0.011 for bursting, 0.015 for constant stimulation)

ECS_START       = 30000;      % Start of ECS K+ input
ECS_END         = 32000;      % End of ECS K+ input

J_PLC           = 0.11;      % Jplc value in EC: 0.11 for steady state, 0.3 for oscillations

GLU_SWITCH      = 1;        % Turn on glutamate input (for NO and Ca2+ pathways)
NO_PROD_SWITCH  = 1;        % Turn on Nitric Oxide production 
TRPV_SWITCH     = 1;        % Turn on TRPV4 Ca2+ channel from AC to PVS
RK_SWITCH       = 0;        % Make R_k variable (1) or constant (0)
O2SWITCH        = 1;        % 0: ATP is plentiful, 1: ATP is limited (oxygen-limited regime, default)

% Output start time and estimated finish time based on how many sec of
% stimulation there will be and speed of computer 
%(home: 0.38, work: 0.88)
estimatedMinutes = (NEURONAL_END - NEURONAL_START)*0.88;
timeStart = datetime('now');
estimatedTimeEnd = timeStart + minutes(estimatedMinutes); 
fprintf('Start time is %s\nEstimated end time is %s\n', char(timeStart), char(estimatedTimeEnd));

nv = NVU(Neuron('SC_coup', 10, 'O2switch', O2SWITCH, 'startpulse', NEURONAL_START, 'lengthpulse', NEURONAL_END - NEURONAL_START, 'Istrength', CURRENT_STRENGTH, 'GluSwitch', GLU_SWITCH, 'NOswitch', NO_PROD_SWITCH, 't0_ECS', ECS_START, 'ECS_input', 9), ...
    Astrocyte('R_decay', 0.15, 'Rk_switch', RK_SWITCH, 'trpv_switch', TRPV_SWITCH, 'startpulse', NEURONAL_START, 'lengthpulse', NEURONAL_END - NEURONAL_START, 't0_ECS', ECS_START, 'tend_ECS', ECS_END), ...
    WallMechanics('wallMech', 3), ...
    SMCEC('J_PLC', J_PLC, 'NOswitch', NO_PROD_SWITCH), 'odeopts', odeopts);

dt = 0.001;
nv.T = 0:dt:XLIM2;
numTimeSteps = length(nv.T);
nv.simulate()

% Find a point in time 50 sec before neuronal stimulation has begun (preNeuronalStimTime1)
% and the point in time when stimulation begins (preNeuronalStimTime2) and get the values for 
% CBF, CBV, DHG at that time, then find the midpoint between the min and max and normalise with that value. 
% Done in this way so that it also works when the variables are oscillatory (i.e. J_PLC = 0.3)
preNeuronalStimTime1 = floor((NEURONAL_START-5)*numTimeSteps/XLIM2);
preNeuronalStimTime2 = floor((NEURONAL_START)*numTimeSteps/XLIM2);
CBF = nv.out('CBF'); CBV = nv.out('CBV'); DHG = nv.out('DHG'); CMRO2 = nv.out('CMRO2');
CBF_0 = 0.5*( max(CBF(preNeuronalStimTime1:preNeuronalStimTime2)) + min(CBF(preNeuronalStimTime1:preNeuronalStimTime2)) );
CBV_0 = 0.5*( max(CBV(preNeuronalStimTime1:preNeuronalStimTime2)) + min(CBV(preNeuronalStimTime1:preNeuronalStimTime2)) );
DHG_0 = 0.5*( max(DHG(preNeuronalStimTime1:preNeuronalStimTime2)) + min(DHG(preNeuronalStimTime1:preNeuronalStimTime2)) );
CMRO2_0 = 0.5*( max(CMRO2(preNeuronalStimTime1:preNeuronalStimTime2)) + min(CMRO2(preNeuronalStimTime1:preNeuronalStimTime2)) );

np = nv.neuron.params; % Shortcut for neuron parameters

CBF_N = CBF./CBF_0;
CBV_N = (CBV./CBV_0)';
DHG_N = (DHG./DHG_0)';
CMRO2_N = CMRO2./CMRO2_0;

HBT_N = CBF_N.*DHG_N./CMRO2_N;
HBO_N = (HBT_N - 1) - (DHG_N - 1) + 1;

figure;
plot(nv.T, HBT_N, nv.T, HBO_N, nv.T, DHG_N);
legend('HBT','HBO','DHG');

% % Run this to take the ICs as the end of the last simulation run i.e. steady state ICs
% ICs = (nv.U(end, :))';
% nv.u0 = ICs;
% nv.simulateManualICs() 

% Plot figures - whatever you want

%% Plot for paper
% figure(10);
% subplot(4,2,1);
%     hold all;
%     plot(nv.T, nv.out('v_sa'), 'k', 'LineWidth', 1);
%     ylabel('v_{sa} [mV]');
%     xlim([XLIM1 XLIM2])
%     ylim([-100 20])
%     title('A. Soma membrane potential');
%     xlabel('time (s)');
%     rectangle('Position',[300, -100, 16, 5],'FaceColor',[0 0 0])
%     rectangle('Position',[320, -100, 2, 5],'FaceColor',[0 0 0])
% subplot(4,2,2);
%     hold all;
%     plot(nv.T, nv.out('K_e'), 'k', 'LineWidth', 1);
%     ylabel('K_e [mM]');
%     xlim([XLIM1 XLIM2])
%     title('B. Extracellular potassium');
%     xlabel('time (s)');
%     rectangle('Position',[300, 2, 16, 0.2],'FaceColor',[0 0 0])
%     rectangle('Position',[320, 2, 2, 0.2],'FaceColor',[0 0 0])
% subplot(4,2,3);
%     hold all;
%     plot(nv.T, nv.out('R')*1e6, 'k', 'LineWidth', 1);
%     ylabel('Radius [\mum]');
%     ylim([22.4 25])
%     xlim([XLIM1 XLIM2])
%     title('C. Radius');
%     xlabel('time (s)');
%     rectangle('Position',[300, 22.4, 16, 0.1],'FaceColor',[0 0 0])
%     rectangle('Position',[320, 22.4, 2, 0.1],'FaceColor',[0 0 0])
% subplot(4,2,4);
%     hold all;
%     plot(nv.T, (nv.out('CBF')-CBF_0)./CBF_0, 'k', 'LineWidth', 1);
%     ylabel('\Delta CBF / CBF_0 [-]');
%     xlim([XLIM1 XLIM2])
%     ylim([-0.1 0.4])
%     title('D. Normalised \DeltaCBF');
%     xlabel('time (s)');
%     rectangle('Position',[300, -0.1, 16, 0.02],'FaceColor',[0 0 0])
%     rectangle('Position',[320, -0.1, 2, 0.02],'FaceColor',[0 0 0])
% subplot(4,2,5);
%     hold all;
%     plot(nv.T, nv.out('O2'), 'k', 'LineWidth', 1);
%     ylabel('O_2 [mM]');
%     xlim([XLIM1 XLIM2])
%     ylim([0.0265 0.03])
%     title('E. Oxygen concentration');
%     xlabel('time (s)');
%     rectangle('Position',[300, 0.0265, 16, 0.00015],'FaceColor',[0 0 0])
%     rectangle('Position',[320, 0.0265, 2, 0.00015],'FaceColor',[0 0 0])
% subplot(4,2,6);
%     hold all;
%     plot(nv.T, nv.out('CBV')./CBV_0, 'k', 'LineWidth', 1);
%     ylabel('CBV [-]');
%     xlim([XLIM1 XLIM2])
%     ylim([0.97 1.1])
%     title('F. Normalised CBV');
%     xlabel('time (s)');
%     rectangle('Position',[300, 0.97, 16, 0.005],'FaceColor',[0 0 0])
%     rectangle('Position',[320, 0.97, 2, 0.005],'FaceColor',[0 0 0])
% subplot(4,2,7);
%     hold all;
%     plot(nv.T, nv.out('DHG')./DHG_0, 'k', 'LineWidth', 1);
%     ylabel('DHG [-]');
%     xlim([XLIM1 XLIM2])
%     ylim([0.86 1.05])
%     title('G. Normalised deoxyhemoglobin');
%     xlabel('time (s)');
%     rectangle('Position',[300, 0.86, 16, 0.008],'FaceColor',[0 0 0])
%     rectangle('Position',[320, 0.86, 2, 0.008],'FaceColor',[0 0 0])
% subplot(4,2,8);
%     hold all;
%     plot(nv.T, 100 * np.V_0 * ( np.a_1 * (1 - nv.out('DHG')/DHG_0) - np.a_2 * (1 - nv.out('CBV')/CBV_0) ), 'k', 'LineWidth', 1);
%     ylabel('\Delta BOLD (%)');
%     xlim([XLIM1 XLIM2])
%     ylim([-0.9 1.5])
%     title('H. BOLD response');
%     xlabel('time (s)');
%     rectangle('Position',[300, -0.9, 16, 0.09],'FaceColor',[0 0 0])
%     rectangle('Position',[320, -0.9, 2, 0.09],'FaceColor',[0 0 0])

%% Plot BOLD variables

figure(FIG_NUM);
subplot(3,3,1);
    hold all;
    plot(nv.T, nv.out('v_sa'), 'LineWidth', 1);
    ylabel('v_{sa} [mV]');
    xlim([XLIM1 XLIM2])
subplot(3,3,2);
    hold all;
    plot(nv.T, nv.out('Ca_i'), 'LineWidth', 1);
    ylabel('SMC Ca2+');
    xlim([XLIM1 XLIM2])
subplot(3,3,3);
    hold all;
    plot(nv.T, nv.out('Na_e'), 'LineWidth', 1);
    ylabel('Na_e [mM]');
    xlim([XLIM1 XLIM2])
subplot(3,3,4);
    hold all;
    plot(nv.T, (nv.out('CBF')-CBF_0)./CBF_0, 'LineWidth', 1);
    ylabel('\Delta CBF / CBF_0');
    xlim([XLIM1 XLIM2])
subplot(3,3,5);
    hold all;
    plot(nv.T, nv.out('R')*1e6, 'LineWidth', 1);
    ylabel('Radius [\mum]');
    xlim([XLIM1 XLIM2])
subplot(3,3,6);
    hold all;
    plot(nv.T, nv.out('K_s')/1e3, 'LineWidth', 1);
    ylabel('K_s');
    xlim([XLIM1 XLIM2])
subplot(3,3,7);
    hold all;
    plot(nv.T, nv.out('CBV')./CBV_0, 'LineWidth', 1);
    ylabel('CBV');
    xlim([XLIM1 XLIM2])
subplot(3,3,8);
    hold all;
    plot(nv.T, nv.out('DHG')./DHG_0, 'LineWidth', 1);
    ylabel('DHG');
    xlim([XLIM1 XLIM2])
subplot(3,3,9);
    hold all;
    % BOLD using normalised variables
    plot(nv.T, 100 * np.V_0 * ( np.a_1 * (1 - nv.out('DHG')/DHG_0) - np.a_2 * (1 - nv.out('CBV')/CBV_0) ), 'LineWidth', 1);  
    ylabel('\Delta BOLD (%)');
    xlim([XLIM1 XLIM2])
    
%     % 
% figure(FIG_NUM + 1);
% 
% subplot(2,3,1)
%     hold all;
%     plot(nv.T, nv.out('Ca_k'), 'LineWidth', 1);
%     ylabel('Ca_k [\muM]');
%     xlim([XLIM1 XLIM2])
%     xlabel('Time [s]'); 
% subplot(2,3,2)
%     hold all;
%     plot(nv.T, nv.out('eet_k'), 'LineWidth', 1);
%     ylabel('eet_k [\muM]');
%     xlim([XLIM1 XLIM2])
%     xlabel('Time [s]'); 
% subplot(2,3,3)
%     hold all;
%     plot(nv.T, nv.out('w_k'), 'LineWidth', 1);
%     ylabel('w_k [-]');
%     xlim([XLIM1 XLIM2])
%     xlabel('Time [s]'); 
% subplot(2,3,4)
%     hold all;
%     plot(nv.T, nv.out('R')*1e6, 'LineWidth', 1);
%     ylabel('Radius [\mum]');
%     xlim([XLIM1 XLIM2])
%     xlabel('Time [s]'); 
% subplot(2,3,5)
%     hold all;
%     plot(nv.T, nv.out('K_p')/1e3, 'LineWidth', 1);
%     xlabel('Time [s]'); 
%     ylabel('K_p [mM]');
%     xlim([XLIM1 XLIM2])
% subplot(2,3,6)
%     hold all;
%     plot(nv.T, nv.out('R')*1e6, 'LineWidth', 1)
%     xlabel('Time [s]'); 
%     ylabel('Radius [\mum]');
%     xlim([XLIM1 XLIM2])

timeEnd = datetime('now');
fprintf('End time is %s\n', char(timeEnd));
