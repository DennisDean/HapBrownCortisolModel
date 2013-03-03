function testBrownCortisolModel2001()
%testCortModel2001 Generate simulated cortisol time series
%   Function generates simuated time series for the cortisol time series
%   described in:
%
%   Brown EN, Meehan PM, Demster AP (2001) A stochastic differential 
%   equation model of diurnal cortisol patterns. American Journal of 
%   Physiology - Endocrinology and Metabolism 280: E450-E461.
%
% Simulated cortisol time series are being generated as requested for the
% hierarhical analysis paper that is currently under revision at PLoS One.
%
% ---------------------------------------------
% Dennis A. Dean, II, Ph.D
%
% Program for Sleep and Cardiovascular Medicine
% Brigam and Women's Hospital
% Harvard Medical School
% 221 Longwood Ave
% Boston, MA  02149
%
% File created: February 24, 2013
% Last updated: March 2, 2013
%    
% Copyright © [2013] The Brigham and Women's Hospital, Inc. THE BRIGHAM AND 
% WOMEN'S HOSPITAL, INC. AND ITS AGENTS RETAIN ALL RIGHTS TO THIS SOFTWARE 
% AND ARE MAKING THE SOFTWARE AVAILABLE ONLY FOR SCIENTIFIC RESEARCH 
% PURPOSES. THE SOFTWARE SHALL NOT BE USED FOR ANY OTHER PURPOSES, AND IS
% BEING MADE AVAILABLE WITHOUT WARRANTY OF ANY KIND, EXPRESSED OR IMPLIED, 
% INCLUDING BUT NOT LIMITED TO IMPLIED WARRANTIES OF MERCHANTABILITY AND 
% FITNESS FOR A PARTICULAR PURPOSE. THE BRIGHAM AND WOMEN'S HOSPITAL, INC. 
% AND ITS AGENTS SHALL NOT BE LIABLE FOR ANY CLAIMS, LIABILITIES, OR LOSSES 
% RELATING TO OR ARISING FROM ANY USE OF THIS SOFTWARE.
%



% Test Flags
Test1 = 0;     % Plot randomly generated cortisol time seris
Test2 = 0;     % Get simulated data and plot 
Test3 = 0;     % Get time+simulated data and plot
Test4 = 0;     % Get time+simulated data and cortisol parameters; plot
Test5 = 0;     % Pass parameters 
Test6 = 0;     % Pass parameters and set title
Test7 = 1;     % Set amplitudes and interpulse intervals

%------------------------------------------------------------------- Test 1
if Test1 == 1;
    % Create randomly generated cortisol time series
    BrownCortisolModel2001;
end
%------------------------------------------------------------------- Test 2
if Test2 == 1;
    % Create randomly generated cortisol time series
    Y = BrownCortisolModel2001;
    
    %------------------ Display Simulated Cortisol
    fid = figure('InvertHardcopy','off','Color',[1 1 1]);
    plot ([1:1:length(Y)]-1, Y,'LineWidth',2);
    
    % Reformat Axis
    set(gca, 'LineWidth',2);
    set(gca, 'FontWeight','bold');
    set(gca, 'FontSize',14);
    
    % Label
    xlabel('Time (min)','FontWeight','bold','FontSize',14);
    ylabel('Cortisol Amplitude (ug/dL)','FontWeight','bold','FontSize',14);
end
%------------------------------------------------------------------- Test 2
if Test3 == 1;
    % Create randomly generated cortisol time series
    [t Y] = BrownCortisolModel2001;
    
    %------------------ Display Simulated Cortisol
    fid = figure('InvertHardcopy','off','Color',[1 1 1]);
    plot (t, Y,'LineWidth',2);
    
    % Reformat Axis
    set(gca, 'LineWidth',2);
    set(gca, 'FontWeight','bold');
    set(gca, 'FontSize',14);
    
    % Label
    xlabel('Time (min)','FontWeight','bold','FontSize',14);
    ylabel('Cortisol Amplitude (ug/dL)','FontWeight','bold','FontSize',14);
end
%------------------------------------------------------------------- Test 2
if Test4 == 1;
    % Get time+simulated data and circadian parameter; plot
    [t Y cortStruct] = BrownCortisolModel2001;
    
    % Show return strucut
    fprintf('Test 4: cortStruct');
    disp()
    %------------------ Display Simulated Cortisol
    fid = figure('InvertHardcopy','off','Color',[1 1 1]);
    plot (t, Y,'LineWidth',2);
    
    % Reformat Axis
    set(gca, 'LineWidth',2);
    set(gca, 'FontWeight','bold');
    set(gca, 'FontSize',14);
    
    % Label
    xlabel('Time (min)','FontWeight','bold','FontSize',14);
    ylabel('Cortisol Amplitude (ug/dL)','FontWeight','bold','FontSize',14);
end
%------------------------------------------------------------------- Test 2
if Test5 == 1;
    % Cortisol Model Parameters
    CortPar.id      = 1;              % Simulated dataset id
    CortPar.gamma   = 0.5;            % Clearance parameter
    CortPar.sigma_A = 1;              % Amplitude sigma
    CortPar.mu_A    = 5;              % Amplitude Mean
    CortPar.sigma_w = 12;             % Interpulse interval sigma
    CortPar.mu_w    = 70;             % Interpulse interval mu
    CortPar.sigma_e = 0.5000;         % Error value (not used)
    CortPar.dt      = 1;              % Time increment
    CortPar.T       = 1440;           % Simulation period
    CortPar.alpha_w = 34.0278;        % Interpulse gamma distribution param
    CortPar.beta_w  = 2.0571;         % Interpulse gamma distribution param
    CortPar.w       = [];             % Interpulse intervals (not used)
    CortPar.A       = [];             % Interpulse amplitdues (not used)
    
    % Get time+simulated data and circadian parameter; plot
    BrownCortisolModel2001(CortPar);
end
%------------------------------------------------------------------- Test 6
if Test6 == 1;
    % Get time+simulated data and circadian parameter; plot
    
    % Run constants
    PLOT_SIMULATION = 1;
    
    % Cortisol Model Parameters
    CortPar.id      = 1;              % Simulated dataset id
    CortPar.gamma   = 0.25;            % Clearance parameter
    CortPar.sigma_A = 1;              % Amplitude sigma
    CortPar.mu_A    = 5;              % Amplitude Mean
    CortPar.sigma_w = 12;             % Interpulse interval sigma
    CortPar.mu_w    = 70;             % Interpulse interval mu
    CortPar.sigma_e = 0.5000;         % Error value (not used)
    CortPar.dt      = 1;              % Time increment
    CortPar.T       = 1440;           % Simulation period
    CortPar.alpha_w = 34.0278;        % Interpulse gamma distribution param
    CortPar.beta_w  = 2.0571;         % Interpulse gamma distribution param
    CortPar.w       = [];             % Interpulse intervals (not used)
    CortPar.A       = [];             % Interpulse amplitdues (not used)
    
    % Run simulations for a series of Gamma
    CortPar.gamma = 0.25;
    [t025 Y025 cortStruct025] = ...
        BrownCortisolModel2001(CortPar, 'Gamma = 0.25', PLOT_SIMULATION);
    
    CortPar.gamma = 0.5;
    [t050 Y050 cortStruct050] = ...
        BrownCortisolModel2001(CortPar, 'Gamma = 0.5', PLOT_SIMULATION);
    
    CortPar.gamma = 0.75;
    [t075 Y075 cortStruct075] = ...
        BrownCortisolModel2001(CortPar, 'Gamma = 0.75', PLOT_SIMULATION);
    
    CortPar.gamma = 1.0;
    [t100 Y100 cortStruct100] = ...
        BrownCortisolModel2001(CortPar, 'Gamma = 1.0', PLOT_SIMULATION);
    
    CortPar.gamma = 1.5;
    [t150 Y150 cortStruct150] = ...
        BrownCortisolModel2001(CortPar, 'Gamma = 1.5', PLOT_SIMULATION);
    
    CortPar.gamma = 2.5;
    [t250 Y250 cortStruct250] = ...
        BrownCortisolModel2001(CortPar, 'Gamma = 2.5', PLOT_SIMULATION);
    
    CortPar.gamma = 5.0;
    [t500 Y500 cortStruct500] = ...
        BrownCortisolModel2001(CortPar, 'Gamma = 5.0', PLOT_SIMULATION);
    
    % Save structure for future analysis
    save 'randomCortisolAmplitudeSimulation.mat' ...
        t025 Y025 cortStruct025 t050 Y050 cortStruct050 t075 Y075 cortStruct075...
        t100 Y100 cortStruct100 t150 Y150 cortStruct150 t250 Y250 cortStruct250 ...
        t500 Y500 cortStruct500;
    
end
%------------------------------------------------------------------- Test 7
if Test7 == 1;
    % Get time+simulated data and circadian parameter; plot
    
    % Run constants
    PLOT_SIMULATION = 1;
    
    % Cortisol Model Parameters
    CortPar.id      = 1;              % Simulated dataset id
    CortPar.gamma   = 0.25;            % Clearance parameter
    CortPar.sigma_A = 1;              % Amplitude sigma
    CortPar.mu_A    = 5;              % Amplitude Mean
    CortPar.sigma_w = 12;             % Interpulse interval sigma
    CortPar.mu_w    = 70;             % Interpulse interval mu
    CortPar.sigma_e = 0.5000;         % Error value (not used)
    CortPar.dt      = 1;              % Time increment
    CortPar.T       = 1440;           % Simulation period
    CortPar.alpha_w = 34.0278;        % Interpulse gamma distribution param
    CortPar.beta_w  = 2.0571;         % Interpulse gamma distribution param
    CortPar.w       = [];             % Interpulse intervals (not used)
    CortPar.A       = [];             % Interpulse amplitdues (not used)
    
    % Get randomly generated amplitude and interpulse intervals
    [t Y cortStruct] = BrownCortisolModel2001;
    CortPar.w = cortStruct.w_j;
    CortPar.A = cortStruct.A;
    
    % Run simulations for a series of Gamma
    CortPar.gamma = 0.25;
    [t025 Y025 cortStruct025] = ...
        BrownCortisolModel2001(CortPar, 'Gamma = 0.25', PLOT_SIMULATION);
    
    CortPar.gamma = 0.5;
    [t050 Y050 cortStruct050] = ...
        BrownCortisolModel2001(CortPar, 'Gamma = 0.5', PLOT_SIMULATION);
    
    CortPar.gamma = 0.75;
    [t075 Y075 cortStruct075] = ...
        BrownCortisolModel2001(CortPar, 'Gamma = 0.75', PLOT_SIMULATION);
    
    CortPar.gamma = 1.0;
    [t100 Y100 cortStruct100] = ...
        BrownCortisolModel2001(CortPar, 'Gamma = 1.0', PLOT_SIMULATION);
    
    CortPar.gamma = 1.5;
    [t150 Y150 cortStruct150] = ...
        BrownCortisolModel2001(CortPar, 'Gamma = 1.5', PLOT_SIMULATION);
    
    CortPar.gamma = 2.5;
    [t250 Y250 cortStruct250] = ...
        BrownCortisolModel2001(CortPar, 'Gamma = 2.5', PLOT_SIMULATION);
    
    CortPar.gamma = 5.0;
    [t500 Y500 cortStruct500] = ...
        BrownCortisolModel2001(CortPar, 'Gamma = 5.0', PLOT_SIMULATION);
    
    % Save structure for future analysis
    save 'constantCortisolAmplitudeSimulation.mat' ...
        t025 Y025 cortStruct025 t050 Y050 cortStruct050 t075 Y075 cortStruct075...
        t100 Y100 cortStruct100 t150 Y150 cortStruct150 t250 Y250 cortStruct250 ...
        t500 Y500 cortStruct500;
end


end

