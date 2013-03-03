function varargout = BrownCortisolModel2001(varargin)
% BrownCortisolModel2001 A function for implementing the cortisol model 
% published in 2001.  The function is adapted from code provided by David 
% Nguyen for a related reversible jump Monte Carlo project.
%
% The cortisol model details can be found in:
%
%   Brown EN, Meehan PM, Demster AP (2001) A stochastic differential 
%   equation model of diurnal cortisol patterns. American Journal of 
%   Physiology - Endocrinology and Metabolism 280: E450-E461.
% 
% Additional information about the RJMC project can be found in the
% following two references.
%
%   Nguyen DP, Frank LM, Brown EN (2003) An application of reversible-jump
%   Markov chain Monte Carlo to spike classification of multi-unit 
%   extracellular recordings. Network: Computation in Neural Systems 
%   14: 61-82.
%  
%   Dean II DA (2011) Integrating Formal Language Theory with 
%   Mathematical Modeling to Solve Computational Issues in Sleep and 
%   Circadian Applications: University of Massachusetts. 1-239 p.
%
%
% Function Protypes:
%     BrownCortisolModel2001;
%     BrownCortisolModel2001(CortPar);
%     BrownCortisolModel2001(CortPar,titleStr);
%     Y = BrownCortisolModel2001; 
%     [t Y] = BrownCortisolModel2001;
%     [t Y cortStruct] = BrownCortisolModel2001;
%     [t Y cortStruct] = BrownCortisolModel2001;
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

%----------------------------------- Program Constants and Model Parameters
% Program Constants
RAND_GEN_FACTOR = 50;  % An adhoc constant to insure that there is at least
                       % one secretion time beyond the simulation period
PLOT_SIMULATION = 0;   % If set to true, ampltidues and simulation are 
                       % plotted     
titleStr = '';         % Graph title string                       

% Cortisol Model Parameters
CortPar.id      = 1;                 % Simulated dataset id
CortPar.gamma   = 1.0;               % Clearance parameter
CortPar.sigma_A = 1;                 % Amplitude sigma
CortPar.mu_A    = 5;                 % Amplitude Mean
CortPar.sigma_w = 12;                % Interpulse interval sigma
CortPar.mu_w    = 70;                % Interpulse interval mu
CortPar.sigma_e = 0.5000;            % Error value (not used)
CortPar.dt      = 1;                 % Time increment
CortPar.T       = 1440;              % Simulation period
CortPar.alpha_w = 34.0278;           % Interpulse gamma distribution param
CortPar.beta_w  = 2.0571;            % Interpulse gamma distribution param
CortPar.w       = [];                % Interpulse intervals (not used)
CortPar.A       = [];                % Interpulse amplitdues (not used)
                
% Circadian Variation Parameters
mu_zeta      = [6.10 -4.75 3.93 -3.76 -2.53]';
sigma_zeta   = [ [ 1.937  0.423  -0.095   0.573  -0.214 ]; ...
                 [ 0.423  1.253  -0.197   0.514   0.296 ]; ...
                 [-0.095 -0.197   0.631  -0.254  -0.125 ]; ...
                 [ 0.573  0.514  -0.254   1.276  -0.123 ]; ...
                 [-0.214  0.296  -0.125  -0.123   0.509 ] ];      
%------------------------------------------------------------ Process Input
if nargin == 1
    CortPar = varargin{1};
elseif nargin == 2
    CortPar = varargin{1};
    titleStr = varargin{2};
elseif nargin == 3
    CortPar = varargin{1};
    titleStr = varargin{2};
    PLOT_SIMULATION = varargin{3};
end          
%---------------------------------------------- Initialize Output Variables
% Initialize return varaibles
t = [0:CortPar.dt:CortPar.T]';
Y = zeros(size(t));

%----------------------------------------------------- Draw secretion times
if or(isempty(CortPar.w), isempty(CortPar.A))
    % Randomly draw secretion times
    alpha_w = (CortPar.mu_w/CortPar.sigma_w)^2;
    beta_w = CortPar.sigma_w^2/CortPar.mu_w;
    w_j = gamrnd(alpha_w, beta_w, floor(CortPar.T/RAND_GEN_FACTOR),1);
    w_j_start = cumsum(w_j);

    % Prune list to fit in time 
    indexes = find(w_j_start < CortPar.T);    % Inexes within tiem frame
    w_j = w_j(indexes);                       % Interpulse intervals
    w_j_start = floor(w_j_start(indexes));    % Align secretions with 
                                              %    time incements
    
    % Initialize amplitude array                                          
    A = ones(size(w_j_start));                % Define amplitude size
else
    % Copy interpulse intervals
    w_j = CortPar.w;                          % Interpulse intervals
    w_j_start = cumsum(w_j);                  % Align secretions with 
                                              %    time incements
end
%------------------------------------------------ Draw circadian parameters
if or(isempty(CortPar.w), isempty(CortPar.A))
    % Initialize amplitdue                                           
    A = ones(size(w_j_start));                % Define amplitude size   

    % Draw random circadian amplitude 
    obj_zeta   = gmdistribution(mu_zeta',sigma_zeta);
    zeta_rnd   = random(obj_zeta);
    [mu twoHarmStruct] = two_harmonic_mean(zeta_rnd, t/60);
    A = mu(w_j_start+1);
else
    % Copy amplitude to local variable
    A = CortPar.A;
    
    % Define harmonic structure to nill
    twoHarmStruct = {};
end
%----------------------------------------- Simulate secretion and clearance
% Determine maximum secretion width to simulate
num_pulses = length(w_j);

A_t = zeros();
dt = CortPar.dt;
gamma = CortPar.gamma;
rise(1) = 0.02;

% Identify initial secretion
J = 1;

% Simulate each point
for r = 2:length(t)
    % Detect a secretion
    if J <= num_pulses
        % Process stochastic secretions
        if w_j_start(J) <= t(r)
            Y(r) =Y(r-1)*exp(-gamma/60*dt) + A(J);
            J = J + 1;
        else
            % Not a secretion time
            Y(r) =Y(r-1)*exp(-gamma/60*dt);
        end
    else
        % No more secretions to proces
        Y(r) =Y(r-1)*exp(-gamma/60*dt);
    end    
end
Y = Y;

%------------------------------------------------------------- Plot Results
if or(PLOT_SIMULATION == 1, nargout == 0)
    %------------------ Disply Amplitude
    fid = figure('InvertHardcopy','off','Color',[1 1 1]);
    stem1 = stem(w_j_start,A,'LineWidth',2);
    baseline1 = get(stem1,'BaseLine');
    set(baseline1,'LineWidth',2);

    % Set Title
    title(titleStr, 'FontWeight','bold', 'FontSize',14);
    
    % Adjust time axis
    v = axis;
    v(2) = CortPar.T; 
    axis(v);

    % Reformat Axis
    set(gca, 'LineWidth',2);
    set(gca, 'FontWeight','bold');
    set(gca, 'FontSize',14);
    
    % Label
    xlabel('Time (min)','FontWeight','bold','FontSize',14);
    ylabel('Cortisol Amplitude (ug/dL)','FontWeight','bold','FontSize',14);

    %------------------ Display Simulated Cortisol
    fid = figure('InvertHardcopy','off','Color',[1 1 1]);
    plot (t, Y,'LineWidth',2);

    % Set Title
    title(titleStr, 'FontWeight','bold', 'FontSize',14);
    
    % Adjust time axis
    v = axis;
    v(2) = CortPar.T; 
    axis(v);
    
    % Reformat Axis
    set(gca, 'LineWidth',2);
    set(gca, 'FontWeight','bold');
    set(gca, 'FontSize',14);
    
    % Label
    xlabel('Time (min)','FontWeight','bold','FontSize',14);
    ylabel('Cortisol Amplitude (ug/dL)','FontWeight','bold','FontSize',14);
end
%------------------------------------------------ Create Cortisol Structure
% Store information tp reconstruct signal
cortStruct.twoHarmStruct = twoHarmStruct;
cortStruct.w_j = w_j;
cortStruct.w_j_start = w_j_start;
cortStruct.A = A ;
cortStruct.CortPar = CortPar ;
cortStruct.mu_zeta = mu_zeta ;
cortStruct.sigma_zeta = sigma_zeta ;
%------------------------------------------------ Generate Output Structure
if nargout == 0
    varargout = {};
elseif nargout == 1
    varargout{1} = Y;
elseif nargout == 2
    varargout{1} = t;
    varargout{2} = Y;
elseif nargout == 3
    varargout{1} = t;
    varargout{2} = Y;
    varargout{3} = cortStruct;
end

end
%-------------------------------------------------------- Two Harmonic Mean
function varargout = two_harmonic_mean(zeta_coeff, t)
% two_harmonic_mean create two harmonic circadian amplitude
%
% Input:
%     zeta_coeff : Array of 5 variables (c0, c1, d1, c2, d2)
%              t : Time array in hours
%
% Ouput:
%             mu : Circadian amplitude shifted so min(mu) = 0
%  twoHarmStruct : Structure that contains variables and equation string
%
%
% 
% Program Constants
DEBUG = 1;
NUMBER_OF_HARMONICS = 2;
SHIFT_POSITIVE = 1;
c0 = 1;
c1 = 2;
d1 = 3;
c2 = 4;
d2 = 5;

% Function constant
zeta_coeff(2:5);
w = 2*pi/24;
mu = zeta_coeff(1)+ ...
     [ cos(t*w)  sin(t*w)  cos(t*w*2)  sin(t*w*2)]*zeta_coeff(2:5)';
 
% Shift MU to be positive
mu = mu - min(mu);

% Create output structure
if nargout == 1
   varargout{1} = mu;
elseif nargout == 2
   % Return Amplitudes 
   varargout{1} = mu;
  
   % Create structure to store randomly drawn varaibles
   twoHarmStruct.c0 = 1;
   twoHarmStruct.c1 = 2;
   twoHarmStruct.d1 = 3;
   twoHarmStruct.c2 = 4;
   twoHarmStruct.d2 = 5;
   
   % Return two harmonic parameters
   varargout{2} = twoHarmStruct;
else
   varargout = {};
end
    
end


