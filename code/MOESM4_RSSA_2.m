%% example usage: MOESM4_RSSA_2(4000, 40, 40, -50, 0.1, 0.0001)

function [T, Dynamics] = MOESM4_RSSA_2(tMax, Mtot, Ntot, V, delta, dT)
% - tMax: Maximum simulation time
% - Mtot, Ntot: Total number of calcium and potassium channels
% - V: Voltage
% - delta: Fluctuation rate (0-1)
% - dT: Time step

%% STARTUP

tic % measure simulation time

% startup simulation parameters
T = (0:dT:tMax)'; % Time vector
initialState = [Mtot, ceil(Ntot/2), 0, ceil(Ntot/2)]; % start: all Ca closed, K half and half
initialV = V;     % initial voltage

Dynamics = nan(length(T),length(initialState)+1); % vector storing the dynamics
maxAllowedSteps = tMax*500000;

% pre-generation of some random numbers
initialLenght = 1000;
randV = rand(1,initialLenght);
nRandVResets = 1;
usedRandomNumbers = 0;

% Define reaction matrices: vMinus, vPlus, v
% state [Ca_closed, K_closed, Ca_open, K_open]
vMinus = [1, 0, 0, 0;  % Ca opening
          0, 0, 1, 0;  % Ca closing
          0, 1, 0, 0;  % K opening
          0, 0, 0, 1]; % K closing
vPlus =  [0, 0, 1, 0;
          1, 0, 0, 0;
          0, 0, 0, 1;
          0, 1, 0, 0];

v = vPlus - vMinus;

% Model Parameters
phi_m = 0.4;  % scaling factor (Ca)
phi_n = 0.04; % scaling factor (K)
va = -1.2; vb = 18; % Ca half activation point and slope
vc = 2; vd = 30;    % K half activation point and slope

% Functions for gating variables
global Iapp; Iapp=@(t)100; % applied current
    % Calcium
global xi_m; xi_m=@(v)(v-va)/vb; % scaled argument for m-gate input voltage
global minf; minf=@(v)0.5*(1+tanh(xi_m(v))); % m-gate activation function
global tau_m; tau_m=@(v)1./(phi_m*cosh(xi_m(v)/2)); % m-gate time constant
global alpha_m; alpha_m=@(v)(minf(v)./tau_m(v));   % alpha Ca opening rate
global beta_m; beta_m=@(v)((1-minf(v))./tau_m(v)); % beta Ca closing rate
    % Potassium
global xi_n; xi_n=@(v)(v-vc)/vd; % scaled argument for n-gate input
global ninf; ninf=@(v)0.5*(1+tanh(xi_n(v))); % n-gate activation function
global tau_n; tau_n=@(v)1./(phi_n*cosh(xi_n(v)/2));  % n-gate time constant
global alpha_n; alpha_n=@(v)(ninf(v)./tau_n(v));   % alpha K opening rate
global beta_n; beta_n=@(v)((1-ninf(v))./tau_n(v)); % beta K closing rate

% Reaction rate vector
c = [alpha_m(initialV), beta_m(initialV), alpha_n(initialV), beta_n(initialV)];


%% INITIALIZE SIMULATION VARIABLES
i = 1;    % iteration index
T(i) = 0; % T0 = 0
Dynamics(i, :) = [initialState, initialV];
currentTime = 0;
currentState = initialState;
currentV = initialV;

% Fluctuation bounds (state bounds)
currentStateUp = initialState + round(delta * initialState);
currentStateDown = initialState - round(delta * initialState);
currentStateUp(currentStateDown < 0) = currentStateUp(currentStateDown < 0) + abs(currentStateDown(currentStateDown < 0));
currentStateDown(currentStateDown < 0) = currentStateDown(currentStateDown < 0) + abs(currentStateDown(currentStateDown < 0));

% reaction propensities
aUp = zeros(size(c));
aDown = zeros(size(c));

for j = 1:length(c)
        aUp(j) = computeReactionPropensity(vMinus,c,currentStateUp,j);
        aDown(j) = computeReactionPropensity(vMinus,c,currentStateDown,j);
end
a0Up = sum(aUp);

%% RSSA LOOP

% simulation statistics
nSimulationSteps = 0;
nFastAccept = 0;
nSlowAccept = 0;
nRejections = 0;
nFluctuationIntUpdates = 0;

% simulation loop
while currentTime < tMax && nSimulationSteps < maxAllowedSteps % until tMax or max steps
    stateConsistency = true;
        while (stateConsistency && currentTime < tMax && nSimulationSteps < maxAllowedSteps)
            u = 1;
            reactionAccepted = false;
            if (a0Up > 0) % this if clause should avoid searching for a reaction if no reactions can be selected...
                while (~reactionAccepted)
                    % extraction of three unused random numbers from randV
                    if (usedRandomNumbers + 3 > length(randV))
                        % generation of new random numbers if we reached the end of the array...
                        randV = rand(1,initialLenght);
                        usedRandomNumbers = 0;
                        nRandVResets = nRandVResets + 1;
                    end
                    r1 = randV(usedRandomNumbers+1);
                    r2 = randV(usedRandomNumbers+2);
                    r3 = randV(usedRandomNumbers+3);
                    usedRandomNumbers = usedRandomNumbers + 3;

                    % selection of a reaction candidate by using the upper bound of reaction propensities
                    mu = 1;
                    while sum(aUp(1:mu)) < r1*a0Up
                        mu = mu + 1;
                    end

                    % rejection-based strategy to accept/reject the reaction candidate
                    if (r2 <= aDown(mu)/aUp(mu))
                        reactionAccepted = true; % fast acceptance
                        nFastAccept = nFastAccept + 1;
                    else
                        % computation of the "real" propensity
                        a = computeReactionPropensity(vMinus,c,currentState,mu);
                        if (r2 <= a/aUp(mu))
                            reactionAccepted = true; % slow acceptance
                            nSlowAccept = nSlowAccept + 1;
                        else
                            nRejections = nRejections + 1; % if this is executed, the candidate reaction is rejected
                        end
                    end
                    u = u*r3; % in any case I update u to update the final computation of tau
                end

                % computation of the next tau
                tau = (1/a0Up)*log(1/u);
                
                % update voltage; continuous variable
                % takes in input the time interval [t, t+tau]
                % and the current voltage, state and channel numbers
                [~, uOut] = ode23(@dudtfunc, [currentTime, currentTime+tau], ...
                           [currentV, currentState(3), currentState(4), Mtot, Ntot]);
                currentV = uOut(end, 1); % extract voltage

                % update rates with new voltage
                c = [alpha_m(currentV), beta_m(currentV), alpha_n(currentV), beta_n(currentV)];

                % update time and states; discrete variables
                currentTime = currentTime + tau; % new time
                currentState = currentState + v(mu,:); % new state, applied reaction mu

                % saving of the current state to the dynamics timeseries if needed
                if currentTime >= T(i)+dT
                    i = i+1;
                    T(i) = currentTime;
                    Dynamics(i,:) = [currentState, currentV];
                end
            else
                % if nothing can be done, we just go ahead unitl tMax
                % without further updating the current state
                % PS: I save two times the state to be sure of capturing
                % the steady state condition of the model
                i = i+1;
                T(i) = currentTime;
                Dynamics(i,:) = [currentState, currentV];
                currentTime = tMax;
                i = i+1;
                T(i) = currentTime;
                Dynamics(i,:) = [currentState, currentV];
            end
            
            % update of the number of simulation steps
            nSimulationSteps = nSimulationSteps + 1;
            
            % check if the updated state is still consistent with its fluctuation interval
            stateConsistency = checkStateConsistency(currentState,currentStateUp,currentStateDown);

        end
        % if state is not consistent

        % update of the fluctation interval of the model state
        currentStateUp = currentState + round(delta*currentState); % state upper bound
        currentStateDown = currentState - round(delta*currentState); % state lower bound

        % check for negative state lower bounds (in case we shift the
        % flucutation interval to have the lower bound set to zero)
        currentStateUp(currentStateDown < 0) = currentStateUp(currentStateDown < 0) + abs(currentStateDown(currentStateDown < 0));
        currentStateDown(currentStateDown < 0) = currentStateDown(currentStateDown < 0) + abs(currentStateDown(currentStateDown < 0));
        
        % update of reaction propensities (up and down)
        aUp = zeros(size(c));
        aDown = zeros(size(c));
        for j = 1:length(c)
            aUp(j) = computeReactionPropensity(vMinus,c,currentStateUp,j);
            aDown(j) = computeReactionPropensity(vMinus,c,currentStateDown,j);
        end

        % update of the sum of the upper bound of reaction propensities
        a0Up = sum(aUp);
        
        nFluctuationIntUpdates = nFluctuationIntUpdates + 1;
    end
    
    % cut of the residual part of timeseries that remained set to NaN
    T = T(~isnan(Dynamics(:,1)));
    Dynamics = Dynamics(~isnan(Dynamics(:,1)),:);
    
    if (nSimulationSteps >= maxAllowedSteps)
        % printing of a warning message to tell to the user that the simulation has been stopped in advance
        disp(' '); % to print an empty line
        disp(['WARNING: the simulation reached the maximum allowed number of simulation steps (' num2str(maxAllowedSteps) ')!']);
        disp(' '); % to print an empty line
    end
    
    disp(['Total number of computed simulation steps: ' num2str(nSimulationSteps)]);
    disp(['Total number of used random numbers: ' num2str((nRandVResets-1)*initialLenght+usedRandomNumbers)]);
    disp(['Number of fluctuation interval updates: ' num2str(nFluctuationIntUpdates) ' (' num2str(nFluctuationIntUpdates/nSimulationSteps*100) '%)']);
    disp(['Number of fast acceptance steps: ' num2str(nFastAccept) ' (' num2str(nFastAccept/(nFastAccept+nSlowAccept+nRejections)*100) '%)']);
    disp(['Number of slow acceptance steps: ' num2str(nSlowAccept) ' (' num2str(nSlowAccept/(nFastAccept+nSlowAccept+nRejections)*100) '%)']);
    disp(['Number of rejection steps: ' num2str(nRejections) ' (' num2str(nRejections/(nFastAccept+nSlowAccept+nRejections)*100) '%)']);
    toc % display simulation runtime

%% PLOTTING

% image 1
figure(1)
subplot(3,1,1),plot(T,Dynamics(:, 3)),ylabel('M','FontSize',20),set(gca,'FontSize',20)
    xlim([0, 4000]) % M(t), Calcium dynamics in time
subplot(3,1,2),plot(T,Dynamics(:, 4)),ylabel('N','FontSize',20),set(gca,'FontSize',20)
    xlim([0, 4000]) % N(t), Potassium dynamics in time
subplot(3,1,3),plot(T,Dynamics(:, 5)),xlabel('Time','FontSize',20)
    ylabel('V','FontSize',20),set(gca,'FontSize',20)
    xlim([0, 4000]) % V(t), Voltage dynamics in time

% image 2
figure(2)
subplot(6,1,1:6),plot3(Dynamics(:, 5),Dynamics(:, 3),Dynamics(:, 4),'.-'),xlabel('V','FontSize',20)
    ylabel('M','FontSize',20),zlabel('N','FontSize',20),set(gca,'FontSize',20)
grid on, rotate3d, shg % 3D plot --> [V,M,N]
end

%% FUNCTIONS

% computes the reaction propensity
function a = computeReactionPropensity(vMinus,c,state,reactionIndex)
    a = c(reactionIndex);
    if (sum(vMinus(reactionIndex,:) > 0))
        for i = 1:length(state)
            % the following if clauses allow to limit the usage of the nchoosek function (needed only when vMinus(ReactionIndex,i) > 1)
            if vMinus(reactionIndex,i) == 1
                a = a*state(i);
            end
            if vMinus(reactionIndex,i) > 1
                if (state(i) >= vMinus(reactionIndex,i))
                    a = a*nchoosek(state(i),vMinus(reactionIndex,i)); % nchoosek(n,k) returns the binomial coefficient (n k)
                else
                    a = a*0; % no available reactants...
                end
            end
        end
    end
end

% checks if we need to recompute the bounds or not, returns boolean value
function test = checkStateConsistency(currentState,currentStateUp,currentStateDown)
    test = true;
    for i = 1:length(currentState)
        test = test && currentState(i) >= currentStateDown(i) && currentState(i) <= currentStateUp(i);
    end
end

% computes the ODE of the voltage during an event, between (t, t+tau)
function dudt=dudtfunc(t,u)
    global Iapp
    vK = -84; vL = -60; vCa = 120; % reversal potentials of K, Leak, Ca
    gK = 8; gL = 2; gCa = 4.4;     % maximum conductance of K, Leak, Ca
    C = 20; % Capacitance

    % u = [voltage, Ca_open, K_open, Ca_tot, K_tot] -> 
    % [currentV, currentState(3), currentState(4), Mtot, Ntot]

    % Follows Kirchhoff’s current conservation law with Ca and K ion channels and leak
    dudt = [(Iapp(t) - gCa * (u(2) / u(4)) * (u(1) - vCa) ...
           - gL * (u(1) - vL) - gK * (u(3) / u(5)) * (u(1) - vK)) / C;
           u(2); % the discrete states remain the same, the voltage changes during the transition
           u(3);
           u(4); % the total number of channels remains the same
           u(5)];
end