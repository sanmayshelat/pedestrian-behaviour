%% Get Settings
% Change settings for network creation, agent properties
% Use this to create GUI later
% 
%%

%% Runtime Parameters
settings.runtime.setRandom = 1;
settings.runtime.plot = 1;

%% Location Capacities
settings.capacity.locationChair = 1;
settings.capacity.locationOutside = Inf;
settings.capacity.locationToilet = 1;
settings.capacity.serviceCoffee = 5;
settings.capacity.wormholeTurnstile = Inf;
settings.capacity.wormholeLift = Inf;
settings.capacity.wormholeStair = Inf;
settings.capacity.wormholeStair = Inf;

%% Service Parameters
% Note: The name of service has to be exactly as it is in the layers
settings.serviceLocations.serviceCoffee = 0.50; %[m]
settings.serviceLocations.servicePrinter = 0.50; %[m]
settings.serviceTime.serviceCoffee = 45; %[s]
settings.serviceTime.servicePrinter = 60; %[s]
%% Wormhole Parameters
% Note: The name of wormhole has to be exactly as it is in the layers
settings.wormholeLocations.wormholeTurnstile = 0.50; %[m]
settings.wormholeLocations.wormholeLift = 1; %[m]
settings.wormholeLocations.wormholeStair = 1; %[m]

settings.serviceTime.wormholeTurnstile = 1; %[s]
settings.serviceTime.wormholeLift = 1; %[s]
settings.serviceTime.wormholeStair = 1; %[s]

%% Navigation Graph Parameters
settings.navGraph.distClearance = 0.5; %[m] clearance distance from obstacle corners
settings.navGraph.minClearance = 0.25; %[m] minimum space required for movemement
settings.navGraph.wormholeClearance = 1; %[m] wormholes are a form of navPoint
settings.navGraph.halfConeAngle = 10; % [degrees]

%% Individual Walking Parameters
settings.indiWalk.numAgents = 26;

% General Parameters (Helbing, 2000)
settings.indiWalk.powerRepulsion = 2000; %% [N] A-A;A-O repulsion
settings.indiWalk.distInteractionAgent = 0.045; %0.08 %Helbing's parameter;% [m] characteristic length of repulsion forces
settings.indiWalk.distInteractionWall = 0.045; %[m] Allows better movements in narrow corridors
settings.indiWalk.relaxationTime = 0.5; % [s] acceleration time
settings.indiWalk.desiredSpeed = 1; % [m/s]
settings.indiWalk.desiredSpeedLow = 1.2; % [m/s]
settings.indiWalk.desiredSpeedHigh = 1.4; % [m/s]
settings.indiWalk.maxSpeed = 2.0; % [m/s] NOT FROM Helbing
% Fluctuation Gaussian Param. (Spaarnaij, 2015)
if settings.runtime.setRandom == 1
    settings.indiWalk.fluctuationSD = 0.01; % [m/s2]
    settings.indiWalk.fluctuationMean = 0;
else % if not random then don't fluctuate otherwise will get biased on one side
    settings.indiWalk.fluctuationSD = 0; % [m/s2]
    settings.indiWalk.fluctuationMean = 0;
end
% Ellipse (Campanella, 2016)
settings.indiWalk.ellipseFront = 0.85;
settings.indiWalk.ellipseBack = 1.25;
% Face Validation
settings.indiWalk.anticipationTime = 0.5; %[s]
settings.indiWalk.wayfindingPointClearance = 0.25; %[m]
settings.indiWalk.lateralAccDistAgent = 1.2; %[m]
settings.indiWalk.lateralAccDistEllipsey = 0.25; %[m]
settings.indiWalk.lateralAccApproachAngle = 10; %[degrees]

%% Walking Model Parameters
settings.walkModel.timeStep = 0.05; %[s]

%% Tactical Level
settings.tactical.timeStep = 1; %[min]
settings.tactical.plannedTimeStep = 15; %[min]
settings.tactical.preferredMeetingGap = 1*...
    settings.tactical.plannedTimeStep; %[min] Agents prefer to keep a gap between meetings
settings.tactical.mcGenTol = 10^-10;