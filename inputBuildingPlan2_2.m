%% Version Change
% 1. Redefined allTypes
% 2. Renamed according to the new type nomenclature
% 3. Added location definitions

%% Notes
% Takes files from "getFiles"

%% Building Plan Input
% Save in AutoCAD as .dxf; import to GIS in  EPSG: 32633; save as .shp
% AutoCAD file should use explode so that all edges are separated. Layers
% should be split according to: (1) Locations, (2) Obstacles, (3) Doors,
% (4) Wormhole, (5) Services
% See AutoCAD file for example on how to create drawings

buildingPlanEdges = shaperead(files.buildingPlan.edges);
% Coordinates
temp_x = vertcat(buildingPlanEdges.X);
temp_y = vertcat(buildingPlanEdges.Y);
% Type
allTypes = unique({buildingPlanEdges.Layer}');
typeObstacle = find(~cellfun('isempty',regexp(allTypes,'obstacle\w*')));
typeDoor = find(~cellfun('isempty',regexp(allTypes,'door\w*')));
typeWormhole = find(~cellfun('isempty',regexp(allTypes,'wormhole\w*')));
typeService = find(~cellfun('isempty',regexp(allTypes,'service\w*')));
temp_types = {buildingPlanEdges.Layer}';

% Doors
[~,temp_doors] = ismember(temp_types,{allTypes{typeDoor}});
for i = 1:numel(typeDoor)
    doors(i).x1 = temp_x(temp_doors==i,1);
    doors(i).y1 = temp_y(temp_doors==i,1);
    doors(i).x2 = temp_x(temp_doors==i,2);
    doors(i).y2 = temp_y(temp_doors==i,2);
    doors(i).type = allTypes{typeDoor(i)};
    allDoors{i} = [doors(i).x1 doors(i).y1 doors(i).x2 doors(i).y2];
    numDoors(i) = numel(doors(i).x1);
end

% Wormhole
[~,temp_wormholes] = ismember(temp_types,{allTypes{typeWormhole}});
for i = 1:numel(typeWormhole)
    wormholes(i).x1 = temp_x(temp_wormholes==i,1);
    wormholes(i).y1 = temp_y(temp_wormholes==i,1);
    wormholes(i).x2 = temp_x(temp_wormholes==i,2);
    wormholes(i).y2 = temp_y(temp_wormholes==i,2);
    wormholes(i).type = allTypes{typeWormhole(i)};
    allWormholes{i} = [wormholes(i).x1 wormholes(i).y1 wormholes(i).x2 wormholes(i).y2];
    numWormholes(i) = numel(wormholes(i).x1);
end

% Services
[~,temp_services] = ismember(temp_types,{allTypes{typeService}});
for i = 1:numel(typeService)
    services(i).x1 = temp_x(temp_services==i,1);
    services(i).y1 = temp_y(temp_services==i,1);
    services(i).x2 = temp_x(temp_services==i,2);
    services(i).y2 = temp_y(temp_services==i,2);
    services(i).type = allTypes{typeService(i)};
    allServices{i} = [services(i).x1 services(i).y1 services(i).x2 services(i).y2];
    numServices(i) = numel(services(i).x1);
end

% Obstacles
temp_obstacles = ismember(temp_types,{allTypes{typeObstacle}});
obstacles.x1 = temp_x(temp_obstacles,1);
obstacles.y1 = temp_y(temp_obstacles,1);
obstacles.x2 = temp_x(temp_obstacles,2);
obstacles.y2 = temp_y(temp_obstacles,2);
% Add Services,Wormholes as Obstacles
% obstacles.x1 = [obstacles.x1;vertcat(services(:).x1);vertcat(wormholes(:).x1)];
% obstacles.x2 = [obstacles.x2;vertcat(services(:).x2);vertcat(wormholes(:).x2)];
% obstacles.y1 = [obstacles.y1;vertcat(services(:).y1);vertcat(wormholes(:).y1)];
% obstacles.y2 = [obstacles.y2;vertcat(services(:).y2);vertcat(wormholes(:).y2)];
allObstacles = [obstacles.x1 obstacles.y1,obstacles.x2 obstacles.y2];
numObstacles = numel(obstacles.x1);

% Abstract Spaces
typeAbstractSpaces = find(~cellfun('isempty',regexp(allTypes,'abstract\w*')));
[~,temp_abstract] = ismember(temp_types,{allTypes{typeAbstractSpaces}});
for i = 1:numel(typeAbstractSpaces)
    abstractSpace(i).x1 = temp_x(temp_abstract==i,1);
    abstractSpace(i).y1 = temp_y(temp_abstract==i,1);
    abstractSpace(i).x2 = temp_x(temp_abstract==i,2);
    abstractSpace(i).y2 = temp_y(temp_abstract==i,2);
    abstractSpace(i).type = allTypes{typeAbstractSpaces(i)};
    allAbstractSpaces{i} = [abstractSpace(i).x1 abstractSpace(i).y1 abstractSpace(i).x2 abstractSpace(i).y2];
end

% Locations
buildingPlanLocations = shaperead(files.buildingPlan.locations);
% Coordinates
temp_x = vertcat(buildingPlanLocations.X);
temp_y = vertcat(buildingPlanLocations.Y);
% Type
allLocationTypes = unique({buildingPlanLocations.Layer}');
typeLocation = [numel(allTypes)+1:numel(allTypes)+numel(allLocationTypes)]';
temp_types = {buildingPlanLocations.Layer}';
for i = 1:numel(typeLocation)
    temp = strcmp(temp_types,allLocationTypes{i});
    locations(i).x = temp_x(temp);
    locations(i).y = temp_y(temp);
    locations(i).type = allLocationTypes{i};
    allLocations{i} = [locations(i).x locations(i).y];
    numLocations(i) = numel(locations(i).x);
end
allTypes = [allTypes;allLocationTypes];

%% Service Locations
for i = 1:numel(typeService)
    temp_x3 = (services(i).x1 + services(i).x2)/2;
    temp_y3 = (services(i).y1 + services(i).y2)/2;
    temp_d = settings.serviceLocations.(services(i).type);
    temp_x = services(i).x2-services(i).x1;
    temp_y = services(i).y2-services(i).y1;
    services(i).locationx = [];
    services(i).locationy = [];
    services(i).location1x = temp_x3 + (-1)*temp_y.*temp_d./sqrt(temp_x.^2+temp_y.^2);
    services(i).location1y = temp_y3 + (-1)*temp_x.*temp_d./sqrt(temp_x.^2+temp_y.^2);
    services(i).location2x = temp_x3 + temp_y.*temp_d./sqrt(temp_x.^2+temp_y.^2);
    services(i).location2y = temp_y3 + temp_x.*temp_d./sqrt(temp_x.^2+temp_y.^2);
end

%% Wormhole Locations
for i = 1:numel(typeWormhole)
    temp_x3 = (wormholes(i).x1 + wormholes(i).x2)/2;
    temp_y3 = (wormholes(i).y1 + wormholes(i).y2)/2;
    temp_d = settings.wormholeLocations.(wormholes(i).type);
    temp_x = wormholes(i).x2-wormholes(i).x1;
    temp_y = wormholes(i).y2-wormholes(i).y1;
    wormholes(i).locationx = [];
    wormholes(i).locationy = [];
    numWormholes(i) = numel(wormholes(i).x1);
    allWormholeLocations{i} = zeros(numWormholes(i)*2,2);
    for j = 1:2
        wormholes(i).locationx{j} = temp_x3 + ((-1)^j)*temp_y.*temp_d./sqrt(temp_x.^2+temp_y.^2);
        wormholes(i).locationy{j} = temp_y3 + ((-1)^j)*temp_x.*temp_d./sqrt(temp_x.^2+temp_y.^2);
        allWormholeLocations{i}(j:2:end) = [wormholes(i).locationx{j} wormholes(i).locationy{j}];
    end
end

%% Plot
%{
% plot check?
close all
figure
hold on
for i = 1:numel(obstacles.x1)
    plot([obstacles.x1(i) obstacles.x2(i)],[obstacles.y1(i) obstacles.y2(i)],'black')
end

for j = 1:numel(typeLocation)
    plot(locations(j).x,locations(j).y,'o','color',[1 0.4 0],'MarkerSize',3);
end

for j = 1:numel(typeService)
    plot(services(j).location1x,services(j).location1y,'+','color',[1 0.4 0],'MarkerSize',3);
    plot(services(j).location2x,services(j).location2y,'+','color',[1 0.4 0],'MarkerSize',3);
end

for j = 1:numel(typeWormhole)
    for i = 1:2
        plot(wormholes(j).locationx{i},wormholes(j).locationy{i},'+','color',[2/3 2/4 2/2],'MarkerSize',3);
    end
end

axis equal
%}

clear temp*