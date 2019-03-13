%% Version Change
% 1. Everything NOT in one graph: 3 Graphs: navigation,location to
% navigation (asymmetric), and wormhole connections
% 2. All graphs have their own indexing

%% Notes
% Takes building plan from "inputBuildingPlan2_2" and onwards
% Takes settings from "getSettings"

%% Find corners
removePrecisionError = 1e-8;
temp_x11 = abs(repmat(obstacles.x1,1,numObstacles)-repmat(obstacles.x1',numObstacles,1))<removePrecisionError;
temp_x22 = abs(repmat(obstacles.x2,1,numObstacles)-repmat(obstacles.x2',numObstacles,1))<removePrecisionError;
temp_x12 = abs(repmat(obstacles.x1,1,numObstacles)-repmat(obstacles.x2',numObstacles,1))<removePrecisionError;
temp_x21 = abs(repmat(obstacles.x2,1,numObstacles)-repmat(obstacles.x1',numObstacles,1))<removePrecisionError;

temp_y11 = abs(repmat(obstacles.y1,1,numObstacles)-repmat(obstacles.y1',numObstacles,1))<removePrecisionError;
temp_y22 = abs(repmat(obstacles.y2,1,numObstacles)-repmat(obstacles.y2',numObstacles,1))<removePrecisionError;
temp_y12 = abs(repmat(obstacles.y1,1,numObstacles)-repmat(obstacles.y2',numObstacles,1))<removePrecisionError;
temp_y21 = abs(repmat(obstacles.y2,1,numObstacles)-repmat(obstacles.y1',numObstacles,1))<removePrecisionError;

temp_11 = triu(and(temp_x11,temp_y11))-eye(numObstacles); %to remove (i,i) = 1
temp_22 = triu(and(temp_x22,temp_y22))-eye(numObstacles);
temp_12 = triu(and(temp_x12,temp_y12));
temp_21 = triu(and(temp_x21,temp_y21));

% Corner struct: p2 is the corner point
[temp_r,temp_c] = find(temp_11);
corner.p1x = obstacles.x2(temp_r); corner.p1y = obstacles.y2(temp_r);
corner.p2x = obstacles.x1(temp_r); corner.p2y = obstacles.y1(temp_r);
corner.p3x = obstacles.x2(temp_c); corner.p3y = obstacles.y2(temp_c);
corner.edge1 = temp_r; corner.edge2 = temp_c;
[temp_r,temp_c] = find(temp_22);
corner.p1x = [corner.p1x; obstacles.x1(temp_r)]; corner.p1y = [corner.p1y; obstacles.y1(temp_r)];
corner.p2x = [corner.p2x; obstacles.x2(temp_r)]; corner.p2y = [corner.p2y; obstacles.y2(temp_r)];
corner.p3x = [corner.p3x; obstacles.x1(temp_c)]; corner.p3y = [corner.p3y; obstacles.y1(temp_c)];
corner.edge1 = [corner.edge1; temp_r]; corner.edge2 = [corner.edge2; temp_c];
[temp_r,temp_c] = find(temp_12);
corner.p1x = [corner.p1x; obstacles.x2(temp_r)]; corner.p1y = [corner.p1y; obstacles.y2(temp_r)];
corner.p2x = [corner.p2x; obstacles.x1(temp_r)]; corner.p2y = [corner.p2y; obstacles.y1(temp_r)];
corner.p3x = [corner.p3x; obstacles.x1(temp_c)]; corner.p3y = [corner.p3y; obstacles.y1(temp_c)];
corner.edge1 = [corner.edge1; temp_r]; corner.edge2 = [corner.edge2; temp_c];
[temp_r,temp_c] = find(temp_21);
corner.p1x = [corner.p1x; obstacles.x1(temp_r)]; corner.p1y = [corner.p1y; obstacles.y1(temp_r)];
corner.p2x = [corner.p2x; obstacles.x2(temp_r)]; corner.p2y = [corner.p2y; obstacles.y2(temp_r)];
corner.p3x = [corner.p3x; obstacles.x2(temp_c)]; corner.p3y = [corner.p3y; obstacles.y2(temp_c)];
corner.edge1 = [corner.edge1; temp_r]; corner.edge2 = [corner.edge2; temp_c];

%% Remove incorrect corners
temp_remove = ... % remove corners formed by collinear obstacles
    (([corner.p1y]' == [corner.p2y]') &...
    ([corner.p3y]' == [corner.p2y]')) |...
    (([corner.p1x]' == [corner.p2x]') &...
    ([corner.p3x]' == [corner.p2x]'));
corner.p1x = corner.p1x(~temp_remove); corner.p1y = corner.p1y(~temp_remove);
corner.p2x = corner.p2x(~temp_remove); corner.p2y = corner.p2y(~temp_remove);
corner.p3x = corner.p3x(~temp_remove); corner.p3y = corner.p3y(~temp_remove);
corner.edge1 = corner.edge1(~temp_remove); corner.edge2 = corner.edge2(~temp_remove);
% corner = corner(~temp_remove);
allCorners = [vertcat(corner.p1x) vertcat(corner.p1y),...
    vertcat(corner.p2x) vertcat(corner.p2y),...
    vertcat(corner.p3x) vertcat(corner.p3y)];
numCorners = numel(corner.p2x);

%% Find angle bisector
temp_x1 = corner.p2x-corner.p1x;
temp_y1 = corner.p2y-corner.p1y;
temp_x2 = corner.p2x-corner.p3x;
temp_y2 = corner.p2y-corner.p3y;
temp_1 = sqrt(temp_x1.^2+temp_y1.^2);
temp_2 = sqrt(temp_x2.^2+temp_y2.^2);
temp_x3 = temp_x1./temp_1 + temp_x2./temp_2;
temp_y3 = temp_y1./temp_1 + temp_y2./temp_2;
temp_3 = sqrt(temp_x3.^2+temp_y3.^2);
temp_x3 = temp_x3./temp_3; % unit angle bisector
temp_y3 = temp_y3./temp_3;

%% Find nav points
corner.navPointx = corner.p2x + temp_x3.*settings.navGraph.distClearance;
corner.navPointy = corner.p2y + temp_y3.*settings.navGraph.distClearance;
% plot check?
plotWalkingArea2

%% Min Clearance & Clear line of sight Initial Calculations
corner.clearance = sqrt(sum((corner.p2x-corner.navPointx).^2 + ...
    (corner.p2y-corner.navPointy).^2,2));
temp_obstacles = [corner.navPointx corner.navPointy,...
    corner.p2x corner.p2y]; % Obstacles b/w nav points and corresponding corners
[temp] = lineSegmentIntersect(allObstacles,temp_obstacles);
temp_intersection = temp.intAdjacencyMatrix;
temp_x = temp.intMatrixX;
temp_y = temp.intMatrixY;

%% Min Clearance
% If min clearance to other obstacles not available set distance between nav
% point and corresponding corner to min clearance
temp_distNavPoint = -temp.intNormalizedDistance2To1*settings.navGraph.distClearance;

temp_distNavPoint(temp_distNavPoint<0) = Inf; % do not consider nav point/obst edge where the edge (or its extension) is not on in this order: intersection pt., nav pt., corresponding corner
temp_distNavPoint(0>temp.intNormalizedDistance1To2 | temp.intNormalizedDistance1To2>1) = Inf; % do not consider nav pt./obst edge if the intersection is with extension of edge
temp_distNavPoint(corner.edge1 + ([1:numCorners]'-1).*numObstacles) = Inf; % intersection with own edge excluded
temp_distNavPoint(corner.edge2 + ([1:numCorners]'-1).*numObstacles) = Inf; % intersection with own edge excluded

temp_indexNavPoint = temp_distNavPoint<settings.navGraph.minClearance;
temp_indexNavPoint = sum(temp_indexNavPoint)>0;
corner.navPointx = corner.p2x + ...
    temp_x3.*(temp_indexNavPoint'.*settings.navGraph.minClearance + ...
    ~temp_indexNavPoint'.*settings.navGraph.distClearance);
corner.navPointy = corner.p2y + ...
    temp_y3.*(temp_indexNavPoint'.*settings.navGraph.minClearance + ...
    ~temp_indexNavPoint'.*settings.navGraph.distClearance);
% plot check?
plotWalkingArea2

%% Clear line of sight between nav point and corresponding corner
% If an edge lies between a nav point and its corresponding corner, shift
% the nav point to the halfway point between the nearest edge and
% corresponding corner
temp_distx = temp_x - repmat(corner.p2x',size(allObstacles,1),1);
temp_disty = temp_y - repmat(corner.p2y',size(allObstacles,1),1);
temp_intersection(corner.edge1 + ([1:numCorners]'-1).*numObstacles) = 0; % intersection with own edge excluded
temp_intersection(corner.edge2 + ([1:numCorners]'-1).*numObstacles) = 0; % intersection with own edge excluded
temp_distx(~temp_intersection) = Inf;
temp_disty(~temp_intersection) = Inf;
temp_dist = (temp_distx.*temp_intersection).^2 + (temp_disty.*temp_intersection).^2;
[~,temp_index] = min(temp_dist);
temp_navPointx = (corner.p2x + ...
    temp_x([0:size(temp_x,1):numel(temp_x)-size(temp_x,1)]...
    +temp_index)')/2;
temp_navPointy = (corner.p2y + ...
    temp_y([0:size(temp_y,1):numel(temp_y)-size(temp_y,1)]...
    +temp_index)')/2;
corner.navPointx(~isnan(min(temp_dist))') = ...
    temp_navPointx(~isnan(min(temp_dist))');
corner.navPointy(~isnan(min(temp_dist))') = ...
    temp_navPointy(~isnan(min(temp_dist))');
corner.clearance = sqrt(sum((corner.p2x-corner.navPointx).^2 + ...
    (corner.p2y-corner.navPointy).^2,2));
% plot check?
plotWalkingArea2

%% Delete nav points without min clearance
corner.navPointx(corner.clearance<settings.navGraph.minClearance) = NaN;
corner.navPointy(corner.clearance<settings.navGraph.minClearance) = NaN;
% plot check?
plotWalkingArea2

%% Merge Nav Points
allNavPoints = [corner.navPointx corner.navPointy];
allClearance = corner.clearance;

temp_x = repmat(corner.p2x,1,numCorners);
temp_y = repmat(corner.p2y,1,numCorners);
temp_dist = sqrt((temp_x-corner.navPointx').^2 + ...
    (temp_y-corner.navPointy').^2); % dist from (row) corner points to (col) nav points (19/7/2017)
temp = temp_dist<=corner.clearance'; %USED: if a corner point is nearer to a nav point than its (the nav point's (19/8/2017)) corresponding corner a merge takes place NOT USED:"...if another corner point lies closer to the corresponding corner than the navigation point itselft the corresponding corner is added to this point..."
temp_points = allNavPoints;
for i = 1:numCorners
    temp_index = ~any(temp-temp(:,i)); %"...if...true vice versa...points are merged..."
    if sum(temp_index)>1
        temp_points(i,:) = sum(allNavPoints(temp_index,:))/2;
    end
end
allNavPoints = temp_points;

%% Min Clearance between Nav Points
temp_x = repmat(allNavPoints(:,1),1,size(allNavPoints,1));
temp_y = repmat(allNavPoints(:,2),1,size(allNavPoints,1));
temp_dist = sqrt((temp_x-allNavPoints(:,1)').^2 + ...
    (temp_y-allNavPoints(:,2)').^2); % dist from (row) corner pts to (col) nav pts (19/7/2017)
temp = temp_dist<=allClearance'; % minimum clearance between points: (if a nav point is closer to another nav point than that nav point is to its corresponding corner then merge (19/8/2017)) (the paper says not to use this but let's see
temp_points = allNavPoints;
for i = 1:numCorners
    temp_index = ~any(temp-temp(:,i));
    if sum(temp_index)>1
        temp_points(i,:) = sum(allNavPoints(temp_index,:))/2;
    end
end
allNavPoints = temp_points;

%% Remove NavPoints near Wormhole Locations
% Wormhole Locations are NavPoints for those who make use of the wormhole
temp_points = allNavPoints;
for k = 1:numel(typeWormhole)
    temp_x = repmat(allWormholeLocations{k}(:,1),1,size(allNavPoints,1));
    temp_y = repmat(allWormholeLocations{k}(:,2),1,size(allNavPoints,1));
    temp_distx = allNavPoints(:,1)' - temp_x;
    temp_disty = allNavPoints(:,2)' - temp_y;
    temp_dist = sqrt((temp_distx).^2 + (temp_disty).^2);
    temp_points(sum(temp_dist<settings.navGraph.wormholeClearance)>0,:) = NaN;
end
allNavPoints = temp_points;

corner.navPointx = temp_points(:,1);
corner.navPointy = temp_points(:,2);

plotWalkingArea2

allNavPoints = allNavPoints(isnan(allNavPoints(:,1))==0,:);
allNavPoints = unique(allNavPoints,'rows','stable');

%% Deciding side of service locations
% So that user does not have to input service location but only the service
% edge which will already be part of the building plan
for i = 1:numel(typeService)
    for k = 1:numel(services(i).x1)
        temp = lineSegmentIntersect...
            (allObstacles,[repmat([services(i).location1x(k) services(i).location1y(k)],...
            size(allNavPoints,1),1) allNavPoints]);
        temp_intersection = temp.intAdjacencyMatrix;
        temp_visible = sum(temp_intersection)==0;
        if sum(temp_visible)==0
            services(i).location1x(k) = services(i).location2x(k);
            services(i).location1y(k) = services(i).location2y(k);
        end
    end
    services(i).locationx = services(i).location1x;
    services(i).locationy = services(i).location1y;
    allServiceLocations{i} = [services(i).locationx services(i).locationy];
end

for j = 1:numel(typeService)
    plot(services(j).locationx,services(j).locationy,'+','color',[1 0.4 0],'MarkerSize',3); 
end

%% Connecting All NavPoints
% Uni-directional edges
temp_x = repmat(allNavPoints(:,1),1,size(allNavPoints,1));
temp_y = repmat(allNavPoints(:,2),1,size(allNavPoints,1));
navGraph.distx = allNavPoints(:,1)' - temp_x;
navGraph.disty = allNavPoints(:,2)' - temp_y;
navGraph.dist = sqrt((navGraph.distx).^2 + ...
    (navGraph.disty).^2);
navGraph.angle = atan2d(navGraph.disty,-navGraph.distx); % negative to put 0 degrees on West, 90 on North
[~,navGraph.distAscending] = sort(navGraph.dist,2);
navGraph.visible = zeros(size(allNavPoints,1),size(allNavPoints,1));
navGraph.Lspace = zeros(size(allNavPoints,1),size(allNavPoints,1));
temp_available = ones(size(allNavPoints,1),size(allNavPoints,1));
for i = 1:size(allNavPoints,1)
    temp = lineSegmentIntersect...
        (allObstacles,[repmat(allNavPoints(i,:),size(allNavPoints,1),1) ...
        allNavPoints]);
    temp_intersection = temp.intAdjacencyMatrix;
    navGraph.visible(i,:) = sum(temp_intersection)==0;
    
    for j = 1:size(allNavPoints,1)
        if navGraph.dist(i,navGraph.distAscending(i,j))~=0 && ... %not the same point
                navGraph.visible(i,navGraph.distAscending(i,j))~=0 && ... %visible
                temp_available(i,navGraph.distAscending(i,j))==1 %not excluded by cone
            navGraph.Lspace(i,navGraph.distAscending(i,j)) = ...
                navGraph.dist(i,navGraph.distAscending(i,j));
            % Cone check
            if navGraph.angle(i,navGraph.distAscending(i,j)) > -180+settings.navGraph.halfConeAngle && ...
                    navGraph.angle(i,navGraph.distAscending(i,j)) < 180-settings.navGraph.halfConeAngle
                temp_available(i,:) = temp_available(i,:).*...
                    [navGraph.angle(i,:)<=navGraph.angle(i,navGraph.distAscending(i,j))-settings.navGraph.halfConeAngle | ...
                    navGraph.angle(i,:)>=navGraph.angle(i,navGraph.distAscending(i,j))+settings.navGraph.halfConeAngle];
            elseif navGraph.angle(i,navGraph.distAscending(i,j)) <= -180+settings.navGraph.halfConeAngle
                temp_available(i,:) = temp_available(i,:).*...
                    [navGraph.angle(i,:)>navGraph.angle(i,navGraph.distAscending(i,j))+settings.navGraph.halfConeAngle & ... % added +halfConeAngle(20/7/17)
                    navGraph.angle(i,:)<360+navGraph.angle(i,navGraph.distAscending(i,j))-settings.navGraph.halfConeAngle];
            elseif navGraph.angle(i,navGraph.distAscending(i,j)) >= 180-settings.navGraph.halfConeAngle
                temp_available(i,:) = temp_available(i,:).*...
                    [navGraph.angle(i,:)<navGraph.angle(i,navGraph.distAscending(i,j))-settings.navGraph.halfConeAngle & ... % added -halfConeAngle(20/7/17)
                    navGraph.angle(i,:)>-360+navGraph.angle(i,navGraph.distAscending(i,j))+settings.navGraph.halfConeAngle];
            end
            
        end
    end
end

%% Connecting All Nodes to NavPoints
% Uni-directional edges

allNodes = [vertcat(allLocations{:});vertcat(allServiceLocations{:});vertcat(allWormholeLocations{:})];
numNodes = [numLocations(:);numServices(:);numWormholes(:)*2];
typeNode = [typeLocation;typeService;typeWormhole];

% Get node type for allNodes
temp_a = cumsum(numNodes(numNodes>0));
temp_b = zeros(1,temp_a(end));
temp_b(temp_a-numNodes(numNodes>0)+1) = 1;
temp = typeNode(numNodes>0);
allNodeTypes = temp(cumsum(temp_b));

temp_x = repmat(allNodes(:,1),1,size(allNavPoints,1));
temp_y = repmat(allNodes(:,2),1,size(allNavPoints,1));
locationGraph.distx = allNavPoints(:,1)' - temp_x;
locationGraph.disty = allNavPoints(:,2)' - temp_y;
locationGraph.dist = sqrt((locationGraph.distx).^2 + ...
    (locationGraph.disty).^2);
locationGraph.angle = atan2d(locationGraph.disty,-locationGraph.distx); % negative to put 0 degrees on West, 90 on North
[~,locationGraph.distAscending] = sort(locationGraph.dist,2);
locationGraph.visible = zeros(size(allNodes,1),size(allNavPoints,1));
locationGraph.Lspace = zeros(size(allNodes,1),size(allNavPoints,1));
temp_available = ones(size(allNodes,1),size(allNavPoints,1));
for i = 1:size(allNodes,1)
    temp = lineSegmentIntersect...
        (allObstacles,[repmat(allNodes(i,:),size(allNavPoints,1),1) ...
        allNavPoints]);
    temp_intersection = temp.intAdjacencyMatrix;
    locationGraph.visible(i,:) = sum(temp_intersection)==0;
    
    for j = 1:size(allNavPoints,1)
        if locationGraph.dist(i,locationGraph.distAscending(i,j))~=0 && ... %not the same point
                locationGraph.visible(i,locationGraph.distAscending(i,j))~=0 && ... %visible
                temp_available(i,locationGraph.distAscending(i,j))==1 %not excluded by cone
            locationGraph.Lspace(i,locationGraph.distAscending(i,j)) = ...
                locationGraph.dist(i,locationGraph.distAscending(i,j));
            % Cone check
            if locationGraph.angle(i,locationGraph.distAscending(i,j)) > -180+settings.navGraph.halfConeAngle && ...
                    locationGraph.angle(i,locationGraph.distAscending(i,j)) < 180-settings.navGraph.halfConeAngle
                temp_available(i,:) = temp_available(i,:).*...
                    [locationGraph.angle(i,:)<=locationGraph.angle(i,locationGraph.distAscending(i,j))-settings.navGraph.halfConeAngle | ...
                    locationGraph.angle(i,:)>=locationGraph.angle(i,locationGraph.distAscending(i,j))+settings.navGraph.halfConeAngle];
            elseif locationGraph.angle(i,locationGraph.distAscending(i,j)) <= -180+settings.navGraph.halfConeAngle
                temp_available(i,:) = temp_available(i,:).*...
                    [locationGraph.angle(i,:)>locationGraph.angle(i,locationGraph.distAscending(i,j))+settings.navGraph.halfConeAngle & ... % added +halfConeAngle(20/7/17)
                    locationGraph.angle(i,:)<360+locationGraph.angle(i,locationGraph.distAscending(i,j))-settings.navGraph.halfConeAngle];
            elseif locationGraph.angle(i,locationGraph.distAscending(i,j)) >= 180-settings.navGraph.halfConeAngle
                temp_available(i,:) = temp_available(i,:).*...
                    [locationGraph.angle(i,:)<locationGraph.angle(i,locationGraph.distAscending(i,j))-settings.navGraph.halfConeAngle & ... % added -halfConeAngle(20/7/17)
                    locationGraph.angle(i,:)>-360+locationGraph.angle(i,locationGraph.distAscending(i,j))+settings.navGraph.halfConeAngle];
            end
            
        end
    end
end
%}

%% Connecting Wormhole Locations
allWormholeNodes = vertcat(allWormholeLocations{:});
numWormholeNodes = numWormholes*2;
temp_a = cumsum(numWormholeNodes(numWormholeNodes>0));
temp_b = zeros(1,temp_a(end));
temp_b(temp_a-numWormholeNodes(numWormholeNodes>0)+1) = 1;
temp = typeWormhole(numWormholeNodes>0);
allWormholeTypes = temp(cumsum(temp_b));

wormholeGraph.Lspace = zeros(sum(numWormholes));
for i = 1:numel(typeWormhole)
    temp = find(allWormholeTypes==typeWormhole(i));
    for j = 1:numWormholes(i)
        wormholeGraph.Lspace(temp(2*j-1),temp(2*j)) = settings.serviceTime.(wormholes(i).type);
        wormholeGraph.Lspace(temp(2*j),temp(2*j-1)) = settings.serviceTime.(wormholes(i).type);
    end    
end
%}

%% Plot

% NavPoints-NavPoints
for i = 1:size(allNavPoints,1)
    for j = 1:size(allNavPoints,1)
        if navGraph.Lspace(i,j)~=0
            plot([allNavPoints(i,1) allNavPoints(j,1)],...
                [allNavPoints(i,2) allNavPoints(j,2)],'-','color',[0 1 0])
        end
    end
end
text(allNavPoints(:,1),allNavPoints(:,2),(num2str([1:size(allNavPoints,1)]')),...
    'Color',[0 0.7 0],'FontSize',5);

% Nodes-NavPoints
for i = 1:size(allNodes,1)
    for j = 1:size(allNavPoints,1)
        if locationGraph.Lspace(i,j)~=0
            plot([allNodes(i,1) allNavPoints(j,1)],...
                [allNodes(i,2) allNavPoints(j,2)],'-','color',[0 1 1])
        end
    end
end
text(allNodes(:,1),allNodes(:,2),(num2str([1:size(allNodes,1)]')),...
    'Color','blue','FontSize',5);

% Wormhole-Wormhole
for i = 1:size(allWormholeNodes,1)
    for j = 1:size(allWormholeNodes,1)
        if wormholeGraph.Lspace(i,j)~=0
            plot([allWormholeNodes(i,1) allWormholeNodes(j,1)],...
                [allWormholeNodes(i,2) allWormholeNodes(j,2)],'-','color',[0.5 0 0])
        end
    end
end
text(allWormholeNodes(:,1)-0.25,allWormholeNodes(:,2)-0.25,...
    (num2str([1:size(allWormholeNodes,1)]')),...
    'Color',[0.8 0 0],'FontSize',5);

%% Clear Temp Vars
clear temp*
