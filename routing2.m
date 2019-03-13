%% Routing
% Routing according to shortest available path, includes consideration of
% queues at wormholes and services
%
% Uses Dijkstra's Algorithm: MATLAB code from File Exchange:
% https://nl.mathworks.com/matlabcentral/fileexchange/20025-dijkstra-s-minimum-cost-path-algorithm
%
% Three graphs are present: navigation, location-to-navigation, wormhole
% connections. NavPoints nearest to each location is found; Dijkstra's is
% used to get the shortest paths between NavPoints. Base paths are formed
% by location_or-{navPoints}-location_dest all as {Points}. Full path are formed by the
% addition of wormhole connections and are of the form
% location_or-{wormholeNodes}-location_dest all as{Nodes}. The complete path can be
% derived by using information from the base paths in the full paths during
% runtime.

%% Version Change
% 1. Readied for integration

%% Notes
% 1. Takes "createNetwork3_2" and onwards
% 2. Uses "dijkstra.m"
% 3. ONLY ALLOWS ONE WORMHOLE TRANSITION IN A PATH

%% Base Paths
[navGraph.cost,navGraph.path] = dijkstra(navGraph.Lspace,navGraph.Lspace);
% Find connected navPoints to Nodes
nodeNavNbrs = cell(size(allNodes,1),1);
nodeNavNbrsCost = cell(size(allNodes,1),1);
navNodeNbrs = cell(size(allNavPoints,1),1);
navNodeNbrsCost = cell(size(allNavPoints,1),1);
for i = 1:size(allNodes,1)
    nodeNavNbrs{i} = find(locationGraph.Lspace(i,:)>0);
    nodeNavNbrsCost{i} = locationGraph.Lspace(i,nodeNavNbrs{i});
    
end

routing.baseCost = ones(size(allNodes,1))*Inf;
routing.basePath = cell(size(allNodes,1)); %no wormhole travel
for i = 1:size(allNodes,1)
    for j = 1:size(allNodes,1)
        if i==j; routing.baseCost(i,j)=0;continue; end
        
        for k = 1:numel(nodeNavNbrs{i})
            for m = 1:numel(nodeNavNbrs{j})
                temp_navO = nodeNavNbrs{i}(k);
                temp_navD = nodeNavNbrs{j}(m);
                temp_cost = nodeNavNbrsCost{i}(k) + ...
                    navGraph.cost(temp_navO,temp_navD) + ...
                    nodeNavNbrsCost{j}(m);
                if temp_cost<routing.baseCost(i,j)
                    routing.baseCost(i,j) = temp_cost;
                    routing.basePath{i,j} = ...
                        {[node2point(i) ...
                        navPoint2point(navGraph.path{temp_navO,temp_navD})'...
                        node2point(j)]};
                end
            end
        end
        
    end
end

%% Full Paths
% Find all wormhole nodes connected to a location; find the wormhole nodes
% corresponding to these, i.e., use the wormhole graph to find all wormhole
% nodes connected to the nodes that are connected to the location node; for
% all the location nodes unconnected to a location
routing.fullCost = num2cell(routing.baseCost);
routing.fullPath = cell(size(allNodes,1));
for i = 1:size(allNodes,1)
    for j = 1:size(allNodes,1)
        routing.fullPath{i,j} = {[i,j]};
    end
end
numFullPaths = ones(size(allNodes,1));


temp_unconnected = cell(size(allNodes,1),1);
temp_wormholeConnections = false(size(allNodes,1),size(allWormholeNodes,1));
temp_correspondingWormholes = zeros(size(allNodes,1),size(allWormholeNodes,1));
for i = 1:size(allNodes,1)
    temp_unconnected{i} = find(routing.baseCost(i,:)==Inf);
    temp_wormholeConnections(i,:) = ...
        routing.baseCost(i,ismember(allNodeTypes,typeWormhole))<Inf;
    temp_correspondingWormholes(i,:) = ...
        sum(wormholeGraph.Lspace(temp_wormholeConnections(i,:)',:));
end

for i = 1:size(allNodes,1) % origin
    for j = 1:numel(temp_unconnected{i}) % destination
        temp_d = temp_unconnected{i}(j);
        if i==temp_d; continue; end
        temp = temp_correspondingWormholes(temp_d,:).*...
            temp_wormholeConnections(i,:); % can they be connected through any wormhole?
        if sum(temp) > 0
            temp_fromWormholes = find(temp); % all 'from' wormholes that have a corresponding wormhole the dest is connected (which wormhole the dest is connected to is unknown)
            routing.fullCost{i,temp_d} = []; % to be able to assign more than one cost as a nested cell array
            for k = 1:numel(temp_fromWormholes)
                temp_toWormholes = ...
                    wormholeGraph.Lspace(temp_fromWormholes(k),:); % all (no indices) the corresponding 'to' wormholes (think of a staircase)
                temp_toWormholes = ...
                    temp_toWormholes.*temp_wormholeConnections(temp_d,:); % keep only (no indices) those that are connected to the destination
                temp_toWormholes = find(temp_toWormholes); % actually 'find' indices
                temp_fromNode = wormhole2node(temp_fromWormholes(k)); % find indice within all nodes
                for m = 1:numel(temp_toWormholes)
                    temp_toNode = wormhole2node(temp_toWormholes(m)); % find indice within all nodes
                    routing.fullCost{i,temp_d}(numFullPaths(i,temp_d)) = ...
                        routing.baseCost(i,temp_fromNode) + ...
                        wormholeGraph.cost(temp_fromWormholes(k),...
                        temp_toWormholes(m)) + ...
                        routing.baseCost(temp_toNode,temp_d);
                        
                    routing.fullPath{i,temp_d}{numFullPaths(i,temp_d)} = ...
                        [i temp_fromNode temp_toNode temp_d];
                    numFullPaths(i,temp_d) = numFullPaths(i,temp_d) + 1;
                end
            end
        end
        numFullPaths(i,temp_d) = numFullPaths(i,temp_d) - 1;
    end
end

%% Points to Nodes
navToNodeCost = zeros(size(allNavPoints,1),size(allNodes,1));
navToNodeNbr = zeros(size(allNavPoints,1),size(allNodes,1));
for i = 1:size(allNavPoints,1)
    navNodeNbrs{i} = find(locationGraph.Lspace(:,i)>0);
    navNodeNbrsCost{i} = locationGraph.Lspace(navNodeNbrs{i},i);
    for k = 1:size(allNodes,1)
        temp = Inf;
        for j = 1:numel(navNodeNbrs{i})
            temp = min(temp,navNodeNbrsCost{i}(j)+...
                min(routing.fullCost{navNodeNbrs{i}(j),k}));
        end
        if ~isempty(temp) && ~isempty(j)
            navToNodeCost(i,k) = temp;
            navToNodeNbr(i,k) = navNodeNbrs{i}(j);
        end
    end
end