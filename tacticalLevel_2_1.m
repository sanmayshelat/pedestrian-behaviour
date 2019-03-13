%% Citation Information
%Please cite: Shelat, S. (2017). Developing an Integrated Pedestrian Behaviour Model for Office Buildings. (MSc), Delft University of Technology.   

%% Tactical Level
% Uses:
% 1. checkAvailability
% 2. mcGenerator_scheduler
% 3. mcSimulator_scheduler

%% Settings
if settings.runtime.setRandom == 0
        rng(1);
else
    rng(day);
end

timeStep = settings.tactical.timeStep;
plannedTimeStep = settings.tactical.plannedTimeStep;
meetingGap = settings.tactical.plannedTimeStep;
agentAvail = ones(numAgents,24*60); % will include preferred meeting gaps
agentLocation = zeros(numAgents,24*60);
agentActivity = zeros(numAgents,24*60);
%}

%% Event Scheduler
% Only activity 'Meeting' considered
tempActivityType = find(strcmp('Meeting',allActivities));
tempLocationType = allTypes(activityLocationType{tempActivityType});
temp = ismember({abstractSpace(:).type},tempLocationType);

meetingRooms = find(temp);
numMeetingRooms = numel(meetingRooms);
meetingRoomAvail = ones(numMeetingRooms,24*60);
meetingRoomCapacity = sum(abstractSpace(meetingRooms).node);

teamMeetingSchedule = zeros(numTeamActivities,2);

tempAgentAvail = agentAvail;
tempArrive = startHabitual(:,1).Hour*60+startHabitual(:,1).Minute;
tempDep = startHabitual(:,2).Hour*60+startHabitual(:,2).Minute;
for  i = 1:numAgents
    tempAgentAvail(i,1:tempArrive(i)-1)=0;
    tempAgentAvail(i,tempDep(i)+1:end)=0;
end
temp_allowedTimes = zeros(1,24*60);
temp_allowedTimes(mod(1:24*60,plannedTimeStep)==0) = 1;
%}

%% Team Activities

for i = 1:numTeamActivities
    % Check if meeting occurs (meeting doesn't occur if room cannot be found)
    if rand>teamMeetingProbability(i)
        continue;
    end
    % Prefer room with capacity closer to number of team meeting members
    [temp_meetingCapacity,temp_index] = sort(meetingRoomCapacity);
    for j = 1:numMeetingRooms
        if temp_meetingCapacity<numel(teamMeeting{i})
            continue;
        end
        
        % Combined team members and room availability check
        tempAgentRoomMatch = ...
            checkAvailability([tempAgentAvail(teamMeeting{i},:);...
            meetingRoomAvail(temp_index(j),:)]);
        temp_possibilities = tempAgentRoomMatch>teamMeetingDuration(i);
        temp_possibilities = temp_possibilities.*temp_allowedTimes;
        temp_possibilities = find(temp_possibilities);
        
        % Check if a meeting is possible in this room
        if isempty(temp_possibilities)
            continue;
        end
        teamMeetingSchedule(i,1) = temp_possibilities(...
            randi([1 numel(temp_possibilities)]));
        teamMeetingSchedule(i,2) = teamMeetingSchedule(i,1) + ...
            teamMeetingDuration(i);
        
        % Agents and room no longer available for this duration
        tempAgentAvail(teamMeeting{i},...
            teamMeetingSchedule(i,1)-meetingGap:...
            teamMeetingSchedule(i,2)+meetingGap) = 0;
        agentAvail(teamMeeting{i},...
            teamMeetingSchedule(i,1)-meetingGap:...
            teamMeetingSchedule(i,2)+meetingGap) = 0;
        meetingRoomAvail(temp_index(j),...
            teamMeetingSchedule(i,1):teamMeetingSchedule(i,2)) = 0;
        
        % Assign locations to agents
        tempNodes = find(abstractSpace(meetingRooms(temp_index(j))).node);
        tempNodes = tempNodes(randperm(numel(tempNodes),...
            numel(teamMeeting{i}))); % Randomize where people sit
        for k = 1:numel(teamMeeting{i})
            agentLocation(teamMeeting{i}(k),...
                teamMeetingSchedule(i,1):...
                teamMeetingSchedule(i,2)) = tempNodes(k);
            agentActivity(teamMeeting{i}(k),...
                teamMeetingSchedule(i,1):...
                teamMeetingSchedule(i,2)) = ...
                find(strcmp('Meeting',allActivities));
        end
        
        % Room found so done
        break;
    end
end
%}

%% Habitual Activities

% Arrival & Departure
agentArrival = zeros(numAgents,1);
agentDep = zeros(numAgents,1);
tempNodes = ...
    find(allNodeTypes==find(strcmp('locationOutside',allTypes)));
for i = 1:numAgents
    % Keep minimum gap with first/last meeting
    temp = find(agentLocation(i,:));
    if isempty(temp)==0
        temp1 = min(startHIndex(i,1)+startHVarIndex(i,1),...
            temp(1)-meetingGap);
        tempArrival = randi([startHIndex(i,1)-startHVarIndex(i,1) temp1]);
        temp1 = max(startHIndex(i,2)-startHVarIndex(i,2),...
            temp(end)+meetingGap);
        tempDep = randi([temp1 startHIndex(i,2)+startHVarIndex(i,2)]);
    else
        tempArrival = randi([startHIndex(i,1)-startHVarIndex(i,1) ...
            startHIndex(i,1)+startHVarIndex(i,1)]);
        tempDep = randi([startHIndex(i,2)-startHVarIndex(i,2) ...
            startHIndex(i,2)+startHVarIndex(i,2)]);
    end
    agentLocation(i,1:tempArrival) = ...
                tempNodes(randperm(numel(tempNodes),1));
    agentLocation(i,tempDep:end) = ...
                tempNodes(randperm(numel(tempNodes),1));
        
    agentActivity(i,1:tempArrival) = ...
        find(strcmp('Arrival',allActivities));
    agentActivity(i,tempDep:end) = ...
        find(strcmp('Departure',allActivities));
    
    agentAvail(i,1:tempArrival) = 0;
    agentAvail(i,tempDep:end) = 0;
    
    agentArrival(i) = tempArrival;
    agentDep(i) = tempDep;
end

% Others
for i = 3:activity.habitual.num
    tempNodes = ...
        activityLocation{strcmp(activity.habitual.type(i),allActivities)};
    for j = 1:numAgents
        % Check availability of agent
        tempAvailability = checkAvailability(agentAvail(j,:));
        % Check availability within preferred time
        tempAvailTimeWindow = find(tempAvailability(...
            startHIndex(j,i)-startHVarIndex(j,i):...
            startHIndex(j,i)+startHVarIndex(j,i))>=durationHIndex(j,i));
        tempAvailpreTimeWindow = find(tempAvailability(1:...
            startHIndex(j,i)-startHVarIndex(j,i)-1)>=durationHIndex(j,i));
        tempAvailpostTimeWindow = find(tempAvailability(...
            startHIndex(j,i)+startHVarIndex(j,i)+1:...
            end)>=durationHIndex(j,i));
        if isempty(tempAvailTimeWindow)==0
            tempH = startHIndex(j,i)-startHVarIndex(j,i) + ...
                tempAvailTimeWindow(...
                randperm(numel(tempAvailTimeWindow),1)) - 1;
            tempH = tempH:tempH+durationHIndex(j,i);
            agentAvail(j,tempH) = 0;
            agentActivity(j,tempH) = ...
                find(strcmp(activity.habitual.type(3),allActivities));
            agentLocation(j,tempH) = ...
                tempNodes(randperm(numel(tempNodes),1));
        % Assuming activities are conducted even if outside time window
        elseif isempty(tempAvailpostTimeWindow)==0
            % As early as possible post time window
            tempH = startHIndex(j,i)+startHVarIndex(j,i)+1 + ...
                tempAvailpostTimeWindow(1) - 1;
            tempH = tempH:tempH+durationHIndex(j,i);
            agentAvail(j,tempH) = 0;
            agentActivity(j,tempH) = ...
                find(strcmp(activity.habitual.type(3),allActivities));
            agentLocation(j,tempH) = ...
                tempNodes(randperm(numel(tempNodes),1));
        elseif isempty(tempAvailpreTimeWindow)==0
            % As late as possible pre time window
            tempH = startHIndex(j,i)-startHVarIndex(j,i)-1 - ...
                tempAvailpreTimeWindow(1) + 1;
            tempH = tempH:tempH+durationHIndex(j,i);
            agentAvail(j,tempH) = 0;
            agentActivity(j,tempH) = ...
                find(strcmp(activity.habitual.type(3),allActivities));
            agentLocation(j,tempH) = ...
                tempNodes(randperm(numel(tempNodes),1));
        % Not scheduled if no time available    
        end
    end
end
%}

%% Event Schedule
agentEventSchedule = cell(numAgents,1);
agentEventLocation = cell(numAgents,1);
agentEventActivity = cell(numAgents,1);
for i = 1:numAgents
    agentEventSchedule{i} = find(diff(agentLocation(i,:)));
    agentEventSchedule{i}(2:2:end) = agentEventSchedule{i}(2:2:end)+1;
    agentEventLocation{i} = agentLocation(i,agentEventSchedule{i});
    agentEventActivity{i} = agentActivity(i,agentEventSchedule{i});
end

%% Assign flex spaces for the day
% Randomly assigns a chair in the flex space to agents with
% first-come-first-serve
tempAgentBaseSpaces = find(baseFlexSpaces);
tempSpaces = baseFlexSpaces(baseFlexSpaces~=0);
tempUniqueSpace = unique(tempSpaces,'stable');
for i = 1:numel(tempUniqueSpace)
    tempNodes = find(abstractSpace...
        (typeAbstractSpaces==tempUniqueSpace(i)).node);
    tempFlexAvail = ones(numel(tempNodes),24*60);
    tempAgentsInThisSpace = find(baseFlexSpaces==tempUniqueSpace(i));
    [~,tempArrivalOrder] = sort(agentArrival(tempAgentsInThisSpace));
    for j = 1:numel(tempArrivalOrder)
        tempAgent = tempAgentsInThisSpace(tempArrivalOrder(j));
        tempAvail = sum(tempFlexAvail(:,...
            agentArrival(tempAgent):agentDep(tempAgent)),2)==...
            agentDep(tempAgent)-agentArrival(tempAgent)+1;
        tempAvail = find(tempAvail);
        baseLocations(tempAgent) = tempAvail(randperm(numel(tempAvail),1));
        tempFlexAvail(baseLocations(tempAgent),...
            agentArrival(tempAgent):agentDep(tempAgent)) = 0;
    end
    baseLocations(baseFlexSpaces==tempUniqueSpace(i)) = ...
        tempNodes(randperm(sum(baseFlexSpaces==tempUniqueSpace(i))));
end
%}


%% Settings
tol = settings.tactical.mcGenTol;

%% Fine-grain Schedule
% Convert minutes schedule to seconds schedule to consider movement times
agentFineActivity = zeros(numAgents,24*60*60);
agentFineLocation = zeros(numAgents,24*60*60);
for i = 1:numAgents
    temp = repmat(agentActivity(i,:),60,1);
    agentFineActivity(i,:) = temp(:)';
    temp = repmat(agentLocation(i,:),60,1);
    agentFineLocation(i,:) = temp(:)';
end

%% Fine-grain Event Schedule
agentFineEventSchedule = cell(numAgents,1);
agentFineEventLocation = cell(numAgents,1);
agentFineEventActivity = cell(numAgents,1);
for i = 1:numAgents
    agentFineEventSchedule{i} = find(diff(agentFineLocation(i,:)));
    agentFineEventSchedule{i}(2:2:end) = ...
        agentFineEventSchedule{i}(2:2:end)+1;
    agentFineEventLocation{i} = ...
        agentFineLocation(i,agentFineEventSchedule{i});
    agentFineEventActivity{i} = ...
        agentFineActivity(i,agentFineEventSchedule{i});
end

%% MC Scheduler
P = zeros(numMCActivities,numMCActivities,3);
mcActivitySchedule = cell(numAgents,1);
mcFineActivitySchedule = cell(numAgents,1);

%% Specific Rule Assumed
% After getting coffee 1/3rd probability of taking a break
% Ignored for receptionist because doesn't take a break
tempNumVars = numMCActivities^2-numMCActivities;
A = eye(tempNumVars);
b = ones(tempNumVars,1);
lb = zeros(tempNumVars,1);
b(9) = 0.33*(1-0.5);
lb(9) = 0.32*(1-0.5);
% b([8 6 12]) = 0.01;

%% Spontaneous Activities Rules
% Change tempNotInduced

for i = 1:activity.spontaneous.num
    temp = activity.spontaneous.index(i);
    tempNotInduced = 2:4;
    for j = 1:numel(tempNotInduced)
        b(sub2ind([numMCActivities numMCActivities],temp,...
            tempNotInduced(j)) - ...
            tempNotInduced(j)*(tempNotInduced(j)<temp)) = 0.0001;
    end
end

%% Merging Schedules
for i = 1:numAgents
    %% Generate Transition Matrices
    %--- Ignore constraints for receptionist ---
    if i ==18 || i==1
        [P(:,:,i),tempf] = mcGenerator_scheduler...
            (ltpAgents(i,:),stAgents(i,:),timeStep,'LsqLin');
        
    else
    
    [P(:,:,i),tempf] = mcGenerator_scheduler...
        (ltpAgents(i,:),stAgents(i,:),timeStep,'LsqLin',...
        {A b [] [] lb []});
    end
    if tempf>tol
        error('LTP and ST combination invalid for Agent %d',i)
    end
    
    %% Simulate MC activities for all remaining available time contiguously
    temp = agentEventSchedule{i};
    tempTotalTime = sum(temp(2:2:end)-temp(1:2:end));%sum(agentAvail(i,:)~=0);
    [mcActivitySchedule{i}] = ...
        mcSimulator_scheduler(P(:,:,i),2*tempTotalTime,ltpAgents(i,:));
    % Simulate extra activities as duration of service activities will be
    % automatically reduced to the service time and therefore the ends
    % cannot be met exactly
    temp = repmat(mcActivitySchedule{i},1,60)';
    mcFineActivitySchedule{i} = temp(:);
%     tempMCActivityChange = find(diff(mcFineActivitySchedule{i})); % change occurs at the end of the time step
    % Smooth agenda (remove activities without min. duration)?
    
    %% Arrange Schedule
    % Insert MC activities and movement times between event times
    tempSpeed = settings.indiWalk.desiredSpeedHigh;
    for j = 1:2:numel(agentFineEventSchedule{i}) % Only upto activity before dep.
        t = agentFineEventSchedule{i}(j); % Set time to arrival time
        tempCurrentPos = agentFineLocation(i,t); % Current Location
        tempNextTime = agentFineEventSchedule{i}(j+1); % Next planned/habitual time
        tempNextPos = agentFineEventLocation{i}(j+1); % Next planned/habitual location
        tempNextCost = min(routing.fullCost{tempCurrentPos,tempNextPos})/...
            tempSpeed; % Time to next p/h location
%         tempNextCost = tempNextCost/(60*timeStep); % movement time is in [s]
        tempTimeC2E = tempNextCost;
        
        temp_index = 1;
        tempMCActivityChange = find(diff(mcFineActivitySchedule{i})); % change occurs at the end of the time step
        tempMCActivityChange = ...
            [tempMCActivityChange; numel(mcFineActivitySchedule{i})]; % add last activity OR BETTER just generate longer mc schedule by doubling totaltime
        tempTotalMCDuration = 0;
        tempTotalMovement = 0;
        tempActivity = [];
        tempLocation = [];
        % Check if MC activity should be inserted
        while tempNextTime - t > tempNextCost
            tempNextMCTime = tempMCActivityChange(temp_index);
            if temp_index==1
                tempMCActivityDuration = tempMCActivityChange(temp_index);
            else
                tempMCActivityDuration = ...
                    tempMCActivityChange(temp_index) - ...
                    tempMCActivityChange(temp_index-1);
            end
            tempMCActivity = mcFineActivitySchedule{i}(tempNextMCTime);
            if tempMCActivity == 1 % base, continuous activity
                tempLocationOptions = baseLocations(i);
            else % other activities
                tempLocationOptions = activityLocation{tempMCActivity};
            end
            tempLocationCosts = zeros(numel(tempLocationOptions),1);
            for k = 1:numel(tempLocationOptions)
                tempLocationCosts(k) = ...
                    min(routing.fullCost{tempCurrentPos,...
                    tempLocationOptions(k)})/tempSpeed;
            end
            [tempTimeC2MC,tempNearestPos] = min(tempLocationCosts);
%             tempTimeC2MC = tempTimeC2MC/(60*timeStep);
            tempCurrentPos = tempLocationOptions(tempNearestPos);
            % If location is service then fixed activity duration
            if regexp(allTypes{allNodeTypes(tempCurrentPos)},'service*\w')
                tempMCActivityDuration = ...
                    settings.serviceTime.(allTypes{...
                    allNodeTypes(tempCurrentPos)});
            end
            tempTimeMC2E = ...
                min(routing.fullCost{tempCurrentPos,tempNextPos})/tempSpeed;
%             tempTimeMC2E = tempTimeMC2E/(60*timeStep);
            tempNextCost = tempNextCost - floor(tempTimeC2E) + ...
                floor(tempTimeC2MC) + tempMCActivityDuration + ...
                floor(tempTimeMC2E);
            
            tempTotalMovement = tempTotalMovement + ...
                tempTimeC2MC + tempTimeMC2E - tempTimeC2E;
            
            temp_index = temp_index + 1;
            tempTimeC2E = tempTimeMC2E;
            tempActivity = [tempActivity ...
                ones(1,floor(tempTimeC2MC)).*(-1) ...
                ones(1,tempMCActivityDuration).*(tempMCActivity)];
            tempLocation = [tempLocation ...
                ones(1,floor(tempTimeC2MC)).*(-1) ...
                ones(1,tempMCActivityDuration).*(tempCurrentPos)];
            tempTotalMCDuration = tempTotalMCDuration + ...
                tempMCActivityDuration;
            
        end
        
        % Add connection to end Event
        if temp_index==1 % if no MC activities in the middle
            tempActivity = ...
                [tempActivity ones(1,floor(tempTimeC2E)).*(-1)];
            tempLocation = ...
                [tempLocation ones(1,floor(tempTimeC2E)).*(-1)];
        else % if there are MC activities in the middle
            % Reduce MC Activity Duration to be in time for next Event
            % numel(tempActivity) used to calculate time to ensure that
            % agents never reach an event early as then they will have
            % nowhere to go
            tempMCActivityReduction = ...
                numel(tempActivity)+floor(tempTimeC2E)-(tempNextTime-t);
            if tempMCActivityReduction >= tempMCActivityDuration
                tempMCActivityReduction = tempMCActivityDuration - 60;
            end
            tempActivity(end-tempMCActivityReduction:end) = [];
            tempActivity = ...
                [tempActivity ones(1,floor(tempTimeC2E)).*(-1)];
            tempLocation(end-tempMCActivityReduction:end) = [];
            tempLocation = ...
                [tempLocation ones(1,floor(tempTimeC2E)).*(-1)];
            % Delete mcFineActivityScheduler upto what is already performed
            tempTotalMCDuration = tempTotalMCDuration - ...
                tempMCActivityReduction;
            mcFineActivitySchedule{i}(1:tempTotalMCDuration) = [];
        end
        
        % Check if any activities are getting overwritten
        if any(agentFineActivity(i,t+1:tempNextTime-1)) || ...
                any(agentFineLocation(i,t+1:tempNextTime-1))
            keyboard;
        end
        agentFineActivity(i,t+1:t+numel(tempActivity)) = tempActivity;
        agentFineLocation(i,t+1:t+numel(tempActivity)) = tempLocation;
    end
    
end

%% Final Schedule
agentFinalSchedule = cell(numAgents,1);
agentFinalLocation = cell(numAgents,1);
agentFinalActivity = cell(numAgents,1);
movementSchedule = []; % Input for Operational Level
for i = 1:numAgents
    agentFinalSchedule{i} = find(diff(agentFineLocation(i,:)));
    agentFinalSchedule{i}(2:2:end) = agentFinalSchedule{i}(2:2:end)+1;
    agentFinalLocation{i} = agentFineLocation(i,agentFinalSchedule{i});
    agentFinalActivity{i} = agentFineActivity(i,agentFinalSchedule{i});
    movementSchedule = [movementSchedule; ...
        [ones(numel(agentFinalSchedule{i})/2,1).*i,...
        agentFinalSchedule{i}(1:2:end)',...
        agentFinalLocation{i}(1:2:end)',...
        agentFinalLocation{i}(2:2:end)',...
        agentFinalActivity{i}(2:2:end)',...
        agentFinalSchedule{i}(2:2:end)']];
end

%% Plot
%{
close all
i = 12;
plot(seconds(1:24*60*60),agentFineActivity(i,:),...
    'DurationTickFormat','hh:mm:ss','Color','red');
xlim([seconds(agentArrival(i)*60-600) ...
    seconds(agentDep(i)*60+600)]);
ylabel('Activity Type',...
    'FontName','Times New Roman','FontSize',9);
xlabel('Time',...
    'FontName','Times New Roman','FontSize',9);
ylim([-1 numActivities])
yticks([-1 1:numActivities])
yticklabels({'Moving' allActivities{:}})
set(gca,'FontName','Times New Roman','FontSize',9);
set(gcf, 'Color', 'w');
%}
%% Clear temp
clear temp*