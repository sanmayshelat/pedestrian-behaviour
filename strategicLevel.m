%% Strategic Level
% Uses:
% 1. GetExcelRange

%% Settings
if settings.runtime.setRandom == 0
    rng(1);
else
    rng('shuffle');
end

%% Organization

temp = xlsread(files.st.input,'OrganizationUnit','B1:B4');
numAgents = temp(2);
numRoles = temp(1);
numTeams = temp(3);
numTeamActivities = temp(4);

[~,roles,~] = xlsread(files.st.input,'Roles',...
    GetExcelRange(2,numRoles+1,2,2));
[~,roleAgents,~] = xlsread(files.st.input,'MCInput',...
    GetExcelRange(2,numAgents+1,2,2));
tempRoles = zeros(numAgents,1);
for i = 1:numRoles
    temp = strcmp(roleAgents,roles{i});
    tempRoles(temp) = i;
end
roleAgents = tempRoles;
%}

%% Teams

teamComposition = xlsread(files.st.input,'Teams',...
    GetExcelRange(2,numTeams+1,2,numRoles+1)); % r: id; c: role
team = cell(numTeams,1);
temp_available = ones(numAgents,1);
% Even distribution of agents over teams according to composition
for i = 1:numTeams
    for j = 1:numRoles
        if teamComposition(i,j)==0
            continue;
        end
        temp = find(temp_available.*(roleAgents==j),teamComposition(i,j));
        if numel(temp)<teamComposition(i,j)
            temp_available(roleAgents==j) = 1; % reset availability
            temp = find(temp_available.*(roleAgents==j),teamComposition(i,j));
        end
        team{i} = [team{i} temp'];
        temp_available(temp) = 0;
    end
end

teamActivities = xlsread(files.st.input,'TeamActivities',...
    GetExcelRange(2,numTeamActivities+1,2,numRoles+4)); % r: id; c: role
teamMeetingTeamId = teamActivities(:,1);
teamMeetingFrequency = teamActivities(:,2);
teamMeetingProbability = teamMeetingFrequency/5; % 5 working days in a week
teamMeetingDuration = teamActivities(:,3);
teamMeetingComposition = teamActivities(:,4:end);
teamMeeting = cell(numTeamActivities,1);
% Assign team members to team activities
% Even distribution of agents over teams according to composition
for i = 1:numTeamActivities
    if sum(teamMeetingComposition(i,:))==...
            sum(teamComposition(teamMeetingTeamId(i),:))
        teamMeeting{i} = team{teamMeetingTeamId(i)};
        continue;
    end
    for j = 1:numRoles
        temp1 = roleAgents(team{teamMeetingTeamId(i)})==j;
        if teamMeetingComposition(i,j)==...
                teamComposition(teamMeetingTeamId(i),j)
            teamMeeting{i} = [teamMeeting{i} ...
                team{teamMeetingTeamId(i)}(temp1)];
            continue;
        end
        temp2 = randperm(sum(temp1),teamMeetingComposition(i,j));
        temp1 = find(temp1);
        teamMeeting{i} = [teamMeeting{i} ...
            team{teamMeetingTeamId(i)}(temp1(temp2))];
    end
end
%}

%% Activities
% 5 categories of activities

temp = xlsread(files.st.input,'Activities','H1:H5');
activity.continuous.num = temp(1);
activity.intermediate.num = temp(2);
activity.spontaneous.num = temp(3);
activity.habitual.num = temp(4);
activity.planned.num = temp(5);
numActivities = sum(temp);

activity.continuous.index = 1:temp(1);
activity.intermediate.index = temp(1)+(1:temp(2));
activity.spontaneous.index = sum(temp(1:2))+(1:temp(3));
activity.habitual.index = sum(temp(1:3))+(1:temp(4));
activity.planned.index = sum(temp(1:4))+(1:temp(5));

[~,temp,~] = xlsread(files.st.input,'Activities',...
    GetExcelRange(2,numActivities+1,1,5));
activity.continuous.type = temp(1:activity.continuous.num,1);
activity.intermediate.type = temp(1:activity.intermediate.num,2);
activity.spontaneous.type = temp(1:activity.spontaneous.num,3);
activity.habitual.type = temp(1:activity.habitual.num,4);
activity.planned.type = temp(1:activity.planned.num,5);
allActivities = [activity.continuous.type;activity.intermediate.type;...
    activity.spontaneous.type;activity.habitual.type;activity.planned.type];
%}

%% Base Locations
baseLocations = xlsread(files.st.input,'BaseLocations',...
    GetExcelRange(2,numAgents+1,3,3));
baseFlexSpaces = xlsread(files.st.input,'BaseLocations',...
    GetExcelRange(2,numAgents+1,4,4));
%}

%% Activity Locations
activityLocationType = cell(numActivities,1);
activityLocation = cell(numActivities,1);
%-1 to remove base location continuous activity
for i = 1:numActivities-1
    activityLocationType{i+1} = xlsread(files.st.input,...
        'ActivityLocations',GetExcelRange(i,i,2,9));
    for j = 1:numel(activityLocationType{i+1})
        % If location type is a space then select positions there
        if any(typeAbstractSpaces==activityLocationType{i+1}(j))
            tempNodes = find(abstractSpace...
                (typeAbstractSpaces==activityLocationType{i+1}).node);
            activityLocation{i+1} = [activityLocation{i+1} ...
                tempNodes];
        end
        tempNodes = find(allNodeTypes==activityLocationType{i+1});
        activityLocation{i+1} = [activityLocation{i+1} tempNodes];
    end
end


%% MC Input
numMCActivities = activity.continuous.num + activity.intermediate.num + ...
    activity.spontaneous.num;
ltpAgents = xlsread(files.st.input,'MCInput',...
    GetExcelRange(2,numAgents+1,3,2+numMCActivities));
stAgents = xlsread(files.st.input,'MCInput',...
    GetExcelRange(2,numAgents+1,3+numMCActivities,2+2*numMCActivities));

%% Habitual Input
startHabitual = xlsread(files.st.input,'Habitual',...
    GetExcelRange(2,numAgents+1,3,2+activity.habitual.num));
startHabitual = datetime(startHabitual,'ConvertFrom','excel');
startHIndex = startHabitual.Hour.*60 + startHabitual.Minute;

startHabitualVar = xlsread(files.st.input,'Habitual',...
    GetExcelRange(2,numAgents+1,3+activity.habitual.num,...
    2+2*activity.habitual.num));
startHabitualVar = datetime(startHabitualVar,'ConvertFrom','excel');
startHVarIndex = startHabitualVar.Hour.*60 + startHabitualVar.Minute;

durationHabitual = xlsread(files.st.input,'Habitual',...
    GetExcelRange(2,numAgents+1,3+2*activity.habitual.num,...
    2+3*activity.habitual.num));
durationHabitual = datetime(durationHabitual,'ConvertFrom','excel');
durationHIndex = durationHabitual.Hour.*60 + durationHabitual.Minute;

durationHabitualVar = xlsread(files.st.input,'Habitual',...
    GetExcelRange(2,numAgents+1,3+3*activity.habitual.num,...
    2+4*activity.habitual.num));
durationHabitualVar = datetime(durationHabitualVar,'ConvertFrom','excel');
durationHVarIndex = durationHabitualVar.Hour.*60 + durationHabitualVar.Minute;

%% Clear temp
clear temp*





