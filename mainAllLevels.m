%% Main
% Store all the programs below and (1) nearestPointOnLine.m, (2)
% lineSegmentIntersect.m, (3) dijkstra.m, (4) location and edges files,

getFiles
getSettings
inputBuildingPlan2_2
createNetwork3_2
routing2
strategicLevel
%}
% numDays = 15;
% activity_day = cell(numDays,1);
% location_day = cell(numDays,1);
% schedule_day = cell(numDays,1);
for day = 1%[1:5 7:10 13]%1:numDays
    tacticalLevel_2_1
    walking2_1
    scheduleExecuted
%     tempDay = sprintf('day%d',day);
%     csvwrite([pwd '/Output/Activity/' tempDay '.dat'],agentExecutedActivity);
%     csvwrite([pwd '/Output/Location/' tempDay '.dat'],agentExecutedLocation);
%     csvwrite([pwd '/Output/Schedule/' tempDay '.dat'],agentExecutedSchedule);
%     dlmwrite([pwd '/Output/Trajectory/' tempDay '.dat'],trajectoryAll,'precision',8);
    
%     location_day{day} = agentExecutedLocation;
%     schedule_day{day} = agentExecutedSchedule;
end

