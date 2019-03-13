%% Citation Information
%Please cite: Shelat, S. (2017). Developing an Integrated Pedestrian Behaviour Model for Office Buildings. (MSc), Delft University of Technology.   

%% Markov Chain Simulator
function [mc,st,rt,ta,ltp,tm] = mcSimulator_scheduler(P,totalTime,ltpGiven)

numStates = size(P,1);

mc = zeros(totalTime,1); % Markov Chain
rt = cell(numStates,1); % Recurrence Time
st = cell(numStates,1); % Sojourn Time
ta = cell(numStates,1); % Time away
tm = zeros(numStates); % Transition Frequency Matrix

t = 1;
startT = zeros(numStates,1); % Sojourn start time of each state
endT = zeros(numStates,1); % Sojourn end time of each state
ltpNz = find(~ltpGiven==0);
mc(t) = ltpNz(randperm(numel(ltpNz),1));
% mc(t) = round((numStates-1)*rand(1)+1); % Start at random point
startT(mc(t)) = t;

for t = 2:totalTime
    temp = cumsum(P(mc(t-1),:));
    temp(end) = 1; % make it 1 instead of almost 1
    mc(t) = find(temp>=rand(1),1);
    if mc(t) == mc(t-1)
        rt{mc(t)} = [rt{mc(t)} 1];
    else
        st{mc(t-1)} = [st{mc(t-1)} t-startT(mc(t-1))];
        rt{mc(t)} = [rt{mc(t)} t-endT(mc(t))];
        ta{mc(t)} = [ta{mc(t)} t-endT(mc(t))];
        startT(mc(t)) = t;
        endT(mc(t-1)) = t-1;
    end
    tm(mc(t-1),mc(t)) = tm(mc(t-1),mc(t)) + 1;
end

if mc(t) == mc(t-1) % Add last sojourn time
    st{mc(t)} = [st{mc(t)} (t)-startT(mc(t-1))];
end

ltp = zeros(numStates,1);
for i = 1:numStates
    ltp(i) = sum(st{i});
end
ltp = ltp/sum(ltp);
