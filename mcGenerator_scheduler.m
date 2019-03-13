%% Markov Chain Generator
function [P,fvalf] = mcGenerator_scheduler(ltp,st,timeStep,method,constraints)

ltpzeros = find(ltp==0);
ltpnz = find(ltp>0);
ltp(ltpzeros) = [];
st(ltpzeros) = [];

numStates = numel(ltp);
st = st/timeStep;

probSum = zeros(numStates,numStates*(numStates-1));
temp = repmat(1:numStates:numStates*(numStates-1),numStates,1)+...
    repmat((0:numStates-1)'*(numStates*(numStates-1)+1),1,numStates-1);
probSum(temp) = 1;
probSteady = zeros(numStates,numStates*(numStates-1));
for n = 1:numStates
    for j = 0:numStates-1
        if n~=j+1
            if j+1>n
                probSteady(n,j*(numStates-1)+n) = ltp(j+1);
            elseif j+1<n
                probSteady(n,j*(numStates-1)+(n-1)) = ltp(j+1);
            end
        end
    end
end
temp_C = [probSteady;probSum];
temp_pii = 1 - 1./st(:)';
temp_pii(isinf(temp_pii)) = 0; % [7/9/2017] when st,ltp for a state = 0;
temp_d = [ltp(:).*(1-temp_pii)';(1-temp_pii)'];

switch method
    case 'LsqLin'
        m = 10000;
        if nargin == 5
            A = constraints{1};
            b = constraints{2};
            Aeq = constraints{3};
            beq = constraints{4};
            lb = constraints{5};
            ub = constraints{6};
        else
            A = [];
            b = [];
            Aeq = [];
            beq = [];
            lb = zeros(size(temp_C,2),1);
            ub = [];
        end
        temp_options = optimoptions('lsqlin','Algorithm','interior-point',...
            'Display','off','MaxIterations',10000);
        [temp_sol,fvalf] = lsqlin(temp_C*m,temp_d*m,A,b,Aeq,beq,...
            lb,ub,[],temp_options);
        P = zeros(numStates);
        P(1,1) = temp_pii(1);
        P(1,2:numStates) = temp_sol(1:numStates-1);
        for n = 2:numStates
            P(n,n) = temp_pii(n);
            P(n,1:n-1) = temp_sol((n-1)*(numStates-1)+...
                1:(n-1)*(numStates-1)+n-1);
            P(n,n+1:numStates) = temp_sol((n-1)*(numStates-1)+...
                n:(n-1)*(numStates-1)+numStates-1);
        end
        
    case 'LsqNonNeg'
        [temp_sol,fvalf,~] = lsqnonneg(temp_C,temp_d);
        P = zeros(numStates);
        P(1,1) = temp_pii(1);
        P(1,2:numStates) = temp_sol(1:numStates-1);
        for n = 2:numStates
            P(n,n) = temp_pii(n);
            P(n,1:n-1) = temp_sol((n-1)*(numStates-1)+1:(n-1)*(numStates-1)+n-1);
            P(n,n+1:numStates) = temp_sol((n-1)*(numStates-1)+n:(n-1)*(numStates-1)+numStates-1);
        end
        
    case 'FminconIP'
        probSum = zeros(numStates,numStates^2); % Sum of probs = 1
        probKnown = zeros(numStates,numStates^2); % p_ii's known
        %see matrix arguments for fmincon MATLAB documentation
        temp = repmat(1:numStates^2:numStates^3,numStates,1)+...
            repmat((0:numStates-1)'*(numStates+1),1,numStates);
        probSum(temp) = 1;
        probKnown(diag(temp)) = 1;
        Aeq = [probSum;probKnown];
        beq = [ones(numStates,1);(1-1./st(:))];
        lb = zeros(numStates);
        ub = ones(numStates); % All probs between 0 & 1
        
        fvalf = 10;
        while fvalf>1e-6
            x0=2*rand(numStates)*1-1;%rand(numStates);%zeros(numStates)
            options = optimoptions('fmincon','Algorithm','interior-point',...
                'Display','off','MaxFunctionEvaluations', 100);
            [P,fvalf] = ...
                fmincon(@(x)ltpmin(x,numStates,ltp(:)'),...
                x0,[],[],Aeq,beq,lb,ub,[],options);
        end
                
end

if isempty(ltpzeros)==0
    temp = zeros(numel(ltpnz)+numel(ltpzeros));
    for i = 1:numel(ltpnz)
        for j = 1:numel(ltpnz)
            temp(ltpnz(i),ltpnz(j)) = P(i,j);
        end
    end
    P = temp;
end

end