function [ J_opt, u_opt_ind ] = PolicyIteration( P, G )
%POLICYITERATION Value iteration
%   Solve a stochastic shortest path problem by Policy Iteration.
%
%   [J_opt, u_opt_ind] = PolicyIteration(P, G) computes the optimal cost and
%   the optimal control input for each state of the state space.
%
%   Input arguments:
%
%       P:
%           A (K x K x L)-matrix containing the transition probabilities
%           between all states in the state space for all control inputs.
%           The entry P(i, j, l) represents the transition probability
%           from state i to state j if control input l is applied.
%
%       G:
%           A (K x L)-matrix containing the stage costs of all states in
%           the state space for all control inputs. The entry G(i, l)
%           represents the cost if we are in state i and apply control
%           input l.
%
%   Output arguments:
%
%       J_opt:
%       	A (K x 1)-matrix containing the optimal cost-to-go for each
%       	element of the state space.
%
%       u_opt_ind:
%       	A (K x 1)-matrix containing the index of the optimal control
%       	input for each element of the state space.

% put your code here

K = size(P,1);
err = 1e-10;

%% Policy Iteration

% Initialize with random proper policy
policy = ones(1,K); % Go up for all states

J_i = sym('j', [1 K]);
eqn = []; 
J_curr = zeros(1,K);
J_prev = zeros(1,K);

while true
    % Policy Evaluation 
    cost = 0;
    for curr_state = 1:K % for every state i
        for next_state = 1:K % for every state j
           cost =+ P(curr_state, next_state, policy(curr_state))*J_i(next_state);
        end
        eqn = [eqn, G(curr_state, policy(curr_state)) + cost == J_i(curr_state)];
    end
    soln = struct2cell(solve(eqn, J_i));
    J_curr = double([soln{:}]);        
        
    % Policy Improvement
    new_cost = 0;
    net_cost = zeros(1,5);
    for curr_state = 1:K % for every state i
        for u = 1:5
            for next_state = 1:K % for every state j                
                new_cost =+ P(curr_state, next_state, u)*J_curr(next_state);
            end
            net_cost(u) = G(curr_state, u) + new_cost;
        end
        [~, policy(curr_state)] = min(net_cost);
    end    
    if all(abs(J_curr - J_prev) <= err)        
        break;
    end
    J_prev = J_curr;
    J_curr = zeros(1,K);
    eqn = [];
end

J_opt = J_curr;
u_opt_ind = policy;

end

