function G = ComputeStageCosts( stateSpace, controlSpace, map, gate, mansion, cameras )
%COMPUTESTAGECOSTS Compute stage costs.
% 	Compute the stage costs for all states in the state space for all
%   control inputs.
%
%   G = ComputeStageCosts(stateSpace, controlSpace, map, gate, mansion,
%   cameras) computes the stage costs for all states in the state space
%   for all control inputs.
%
%   Input arguments:
%
%       stateSpace:
%           A (K x 2)-matrix, where the i-th row represents the i-th
%           element of the state space.
%
%       controlSpace:
%           A (L x 1)-matrix, where the l-th row represents the l-th
%           element of the control space.
%
%       map:
%           A (M x N)-matrix describing the terrain of the estate map.
%           Positive values indicate cells that are inaccessible (e.g.
%           trees, bushes or the mansion) and negative values indicate
%           ponds or pools.
%
%   	gate:
%          	A (1 x 2)-matrix describing the position of the gate.
%
%    	mansion:
%          	A (F x 2)-matrix indicating the position of the cells of the
%           mansion.
%
%    	cameras:
%          	A (H x 3)-matrix indicating the positions and quality of the 
%           cameras.
%
%   Output arguments:
%
%       G:
%           A (K x L)-matrix containing the stage costs of all states in
%           the state space for all control inputs. The entry G(i, l)
%           represents the expected stage cost if we are in state i and 
%           apply control input l.
map = map';

K = size(stateSpace,1);
P = ComputeTransitionProbabilities( stateSpace, controlSpace, map', gate, mansion, cameras );
G = zeros(K,5);
x_gate = gate(1);
y_gate = gate(2);
map_camera_prob = computecameraprob( stateSpace, controlSpace, map, gate, mansion, cameras );
for i = 1:K
    for l = 1:5
        G(i, l) = inf;
        cost = 0.0;
        cumul_proba = 0.0;
        for j = 1:K
            proba = P(i,j,l);
%             cumul_proba = cumul_proba + proba;
            if(proba>0.0)
                x = stateSpace(i,1);
                y = stateSpace(i,2);
                x_dest = stateSpace(j,1);
                y_dest = stateSpace(j,2);

                if(x_dest == x_gate && y_dest == y_gate && x == x_gate && y == y_gate)
                    caught = (1 - 0.001) * map_camera_prob(x_gate, y_gate);
                    cost = cost + caught*7 + (1 - caught)*(1-0.001);
                    cumul_proba = caught + (1 - caught)*(1-0.001);
                elseif(x_dest == x_gate && y_dest == y_gate && x == x_gate && y == y_gate - 1 && controlSpace(l) == 'n')
                    caught_3 = proba * map_camera_prob(x_gate, y_gate);
                    cost = cost + caught_3*7 + (1-caught_3);
                    cumul_proba = cumul_proba + caught_3 + (1-caught_3);
                elseif(x_dest == x_gate && y_dest == y_gate && x == x_gate && y == y_gate + 1 && controlSpace(l) == 's')
                    caught_2 = proba * map_camera_prob(x_gate, y_gate);
                    cost = cost + caught_2*7 + (1-caught_2);
                    cumul_proba = cumul_proba + caught_2 + (1-caught_2);
                elseif(x_dest == x_gate && y_dest == y_gate && x == x_gate - 1 && y == y_gate && controlSpace(l) == 'w')
                    caught_3 = proba * map_camera_prob(x_gate, y_gate);
                    cost = cost + caught_3*7 + (1-caught_3);
                    cumul_proba = cumul_proba + caught_3 + (1-caught_3);
                elseif(x_dest == x_gate && y_dest == y_gate && x == x_gate + 1 && y == y_gate && controlSpace(l) == 'e')
                    caught_2 = proba * map_camera_prob(x_gate, y_gate);
                    cost = cost + caught_2*7 + (1-caught_2);
                    cumul_proba = cumul_proba + caught_2 + (1-caught_2);
                elseif(x_dest == x_gate && y_dest == y_gate && (x ~= x_gate || y ~= y_gate))
                    cost = cost + proba * 7;
                    cumul_proba = cumul_proba + proba;
                elseif(map(x_dest, y_dest) < 0.0 && (x_dest ~= x || y_dest ~=y))
                    cost = cost + proba * 4;
                    cumul_proba = cumul_proba + proba;
                elseif(map(x_dest, y_dest) <= 0.0)
                    cost = cost + proba * 1;
                    cumul_proba = cumul_proba + proba;

                end
            end
        end
        if(cost > 0.0)
            G(i,l) = cost + (1 - cumul_proba);
        end
    end
end
end

function map_camera_prob = computecameraprob(stateSpace, controlSpace, map, gate, mansion, cameras )
M = size(map,1);
N = size(map,2);
H = size(cameras,1);

%% Computing a matrix which gives at each spot the probability of being seen by a camera
map_camera_prob = zeros(M,N);
for i = 1:H
    C = cameras(i,1:3);
    x = C(1);
    y = C(2);
    z = C(3);
    for m = x+1:M
        if(map(m,y)>0.0)
            break;
        else
            prob = z/(m-x);
            if(map(m,y)<0.0)
                prob = 1 - (1 - prob)^4;
            end
            map_camera_prob(m,y) = 1 - (1 - map_camera_prob(m,y))*(1 - prob);
        end
    end
    for m = 1:x-1
        invert_m = x - m;
        if(map(invert_m,y)>0.0)
            break;
        else
            prob = z/(x-invert_m);
            if(map(invert_m,y)<0.0)
                prob = 1 - (1 - prob)^4;
            end
            map_camera_prob(invert_m,y) = 1 - (1 - map_camera_prob(invert_m,y))*(1 - prob);
        end
    end    
    for n = y+1:N
        if(map(x,n)>0.0)
            break;
        else
            prob = z/(n-y);
            if(map(x,n)<0.0)
                prob = 1 - (1 - prob)^4;
            end
            map_camera_prob(x,n) = 1 - (1 - map_camera_prob(x,n))*(1 - prob);
        end
    end
    for n = 1:y-1
        invert_n = y - n;
        if(map(x,invert_n)>0.0)
            break;
        else
            prob = z/(y-invert_n);
            if(map(x,invert_n)<0.0)
                prob = 1 - (1 - prob)^4;
            end
            map_camera_prob(x,invert_n) = 1 - (1 - map_camera_prob(x,invert_n))*(1 - prob);
        end
    end 
end
end
