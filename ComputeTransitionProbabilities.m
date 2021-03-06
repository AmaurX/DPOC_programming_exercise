function P = ComputeTransitionProbabilities( stateSpace, controlSpace, map, gate, mansion, cameras )
%COMPUTETRANSITIONPROBABILITIES Compute transition probabilities.
% 	Compute the transition probabilities between all states in the state
%   space for all control inputs.
%
%   P = ComputeTransitionProbabilities(stateSpace, controlSpace,
%   map, gate, mansion, cameras) computes the transition probabilities
%   between all states in the state space for all control inputs.
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
%       P:
%           A (K x K x L)-matrix containing the transition probabilities
%           between all states in the state space for all control inputs.
%           The entry P(i, j, l) represents the transition probability
%           from state i to state j if control input l is applied.

% put your code here
K = size(stateSpace,1);
L = size(controlSpace,1);
P = zeros(K,K,L);
pc=0.001;

%for some reason, map is transverse, so I untransverse it.
map = map';


M = size(map,1);
N = size(map,2);

H = size(cameras,1);
F = size(mansion, 1);
x_gate = gate(1);
y_gate = gate(2);

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

map_camera_prob

%% Computing a matrix that gives at each spot the probability of successfully taking a picture
map_pic_prob = zeros(M,N);
for i = 1:M
    for j = 1:N
        if(map(i,j)<= 0.0)
            map_pic_prob(i,j) = pc;
        end
    end
end
z = 0.5;
for i = 1:F
    C = mansion(i,1:2);
    x = C(1);
    y = C(2);
    
    for m = x+1:M
        if(map(m,y)>0.0)
            break;
        else
            prob = max(pc, z/(m-x));
            map_pic_prob(m,y) = max(map_pic_prob(m,y), prob);
        end
    end
    for m = 1:x-1
        invert_m = x - m;
        if(map(invert_m,y)>0.0)
            break;
        else
            prob = max(pc, z/(x-invert_m));

            map_pic_prob(invert_m,y) = max(map_pic_prob(invert_m,y), prob);
        end
    end    
    for n = y+1:N
        if(map(x,n)>0.0)
            break;
        else
            prob = max(pc, z/(n-y));

            map_pic_prob(x,n) = max(map_pic_prob(x,n), prob);
        end
    end
    for n = 1:y-1
        invert_n = y - n;
        if(map(x,invert_n)>0.0)
            break;
        else
            prob = max(pc, z/(y-invert_n));
            map_pic_prob(x,invert_n) = max(map_pic_prob(x,invert_n), prob);
        end
    end 
end

%% computing the transitionprobmatrix
for i = 1:K
    for j = 1:K
        for l = 1:L
            % to compute
            S = stateSpace(i,:);
            x = S(1);
            y = S(2);
            
            S_dest = stateSpace(j,:);
            x_dest = S_dest(1);
            y_dest = S_dest(2);
            
            control = controlSpace(l); 
            if(control== 'p')
                if(x_dest == x && x == x_gate && y_dest == y && y == y_gate)
                    P(i,j,l) = 1 - map_pic_prob(x,y);
                elseif(x_dest == x_gate && y_dest == y_gate)
                    P(i,j,l) = (1 - map_pic_prob(x,y)) * map_camera_prob(x,y);
                elseif(x_dest==x && y_dest==y)
                    P(i,j,l) = (1 - map_pic_prob(x,y)) * (1 - map_camera_prob(x,y));
                end
            elseif(control == 'e')
                if(x < M && map(x+1, y) <= 0.0)
                    if(x_dest == x +1 && x +1 == x_gate && y_dest == y && y == y_gate)
                        P(i,j,l) = 1;
                    elseif(x_dest == x+1 && y_dest == y)
                        P(i,j,l) = 1 - map_camera_prob(x+1,y);
                    elseif(x_dest == x_gate && y_dest == y_gate)
                        P(i,j,l) = map_camera_prob(x+1,y);
                    end
                end
            elseif(control == 'w')
                if(x > 1 && map(x-1, y) <= 0.0)
                    if(x_dest == x -1 && x -1 == x_gate && y_dest == y && y == y_gate)
                        P(i,j,l) = 1;
                    elseif(x_dest == x-1 && y_dest == y)
                        P(i,j,l) = 1 - map_camera_prob(x-1,y);
                    elseif(x_dest == x_gate && y_dest == y_gate)
                    P(i,j,l) = map_camera_prob(x-1,y);
                    end
                end
            elseif(control == 'n')
                if(y<N && map(x, y+1) <= 0.0)
                    if(x_dest == x && x == x_gate && y_dest == y+1 && y+1 == y_gate)
                        P(i,j,l) = 1;
                    elseif(x_dest == x && y_dest == y+1)
                        P(i,j,l) = 1 - map_camera_prob(x,y+1);
                    elseif(x_dest == x_gate && y_dest == y_gate)
                        P(i,j,l) = map_camera_prob(x,y+1);
                    end
                end
            else
                
                if(y>1 && map(x, y-1) <= 0.0)
                    if(x_dest == x && x == x_gate && y_dest == y-1 && y-1 == y_gate)
                        P(i,j,l) = 1;
                    elseif(x_dest == x && y_dest == y-1)
                        P(i,j,l) = 1 - map_camera_prob(x,y-1);
                    elseif(x_dest == x_gate && y_dest == y_gate)
                        P(i,j,l) = map_camera_prob(x,y-1);
                    end
                end
            end
        end
    end
end

end
