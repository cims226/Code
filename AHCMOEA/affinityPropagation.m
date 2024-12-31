function [idx, netsim, dpsim, expref] = affinityPropagation(S, max_iter, conv_iter, damping)
%--------------------------------------------------------------------------%
    if nargin < 4
        damping = 0.5; 
    end
    if nargin < 3
        conv_iter = 100; 
    end
    if nargin < 2
        max_iter = 1000; 
    end

    N = size(S, 1); 

    
    R = zeros(N, N);
    A = zeros(N, N);

    for iter = 1:max_iter
        
        old_R = R;
        AS = A + S;
        [max_AS, max_ind] = max(AS, [], 2);
        for i = 1:N
            AS(i, max_ind(i)) = -Inf;
        end
        [second_max_AS, ~] = max(AS, [], 2);
        R = S - repmat(max_AS, 1, N);
        for i = 1:N
            R(i, max_ind(i)) = S(i, max_ind(i)) - second_max_AS(i);
        end
        R = (1 - damping) * R + damping * old_R;

        
        old_A = A;
        Rp = max(R, 0);
        for k = 1:N
            Rp(k, k) = R(k, k);
        end
        A = repmat(sum(Rp, 1), N, 1) - Rp;
        dA = diag(A);
        A = min(A, 0);
        for k = 1:N
            A(k, k) = dA(k);
        end
        A = (1 - damping) * A + damping * old_A;

       
        E = (R + A) > 0;
        if iter >= conv_iter
            se = sum(E, 2);
            unconverged = length(find(sum(se, 1) ~= conv_iter));
            if unconverged == 0
                break;
            end
        end
    end

    
    [~, c] = max(R + A, [], 2);
    idx = unique(c);
    cluster_assignments = zeros(N, 1);
    for i = 1:length(idx)
        cluster_assignments(c == idx(i)) = i;
    end

    netsim = sum(S(sub2ind(size(S), (1:N)', c)));
    dpsim = sum(S(sub2ind(size(S), c, c)));
    expref = sum(S(sub2ind(size(S), (1:N)', c)) - S(sub2ind(size(S), c, c)));
    idx = cluster_assignments;
end



