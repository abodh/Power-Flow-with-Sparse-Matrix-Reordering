function [Volt, Angle, error_avg] = ...
    NewtonRhapson(tolerance, n_bus, n_pv, n_pq, pq_bus_id, V, delta, ...
    Y_next, Y_first, Ps, Qs, Scheme)
    
    % Power flow iteration
    iter = 0;
    mismatch = ones(2 * n_bus - n_pv - 2, 1);
    
    while any(abs(mismatch(:)) > tolerance)
        iter = iter + 1;
        
        % accumulate V and delta for all iteration
        Volt(:,iter) = V;
        Angle(:,iter) = delta;
        
        % calculate Jacobian
        [J_next, J_first] = Jacobian_storage(V, delta, n_bus, n_pq, ...
            pq_bus_id, Y_next, Y_first);
        
        % reorder Jacobian and return the ordered index given the Scheme
        [ordering_idx, original_idx] = ...
            Sparse_reordering(J_next, J_first, Scheme);
        
        % calculate power mismatch
        mismatch = power_mismatch(Ps, Qs, Y_next,Y_first, V, delta,...
            n_bus, pq_bus_id);
        
        % orders the mismatches according to the ordering scheme
        mismatch = mismatch(ordering_idx);
                
        % find the error [del_delta | del_V]
        error = Sparse_LU(J_next,J_first,mismatch,...
            ordering_idx);
        
        % error must be in original order for proper calculation 
        error = error(original_idx);
        
        % update V and delta
        delta(2:end) = delta(2:end) + error(1 : n_bus-1);
        V(pq_bus_id) = V(pq_bus_id) + error(n_bus : end);
        
        % average error
        error_avg(iter) = abs(mean(error(:)));
    end
end