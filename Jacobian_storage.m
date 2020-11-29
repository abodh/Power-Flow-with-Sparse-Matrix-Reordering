function Jacob_storage = Jacobian_storage(V, delta, n_bus, n_pq,...
    pq_bus_id, Y_next)

    % initialize P and Q
    P = zeros(n_bus,1);
    Q = zeros(n_bus,1);

    % calculating the active and reactive power at each bus
    %{
        P(i) = sum(j=1->n) |Vi||Vj|(Gij * cos(delta_i - delta_j) + 
                                    Bij * sin(delta_i - delta_j)
        Q(i) = sum(j=1->n) |Vi||Vj|(Gij * sin(delta_i - delta_j) - 
                                    Bij * cos(delta_i - delta_j)
        Note: 
             Gij and Bij are the real and imaginary part of Yij and hence
             are extracted accordingly from the sparse_table class
    %}
    for i = 1 : n_bus
        for j = 1 : n_bus
            
            P(i) = P(i) + V(i)*V(j)*(real(sparse_table.retrieve(Y_next,...
                i,j))* cos(delta(i)-delta(j))+imag(sparse_table.retrieve...
                (Y_next,i,j))*sin(delta(i)-delta(j)));
            
            Q(i) = Q(i) + V(i)*V(j)*(real(sparse_table.retrieve(Y_next,...
                i,j))*sin(delta(i)-delta(j))- imag(sparse_table.retrieve...
                (Y_next,i,j))*cos(delta(i)-delta(j)));
        end
    end

    %{
        initialize the sparse table for Jacobian
    Note:
    
    * The concept of sub-Jacobian does not hold here as we would have a
      single storage table instead of sub-matrices.
    
    * REMEMBER: 
        Jacob_storage = [idx value row col NIR NIC]
    
        store (Jacob_storage, [], values, row, col)
    
    Here, ~/[] means that we are not interested on that argument as we do 
          not need to know about the FIR or FIC values 
    
    * refer to the link below for sub-matrices equations:
        
        https://bit.ly/2GU1Xx0    
        B. Sereeter, C. Vuik, and C. Witteveen
        REPORT 17-07 
        On a comparison of Newton-Raphson solvers for power flow problems    
    %}
    
    [Jacob_storage, ~] = sparse_table.Blank_array();
    
    %{
      although we will not calculate each of the sub-matrices as mentioned
      above, we will still name such sub-matrices in the comment in order
      to understand the calculations better.
    %}
    
    % Calculating J11    
    for i = 2 : n_bus
        for k = 2 : n_bus
            if (i == k)
                [Jacob_storage, ~] = sparse_table.store(Jacob_storage,...
                    [], - Q(i) - (V(i)^2 * imag(sparse_table.retrieve...
                    (Y_next,i,i))), i-1, k-1);
            else
                [Jacob_storage, ~] = sparse_table.store(Jacob_storage,...
                    [],V(i)*V(k)*(real(sparse_table.retrieve...
                    (Y_next,i,k))*sin(delta(i)-delta(k))...
                    - imag(sparse_table.retrieve(Y_next,i,k))*...
                    cos(delta(i)-delta(k))), i-1, k-1);
            end
        end
    end

    % Calculating J21
    for i = 2 : n_pq + 1 
        j = pq_bus_id(i-1);
        for k = 2 : n_bus
            if (j == k)
                [Jacob_storage, ~] = sparse_table.store(Jacob_storage,...
                    [], P(j) - (V(j)^2 * real(sparse_table.retrieve...
                    (Y_next,j,j))), i-1+ n_bus-1, k-1);
            else
                [Jacob_storage, ~] = sparse_table.store(Jacob_storage,...
                    [],-V(j)*V(k)*(real(sparse_table.retrieve...
                    (Y_next,j,k))*cos(delta(j)-delta(k))...
                    + imag(sparse_table.retrieve(Y_next,j,k))*...
                    sin(delta(j)-delta(k))), i-1 + n_bus-1, k-1);                
            end
        end
    end
    
    % Calculating J12
    for i = 2 : n_bus  
        for k = 2 : n_pq + 1
            j = pq_bus_id(k-1);
            if (i == j)
                [Jacob_storage, ~] = sparse_table.store(Jacob_storage,...
                    [], P(j) + (V(j)^2 * real(sparse_table.retrieve...
                    (Y_next,j,j))), i-1, k-1 + n_bus-1);
                %J12(i-1,k-1) = P(j) + (V(j)^2 * G(j,j));
            else
                [Jacob_storage, ~] = sparse_table.store(Jacob_storage,...
                    [],V(i)*(real(sparse_table.retrieve...
                    (Y_next,i,j))*cos(delta(i)-delta(j))...
                    + imag(sparse_table.retrieve(Y_next,i,j))*...
                    sin(delta(i)-delta(j))), i-1, k-1 + n_bus-1);
                
                %J12(i-1,k-1) = V(i)*(G(i,j)*cos(delta(i)-delta(j))...
                    %+ B(i,j)*sin(delta(i)-delta(j)));
            end
        end
    end

    % Calculating J22
    for i = 2 : n_pq + 1
        j = pq_bus_id(i-1);
        for k = 2 : n_pq + 1
            l = pq_bus_id(k-1);
            if (j == l)
                [Jacob_storage, ~] = sparse_table.store(Jacob_storage,...
                    [], Q(j) - (V(j)^2 * imag(sparse_table.retrieve...
                    (Y_next,j,j))), i-1 + n_bus-1, k-1 + n_bus-1);
                
                %J22(i-1,k-1) = Q(j) - (V(j)^2 * B(j,j));
            else
                [Jacob_storage, ~] = sparse_table.store(Jacob_storage,...
                    [],V(j)*(real(sparse_table.retrieve...
                    (Y_next,j,l))*sin(delta(j)-delta(l))...
                    - imag(sparse_table.retrieve(Y_next,j,l))*...
                    cos(delta(j)-delta(l))), i-1 + n_bus-1, k-1 + n_bus-1);
                
                %J22(i-1,k-1) = V(j)*(G(j,l)*sin(delta(j)-delta(l))...
                %    - B(j,l)*cos(delta(j)-delta(l)));
            end
        end
    end
end