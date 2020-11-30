function error_VD = Sparse_LU(A, b)

    % calculate the number of rows/cols in the Jacobian
    N = max(A(:,3));
    
    % create storage for L and U
    [L_next,L_first] = sparse_table.Blank_array();
    [U_next,U_first] = sparse_table.Blank_array();
    
    for j = 1:N
        for k = j:N
            multiple = 0;
            for m = 1:N
                %{
                    a separate loop, different from the previous work, was
                    created as the class function did not return any value
                    for multiple slicings for eg. 1:n argument did not work
                %}
                multiple = multiple + sparse_table.retrieve(L_next,k,m)*...
                    sparse_table.retrieve(U_next,m,j);
            end
            [L_next,L_first] = sparse_table.store(L_next,L_first,...
                sparse_table.retrieve(A,k,j) - multiple,k,j);
        end
        if j == N
            % for U the final value is 1 so no need to calculate
            break
        end
        for k = j+1 : N
            multiple = 0;
            for m = 1:j-1
                multiple = multiple + sparse_table.retrieve(L_next,j,m)...
                    *sparse_table.retrieve(U_next,m,k);
            end
            [U_next,U_first] = sparse_table.store(U_next,U_first,...
                (sparse_table.retrieve(A,j,k) - multiple)/...
                sparse_table.retrieve(L_next,j,j),j,k);
        end
    end
    for i=1:N
        [U_next,U_first]=sparse_table.store(U_next,U_first,...
            sparse_table.retrieve(U_next,i,i)+1,i,i);
    end

    % forward substitution
    Nb = length(b); % total length of b (A is not a matrix so b is used)
    y = zeros(Nb,1);
    for i = 1:Nb
            if i == 1
                y(1) = b(1)/sparse_table.retrieve(L_next,1,1);
            else
                multiple = 0;
                for j = 1:i-1
                    multiple = multiple + sparse_table.retrieve...
                        (L_next,i,j) * y(j);
                end
                y(i) = (b(i) - multiple)/sparse_table.retrieve(L_next,i,i);
            end
    end
    
    % backward substitution
    x = zeros(Nb,1);
    for i = Nb:-1:1   
        if i == Nb
            x(Nb) = y(Nb)/sparse_table.retrieve(U_next,Nb,Nb);
        else
            multiple = 0;
            for j = i+1:Nb
                multiple = multiple + sparse_table.retrieve...
                        (U_next,i,j) * x(j);
            end
            x(i) = (y(i)-multiple)/sparse_table.retrieve(U_next,i,i);
        end
    end
    error_VD = x;
end