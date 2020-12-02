function error_VD = Sparse_LU(A_next,A_first, b, orders)

    % calculate the number of rows/cols in the Jacobian
    N = max(A_next(:,3));
    
    % create storage for L and U
    [L_next,L_first] = sparse_table.Blank_array();
    [U_next,U_first] = sparse_table.Blank_array();
    for j = 1:N
        for k = j:N
%             multiple = 0;
%             for m = 1:N
%                 %{
%                  a separate loop, different from the previous work, was
%                  created as the class function did not return any value
%                  for multiple slicings for eg. 1:n argument did not work
%                 %}
%                 multiple = multiple + sparse_table.retrieve(L_next,...
%                 k,m)* sparse_table.retrieve(U_next,m,j);
%             end

            % Only the non-zero values are to be selected and hence we look
            % at either the FIR/FIC column or NIR/NIC column to get the
            % values that are necessary for the calculation
            
            [L_next,L_first] = sparse_table.store(L_next,L_first,...
                sparse_table.retrieve(A_next,A_first,orders(k),...
                orders(j)) - Product(L_next,L_first, U_next, U_first,...
                k,j), k, j);
        end
        if j == N
            % for U the final value is 1 so no need to calculate
            break
        end
        for k = j+1 : N
%             multiple = 0;
%             for m = 1:j-1
%                 multiple = multiple + sparse_table.retrieve(L_next,...
%                 j,m)*sparse_table.retrieve(U_next,m,k);
%             end
            [U_next,U_first] = sparse_table.store(U_next,U_first,...
                (sparse_table.retrieve(A_next,A_first,orders(j),...
                orders(k)) - Product(L_next,L_first, U_next, U_first,...
                j,k))/sparse_table.retrieve(L_next,L_first,j,j),j,k);
        end
    end
    for i=1:N
        % the diagonal element of U is just 1
        [U_next,U_first]=sparse_table.store(U_next,U_first,...
            sparse_table.retrieve(U_next,U_first,i,i)+1,i,i);
    end

    % forward substitution
    Nb = length(b); % total length of b (A is not a matrix so b is used)
    y = zeros(Nb,1);
    for i = 1:Nb
            if i == 1
                y(1) = b(1)/sparse_table.retrieve(L_next,L_first,1,1);
            else
                multiple = 0;
                for j = 1:i-1
                    multiple = multiple + sparse_table.retrieve...
                        (L_next,L_first,i,j) * y(j);
                end
                y(i) = (b(i) - multiple)/sparse_table.retrieve(L_next,...
                    L_first,i,i);
            end
    end
    
    % backward substitution
    x = zeros(Nb,1);
    for i = Nb:-1:1   
        if i == Nb
            x(Nb) = y(Nb)./sparse_table.retrieve(U_next,U_first,Nb,Nb);
        else
            multiple = 0;
            for j = i+1:Nb
                multiple = multiple + sparse_table.retrieve...
                        (U_next,U_first,i,j) * x(j);
            end
            x(i) = (y(i)-multiple)./sparse_table.retrieve(U_next,U_first,...
                i,i);
        end
    end
    error_VD = x;
end

function multiple = Product(A_next,A_first, B_next, B_first,...
                a,b)
    % here A is Lower-triangular and B is Upper-triangular
            
    if isempty(A_next) || a > size(A_first)*[1;0] || b > size(B_first)...
            *[1;0]
        % checks if the NIR/NIC table is empty or if FIR/FIC table is 
        % smaller than the searched index
        multiple=0;
        return
    end
    multiple = 0;

    % NOTE: table is stored as [ Idx | FIR | FIC ]
    
    % first in row of A
    FIR_A = A_first(a,2);
    
    % first in column of B
    FIC_B = B_first(b,3);      


    while 1
        % Infinite loop until returned
        % i.e. returns only if we reach the end of NIR or NIC 

        if FIR_A == 0 || FIC_B == 0
            % returns if null
            return

        elseif A_next(FIR_A,4) > B_next(FIC_B,3)
            % A_next(FIR,col) > B_next(FIC, row) 
            % if the column of A is greater than row of B then the FIC of B
            % should be its NIC
            
            % NIC
            FIC_B = B_next(FIC_B,6);   

        elseif A_next(FIR_A,4) == B_next(FIC_B,3)
            % calculation is done if the column of A is equal to the 
            % row of B i.e. go one column down for A and same row right for
            % B
            
            multiple = multiple + A_next(FIR_A,2) * B_next(FIC_B,2);
            
            % now go for NIR of A
            FIR_A = A_next(FIR_A,5);
            
            % go for NIC of B
            FIC_B = B_next(FIC_B,6);            
        else
            % same as A_next(FIR,col) < B_next(FIC, row)
            % NIR
            FIR_A = A_next(FIR_A,5);
        end
    end               
end