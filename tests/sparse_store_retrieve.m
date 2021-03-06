% [A,B] = Sparse();
% example from the book

A = []; %[index VALUE NROW NCOL NIR NIC] 
B = []; %[index FIR FIC]
[A,B] = store(A,B,-1,1,1);
[A,B] = store(A,B,-2,1,3);
[A,B] = store(A,B,2,2,1);
[A,B] = store(A,B,8,2,2);
[A,B] = store(A,B,1,2,4);
[A,B] = store(A,B,3,3,3);
[A,B] = store(A,B,-2,3,5);
[A,B] = store(A,B,-3,4,2);
[A,B] = store(A,B,2,4,3);
[A,B] = store(A,B,1,5,1);
[A,B] = store(A,B,2,5,2);
[A,B] = store(A,B,-4,5,5);

% test on adding an extra element
% [A,B] = store(A,B,2,3,4);
 
% function to store the elements in the sparse tabular format
function [A,B] = store(A,B,value,row,col)
    if value == 0
        % does not do anything when the value is 0 
        return

    elseif isempty(A) == 1
        % initialize the storage table for the first entry
        A=[1,value,row,col,0,0];
        B=[1,1,1];
        return

    else
        %identify the element just before the specified element in the row
        prev_in_row = max(A((A(:,3) == row)&(A(:,4) < col),4));
        if isempty(prev_in_row)
            % if no element before the specified element prev_in_row = 0
            prev_in_row = 0;
        end
        
        %identify the element just after the specified element in the row
        next_in_row = min(A((A(:,3) == row)&(A(:,4) > col),4));
        if isempty(next_in_row)
            % if no element before the specified element next_in_row = 0
            next_in_row = 0;
        end
        
        %identify the element just before the specified element in the col
        prev_in_col = max(A((A(:,3) < row) & (A(:,4) == col),3));
        if isempty(prev_in_col)
            % if no element before the specified element prev_in_col = 0
            prev_in_col = 0;
        end
        
        %identify the element just after the specified element in the col
        next_in_col = min(A((A(:,3) > row) & (A(:,4) == col),3));
        if isempty(next_in_col)
             % if no element before the specified element next_in_col = 0
            next_in_col=0;
        end
        
        NIR = A((A(:,3) == row) & (A(:,4) == next_in_row),1);
        if isempty(NIR)
            NIR=0;
        end
        NIC = A((A(:,3) == next_in_col) & (A(:,4) == col),1);
        if isempty(NIC)
            NIC=0;
        end
        
        % increase the counter for the index and insert the elements in A
        A=[A;size(A)*[1;0]+1,value,row,col,NIR,NIC];
        
        % if any prev in row or prev in col elements are identified then
        % NIR and NIC for the previous element is updated accordingly
        A(A(:,3) == row & A(:,4) == prev_in_row,5) = size(A)*[1;0];
        A(A(:,3) == prev_in_col & A(:,4) == col,6) = size(A)*[1;0];
        
        % the snippet below fills the FIR and FIC
        if prev_in_row == 0
            % if there is no prev in row then this is the first element
            B(row,1) = row;
            B(row,2) = size(A)*[1;0];
        end
        
        if prev_in_col == 0
            % if there is no prev in col then this is the first element
            B(col,3) = size(A)*[1;0];
        end
    end
end

% function to retrieve the elements from the sparse tabular format
function value = retrieve(A,row,col)
    if isempty(A)
        % return 0 if the sparse storage is empty
        value=0;
        return
    end
    
    % retreive the value for which row and col match in the sparse storage
    value = A((A(:,3) == row)&(A(:,4) == col),2);
    if isempty(value)
        % return 0 if the value is empty
        value = 0;
    end
end