function [ val ] = check_user_input( answer,i )
%check_user_input Always return a valid user input
%   val is the value returned
%   i is the index of the value currently checked

[val, status] = str2double(answer{i});  
if ~status
    % Handle empty value returned 
    % for unsuccessful conversion
    val = 0;
end
end

