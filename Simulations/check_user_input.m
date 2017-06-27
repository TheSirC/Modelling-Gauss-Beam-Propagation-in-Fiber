function [ val ] = check_user_input( answer,i,defaultans )
%check_user_input Always return a valid user input
%   val is the value returned
%   i is the index of the value currently checked
%   defaultans is the cell containing default value to catch a value if
%   conversion fail

val = str2double(answer{i});  
if isnan(val)
    % Handle empty value returned 
    % for unsuccessful conversion
    val = str2double(defaultans{i});
end
end

