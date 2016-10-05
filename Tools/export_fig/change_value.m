%CHANGE_VALUE  Change the value assigned to a unique variable in a file
%
% Examples:
%   fail = change_value(fullPath, variableName, value)
%
% Function to change the value assigned to a variable in a text file. The
% assignment must exist already, and must be on a line on its own. Only the
% first such assignment is changed.
%
% IN:
%   filePath - Full path of the file to change.
%   variableName - String containing the name of the variable whose value
%                  is to be set.
%   value - The value to be assigned to the variable.
%
% OUT:
%   fail - true if change failed, false otherwise.

function fail = change_value(filePath, variableName, value)
fail = true;
% Read in the file
fh = fopen(filePath, 'rt');
if fh < 0
    return
end
fstrm = fread(fh, '*char')';
fclose(fh);
% Find the path
first_sec = regexp(fstrm, ['[\n\r]+ *' variableName ' *= *'], 'end', 'once');
second_sec = first_sec + regexp(fstrm(first_sec+1:end), ';? *[\n\r]+', 'once');
if isempty(first_sec) || isempty(second_sec)
    return
end
% Create the changed line
if ischar(value)
        str = '''%s''';
else
        str = '%1.50g';
end
str = sprintf(str, value);
% Save the file with the value changed
fh = fopen(filePath, 'wt');
if fh < 0
    return
end
fprintf(fh, '%s%s%s', fstrm(1:first_sec), str, fstrm(second_sec:end));
fclose(fh);
fail = false;
return