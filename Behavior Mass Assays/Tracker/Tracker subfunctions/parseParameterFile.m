
function file_struct = parseParameterFile(filename, keywords)
%-- code modfied from:
%-- http://www.mathworks.com/matlabcentral/newsreader/view_thread/12122
%-- original author:
%-- derrico@kodak.com (John D'Errico)

%-- lines starting with '%' will be ignored
%-- only keys in the provided list of keywords will be included in the file structure
%-- values have to be terminated with ";", everything after the semicolon
%--     will be ignored
%-- values not terminated with a semicolon will raise an error

%-- start by checking that the file exists at all
if exist(filename, 'file') ~= 2
    error(['File not found: ', filename])
end;

%-- open file
fid = fopen(filename, 'r');
currLine = fgetl(fid);

%-- initialize output file structure
file_struct = struct;

checkDoubleKey = {};

i = 1;

%-- walk through each line of the file
while ischar(currLine)
    %-- ignore comment and empty lines
	if isempty(regexpi(currLine,'^%*')) && ~isempty(strtrim(currLine))
        if ~isempty(regexpi(currLine,'='))
            %-- get keyword and keyvalue in current line
            %-- remove leading and trailing blanks
            currKey = strtrim(currLine(1 : strfind(currLine, '=') - 1));
            currRest = strtrim(strrep(currLine(strfind(currLine,'=') + 1 : end),'''',''));

            %-- only include value, if line is terminated properly
            if ~isempty(regexpi(currRest,';'))
                %-- check if there is a value in front of the semicolon
                %-- otherwise escape with empty string
                if regexpi(currRest,';') > 1
                    %-- include only value before terminating semicolon
                    [currValue, ~] = strtok(currRest, ';');
                    currValue = strtrim(currValue);
                    %-- check if current value is a number and convert
                    %-- accordingly
                    if ~isempty(str2num(currValue))
                        currValue = str2num(currValue);
                    %-- check if input parameter is a cell array and add
                    %-- the values accordingly
                    elseif ~isempty(regexpi(currValue,'^{')) && ~isempty(regexpi(currValue,'}')) ...
                            && ~isempty(regexpi(currValue,','))
                        currValue = regexprep(strrep(strrep(strrep(currValue,'{',''),'}',''),'''',''),'\s+\,\s+',',');
                        tmpCell = {};
                        getValues = strfind(currValue, ',');
                        lastPosition = 1;
                        
                        for currCellVal = 1 : length(getValues)
                            tmpCell(1, size(tmpCell, 2) + 1) = {strtrim(currValue(lastPosition : getValues(currCellVal)-1))};
                            lastPosition = getValues(currCellVal)+1;
                        end;
                        tmpCell(1, size(tmpCell, 2) + 1) = {strtrim(currValue(lastPosition : end))};
                        currValue = tmpCell;
                    end;
                else
                    currValue = '';
                end;
            else
                errorMsg = strcat('Error at key ', 32, currKey, 32, ', line ', 32, num2str(i), ...
                    '. Key value pair was not properly terminated using a semicolon.');
                error(errorMsg);
            end;
        else
            errorMsg = strcat('Error at line ', 32, num2str(i), ...
                '. Key - value separator "=" is missing.');
            error(errorMsg);
        end;
        %-- check with keylist if current key is valid
        if sum(ismember(keywords, currKey));
            %-- only add first occurrence of current key to parameter structure
            if ~sum(ismember(checkDoubleKey, currKey))
                file_struct.(currKey) = currValue;
                checkDoubleKey{end+1} = currKey;
            else
                disp(strcat('... WARNING: key ', 32, currKey, 32, ...
                    'appears more often than once. Secondary occurrences will be ignored!'));
            end;
        end;
	end;
    %-- get next line in file
    currLine = fgetl(fid);
    i = i + 1;
end;

%-- close parameter file
fclose(fid);
end

