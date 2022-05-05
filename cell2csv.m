function [Q] = cell2csv(cell_array, filename)

%
% S. HSIANG 
% SMH2137@COLUMBIA.EDU
% 4/10
%
% ----------------------------------
%
% cell2csv(cell_array, filename)
%
% Writes out the elements in a cell array to a comma-seperated values (csv)
% file (eg. for use in STATA or EXCEL). CELL_ARRAY is the cell array to be
% converted and FILENAME is a string to be used and will '.csv'
% appended to it.  FILENAME.csv will be written out in the current
% working directory. Elements in CELL_ARRAY can be strings or scalar-valued
% numeric. Note that if strings have commas in them, this will cause
% individual rows of the output file to be miss-aligned. A header row, for
% example a row that specifies the names of variables in statistical data,
% will be written out properly so long as it is formatted in the original
% cell array.

p = cell_array;
S = size(p);

filename = [filename '.csv'];

fid= fopen(filename,'wt');

for eta = 1:S(1);
    
    if ischar(p{eta,1})
        var_list =p{eta,1};
    else
        var_list =num2str(p{eta,1});
    end

    if length(S) > 1
        for i = 2:S(2)
            if ischar(p{eta,i}) == 1
                var_list = [var_list ',' p{eta,i}];
            else
                var_list = [var_list ',' num2str(p{eta,i})];
            end
        end
    end
    
    fprintf(fid,[var_list '\n']);
    
end


fclose(fid);
Q = true;


