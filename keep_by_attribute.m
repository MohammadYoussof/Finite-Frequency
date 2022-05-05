function [S, A] = keep_by_attribute(s, a, field, value)



%
% S. HSIANG
% SMH2137@COLUMBIA.EDU
% 4.11
%
% ---------------------
%
% [S, A] = keep_by_attribute(s, a, field, value)
%
% Generates a new polygon and attribute structure [S, A] that is identical
% to [s, a] except that any element for which strcmp(FIELD,VALUE) == false
% is dropped. FIELD and VALUE must both be strings.
%
% Useful for keeping sets of polygons, for example, counties in a specific 
% state.
%
% see also combine_attributes, drop_by_attribute



elements = [];

for i = 1:length(a)
    command = ['if strcmp(a(i).' field ',''' value ''')==true; elements = [elements, i]; end'];
    eval(command)    
end


S = s(elements);
A = a(elements);



