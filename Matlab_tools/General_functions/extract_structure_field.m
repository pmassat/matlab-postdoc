function y = extract_structure_field(struct,sfield,x)
% For a structure 'struct' containing a field 'sfield',
% extract the value of field 'x' of 'sfield' as a numeric vector. 
% This is the equivalent for a structure of extracting a column from a
% table, but there is no simpler way of doing it.
% Note: this only works if field 'x' is a numeric array.
    y = cell2mat( arrayfun(@(c) c.(sfield).(x), struct.', 'Uniform', 0) );
end