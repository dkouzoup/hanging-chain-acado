function [ found ] = contains( str, substr )

    if ~isempty(strfind(str, substr))
        found = 1;
    else
        found = 0;
    end

end

