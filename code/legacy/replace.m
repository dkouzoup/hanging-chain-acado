function [ str ] = replace( str, substr, newsubstr )

    len = length(substr);
      
    allpos = strfind(str, substr);
    pos = allpos(1);

    str = [str(1:pos-1) newsubstr str(pos+len:end)];

    
end

