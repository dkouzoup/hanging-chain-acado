function updateProblemDimensions(NI, NX, NU)

cd include

f   = fileread('gdfgm_dimensions.h');
lim = 10; % limit in characters

if nargin == 1
    % only change horizon
    NX = [];
    NU = [];
end
if nargin ~= 1 && nargin ~= 3
    error('Wrong number of inputs')
end

% Set NX
if ~isempty(NX)
    
    NXstr = num2str(NX);
    
    if length(NXstr) > lim
        error('Too many states');
    end
    
    while length(NXstr) < lim
        NXstr = [NXstr ' '];
    end
    
    NXpos = strfind(f,'#define NX');
    
    f(NXpos+11:NXpos+10+lim) = NXstr;
    
end

% Set NU
if ~isempty(NU)
    
    NUstr = num2str(NU);
    
    if length(NUstr) > lim
        error('Too many inputs');
    end
    
    while length(NUstr) < lim
        NUstr = [NUstr ' '];
    end
    
    NUpos = strfind(f,'#define NU');
    
    f(NUpos+11:NUpos+10+lim) = NUstr;
    
end

% Set NI
NIstr = num2str(NI);

if length(NIstr) > lim
    error('Too long horizon');
end

while length(NIstr) < lim
    NIstr = [NIstr ' '];
end

NIpos = strfind(f,'#define NI');

f(NIpos+11:NIpos+10+lim) = NIstr;


% Update NP
NPstr = num2str((NX+NU)*NI+NX);

if length(NPstr) > lim
    error('Too many primal variables');
end

while length(NPstr) < lim
    NPstr = [NPstr ' '];
end

NPpos = strfind(f,'#define NP');

f(NPpos+11:NPpos+10+lim) = NPstr;


% Update ND
NDstr = num2str((NI+1)*NX);

if length(NDstr) > lim
    error('Too many dual variables');
end

while length(NDstr) < lim
    NDstr = [NDstr ' '];
end

NDpos = strfind(f,'#define ND');

f(NDpos+11:NDpos+10+lim) = NDstr;


% Write new file and overwrite

fid = fopen('temp.h','wt');
fprintf(fid,'%s\n',f);
fclose(fid);

copyfile('temp.h', 'gdfgm_dimensions.h');
delete('temp.h');

cd ..

end




