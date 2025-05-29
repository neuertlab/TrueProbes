function m = Match(s)
if (strcmpi(s,'a')==1)
    m = 't';
elseif (strcmpi(s,'t')==1)
    m = 'a';
elseif (strcmpi(s,'g')==1)
    m = 'c';
elseif (strcmpi(s,'c')==1)
    m = 'g';
elseif (strcmpi(s,'n')==1)
    m = 'x';
else 
    m = [];
end

end
