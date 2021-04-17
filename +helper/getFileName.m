function filename=getFileName(picIND, DirName, verbose)
% function x=getFileName(picIND)
% Created: M. Heinz 19Mar2004
% For: CNexps (GE/MH)
%
% Returns filename for given picture number
if ~exist('verbose', 'var')
    verbose= 1;
end

if ~exist('DirName', 'var')
    DirName= ['.' filesep];
else % a directory is passed
    if ~strcmp(DirName(end), filesep)
        DirName= [DirName filesep];
    end
end

d=dir(strcat(DirName, '*',sprintf('p%04d',picIND),'*'));
if length(d)>1
    warning('More than 1 file with this picture number');
    filename='';
elseif isempty(d)
    if verbose
        warning('Picture does not exists');
    end
    filename='';
else
    filename=d.name;
end

return;