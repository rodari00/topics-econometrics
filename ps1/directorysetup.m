% ----------------------------------------------------------------------- %
% Function directorysetup(root) creates a set of child directories to     %
% instantiate a local project                                             %
%                                                                         %
%   Input parameters:                                                     %
%       - root:        Root directory                          
% ----------------------------------------------------------------------- %
%   Version: 1.1                                                          %
%   Author:  Federico Rodari                                              %
%   Date:    22/09/2022                                                   %
%   E-mail:  rodari (at) bc (dot) edu                                     %
% ----------------------------------------------------------------------- %

function [f,d,t] = directorysetup(root,sys)

if strcmp(sys,'Win')
    f = strcat(root,'\figures');
    d = strcat(root,'\data');
    t = strcat(root,'\tables');
else
f = strcat(root,'/figures');
d = strcat(root,'/data');
t = strcat(root,'/tables');
end

if ~exist(f, 'dir')
       mkdir(f)
       %printf('\nCreated folder at %s',f,'...\n')
end

if ~exist(d, 'dir')
       mkdir(d)
end

if ~exist(t, 'dir')
       mkdir(t)

end





end