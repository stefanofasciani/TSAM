function oscsend(u,path,varargin)
% Sends a Open Sound Control (OSC) message through a UDP connection
%
% oscsend(u,path)
% oscsend(u,path,types,arg1,arg2,...)
% oscsedn(u,path,types,[args])
%
% u = UDP object with open connection.
% path = path-string
% types = string with types of arguments,
%    supported:
%       i = integer
%       f = float
%       s = string
%       N = Null (ignores corresponding argument)
%       I = Impulse (ignores corresponding argument)
%       T = True (ignores corresponding argument)
%       F = False (ignores corresponding argument)
%       B = boolean (not official: converts argument to T/F in the type)
%    not supported:
%       b = blob
%
% args = arguments as specified by types.
%
% EXAMPLE
%       u = udp('127.0.0.1',7488);  
%       fopen(u);
%       oscsend(u,'/test','ifsINBTF', 1, 3.14, 'hello',[],[],false,[],[]);
%       fclose(u);
%
% See http://opensoundcontrol.org/ for more information about OSC.

% MARK MARIJNISSEN 10 may 2011 (markmarijnissen@gmail.com)



% The original file was modified for the integration with the
% Timbre Space Analyzer and Mapper (TSAM).
% 
% The TSAM is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% The TSAM is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU Less General Public License for more details.
% 
% You should have received a copy of the GNU Lesser General Public License
% along with TSAM.  If not, see <http://www.gnu.org/licenses/>.
%
% The TSAM can be obtained at http://stefanofasciani.com/tsam.html
% TSAM Copyright (C) 2016 Stefano Fasciani, University of Wollongong
% Inquiries: stefanofasciani@stefanofasciani.com


    
    %figure out little endian for int/float conversion
    [~, ~, endian] = computer;
    littleEndian = endian == 'L';

    % set type
    if nargin >= 3,
        types = oscstr([',' varargin{1}]);
    else
        types = oscstr(',');
    end;
    
    % set args (either a matrix, or varargin)
    if nargin == 3 && length(types) > 2
        args = varargin{2};
    else
        args = varargin(2:end);
    end;

    % convert arguments to the right bytes
    data = [];
    for i=1:length(args)
        switch(types(i+1))
            case 'i'
                data = [data oscint(args{i},littleEndian)];
            case 'f'
                data = [data oscfloat(args{i},littleEndian)];
            case 's'
                data = [data oscstr(args{i})];
            case 'B'
                if args{i}
                    types(i+1) = 'T';
                else
                    types(i+1) = 'F';
                end;
            case {'N','I','T','F'}
                %ignore data
            otherwise
                warning(['Unsupported type: ' types(i+1)]);
        end;
    end;
    
    %write data to UDP
    data = [oscstr(path) types data];
    fwrite(u,data);
end

%Conversion from double to float
function float = oscfloat(float,littleEndian)
   if littleEndian
        float = typecast(swapbytes(single(float)),'uint8');
   else
        float = typecast(single(float),'uint8');
   end;
end

%Conversion to int
function int = oscint(int,littleEndian)
   if littleEndian
        int = typecast(swapbytes(int32(int)),'uint8');
   else
        int = typecast(int32(int),'uint8');
   end;
end

%Conversion to string (null-terminated, in multiples of 4 bytes)
function string = oscstr(string)
    string = [string 0 0 0 0];
    string = string(1:end-mod(length(string),4));
end