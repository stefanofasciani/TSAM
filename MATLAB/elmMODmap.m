function [out] = elmMODmap(in,struct)

    %%%%    Authors:    MR QIN-YU ZHU AND DR GUANG-BIN HUANG
    %%%%    NANYANG TECHNOLOGICAL UNIVERSITY, SINGAPORE
    %%%%    EMAIL:      EGBHUANG@NTU.EDU.SG; GBHUANG@IEEE.ORG
    %%%%    WEBSITE:    http://www.ntu.edu.sg/eee/icis/cv/egbhuang.htm
    %%%%    DATE:       APRIL 2004


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


tempH=struct.inweight*in';                    %   Extend the bias matrix BiasofHiddenNeurons to match the demention of H
tempH=tempH + struct.bias;
switch lower(struct.actfunct)
    case {'sig','sigmoid'}
        %%%%%%%% Sigmoid 
        H_test = 1 ./ (1 + exp(-tempH));
    case {'sin','sine'}
        %%%%%%%% Sine
        H_test = sin(tempH);        
    case {'hardlim'}
        %%%%%%%% Hard Limit
        %H_test = hardlim(tempH);
        H_test = double(tempH >= 0);
    case {'tribas'}
        %%%%%%%% Triangular basis function
        %H_test = tribas(tempH);
        H_test = max(0,1-abs(tempH));
    case {'radbas'}
        %%%%%%%% Radial basis function
        %H_test = radbas(tempH);
        H_test = exp(-(tempH.*tempH));
        %%%%%%%% More activation functions can be added here        
end
out=(H_test' * struct.outweight)';                       %   TY: the actual output of the testing data
