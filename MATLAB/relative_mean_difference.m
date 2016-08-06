function[relmeandiff]=relative_mean_difference(in)

% This file is part of the Timbre Space Analyzer & Mapper (TSAM)
% 
% The TSAM can be obtained at http://stefanofasciani.com/tsam.html
% TSAM Copyright (C) 2016 Stefano Fasciani, University of Wollongong
% Inquiries: stefanofasciani@stefanofasciani.com
% 
% The TSAM is free software: you can redistribute it and/or modify it under the terms 
% of the GNU Lesser General Public License as published by the Free Software Foundation, 
% either version 3 of the License, or (at your option) any later version.
% 
% The TSAM is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
% without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
% See the GNU Less General Public License for more details.
% 
% You should have received a copy of the GNU Lesser General Public License along with TSAM. 
% If not, see <http://www.gnu.org/licenses/>.
% 
% If you use the TSAM or any part of it in any program or publication, please acknowledge 
% its authors by adding a reference to this pubblication:
% 
% S. Fasciani, 2016, "TSAM: a tool for analyzing, modeling, and mapping the timbre of sound 
% synthesizers" in proceedings of the 13th Sound and Music Computing Conference, Hamburg, Germany.


relmeandiff=zeros(size(in,1),1);
n=size(in,2);

for i=1:n
   for j=1:n
        relmeandiff=relmeandiff+(abs(in(:,i)-in(:,j)));
   end
end

relmeandiff=(relmeandiff./(n*(n-1).*abs(mean(in,2))));






