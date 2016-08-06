function [period_hz] = period_detection(window,hop,d,sr)

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


period_hz=zeros(1,size(d,2));

for i=1:size(d,2)
       
        temp=d(:,i);
        temp=temp-max(temp);
        fil_out=filter(((window)/sr), [1 ((window)/sr)-1], temp);
        autocx=xcorr(fil_out,size(d,1),'coeff');
        [pk_loc,pk_val]=findpeaks(autocx(1:size(d,1)));
        pk_val(pk_val==1)=0;
        [~,idx]=max(pk_val);
        pk_loc=pk_loc(idx);
        if ~isempty(pk_loc)
            period_smp=((size(d,1)-pk_loc+2)*hop);
            period_time=(period_smp/sr);
            period_hz(i)=(1.0/period_time);
        else
            period_hz(i)=0;
        end

end










