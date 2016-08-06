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


%global variables
matlab_osc_port=9001;
max_osc_port=9002;
max_ip_address='127.0.0.1';
quit=0;
function_id=0;

fprintf('Starting TSAM Engine\n'); 


%OSC client communication initialization
MAXoscclient=udp(max_ip_address,max_osc_port);
fopen(MAXoscclient);
oscsend(MAXoscclient,'/status/run','i',1);

%OSC server communication initialization
MATLABoscserver=osc_new_server(matlab_osc_port);

%{1}anamode {2}range {3}period {4}dim {5}merge {6}name {7}path
cmpt_args=[];
rt_args=[];

a=[];
while quit==0
   
    
    
    %receive OSC message
    oscrm=osc_recv(MATLABoscserver); 
    %check message validity
    if ~isempty(oscrm)
        % process all messages in cue
        for m=1:length(oscrm)
            %route OSC data
            oscpath=oscrm{m}.path;
            oscdata=cell2mat(oscrm{m}.data);  
            if (isequal(oscpath,'/cmd/quit'))
                quit=oscdata;
            elseif (isequal(oscpath,'/cmpt/funct'))
                function_id=oscdata; %id=1 mapping; id=2 score
            elseif (isequal(oscpath,'/cmpt/anamode'))
                cmpt_args{1}=oscdata;
            elseif (isequal(oscpath,'/cmpt/range'))
                cmpt_args{2}=oscdata;
            elseif (isequal(oscpath,'/cmpt/period'))
                cmpt_args{3}=oscdata;
            elseif (isequal(oscpath,'/cmpt/dim'))
                cmpt_args{4}=oscdata;
            elseif (isequal(oscpath,'/cmpt/merge'))
                cmpt_args{5}=oscdata;
            elseif (isequal(oscpath,'/cmpt/vcprst'))
                cmpt_args{6}=oscdata;
            elseif (isequal(oscpath,'/cmpt/name'))
                cmpt_args{7}=oscdata;
            elseif (isequal(oscpath,'/cmpt/path'))
                cmpt_args{8}=oscdata;
            elseif (isequal(oscpath,'/cmpt/manualdesc'))
                cmpt_args{9}=oscdata;
            elseif (isequal(oscpath,'/cmpt/mskvect'))
                cmpt_args{10}=oscdata;
            elseif (isequal(oscpath,'/cmpt/actfunc'))
                cmpt_args{11}=oscdata;
            elseif (isequal(oscpath,'/cmpt/dimred'))
                cmpt_args{12}=oscdata;
            elseif (isequal(oscpath,'/cmpt/descsel'))
                cmpt_args{13}=oscdata;
            elseif (isequal(oscpath,'/cmpt/sr'))
                cmpt_args{14}=oscdata;
            elseif (isequal(oscpath,'/cmpt/topdesc'))
                cmpt_args{15}=oscdata;  
            elseif (isequal(oscpath,'/cmpt/run'))
                if function_id==1
                    ComputeMapping(MAXoscclient,MATLABoscserver,cmpt_args);
                    oscsend(MAXoscclient,'/status/cmpt/run','i',0);
                elseif function_id==2
                    ComputeScore(MAXoscclient,MATLABoscserver,cmpt_args);
                    oscsend(MAXoscclient,'/status/cmpt/run','i',0);
                else
                    fprintf('TSAM-MAIN function id = %d not valid\n',function_id);
                end
            else
                fprintf('TSAM-MAIN message not valid: path = %s; data = %f \n',oscpath,oscdata);
            end      
        end
    end
    
    
    
    
end


%OSC server communication termination
osc_free_server(MATLABoscserver);
%OSC client communication termination
oscsend(MAXoscclient,'/status/run','i',0);
fclose(MAXoscclient);


fprintf('Terminating TSAM Engine\n'); 

