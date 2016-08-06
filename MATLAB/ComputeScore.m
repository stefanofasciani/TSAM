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



function [ret]=ComputeScore(MAXoscclient,MATLABoscserver,args)

fignum=1;

old=0;

featfile=sprintf('%s%sFeat.txt',args{8},args{7});
paramfile=sprintf('%s%sParams.txt',args{8},args{7});


oscsend(MAXoscclient,'/status/cmpt/percentage','f',1);

D=load(featfile);
I=load(paramfile);

if old
   I=I';
end

oscsend(MAXoscclient,'/status/cmpt/percentage','f',15);

num_params=std(I);
num_params=size((num_params(num_params~=0)),2);
num_desc=size(D,2);

analysis_mode=args{1};
range_enabled=1;
winsize=args{3}(2);
hopsize=args{3}(3);
dimensionality=args{4};
mergegap=args{5};
vectperstate=args{6};
mask=args{10};
sr=args{14};

%apply mask
D=D.*repmat(mask,size(D,1),1);

param_score=zeros(1,10);
corrmat=zeros(10,108*2);

if ((analysis_mode==0)||(analysis_mode==2))
    
    noisiness_sustain=zeros(num_desc,1);
    variance_sustain=zeros(num_desc,1);
    independence_sustain=zeros(num_desc,1);
    paramcorrelation_sustain=zeros(num_desc,1);

    variance_sustain_range=zeros(num_desc,1);
    independence_sustain_range=zeros(num_desc,1);
    paramcorrelation_sustain_range=zeros(num_desc,1);
    
    num_combs=size(I,1);
    
    cnt=1;
    average=[];
    range=[];
    D_post=[];
    if (analysis_mode==2)
        I=I(1:vectperstate:(size(D,1)),:);
    end
    
    %COMPUTING MEAN AND RANGE FOR EACH SUB MATRIX PER STATE
    if vectperstate > 2
        for i=1:vectperstate:(size(D,1))
            temp=D(i:i+vectperstate-1,:);
            noisiness_sustain=noisiness_sustain+relative_mean_difference(temp');
            average=[average ; mean(temp)];
            if range_enabled
                range=[range ; (max(temp)-min(temp))];
            end
        end
        
        noisiness_sustain=noisiness_sustain./num_combs;
    
        oscsend(MAXoscclient,'/status/cmpt/percentage','f',65);
        
        D_post=[average range];
        D_post(:,all(isnan(D_post),1))=[];
        D_post(isnan(D_post))=0;
        variance_sustain=(relative_mean_difference((D_post(:,1:num_desc))'));
        variance_sustain_range=(relative_mean_difference((D_post(:,1+num_desc:2*num_desc))'));
    
    else

        D_post=D;
        D_post(:,all(isnan(D_post),1))=[];
        D_post(isnan(D_post))=0;
        variance_sustain=(relative_mean_difference(D_post'));
        
    end
    
    oscsend(MAXoscclient,'/status/cmpt/percentage','f',70);
    
    %INDEPENDENCE
    %on average
    temp1=D_post(:,1:num_desc);
    independence_sustain=((nansum(abs(corr(temp1,temp1)))-1)/(num_desc-1))';
    independence_sustain=1-independence_sustain;
    
    oscsend(MAXoscclient,'/status/cmpt/percentage','f',75);
    
    %on on range
    temp1=D_post(:,1+num_desc:2*num_desc);
    independence_sustain_range=((nansum(abs(corr(temp1,temp1)))-1)/(num_desc-1))';
    independence_sustain_range=1-independence_sustain_range;
    
    
    oscsend(MAXoscclient,'/status/cmpt/percentage','f',80);    
    
    %CORRELATION WITH PARAM
    %on average
    temp1=D_post(:,1:num_desc);
    for i=1:(size(temp1,2))
        cnt=0;
        temp4=0;
        temp2=temp1(:,i);
        for j=1:size(I,2)
            if sum(I(:,j))~=0
                temp3=I(:,j);
                temp5=nansum(abs(corr(temp2,temp3)));
                paramcorrelation_sustain(i)=paramcorrelation_sustain(i)+temp5;
                param_score(j)=param_score(j)+temp5;
                corrmat(j,i)=temp5;
                cnt=cnt+1;
            end
        end
        paramcorrelation_sustain(i)=paramcorrelation_sustain(i)/cnt;
    end
    
    oscsend(MAXoscclient,'/status/cmpt/percentage','f',85);
    
    %on on range
    temp1=D_post(:,1+num_desc:2*num_desc);
    for i=1:(size(temp1,2))
        cnt=0;
        temp4=0;
        temp2=temp1(:,i);
        for j=1:size(I,2)
            if sum(I(:,j))~=0
                temp3=I(:,j);
                temp5=nansum(abs(corr(temp2,temp3)));
                paramcorrelation_sustain_range(i)=paramcorrelation_sustain_range(i)+temp5;
                param_score(j)=param_score(j)+temp5;
                corrmat(j,i+108)=temp5;
                cnt=cnt+1;
            end
        end
        paramcorrelation_sustain_range(i)=paramcorrelation_sustain_range(i)/cnt;
    end
    
    oscsend(MAXoscclient,'/status/cmpt/percentage','f',95);

    outmat=[noisiness_sustain' ; variance_sustain' ; independence_sustain' ; paramcorrelation_sustain' ; zeros(1,num_desc) ; variance_sustain_range' ; independence_sustain_range' ; paramcorrelation_sustain_range'];
    outfile=sprintf('%s%sScoreDsc.txt',args{8},args{7});
    if exist(outfile, 'file')
        delete(outfile);
    end
    save(outfile,'outmat','-ASCII','-double');
    
    outmat=param_score./(abs(max(param_score)));
    outfile=sprintf('%s%sScorePrm.txt',args{8},args{7});
    if exist(outfile, 'file')
        delete(outfile);
    end
    save(outfile,'outmat','-ASCII','-double');
    
    outmat=corrmat;
    outfile=sprintf('%s%sCorrMat.txt',args{8},args{7});
    if exist(outfile, 'file')
        delete(outfile);
    end
    save(outfile,'outmat','-ASCII','-double');
    
    oscsend(MAXoscclient,'/score','i',1);
    
elseif ((analysis_mode==1)||(analysis_mode==3))
    
    noisiness_envelope=zeros(num_desc,1);
    variance_envelope=zeros(num_desc,1);
    independence_envelope=zeros(num_desc,1);
    paramcorrelation_envelope=zeros(num_desc,1);

    num_combs=size(I,1);
    
    cnt=1;
    D_post=[];
    range=[];
    if (analysis_mode==3)
        I=I(1:vectperstate:(size(D,1)),:);
    end
    %PUTTING TEMPORAL ENVELOPE OF FEATURES IN A SINGLE VECTOR
    
    for i=1:vectperstate:(size(D,1))
        temp=D(i:i+vectperstate-1,:);
        for j=1:num_desc
            noisiness_envelope(j)=noisiness_envelope(j)+sum(abs(diff(diff(temp(:,j))>0)))/(length(diff(temp(:,j)))-1);
            variance_envelope(j)=variance_envelope(j)+relative_mean_difference(temp(:,j)');
        end
        independence_envelope=independence_envelope+((nansum(abs(corr(temp,temp)))-1)/(num_desc-1))';
    end
    
    noisiness_envelope=noisiness_envelope./num_combs;
    variance_envelope=variance_envelope./num_combs;
    independence_envelope=independence_envelope./num_combs;
    independence_envelope=1-independence_envelope;
    
    oscsend(MAXoscclient,'/status/cmpt/percentage','f',80);
    
    for i=1:vectperstate:(size(D,1))
        temp=[];
        for j=1:vectperstate
            temp=[temp D(i+j-1,:)];
        end
        D_post(cnt,:)=temp;
        cnt=cnt+1;
    end
    %D_post(:,find(sum(abs(D_post))==0))=[];

    oscsend(MAXoscclient,'/status/cmpt/percentage','f',85);
    
    cnt=0;
    for j=1:size(I,2)
        if sum(I(:,j))~=0
            temp1=I(:,j);
            for k=1:num_desc
                temp2=D_post(:,k:num_desc:num_desc*vectperstate);
                temp5=(nansum(abs(corr(temp1,temp2)))/vectperstate);
                paramcorrelation_envelope(k)=paramcorrelation_envelope(k)+temp5;
                param_score(j)=param_score(j)+temp5;
                corrmat(j,k)=temp5;
            end
            cnt=cnt+1;
        end
    end
    paramcorrelation_envelope=paramcorrelation_envelope/cnt;
    
    oscsend(MAXoscclient,'/status/cmpt/percentage','f',95);
    
    outmat=[noisiness_envelope' ; variance_envelope' ; independence_envelope' ; paramcorrelation_envelope' ; zeros(1,num_desc) ; zeros(1,num_desc) ; zeros(1,num_desc) ; zeros(1,num_desc)];
    outfile=sprintf('%s%sScoreDsc.txt',args{8},args{7});
    if exist(outfile, 'file')
        delete(outfile);
    end
    save(outfile,'outmat','-ASCII','-double');
    
    outmat=param_score./(abs(max(param_score)));
    outfile=sprintf('%s%sScorePrm.txt',args{8},args{7});
    if exist(outfile, 'file')
        delete(outfile);
    end
    save(outfile,'outmat','-ASCII','-double');
    
    outmat=corrmat;
    outfile=sprintf('%s%sCorrMat.txt',args{8},args{7});
    if exist(outfile, 'file')
        delete(outfile);
    end
    save(outfile,'outmat','-ASCII','-double');
    
    oscsend(MAXoscclient,'/score','i',1);
    
    
elseif (analysis_mode==4)
    
    noisiness_sustain=zeros(num_desc,1);
    variance_sustain=zeros(num_desc,1);
    independence_sustain=zeros(num_desc,1);
    paramcorrelation_sustain=zeros(num_desc,1);

    variance_sustain_range=zeros(num_desc,1);
    independence_sustain_range=zeros(num_desc,1);
    paramcorrelation_sustain_range=zeros(num_desc,1);
    
    cnt=0;
    average=[];
    range=[];
    D_post=[];
    
    for j=1:numel(mergegap)
        ranges=0:mergegap(j):(1+mergegap(j));
        for i=1:(numel(ranges)-1)
            I(find((I(:,j)<=ranges(i+1))&(I(:,j)>ranges(i))),j)=(ranges(i+1)+ranges(i))/2;
        end
    end

    Iuniq=unique(I,'rows');
    
    num_combs=size(I,Iuniq);

    for i=1:size(Iuniq,1)
        idx=find(ismember(I,Iuniq(i,:),'rows'));
        temp=D(idx,:);
        if size(temp,1)>1
            noisiness_sustain=noisiness_sustain+relative_mean_difference(temp');
            cnt=cnt+1;
            average=[average ; mean(temp)];
            if range_enabled
                range=[range ; (max(temp)-min(temp))];
            end
        else
            noisiness_sustain=noisiness_sustain+relative_mean_difference(temp');
            average=[average ; temp];
            if range_enabled
                fprintf('Warning: range on single descriptors vector\n'); 
                range=[range ; temp-temp];
            end
        end
    end
    
    noisiness_sustain=noisiness_sustain./cnt;

    D_post=[average range];
    D_post(:,all(isnan(D_post),1))=[];
    D_post(isnan(D_post))=0;
    %D_post(:,find(sum(abs(D_post))==0))=[];
    
    variance_sustain=(relative_mean_difference((D_post(:,1:num_desc))'));
    variance_sustain_range=(relative_mean_difference((D_post(:,1+num_desc:2*num_desc))'));

    I=Iuniq;
    
    oscsend(MAXoscclient,'/status/cmpt/percentage','f',70);
    
    %INDEPENDENCE
    %on average
    temp1=D_post(:,1:num_desc);
    independence_sustain=((nansum(abs(corr(temp1,temp1)))-1)/(num_desc-1))';
    independence_sustain=1-independence_sustain;
    
    oscsend(MAXoscclient,'/status/cmpt/percentage','f',75);
    
    %on on range
    temp1=D_post(:,1+num_desc:2*num_desc);
    independence_sustain_range=((nansum(abs(corr(temp1,temp1)))-1)/(num_desc-1))';
    independence_sustain_range=1-independence_sustain_range;
    
    
    oscsend(MAXoscclient,'/status/cmpt/percentage','f',80);    
    
    %CORRELATION WITH PARAM
    %on average
    temp1=D_post(:,1:num_desc);
    for i=1:(size(temp1,2))
        cnt=0;
        temp4=0;
        temp2=temp1(:,i);
        for j=1:size(I,2)
            if sum(I(:,j))~=0
                temp3=I(:,j);
                temp5=nansum(abs(corr(temp2,temp3)));
                paramcorrelation_sustain(i)=paramcorrelation_sustain(i)+temp5;
                param_score(j)=param_score(j)+temp5;
                corrmat(j,i)=temp5;
                cnt=cnt+1;
            end
        end
        paramcorrelation_sustain(i)=paramcorrelation_sustain(i)/cnt;
    end
    
    oscsend(MAXoscclient,'/status/cmpt/percentage','f',85);
    
    %on on range
    temp1=D_post(:,1+num_desc:2*num_desc);
    for i=1:(size(temp1,2))
        cnt=0;
        temp4=0;
        temp2=temp1(:,i);
        for j=1:size(I,2)
            if sum(I(:,j))~=0
                temp3=I(:,j);
                temp5=nansum(abs(corr(temp2,temp3)));
                paramcorrelation_sustain_range(i)=paramcorrelation_sustain_range(i)+temp5;
                param_score(j)=param_score(j)+temp5;
                corrmat(j,i+108)=temp5;
                cnt=cnt+1;
            end
        end
        paramcorrelation_sustain_range(i)=paramcorrelation_sustain_range(i)/cnt;
    end
    
    oscsend(MAXoscclient,'/status/cmpt/percentage','f',95);
    
    outmat=[noisiness_sustain' ; variance_sustain' ; independence_sustain' ; paramcorrelation_sustain' ; zeros(1,num_desc) ; variance_sustain_range' ; independence_sustain_range' ; paramcorrelation_sustain_range'];
    outfile=sprintf('%s%sScoreDsc.txt',args{8},args{7});
    if exist(outfile, 'file')
        delete(outfile);
    end
    save(outfile,'outmat','-ASCII','-double');
    
    outmat=param_score./(abs(max(param_score)));
    outfile=sprintf('%s%sScorePrm.txt',args{8},args{7});
    if exist(outfile, 'file')
        delete(outfile);
    end
    save(outfile,'outmat','-ASCII','-double');
    
    outmat=corrmat;
    outfile=sprintf('%s%sCorrMat.txt',args{8},args{7});
    if exist(outfile, 'file')
        delete(outfile);
    end
    save(outfile,'outmat','-ASCII','-double');
    
    oscsend(MAXoscclient,'/score','i',1);
    
    
end


oscsend(MAXoscclient,'/status/cmpt/percentage','f',100);
















