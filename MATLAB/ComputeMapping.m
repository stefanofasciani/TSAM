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



function [ret]=ComputeMapping(MAXoscclient,MATLABoscserver,args)


fignum=1;

old=0;

outfile=sprintf('%s%sMap.mat',args{8},args{7});
featfile=sprintf('%s%sFeat.txt',args{8},args{7});
paramfile=sprintf('%s%sParams.txt',args{8},args{7});
if exist(outfile, 'file')
    delete(outfile);
end

oscsend(MAXoscclient,'/status/cmpt/percentage','f',1);

D=load(featfile);
I=load(paramfile);

if old
   I=I';
end

oscsend(MAXoscclient,'/status/cmpt/percentage','f',15);

num_params=std(I);
num_params=size((num_params(num_params~=0)),2);

num_combs=size(I,1);

analysis_mode=args{1};
range_enabled=args{2};
period_enabled=args{3}(1);
winsize=args{3}(2);
hopsize=args{3}(3);
dimensionality=args{4};
mergegap=args{5};
vectperstate=args{6};
mask=args{10};
actfuncid=args{11};
actfunclist={'sig','sin','hardlim','tribas','radbas'};
actfuncname=actfunclist{actfuncid};
dimredfunct=args{12
sr=args{14};

switch args{13}
    case 0
        featselect=args{9
        featselect=featselect(featselect~=0);
    case 1
        featselect=1:108;
    case 2
        featselect=1:35;
    case 3
        featselect=3:18;
    case 4
        featselect=19:31;
    case 5
        featselect=[19:31 85:108];
    case 6
        featselect=36:47;
    case 7
        featselect=48:60;
    case 8
        featselect=61:84;
    case 9
        featselect=85:108;
    case 10
        featselect=args{15}+1;
    otherwise
        featselect=1:108;
end

%apply mask
D=D.*repmat(mask,size(D,1),1);

D=D(:,featselect);

if ((analysis_mode==0)||(analysis_mode==2))
    
    cnt=1;
    average=[];
    range=[];
    period=[];
    D_post=[];
    if (analysis_mode==2)
        I=I(1:vectperstate:(size(D,1)),:);
    end
    %COMPUTING MEAN RANGE AND PERIOD FOR EACH SUB MATRIX PER STATE
    if vectperstate > 2
        for i=1:vectperstate:(size(D,1))
            temp=D(i:i+vectperstate-1,:);
            average=[average ; mean(temp)];
            if range_enabled
                range=[range ; (max(temp)-min(temp))];
            end
            if period_enabled
                period=[period ; period_detection(winsize,hopsize,temp,sr)];
            end
            cnt=cnt+1;
        end
    
    
        period_post=[];
        if period_enabled
            %removing zeros in period
            if numel(period(period==0)) > (numel(period)/2)
                for i=1:num_combs
                    temp=period(i,:);
                    period_post(i,1)=median(temp(temp~=0));
                end  
            else
                for i=1:num_combs
                    temp=period(i,:);
                    temp(temp==0)=median(temp(temp~=0));
                    period_post(i,:)=temp;
                end       
            end
            %fprintf('feat period detected min = %f max = %f\n',min(min(period_post)),max(max(period_post)));
        end
        if (period_enabled==1 || range_enabled==1)
            average=average./max(max(abs(average)));
            range=range./max(max(abs(range)));
            period_post=period_post./max(max(abs(period_post)));
        end
        D_post=[average range period_post];
    
    else
        
        D_post=D;
        
    end
    
    D_post(:,all(isnan(D_post),1))=[];
    D_post(isnan(D_post))=0;
    D_post(:,find(sum(abs(D_post))==0))=[];
    
    
elseif ((analysis_mode==1)||(analysis_mode==3))
    
    cnt=1;
    D_post=[];
    if (analysis_mode==3)
        I=I(1:vectperstate:(size(D,1)),:);
    end
    %PUTTING TEMPORAL ENVELOPE OF FEATURES IN A SINGLE VECTOR
    for i=1:vectperstate:(size(D,1))
        temp=[];
        for j=1:vectperstate
            temp=[temp D(i+j-1,:)];
        end
        D_post(cnt,:)=temp;
        cnt=cnt+1;
    end
    D_post(:,find(sum(abs(D_post))==0))=[];

elseif (analysis_mode==4)
    
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

    for i=1:size(Iuniq,1)
        idx=find(ismember(I,Iuniq(i,:),'rows'));
        temp=D(idx,:);
        if size(temp,1)>1
            average=[average ; mean(temp)];
            if range_enabled
                range=[range ; (max(temp)-min(temp))];
            end
        else
            average=[average ; temp];
            if range_enabled
                fprintf('Warning: range on single descriptors vector\n'); 
                range=[range ; temp-temp];
            end
        end
    end
    
    if (range_enabled==1)
        average=average./max(max(abs(average)));
        range=range./max(max(abs(range)));
    end
    
    D_post=[average range];
    D_post(:,all(isnan(D_post),1))=[];
    D_post(isnan(D_post))=0;
    D_post(:,find(sum(abs(D_post))==0))=[];

    I=Iuniq;
    
end

num_combs=size(I,1);

if num_combs<8
    fprintf('number of parameter combination is lower than 8');
    oscsend(MAXoscclient,'/status/cmpt/percentage','f',-1);
    return
end

num_feat=size(D_post,2);

oscsend(MAXoscclient,'/status/cmpt/percentage','f',25);

%COMPUTING AND REMOVING THE MEAN FROM EACH FEATURE
mean_vect=mean(D_post);

for i=1:num_combs
    D_post(i,:)=D_post(i,:)-mean_vect;
end

D_post=D_post./(max(max(D_post)));

if dimensionality<num_feat

    if isequal(dimredfunct,'Isomap')
        D_star=1;
        iso_div=20;
        while size(D_post,1)~=size(D_star,1)
            [D_star,~]=compute_mapping(D_post,'Isomap',dimensionality,ceil(size(D_post,1)/iso_div));%PCA, LLE, Laplacian, LaplacianEigenmaps, Isomap
            iso_div=iso_div-2;
            if iso_div<1
                disp('problem in Isomap Dimensionality Reduction, using PCA');
                [D_star,~]=compute_mapping(D_post,'PCA',dimensionality);
                break;
            end
        end
    else
        [D_star,~]=compute_mapping(D_post,dimredfunct,dimensionality);
        if size(D_post,1)~=size(D_star,1)
                disp('problem with selected Dimensionality Reduction technique, using PCA');
                [D_star,~]=compute_mapping(D_post,'PCA',dimensionality);
        end   
    end
    [D_star,~]=compute_mapping(D_star,'PCA',dimensionality);%PCA, LLE, Laplacian, LaplacianEigenmaps, Isomap	

else
    
    D_star=D_post;
    
end

%TO UNIFORM DISTRIBUTION - RANK TRANSFORM
for i=1:num_combs
    for j=1:dimensionality
        ranktransform(i,j)=numel(D_star(D_star(:,j)<D_star(i,j)));
    end
end

oscsend(MAXoscclient,'/status/cmpt/percentage','f',45);

%NORMALIZATION TO RANGE 0 1 and scale down
ranktransform=ranktransform./max(max(max(max(max(ranktransform)))));
ranktransform=0.1+ranktransform.*0.8;


switch dimensionality
    
    case 1
        distmeshout=ranktransform;
        
    case 2
        fd=@(p) drectangle2D_sfa(p,0,1,0,1);
        [distmeshout,ttt]=sfa_distmeshnd_sfast(fd,@huniform_sfa,(1/num_combs)^(1/dimensionality)*0.1,[0,0;1,1],ranktransform);
    case 3
        fd=@(p) drectangle3D_sfa(p,0,1,0,1,0,1);
        [distmeshout,ttt]=sfa_distmeshnd_sfast(fd,@huniform_sfa,(1/num_combs)^(1/dimensionality)*0.1,[0,0,0;1,1,1],ranktransform);

end

oscsend(MAXoscclient,'/status/cmpt/percentage','f',75);

hiddenneurons=10;
accuracy=1;
while ((hiddenneurons<201)&&(accuracy>0.01))
    net=elmMOD(distmeshout',D_star',hiddenneurons,actfuncname);
    net.accuracy=net.accuracy/sqrt(size(D_star,1));
    hiddenneurons=hiddenneurons+10;
    accuracy=net.accuracy;
end


oscsend(MAXoscclient,'/status/cmpt/percentage','f',90);


%MAX AND MIN FOR DIRECT D STAR NAVIGATION
for i=1:size(D_star,2)
    D_star_minmax(1,i)=min(D_star(:,i));
    D_star_minmax(2,i)=max(D_star(:,i));
end

for i=1:size(D_star_minmax,2) 
    D_star_map_min(i)=D_star_minmax(1,i);
    D_star_map_range(i)=D_star_minmax(2,i)-D_star_minmax(1,i);
end

D_star_min_range=[D_star_map_min;D_star_map_range];


outfile=sprintf('%s%sMap-DS.txt',args{8},args{7});
save(outfile,'D_star','-ASCII','-double');

outfile=sprintf('%s%sMap-DU.txt',args{8},args{7});
save(outfile,'distmeshout','-ASCII','-double');

outfile=sprintf('%s%sMap-I.txt',args{8},args{7});
save(outfile,'I','-ASCII','-double');

outfile=sprintf('%s%sMap-Dmnrg.txt',args{8},args{7});
save(outfile,'D_star_min_range','-ASCII','-double');

%minIdist maxIdist accuracy
info=[max(pdist((I))) min(pdist((I)))+0.1*min(pdist((I))) net.accuracy];
outfile=sprintf('%s%sMap-Info.txt',args{8},args{7});
save(outfile,'info','-ASCII','-double');

temp=net.inweight;
outfile=sprintf('%s%sMap-Iw.txt',args{8},args{7});
save(outfile,'temp','-ASCII','-double');

temp=net.outweight;
outfile=sprintf('%s%sMap-Ow.txt',args{8},args{7});
save(outfile,'temp','-ASCII','-double');

temp=net.bias;
outfile=sprintf('%s%sMap-Bs.txt',args{8},args{7});
save(outfile,'temp','-ASCII','-double');

temp=find((strcmp(net.actfunct,actfunclist))==1);
outfile=sprintf('%s%sMap-Act.txt',args{8},args{7});
save(outfile,'temp','-ASCII','-double');

oscsend(MAXoscclient,'/status/cmpt/percentage','f',100);

return




%EXTRA FUNCTIONS

function h=huniform_sfa(p,varargin)

%   Copyright (C) 2004-2012 Per-Olof Persson. See COPYRIGHT.TXT for details.

% The original file was modified for the integration with the
% Timbre Space Mapping for VST Synthesizers (TSM4VSTS).

h=ones(size(p,1),1);

function d=drectangle2D_sfa(p,x1,x2,y1,y2)

%   Copyright (C) 2004-2012 Per-Olof Persson. See COPYRIGHT.TXT for details.

% The original file was modified for the integration with the
% Timbre Space Mapping for VST Synthesizers (TSM4VSTS).

d=-min(min(min(-y1+p(:,2),y2-p(:,2)),-x1+p(:,1)),x2-p(:,1));

function d=drectangle3D_sfa(p,x1,x2,y1,y2,z1,z2)

%   Copyright (C) 2004-2012 Per-Olof Persson. See COPYRIGHT.TXT for details.

% The original file was modified for the integration with the
% Timbre Space Mapping for VST Synthesizers (TSM4VSTS).

d=-min(min(min(min(min(-z1+p(:,3),z2-p(:,3)),-y1+p(:,2)),y2-p(:,2)),-x1+p(:,1)),x2-p(:,1));
