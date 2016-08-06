function [p,t,maxdp]=distmeshnd_sfa(moviefile,fdist,fh,h,box,fix,varargin)
%DISTMESHND N-D Mesh Generator using Distance Functions.
%   [P,T]=DISTMESHND(FDIST,FH,H,BOX,FIX,FDISTPARAMS)
%
%      P:           Node positions (NxNDIM)
%      T:           Triangle indices (NTx(NDIM+1))
%      FDIST:       Distance function
%      FH:          Edge length function
%      H:           Smallest edge length
%      BOX:         Bounding box [xmin,xmax;ymin,ymax; ...] (NDIMx2)
%      FIX:         Fixed node positions (NFIXxNDIM)
%      FDISTPARAMS: Additional parameters passed to FDIST
%
%   Example: Unit ball
%      dim=3;
%      d=inline('sqrt(sum(p.^2,2))-1','p');
%      [p,t]=distmeshnd(d,@huniform,0.2,[-ones(1,dim);ones(1,dim)],[]);
%
%   See also: DISTMESH2D, DELAUNAYN, TRIMESH, MESHDEMOND.

%   Copyright (C) 2004-2012 Per-Olof Persson. See COPYRIGHT.TXT for details.


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




dim=size(box,2);
ptol=.01; ttol=.1; L0mult=1+.4/2^(dim-1); deltat=.1; geps=1e-1*h; deps=sqrt(eps)*h;
min_term=1000000000000000;
term_check=0;
framecount=1;
% 1. Create initial distribution in bounding box
if dim==1
  p=(box(1):h:box(2))';
else
  cbox=cell(1,dim);
  for ii=1:dim
    cbox{ii}=box(1,ii):h:box(2,ii);
  end
  pp=cell(1,dim);
  [pp{:}]=ndgrid(cbox{:});
  p=zeros(prod(size(pp{1})),dim);
  for ii=1:dim
    p(:,ii)=pp{ii}(:);
  end
end

% 2. Remove points outside the region, apply the rejection method
p=p(feval(fdist,p,varargin{:})<geps,:);
r0=feval(fh,p);
p=[fix];
N=size(p,1);

count=0;
p0=inf;
while 1
  % 3. Retriangulation by Delaunay
  if max(sqrt(sum((p-p0).^2,2)))>ttol*h
    p0=p;
    t=delaunayn(p);
    pmid=zeros(size(t,1),dim);
    for ii=1:dim+1
      pmid=pmid+p(t(:,ii),:)/(dim+1);
    end
    t=t(feval(fdist,pmid,varargin{:})<-geps,:);
    % 4. Describe each edge by a unique pair of nodes
    pair=zeros(0,2);
    localpairs=nchoosek(1:dim+1,2);
    for ii=1:size(localpairs,1)
      pair=[pair;t(:,localpairs(ii,:))];
    end
    pair=unique(sort(pair,2),'rows');
    % 5. Graphical output of the current mesh
    if dim==2
     if mod(count,10)==0   
      cla,patch('vertices',p,'faces',t,'edgecol','k','facecol',[0.9,0.9,0.9]);
      hold on
      scatter(p(:,1),p(:,2),10,[1 0.1 0.4]);
      hold off
      axis on
      axis([-0.1 1.1 -0.1 1.1]);
      drawnow
      movie(framecount)=getframe(gcf);
      framecount=framecount+1;
     end 
    elseif dim==3
      if mod(count,20)==0
        simpplot(p,t,'p(:,2)>0');
        title(['Retriangulation #',int2str(count)])
        hold on
        scatter3(p(:,1),p(:,2),p(:,3),10,[1 0.1 0.4]);
        hold off
        axis on
        axis([-0.1 1.1 -0.1 1.1 -0.1 1.1]);
        drawnow
        movie(framecount)=getframe(gcf);
        framecount=framecount+1;
      end
    else
      disp(sprintf('Retriangulation #%d',count))
    end
    count=count+1;
  end

  % 6. Move mesh points based on edge lengths L and forces F
  bars=p(pair(:,1),:)-p(pair(:,2),:);
  L=sqrt(sum(bars.^2,2));
  L0=feval(fh,(p(pair(:,1),:)+p(pair(:,2),:))/2);
  L0=L0*L0mult*(sum(L.^dim)/sum(L0.^dim))^(1/dim);
  F=max(L0-L,0);
  Fbar=[bars,-bars].*repmat(F./L,1,2*dim);
  dp=full(sparse(pair(:,[ones(1,dim),2*ones(1,dim)]), ...
                 ones(size(pair,1),1)*[1:dim,1:dim], ...
                 Fbar,N,dim));
  p=p+deltat*dp;

  % 7. Bring outside points back to the boundary
  d=feval(fdist,p,varargin{:}); ix=d>0;
  gradd=zeros(sum(ix),dim);
  for ii=1:dim
    a=zeros(1,dim);
    a(ii)=deps;
    d1x=feval(fdist,p(ix,:)+ones(sum(ix),1)*a,varargin{:});
    gradd(:,ii)=(d1x-d(ix))/deps;
  end
  p(ix,:)=p(ix,:)-d(ix)*ones(1,dim).*gradd;

  % 8. Termination criterion
  maxdp=max(deltat*sqrt(sum(dp(d<-geps,:).^2,2)));
  treshold=ptol*h;
  
  if maxdp<treshold
      
      break; 
  else
      
      min_term_old=min_term;
      min_term=min(min_term_old,maxdp);
      
      if min_term==min_term_old
        term_check=term_check+1;
      else
        term_check=0;  
      end

      if term_check>500
         break; 
      end
      
  end
  
end

%%save movie to avi file
%movie2avi(movie, moviefile, 'compression', 'None','fps',10);

