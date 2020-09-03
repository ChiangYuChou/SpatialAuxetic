close all; clear all; clc;
%% loading file
[v, l, ~]=obj_input('Input.obj');

%% identifing the points above the surface
isOnGround=(abs(v(:,3))<=1e-2);
plot(v(isOnGround,1),v(isOnGround,2),'r.'); hold on;
plot(v(~isOnGround,1),v(~isOnGround,2),'b.');
axis equal
%% indentifing the ration axes
a=0;
list_RA=zeros(length(l(:,1)),1);
for i=1:length(l(:,1))
    if ~isOnGround(l(i,1)) && ~isOnGround(l(i,2))
        j1=(l(:,2)==l(i,1) & isOnGround(l(:,1)));
        l(j1,:)=l(j1,end:-1:1);
        j1=(l(:,1)==l(i,1) & isOnGround(l(:,2)));
        j2=(l(:,2)==l(i,2) & isOnGround(l(:,1)));
        l(j2,:)=l(j2,end:-1:1);
        j2=(l(:,1)==l(i,2) & isOnGround(l(:,2)));
        if sum(j1)>=3 || sum(j2)>=3
            a=a+1;
            list_RA(a)=i;
        end
    end
end
list_RA=list_RA(1:a);

for i=1:a
    plot(v(l(list_RA(i),:),1),v(l(list_RA(i),:),2),'b')
end
%% determint the direction of the rotation axes
f=list_RA(1);
cand=list_RA(2:end);
while ~isempty(cand)
    ff=[];
    for i=1:length(f)
        j1=l(cand,1)==l(f(i),1);
        j2=l(cand,2)==l(f(i),1);
        l(cand(j2),:)=l(cand(j2),end:-1:1);
        j3=l(cand,2)==l(f(i),2);
        j4=l(cand,1)==l(f(i),2);
        l(cand(j4),:)=l(cand(j4),end:-1:1);
        ff=[ff;cand(j1);cand(j2);cand(j3);cand(j4)];
        cand=cand(~j1 & ~j2 & ~j3 & ~j4);
    end
    f=ff;
end
% plot(v(l(list_RA,1),1),v(l(list_RA,1),2),'r*')
% plot(v(l(list_RA,2),1),v(l(list_RA,2),2),'ro')
%% derive the angle of gap
z=[0 0 -1];
ang_GP=zeros(size(list_RA));
for i=1:length(list_RA)
    ax=v(l(list_RA(i),2),:)-v(l(list_RA(i),1),:);
    ax=ax/sqrt(ax*ax');
    no=cross(ax,z);
    no=no/norm(no);
    
    j1=(l(:,1)==l(list_RA(i),1) & isOnGround(l(:,2)));
    cand1=v(l(j1,2),:)-ones(sum(j1),1)*v(l(list_RA(i),1),:);
    cand1=cand1./(sqrt(diag(cand1*cand1'))*ones(1,3));
    ind=cand1*ax';
    cand1=cand1(abs(ind-max(ind))<=2e-3,:);

    
    j2=(l(:,1)==l(list_RA(i),2) & isOnGround(l(:,2)));
    cand2=v(l(j2,2),:)-ones(sum(j2),1)*v(l(list_RA(i),2),:);
    cand2=cand2./(sqrt(diag(cand2*cand2'))*ones(1,3));
    ind=cand2*ax';
    cand2=cand2(abs(ind-min(ind))<=2e-3,:);
    
    a=cross(cand1(cand1*no'>0,:),cand2(cand2*no'>0,:));
    a=a/norm(a);
    a=acos(a*no');
    b=cross(cand1(cand1*no'<0,:),cand2(cand2*no'<0,:));
    b=b/norm(b);
    b=acos(b*no');
    ang_GP(i)=(a+b);
end
%% list of nodes and connecting nodes
list_P=cell(length(list_RA),2);
list_N=cell(length(list_RA),2);
a=0;
b=0;
c=0;
d=sort(l(list_RA,1));
e=0;
f=sort(l(list_RA,2));
for i=1:length(d)
    if d(i) ~= c
        a=a+1;
        list_P{a,1}(1)=d(i);
    end
    c=d(i);
    if f(i) ~= e
        b=b+1;
        list_N{b,1}(1)=f(i);
    end
    e=f(i);
end
list_P=list_P(1:a,:);
list_N=list_N(1:b,:);
%% sort the coonnecting nodes and corners
flat=eye(3); flat(3,3)=0; % constant matric for make v(:,3)=0
for i=1:length(list_P)
    cand=l(list_RA(l(list_RA,1)==list_P{i}(1)),2);
    [~,b]=sort(atan2(v(cand,2)-v(list_P{i}(1),2),v(cand,1)-v(list_P{i}(1),1)));
    list_P{i,2}=cand(b);        %% coonnecting nodes
    
    cand=l(l(:,1)==list_P{i}(1) & isOnGround(l(:,2)),2);
    [~,b]=sort(atan2(v(cand,2)-v(list_P{i}(1),2),v(cand,1)-v(list_P{i}(1),1)));
    list_P{i,3}=cand(b);        %% connecting corners
    
    j1=[v(list_P{i,2}(1),:);v(list_P{i,3}(1),:)]-v(list_P{i,1}(1),:);
    j2=atan2(j1(:,2),j1(:,1));
    if j2(2)>j2(1)
        list_P{i,3}=list_P{i,3}([end 1:end-1]);
    end
    
    list_P{i,4}=v(list_P{i,1},:)*flat;	%% desired pivot
end
for i=1:length(list_N)
    cand=l(list_RA(l(list_RA,2)==list_N{i}(1)),1);
    [~,b]=sort(atan2(v(cand,2)-v(list_N{i}(1),2),v(cand,1)-v(list_N{i}(1),1)));
    list_N{i,2}=cand(b(end:-1:1));      %% coonnecting nodes
    
    cand=l(l(:,1)==list_N{i}(1) & isOnGround(l(:,2)),2);
    [~,b]=sort(atan2(v(cand,2)-v(list_N{i}(1),2),v(cand,1)-v(list_N{i}(1),1)));
    list_N{i,3}=cand(b(end:-1:1));      %% connecting corners
    
    j1=[v(list_N{i,2}(1),:);v(list_N{i,3}(1),:)]-v(list_N{i,1}(1),:);
    j2=atan2(j1(:,2),j1(:,1));
    if j2(2)<j2(1)
        list_N{i,3}=list_N{i,3}([end 1:end-1]);
    end
    
    list_N{i,4}=v(list_N{i,1},:)*flat;  %% desired pivot
end

%%  Display the index of nodes

for i = 1:length(list_P(:,1))
    text(v(list_P{i,1},1)+.010,v(list_P{i,1},2)+.050,num2str(i))
    plot(v(list_P{i,1},1),v(list_P{i,1},2),'r*')
end
for i = 1:length(list_N(:,1))
    text(v(list_N{i,1},1)+.010,v(list_N{i,1},2)+.050,num2str(i))
    plot(v(list_N{i,1},1),v(list_N{i,1},2),'ro')
end

%% Connection between nodes
list_P_in_N=cell(length(list_P(:,1)),1);
j1=cellfun(@(v)v(1),list_N(:,1));
j2=1:length(j1);
for i=1:length(list_P(:,1))
    
    for k=1:length(list_P{i,2})
        list_P_in_N{i}(k,1)=j2(j1==list_P{i,2}(k));
    end
    
    j3=list_N(list_P_in_N{i}(:,1),2);
    for k=1:length(j3)
        j4=1:length(j3{k});
        list_P_in_N{i}(k,2)=j4(j3{k}==list_P{i,1});
    end
end
%%
list_N_in_P=cell(length(list_N(:,1)),1);
j1=cellfun(@(v)v(1),list_P(:,1));
j2=1:length(j1);
for i=1:length(list_N(:,1))
    
    for k=1:length(list_N{i,2})
        list_N_in_P{i}(k,1)=j2(j1==list_N{i,2}(k));
    end
    
    j3=list_P(list_N_in_P{i}(:,1),2);
    for k=1:length(j3)
        j4=1:length(j3{k});
        list_N_in_P{i}(k,2)=j4(j3{k}==list_N{i,1});
    end
end

%% first rotating connector
tilt=15/180*pi;

i=23;
% generating the normal vector of bisector
pivotN=list_N{i,4}-v(list_N{i,1}(1),:);         pivotN=pivotN./sqrt(pivotN*pivotN');
cornerN=v(list_N{i,3},:)-v(list_N{i,1}(1),:);   cornerN=cornerN./(sqrt(diag(cornerN*cornerN'))*ones(1,3));
biseN=cornerN-pivotN;                            biseN=biseN./(sqrt(diag(biseN*biseN'))*ones(1,3));
% generating the hinges
k=1;
x2=cross(pivotN,biseN(k,:)); x2=x2./sqrt(x2*x2');
x1=cross(biseN(k,:),x2);
hinN=cos(tilt)*x1+sin(tilt)*x2;
list_N{i,5}(k,:)=v(list_N{i,1},:)-hinN/hinN(:,3)*v(list_N{i,1},3);
%%
for k=2:length(list_N{i,2})
    %%
    vec1=v(list_N{i,2}(k-1),:)-v(list_N{i,1}(1),:); x3=vec1./sqrt(vec1*vec1');
    x2=cross(vec1,hinN);                            x2=x2./sqrt(x2*x2');
    x1=cross(x2,vec1);                              x1=x1./sqrt(x1*x1');
    ang=ang_GP(l(list_RA,2)==list_N{i,1} & l(list_RA,1)==list_N{i,2}(k-1));
    pitN=(hinN*x1')*(cos(ang)*x1-sin(ang)*x2)+(hinN*x3')*x3;
    list_N{i,6}(k,:)=v(list_N{i,1},:)-pitN/pitN(:,3)*v(list_N{i,1},3);
    vec2=cos(ang/2)*x2+sin(ang/2)*x1;
    hinN=-cross(biseN(k,:),vec2);
    list_N{i,5}(k,:)=v(list_N{i,1},:)-hinN/hinN(:,3)*v(list_N{i,1},3);
end

vec1=v(list_N{i,2}(end),:)-v(list_N{i,1}(1),:); x3=vec1./sqrt(vec1*vec1');
x2=cross(vec1,hinN);                            x2=x2./sqrt(x2*x2');
x1=cross(x2,vec1);                              x1=x1./sqrt(x1*x1');
ang=ang_GP(l(list_RA,2)==list_N{i,1} & l(list_RA,1)==list_N{i,2}(end));
pitN=(hinN*x1')*(cos(ang)*x1-sin(ang)*x2)+(hinN*x3')*x3;
list_N{i,6}(1,:)=v(list_N{i,1},:)-pitN/pitN(:,3)*v(list_N{i,1},3);
%%

fid=DXF_header('Deg15.dxf');
dxf_a=1;
n=length(list_N{i,5}(:,1));
dxf_a=DXF_lines(fid,list_N{i,5}(:,:),[(1:n)', [(2:n)';1]],dxf_a);
dxf_a=DXF_lines(fid,[list_N{i,5}(:,:);list_N{i,6}(:,:)],[(1:n)', n+(1:n)'],dxf_a);

plot(list_N{i,5}([1:end 1],1),list_N{i,5}([1:end 1],2),'r+-')
plot(list_N{i,6}([1:end 1],1),list_N{i,6}([1:end 1],2),'r+')
%%  second generation rotating connector
for i=[19,27,18,13]
j3=list_N(list_P_in_N{i}(:,1),2);
seed=[];
a=1;
for k=1:length(j3)
    if ~isempty(list_N{list_P_in_N{i}(k,1),5})
        seed(a)=k;
        a=a+1;
    end
end
list_P_in_N{i}(seed,:)
ref=list_N{list_P_in_N{i}(seed,1),5}(list_P_in_N{i}(seed,2),:);

% generating the normal vector of bisector
pivotN=list_P{i,4}-v(list_P{i,1}(1),:);     	pivotN=pivotN./sqrt(pivotN*pivotN');
cornerN=v(list_P{i,3},:)-v(list_P{i,1}(1),:);   cornerN=cornerN./(sqrt(diag(cornerN*cornerN'))*ones(1,3));
biseN=cornerN-pivotN;                          	biseN=biseN./(sqrt(diag(biseN*biseN'))*ones(1,3));
% generating the rotation polygon
k=seed;
vec1=v(list_P{i,2}(k),:)-v(list_P{i,1}(1),:);
vec2=ref-v(list_P{i,1}(1),:);
x2=cross(vec1,vec2); x2=x2./sqrt(x2*x2');
hinN=-cross(biseN(k,:),x2);
list_P{i,5}(k,:)=v(list_P{i,1},:)-hinN/hinN(:,3)*v(list_P{i,1},3);

n=length(list_P{i,2});
for k=(seed+1):(seed+n-1)
    %%
    vec1=v(list_P{i,2}(mod(k-2,n)+1),:)-v(list_P{i,1}(1),:);     x3=vec1./sqrt(vec1*vec1');
    x2=cross(vec1,hinN);    x2=x2./sqrt(x2*x2');
    x1=cross(x2,vec1);      x1=x1./sqrt(x1*x1');
    ang=ang_GP(l(list_RA,1)==list_P{i,1} & l(list_RA,2)==list_P{i,2}(mod(k-2,n)+1));
    pitN=(hinN*x1')*(cos(ang)*x1+sin(ang)*x2)+(hinN*x3')*x3;
    list_P{i,6}(mod(k-1,n)+1,:)=v(list_P{i,1},:)-pitN/pitN(:,3)*v(list_P{i,1},3);
    vec2=cos(ang/2)*x2-sin(ang/2)*x1;
    hinN=-cross(biseN(mod(k-1,n)+1,:),vec2);
    list_P{i,5}(mod(k-1,n)+1,:)=v(list_P{i,1},:)-hinN/hinN(:,3)*v(list_P{i,1},3);
end
k=k+1;
    vec1=v(list_P{i,2}(mod(k-2,n)+1),:)-v(list_P{i,1}(1),:);     x3=vec1./sqrt(vec1*vec1');
    x2=cross(vec1,hinN);    x2=x2./sqrt(x2*x2');
    x1=cross(x2,vec1);      x1=x1./sqrt(x1*x1');
    ang=ang_GP(l(list_RA,1)==list_P{i,1} & l(list_RA,2)==list_P{i,2}(mod(k-2,n)+1));
    pitN=(hinN*x1')*(cos(ang)*x1+sin(ang)*x2)+(hinN*x3')*x3;
    list_P{i,6}(mod(k-1,n)+1,:)=v(list_P{i,1},:)-pitN/pitN(:,3)*v(list_P{i,1},3);


plot(list_P{i,5}([1:end 1],1),list_P{i,5}([1:end 1],2),'r+-')
plot(list_P{i,6}([1:end 1],1),list_P{i,6}([1:end 1],2),'r+')
n=length(list_P{i,5}(:,1));
dxf_a=DXF_lines(fid,list_P{i,5}([1:end 1],:),[(1:n)', [(2:n)';1]],dxf_a);
dxf_a=DXF_lines(fid,[list_P{i,5}(:,:);list_P{i,6}(:,:)],[(1:n)', n+(1:n)'],dxf_a);
end

%%  third generation rotating connector
for i=[29,30,31,32,12,16,17,18,19,24]%]
j3=list_P(list_N_in_P{i}(:,1),2);
seed=[];
a=1;
for k=1:length(j3)
    if ~isempty(list_P{list_N_in_P{i}(k,1),5})
        seed(a)=k;
        a=a+1;
    end
end

if length(seed)==2
    if seed(2)-seed(1)~=1
        seed=seed(end:-1:1);
    end
    n=length(list_P{list_N_in_P{i}(seed(1),1),5});
    ref1=list_P{list_N_in_P{i}(seed(1),1),5}(mod(list_N_in_P{i}(seed(1),2),n)+1,:);
    ref2=list_P{list_N_in_P{i}(seed(2),1),5}(list_N_in_P{i}(seed(2),2),:);
    %%
    % get the first hinge
    vec1=v(list_N{i,2}(seed(1)),:)-v(list_N{i,1}(1),:);
    x1=cross(vec1,ref1-v(list_N{i,1}(1),:));    x1=x1./sqrt(x1*x1');
    vec2=v(list_N{i,2}(seed(2)),:)-v(list_N{i,1}(1),:);
    x2=cross(vec2,ref2-v(list_N{i,1}(1),:));    x2=x2./sqrt(x2*x2');
    hinN=-cross(x1,x2);
    k=seed(2);
    list_N{i,5}(k,:)=v(list_N{i,1},:)-hinN/hinN(:,3)*v(list_N{i,1},3);
    % get the corners
    cornerN=v(list_N{i,3},:)-v(list_N{i,1}(1),:);   cornerN=cornerN./(sqrt(diag(cornerN*cornerN'))*ones(1,3));
    % find the pivot
    hinN=hinN./sqrt(hinN*hinN');
    pivotN=list_N{i,4}-v(list_N{i,1}(1),:);     	pivotN=pivotN./sqrt(pivotN*pivotN');
    vec1=cross(hinN,pivotN);                        vec1=vec1./sqrt(vec1*vec1');
    vec2=cross(vec1,hinN);
    j1=cornerN(k,:)*[hinN',vec1',vec2'];
    vec1=cross(hinN,cornerN(k,:));                  vec1=vec1./sqrt(vec1*vec1');
    vec2=cross(vec1,hinN);
    pivotN=[j1(1) -j1(2) j1(3)]*[hinN;vec1;vec2];
    list_N{i,4}=v(list_N{i,1},:)-pivotN/pivotN(:,3)*v(list_N{i,1},3);
    biseN=cornerN-pivotN;                          	biseN=biseN./(sqrt(diag(biseN*biseN'))*ones(1,3));
else
    ref=list_P{list_N_in_P{i}(seed,1),5}(list_N_in_P{i}(seed,2),:);
    % generating the normal vector of bisector
    pivotN=list_N{i,4}-v(list_N{i,1}(1),:);     	pivotN=pivotN./sqrt(pivotN*pivotN');
    cornerN=v(list_N{i,3},:)-v(list_N{i,1}(1),:);   cornerN=cornerN./(sqrt(diag(cornerN*cornerN'))*ones(1,3));
    biseN=cornerN-pivotN;                          	biseN=biseN./(sqrt(diag(biseN*biseN'))*ones(1,3));
    % generating the rotation polygon
    k=seed;
    vec1=v(list_N{i,2}(k),:)-v(list_N{i,1}(1),:);
    vec2=ref-v(list_N{i,1}(1),:);
    x2=cross(vec1,vec2); x2=x2./sqrt(x2*x2');
    hinN=-cross(biseN(k,:),x2);
    list_N{i,5}(k,:)=v(list_N{i,1},:)-hinN/hinN(:,3)*v(list_N{i,1},3);
end
n=length(list_N{i,2});
for k=(k+1):(k+n)
    %%
    vec1=v(list_N{i,2}(mod(k-2,n)+1),:)-v(list_N{i,1}(1),:);    x3=vec1./sqrt(vec1*vec1');
    x2=cross(vec1,hinN);  x2=x2./sqrt(x2*x2');
    x1=cross(x2,vec1);  x1=x1./sqrt(x1*x1');
    ang=ang_GP(l(list_RA,2)==list_N{i,1} & l(list_RA,1)==list_N{i,2}(mod(k-2,n)+1));
    pitN=(hinN*x1')*(cos(ang)*x1-sin(ang)*x2)+(hinN*x3')*x3;
    list_N{i,6}(mod(k-1,n)+1,:)=v(list_N{i,1},:)-pitN/pitN(:,3)*v(list_N{i,1},3);
    vec2=cos(ang/2)*x2+sin(ang/2)*x1;
    hinN=-cross(biseN(mod(k-1,n)+1,:),vec2);
    list_N{i,5}(mod(k-1,n)+1,:)=v(list_N{i,1},:)-hinN/hinN(:,3)*v(list_N{i,1},3);
end
k=k+1;
    vec1=v(list_N{i,2}(mod(k-2,n)+1),:)-v(list_N{i,1}(1),:);    x3=vec1./sqrt(vec1*vec1');
    x2=cross(vec1,hinN);  x2=x2./sqrt(x2*x2');
    x1=cross(x2,vec1);  x1=x1./sqrt(x1*x1');
    ang=ang_GP(l(list_RA,2)==list_N{i,1} & l(list_RA,1)==list_N{i,2}(mod(k-2,n)+1));
    pitN=(hinN*x1')*(cos(ang)*x1-sin(ang)*x2)+(hinN*x3')*x3;
    list_N{i,6}(mod(k-1,n)+1,:)=v(list_N{i,1},:)-pitN/pitN(:,3)*v(list_N{i,1},3);
    
plot(list_N{i,5}([1:end 1],1),list_N{i,5}([1:end 1],2),'r+-')
plot(list_N{i,6}([1:end 1],1),list_N{i,6}([1:end 1],2),'r+')
n=length(list_N{i,5}(:,1));
dxf_a=DXF_lines(fid,list_N{i,5}(:,:),[(1:n)', [(2:n)';1]],dxf_a);
dxf_a=DXF_lines(fid,[list_N{i,5}(:,:);list_N{i,6}(:,:)],[(1:n)', n+(1:n)'],dxf_a);
end

%%  forth generation rotating connector
for i=[28,30,29,7,9,12,14,15,20,22]
j3=list_N(list_P_in_N{i}(:,1),2);
seed=[];
a=1;
for k=1:length(j3)
    if ~isempty(list_N{list_P_in_N{i}(k,1),5})
        seed(a)=k;
        a=a+1;
    end
end

if length(seed)==2
    if seed(2)-seed(1)~=1
        seed=seed(end:-1:1);
    end
    n=length(list_N{list_P_in_N{i}(seed(1),1),5});
    ref1=list_N{list_P_in_N{i}(seed(1),1),5}(mod(list_P_in_N{i}(seed(1),2),n)+1,:);
    ref2=list_N{list_P_in_N{i}(seed(2),1),5}(list_P_in_N{i}(seed(2),2),:);
    %%
    % get the first hinge
    vec1=v(list_P{i,2}(seed(1)),:)-v(list_P{i,1}(1),:);
    x1=cross(vec1,ref1-v(list_P{i,1}(1),:));    x1=x1./sqrt(x1*x1');
    vec2=v(list_P{i,2}(seed(2)),:)-v(list_P{i,1}(1),:);
    x2=cross(vec2,ref2-v(list_P{i,1}(1),:));    x2=x2./sqrt(x2*x2');
    hinN=-cross(x1,x2);
    k=seed(2);
    list_P{i,5}(k,:)=v(list_P{i,1},:)-hinN/hinN(:,3)*v(list_P{i,1},3);
    % get the corners
    cornerN=v(list_P{i,3},:)-v(list_P{i,1}(1),:);   cornerN=cornerN./(sqrt(diag(cornerN*cornerN'))*ones(1,3));
    % find the pivot
    hinN=hinN./sqrt(hinN*hinN');
    pivotN=list_P{i,4}-v(list_P{i,1}(1),:);     	pivotN=pivotN./sqrt(pivotN*pivotN');
    vec1=cross(hinN,pivotN);                        vec1=vec1./sqrt(vec1*vec1');
    vec2=cross(vec1,hinN);
    j1=cornerN(k,:)*[hinN',vec1',vec2'];
    vec1=cross(hinN,cornerN(k,:));                  vec1=vec1./sqrt(vec1*vec1');
    vec2=cross(vec1,hinN);
    pivotN=[j1(1) -j1(2) j1(3)]*[hinN;vec1;vec2];
    list_P{i,4}=v(list_P{i,1},:)-pivotN/pivotN(:,3)*v(list_P{i,1},3);
    biseN=cornerN-pivotN;                          	biseN=biseN./(sqrt(diag(biseN*biseN'))*ones(1,3));
else
    ref=list_N{list_P_in_N{i}(seed,1),5}(list_P_in_N{i}(seed,2),:);
    % generating the normal vector of bisector
    pivotN=list_P{i,4}-v(list_P{i,1}(1),:);     	pivotN=pivotN./sqrt(pivotN*pivotN');
    cornerN=v(list_P{i,3},:)-v(list_P{i,1}(1),:);   cornerN=cornerN./(sqrt(diag(cornerN*cornerN'))*ones(1,3));
    biseN=cornerN-pivotN;                          	biseN=biseN./(sqrt(diag(biseN*biseN'))*ones(1,3));
    % generating the rotation polygon
    k=seed;
    vec1=v(list_P{i,2}(k),:)-v(list_P{i,1}(1),:);
    vec2=ref-v(list_P{i,1}(1),:);
    x2=cross(vec1,vec2); x2=x2./sqrt(x2*x2');
    hinN=-cross(biseN(k,:),x2);
    list_P{i,5}(k,:)=v(list_P{i,1},:)-hinN/hinN(:,3)*v(list_P{i,1},3);
end
n=length(list_P{i,2});
for k=(k+1):(k+n)
    %%
    vec1=v(list_P{i,2}(mod(k-2,n)+1),:)-v(list_P{i,1}(1),:);        x3=vec1./sqrt(vec1*vec1');
    x2=cross(vec1,hinN);  x2=x2./sqrt(x2*x2');
    x1=cross(x2,vec1);  x1=x1./sqrt(x1*x1');
    ang=ang_GP(l(list_RA,1)==list_P{i,1} & l(list_RA,2)==list_P{i,2}(mod(k-2,n)+1));
    pitN=(hinN*x1')*(cos(ang)*x1+sin(ang)*x2)+(hinN*x3')*x3;
    list_P{i,6}(mod(k-1,n)+1,:)=v(list_P{i,1},:)-pitN/pitN(:,3)*v(list_P{i,1},3);
    vec2=cos(ang/2)*x2-sin(ang/2)*x1;
    hinN=-cross(biseN(mod(k-1,n)+1,:),vec2);
    list_P{i,5}(mod(k-1,n)+1,:)=v(list_P{i,1},:)-hinN/hinN(:,3)*v(list_P{i,1},3);
end
k=k+1;
    vec1=v(list_P{i,2}(mod(k-2,n)+1),:)-v(list_P{i,1}(1),:);        x3=vec1./sqrt(vec1*vec1');
    x2=cross(vec1,hinN);  x2=x2./sqrt(x2*x2');
    x1=cross(x2,vec1);  x1=x1./sqrt(x1*x1');
    ang=ang_GP(l(list_RA,1)==list_P{i,1} & l(list_RA,2)==list_P{i,2}(mod(k-2,n)+1));
    pitN=(hinN*x1')*(cos(ang)*x1+sin(ang)*x2)+(hinN*x3')*x3;
    list_P{i,6}(mod(k-1,n)+1,:)=v(list_P{i,1},:)-pitN/pitN(:,3)*v(list_P{i,1},3);
    
plot(list_P{i,5}([1:end 1],1),list_P{i,5}([1:end 1],2),'r+-')
plot(list_P{i,6}([1:end 1],1),list_P{i,6}([1:end 1],2),'r+')
n=length(list_P{i,5}(:,1));
dxf_a=DXF_lines(fid,list_P{i,5}([1:end 1],:),[(1:n)', [(2:n)';1]],dxf_a);
dxf_a=DXF_lines(fid,[list_P{i,5}(:,:);list_P{i,6}(:,:)],[(1:n)', n+(1:n)'],dxf_a);
end

%%  fifth generation rotating connector
for i=[34,4,33,10,14,27]
j3=list_P(list_N_in_P{i}(:,1),2);
seed=[];
a=1;
for k=1:length(j3)
    if ~isempty(list_P{list_N_in_P{i}(k,1),5})
        seed(a)=k;
        a=a+1;
    end
end

if length(seed)==2
    if seed(2)-seed(1)~=1
        seed=seed(end:-1:1);
    end
    n=length(list_P{list_N_in_P{i}(seed(1),1),5});
    ref1=list_P{list_N_in_P{i}(seed(1),1),5}(mod(list_N_in_P{i}(seed(1),2),n)+1,:);
    ref2=list_P{list_N_in_P{i}(seed(2),1),5}(list_N_in_P{i}(seed(2),2),:);
    %%
    % get the first hinge
    vec1=v(list_N{i,2}(seed(1)),:)-v(list_N{i,1}(1),:);
    x1=cross(vec1,ref1-v(list_N{i,1}(1),:));    x1=x1./sqrt(x1*x1');
    vec2=v(list_N{i,2}(seed(2)),:)-v(list_N{i,1}(1),:);
    x2=cross(vec2,ref2-v(list_N{i,1}(1),:));    x2=x2./sqrt(x2*x2');
    hinN=-cross(x1,x2);
    k=seed(2);
    list_N{i,5}(k,:)=v(list_N{i,1},:)-hinN/hinN(:,3)*v(list_N{i,1},3);
    % get the corners
    cornerN=v(list_N{i,3},:)-v(list_N{i,1}(1),:);   cornerN=cornerN./(sqrt(diag(cornerN*cornerN'))*ones(1,3));
    % find the pivot
    hinN=hinN./sqrt(hinN*hinN');
    pivotN=list_N{i,4}-v(list_N{i,1}(1),:);     	pivotN=pivotN./sqrt(pivotN*pivotN');
    vec1=cross(hinN,pivotN);                        vec1=vec1./sqrt(vec1*vec1');
    vec2=cross(vec1,hinN);
    j1=cornerN(k,:)*[hinN',vec1',vec2'];
    vec1=cross(hinN,cornerN(k,:));                  vec1=vec1./sqrt(vec1*vec1');
    vec2=cross(vec1,hinN);
    pivotN=[j1(1) -j1(2) j1(3)]*[hinN;vec1;vec2];
    list_N{i,4}=v(list_N{i,1},:)-pivotN/pivotN(:,3)*v(list_N{i,1},3);
    biseN=cornerN-pivotN;                          	biseN=biseN./(sqrt(diag(biseN*biseN'))*ones(1,3));
else
    ref=list_P{list_N_in_P{i}(seed,1),5}(list_N_in_P{i}(seed,2),:);
    % generating the normal vector of bisector
    pivotN=list_N{i,4}-v(list_N{i,1}(1),:);     	pivotN=pivotN./sqrt(pivotN*pivotN');
    cornerN=v(list_N{i,3},:)-v(list_N{i,1}(1),:);   cornerN=cornerN./(sqrt(diag(cornerN*cornerN'))*ones(1,3));
    biseN=cornerN-pivotN;                          	biseN=biseN./(sqrt(diag(biseN*biseN'))*ones(1,3));
    % generating the rotation polygon
    k=seed;
    vec1=v(list_N{i,2}(k),:)-v(list_N{i,1}(1),:);
    vec2=ref-v(list_N{i,1}(1),:);
    x2=cross(vec1,vec2); x2=x2./sqrt(x2*x2');
    hinN=-cross(biseN(k,:),x2);
    list_N{i,5}(k,:)=v(list_N{i,1},:)-hinN/hinN(:,3)*v(list_N{i,1},3);
end
n=length(list_N{i,2});
for k=(k+1):(k+n)
    %%
    vec1=v(list_N{i,2}(mod(k-2,n)+1),:)-v(list_N{i,1}(1),:);    x3=vec1./sqrt(vec1*vec1');
    x2=cross(vec1,hinN);  x2=x2./sqrt(x2*x2');
    x1=cross(x2,vec1);  x1=x1./sqrt(x1*x1');
    ang=ang_GP(l(list_RA,2)==list_N{i,1} & l(list_RA,1)==list_N{i,2}(mod(k-2,n)+1));
    pitN=(hinN*x1')*(cos(ang)*x1-sin(ang)*x2)+(hinN*x3')*x3;
    list_N{i,6}(mod(k-1,n)+1,:)=v(list_N{i,1},:)-pitN/pitN(:,3)*v(list_N{i,1},3);
    vec2=cos(ang/2)*x2+sin(ang/2)*x1;
    hinN=-cross(biseN(mod(k-1,n)+1,:),vec2);
    list_N{i,5}(mod(k-1,n)+1,:)=v(list_N{i,1},:)-hinN/hinN(:,3)*v(list_N{i,1},3);
end
k=k+1;
    vec1=v(list_N{i,2}(mod(k-2,n)+1),:)-v(list_N{i,1}(1),:);    x3=vec1./sqrt(vec1*vec1');
    x2=cross(vec1,hinN);  x2=x2./sqrt(x2*x2');
    x1=cross(x2,vec1);  x1=x1./sqrt(x1*x1');
    ang=ang_GP(l(list_RA,2)==list_N{i,1} & l(list_RA,1)==list_N{i,2}(mod(k-2,n)+1));
    pitN=(hinN*x1')*(cos(ang)*x1-sin(ang)*x2)+(hinN*x3')*x3;
    list_N{i,6}(mod(k-1,n)+1,:)=v(list_N{i,1},:)-pitN/pitN(:,3)*v(list_N{i,1},3);
    
plot(list_N{i,5}([1:end 1],1),list_N{i,5}([1:end 1],2),'r+-')
plot(list_N{i,6}([1:end 1],1),list_N{i,6}([1:end 1],2),'r+')
n=length(list_N{i,5}(:,1));
dxf_a=DXF_lines(fid,list_N{i,5}(:,:),[(1:n)', [(2:n)';1]],dxf_a);
dxf_a=DXF_lines(fid,[list_N{i,5}(:,:);list_N{i,6}(:,:)],[(1:n)', n+(1:n)'],dxf_a);
end

%%  6th generation rotating connector
for i=[32,1,2,5]
j3=list_N(list_P_in_N{i}(:,1),2);
seed=[];
a=1;
for k=1:length(j3)
    if ~isempty(list_N{list_P_in_N{i}(k,1),5})
        seed(a)=k;
        a=a+1;
    end
end

if length(seed)==2
    if seed(2)-seed(1)~=1
        seed=seed(end:-1:1);
    end
    n=length(list_N{list_P_in_N{i}(seed(1),1),5});
    ref1=list_N{list_P_in_N{i}(seed(1),1),5}(mod(list_P_in_N{i}(seed(1),2),n)+1,:);
    ref2=list_N{list_P_in_N{i}(seed(2),1),5}(list_P_in_N{i}(seed(2),2),:);
    %%
    % get the first hinge
    vec1=v(list_P{i,2}(seed(1)),:)-v(list_P{i,1}(1),:);
    x1=cross(vec1,ref1-v(list_P{i,1}(1),:));    x1=x1./sqrt(x1*x1');
    vec2=v(list_P{i,2}(seed(2)),:)-v(list_P{i,1}(1),:);
    x2=cross(vec2,ref2-v(list_P{i,1}(1),:));    x2=x2./sqrt(x2*x2');
    hinN=-cross(x1,x2);
    k=seed(2);
    list_P{i,5}(k,:)=v(list_P{i,1},:)-hinN/hinN(:,3)*v(list_P{i,1},3);
    % get the corners
    cornerN=v(list_P{i,3},:)-v(list_P{i,1}(1),:);   cornerN=cornerN./(sqrt(diag(cornerN*cornerN'))*ones(1,3));
    % find the pivot
    hinN=hinN./sqrt(hinN*hinN');
    pivotN=list_P{i,4}-v(list_P{i,1}(1),:);     	pivotN=pivotN./sqrt(pivotN*pivotN');
    vec1=cross(hinN,pivotN);                        vec1=vec1./sqrt(vec1*vec1');
    vec2=cross(vec1,hinN);
    j1=cornerN(k,:)*[hinN',vec1',vec2'];
    vec1=cross(hinN,cornerN(k,:));                  vec1=vec1./sqrt(vec1*vec1');
    vec2=cross(vec1,hinN);
    pivotN=[j1(1) -j1(2) j1(3)]*[hinN;vec1;vec2];
    list_P{i,4}=v(list_P{i,1},:)-pivotN/pivotN(:,3)*v(list_P{i,1},3);
    biseN=cornerN-pivotN;                          	biseN=biseN./(sqrt(diag(biseN*biseN'))*ones(1,3));
else
    ref=list_N{list_P_in_N{i}(seed,1),5}(list_P_in_N{i}(seed,2),:);
    % generating the normal vector of bisector
    pivotN=list_P{i,4}-v(list_P{i,1}(1),:);     	pivotN=pivotN./sqrt(pivotN*pivotN');
    cornerN=v(list_P{i,3},:)-v(list_P{i,1}(1),:);   cornerN=cornerN./(sqrt(diag(cornerN*cornerN'))*ones(1,3));
    biseN=cornerN-pivotN;                          	biseN=biseN./(sqrt(diag(biseN*biseN'))*ones(1,3));
    % generating the rotation polygon
    k=seed;
    vec1=v(list_P{i,2}(k),:)-v(list_P{i,1}(1),:);
    vec2=ref-v(list_P{i,1}(1),:);
    x2=cross(vec1,vec2); x2=x2./sqrt(x2*x2');
    hinN=-cross(biseN(k,:),x2);
    list_P{i,5}(k,:)=v(list_P{i,1},:)-hinN/hinN(:,3)*v(list_P{i,1},3);
end
n=length(list_P{i,2});
for k=(k+1):(k+n)
    %%
    vec1=v(list_P{i,2}(mod(k-2,n)+1),:)-v(list_P{i,1}(1),:);        x3=vec1./sqrt(vec1*vec1');
    x2=cross(vec1,hinN);  x2=x2./sqrt(x2*x2');
    x1=cross(x2,vec1);  x1=x1./sqrt(x1*x1');
    ang=ang_GP(l(list_RA,1)==list_P{i,1} & l(list_RA,2)==list_P{i,2}(mod(k-2,n)+1));
    pitN=(hinN*x1')*(cos(ang)*x1+sin(ang)*x2)+(hinN*x3')*x3;
    list_P{i,6}(mod(k-1,n)+1,:)=v(list_P{i,1},:)-pitN/pitN(:,3)*v(list_P{i,1},3);
    vec2=cos(ang/2)*x2-sin(ang/2)*x1;
    hinN=-cross(biseN(mod(k-1,n)+1,:),vec2);
    list_P{i,5}(mod(k-1,n)+1,:)=v(list_P{i,1},:)-hinN/hinN(:,3)*v(list_P{i,1},3);
end
k=k+1;
    vec1=v(list_P{i,2}(mod(k-2,n)+1),:)-v(list_P{i,1}(1),:);        x3=vec1./sqrt(vec1*vec1');
    x2=cross(vec1,hinN);  x2=x2./sqrt(x2*x2');
    x1=cross(x2,vec1);  x1=x1./sqrt(x1*x1');
    ang=ang_GP(l(list_RA,1)==list_P{i,1} & l(list_RA,2)==list_P{i,2}(mod(k-2,n)+1));
    pitN=(hinN*x1')*(cos(ang)*x1+sin(ang)*x2)+(hinN*x3')*x3;
    list_P{i,6}(mod(k-1,n)+1,:)=v(list_P{i,1},:)-pitN/pitN(:,3)*v(list_P{i,1},3);
    
plot(list_P{i,5}([1:end 1],1),list_P{i,5}([1:end 1],2),'r+-')
plot(list_P{i,6}([1:end 1],1),list_P{i,6}([1:end 1],2),'r+')
n=length(list_P{i,5}(:,1));
dxf_a=DXF_lines(fid,list_P{i,5}([1:end 1],:),[(1:n)', [(2:n)';1]],dxf_a);
dxf_a=DXF_lines(fid,[list_P{i,5}(:,:);list_P{i,6}(:,:)],[(1:n)', n+(1:n)'],dxf_a);
end

%%  7th generation rotating connector
for i=[37,6]
j3=list_P(list_N_in_P{i}(:,1),2);
seed=[];
a=1;
for k=1:length(j3)
    if ~isempty(list_P{list_N_in_P{i}(k,1),5})
        seed(a)=k;
        a=a+1;
    end
end

if length(seed)==2
    if seed(2)-seed(1)~=1
        seed=seed(end:-1:1);
    end
    n=length(list_P{list_N_in_P{i}(seed(1),1),5});
    ref1=list_P{list_N_in_P{i}(seed(1),1),5}(mod(list_N_in_P{i}(seed(1),2),n)+1,:);
    ref2=list_P{list_N_in_P{i}(seed(2),1),5}(list_N_in_P{i}(seed(2),2),:);
    %%
    % get the first hinge
    vec1=v(list_N{i,2}(seed(1)),:)-v(list_N{i,1}(1),:);
    x1=cross(vec1,ref1-v(list_N{i,1}(1),:));    x1=x1./sqrt(x1*x1');
    vec2=v(list_N{i,2}(seed(2)),:)-v(list_N{i,1}(1),:);
    x2=cross(vec2,ref2-v(list_N{i,1}(1),:));    x2=x2./sqrt(x2*x2');
    hinN=-cross(x1,x2);
    k=seed(2);
    list_N{i,5}(k,:)=v(list_N{i,1},:)-hinN/hinN(:,3)*v(list_N{i,1},3);
    % get the corners
    cornerN=v(list_N{i,3},:)-v(list_N{i,1}(1),:);   cornerN=cornerN./(sqrt(diag(cornerN*cornerN'))*ones(1,3));
    % find the pivot
    hinN=hinN./sqrt(hinN*hinN');
    pivotN=list_N{i,4}-v(list_N{i,1}(1),:);     	pivotN=pivotN./sqrt(pivotN*pivotN');
    vec1=cross(hinN,pivotN);                        vec1=vec1./sqrt(vec1*vec1');
    vec2=cross(vec1,hinN);
    j1=cornerN(k,:)*[hinN',vec1',vec2'];
    vec1=cross(hinN,cornerN(k,:));                  vec1=vec1./sqrt(vec1*vec1');
    vec2=cross(vec1,hinN);
    pivotN=[j1(1) -j1(2) j1(3)]*[hinN;vec1;vec2];
    list_N{i,4}=v(list_N{i,1},:)-pivotN/pivotN(:,3)*v(list_N{i,1},3);
    biseN=cornerN-pivotN;                          	biseN=biseN./(sqrt(diag(biseN*biseN'))*ones(1,3));
else
    ref=list_P{list_N_in_P{i}(seed,1),5}(list_N_in_P{i}(seed,2),:);
    % generating the normal vector of bisector
    pivotN=list_N{i,4}-v(list_N{i,1}(1),:);     	pivotN=pivotN./sqrt(pivotN*pivotN');
    cornerN=v(list_N{i,3},:)-v(list_N{i,1}(1),:);   cornerN=cornerN./(sqrt(diag(cornerN*cornerN'))*ones(1,3));
    biseN=cornerN-pivotN;                          	biseN=biseN./(sqrt(diag(biseN*biseN'))*ones(1,3));
    % generating the rotation polygon
    k=seed;
    vec1=v(list_N{i,2}(k),:)-v(list_N{i,1}(1),:);
    vec2=ref-v(list_N{i,1}(1),:);
    x2=cross(vec1,vec2); x2=x2./sqrt(x2*x2');
    hinN=-cross(biseN(k,:),x2);
    list_N{i,5}(k,:)=v(list_N{i,1},:)-hinN/hinN(:,3)*v(list_N{i,1},3);
end
n=length(list_N{i,2});
for k=(k+1):(k+n)
    %%
    vec1=v(list_N{i,2}(mod(k-2,n)+1),:)-v(list_N{i,1}(1),:);    x3=vec1./sqrt(vec1*vec1');
    x2=cross(vec1,hinN);  x2=x2./sqrt(x2*x2');
    x1=cross(x2,vec1);  x1=x1./sqrt(x1*x1');
    ang=ang_GP(l(list_RA,2)==list_N{i,1} & l(list_RA,1)==list_N{i,2}(mod(k-2,n)+1));
    pitN=(hinN*x1')*(cos(ang)*x1-sin(ang)*x2)+(hinN*x3')*x3;
    list_N{i,6}(mod(k-1,n)+1,:)=v(list_N{i,1},:)-pitN/pitN(:,3)*v(list_N{i,1},3);
    vec2=cos(ang/2)*x2+sin(ang/2)*x1;
    hinN=-cross(biseN(mod(k-1,n)+1,:),vec2);
    list_N{i,5}(mod(k-1,n)+1,:)=v(list_N{i,1},:)-hinN/hinN(:,3)*v(list_N{i,1},3);
end
k=k+1;
    vec1=v(list_N{i,2}(mod(k-2,n)+1),:)-v(list_N{i,1}(1),:);    x3=vec1./sqrt(vec1*vec1');
    x2=cross(vec1,hinN);  x2=x2./sqrt(x2*x2');
    x1=cross(x2,vec1);  x1=x1./sqrt(x1*x1');
    ang=ang_GP(l(list_RA,2)==list_N{i,1} & l(list_RA,1)==list_N{i,2}(mod(k-2,n)+1));
    pitN=(hinN*x1')*(cos(ang)*x1-sin(ang)*x2)+(hinN*x3')*x3;
    list_N{i,6}(mod(k-1,n)+1,:)=v(list_N{i,1},:)-pitN/pitN(:,3)*v(list_N{i,1},3);
    
plot(list_N{i,5}([1:end 1],1),list_N{i,5}([1:end 1],2),'r+-')
plot(list_N{i,6}([1:end 1],1),list_N{i,6}([1:end 1],2),'r+')
n=length(list_N{i,5}(:,1));
dxf_a=DXF_lines(fid,list_N{i,5}(:,:),[(1:n)', [(2:n)';1]],dxf_a);
dxf_a=DXF_lines(fid,[list_N{i,5}(:,:);list_N{i,6}(:,:)],[(1:n)', n+(1:n)'],dxf_a);
end

%%  8th generation rotating connector
for i=[3]
j3=list_N(list_P_in_N{i}(:,1),2);
seed=[];
a=1;
for k=1:length(j3)
    if ~isempty(list_N{list_P_in_N{i}(k,1),5})
        seed(a)=k;
        a=a+1;
    end
end

if length(seed)==2
    if seed(2)-seed(1)~=1
        seed=seed(end:-1:1);
    end
    n=length(list_N{list_P_in_N{i}(seed(1),1),5});
    ref1=list_N{list_P_in_N{i}(seed(1),1),5}(mod(list_P_in_N{i}(seed(1),2),n)+1,:);
    ref2=list_N{list_P_in_N{i}(seed(2),1),5}(list_P_in_N{i}(seed(2),2),:);
    %%
    % get the first hinge
    vec1=v(list_P{i,2}(seed(1)),:)-v(list_P{i,1}(1),:);
    x1=cross(vec1,ref1-v(list_P{i,1}(1),:));    x1=x1./sqrt(x1*x1');
    vec2=v(list_P{i,2}(seed(2)),:)-v(list_P{i,1}(1),:);
    x2=cross(vec2,ref2-v(list_P{i,1}(1),:));    x2=x2./sqrt(x2*x2');
    hinN=-cross(x1,x2);
    k=seed(2);
    list_P{i,5}(k,:)=v(list_P{i,1},:)-hinN/hinN(:,3)*v(list_P{i,1},3);
    % get the corners
    cornerN=v(list_P{i,3},:)-v(list_P{i,1}(1),:);   cornerN=cornerN./(sqrt(diag(cornerN*cornerN'))*ones(1,3));
    % find the pivot
    hinN=hinN./sqrt(hinN*hinN');
    pivotN=list_P{i,4}-v(list_P{i,1}(1),:);     	pivotN=pivotN./sqrt(pivotN*pivotN');
    vec1=cross(hinN,pivotN);                        vec1=vec1./sqrt(vec1*vec1');
    vec2=cross(vec1,hinN);
    j1=cornerN(k,:)*[hinN',vec1',vec2'];
    vec1=cross(hinN,cornerN(k,:));                  vec1=vec1./sqrt(vec1*vec1');
    vec2=cross(vec1,hinN);
    pivotN=[j1(1) -j1(2) j1(3)]*[hinN;vec1;vec2];
    list_P{i,4}=v(list_P{i,1},:)-pivotN/pivotN(:,3)*v(list_P{i,1},3);
    biseN=cornerN-pivotN;                          	biseN=biseN./(sqrt(diag(biseN*biseN'))*ones(1,3));
else
    ref=list_N{list_P_in_N{i}(seed,1),5}(list_P_in_N{i}(seed,2),:);
    % generating the normal vector of bisector
    pivotN=list_P{i,4}-v(list_P{i,1}(1),:);     	pivotN=pivotN./sqrt(pivotN*pivotN');
    cornerN=v(list_P{i,3},:)-v(list_P{i,1}(1),:);   cornerN=cornerN./(sqrt(diag(cornerN*cornerN'))*ones(1,3));
    biseN=cornerN-pivotN;                          	biseN=biseN./(sqrt(diag(biseN*biseN'))*ones(1,3));
    % generating the rotation polygon
    k=seed;
    vec1=v(list_P{i,2}(k),:)-v(list_P{i,1}(1),:);
    vec2=ref-v(list_P{i,1}(1),:);
    x2=cross(vec1,vec2); x2=x2./sqrt(x2*x2');
    hinN=-cross(biseN(k,:),x2);
    list_P{i,5}(k,:)=v(list_P{i,1},:)-hinN/hinN(:,3)*v(list_P{i,1},3);
end
n=length(list_P{i,2});
for k=(k+1):(k+n)
    %%
    vec1=v(list_P{i,2}(mod(k-2,n)+1),:)-v(list_P{i,1}(1),:);        x3=vec1./sqrt(vec1*vec1');
    x2=cross(vec1,hinN);  x2=x2./sqrt(x2*x2');
    x1=cross(x2,vec1);  x1=x1./sqrt(x1*x1');
    ang=ang_GP(l(list_RA,1)==list_P{i,1} & l(list_RA,2)==list_P{i,2}(mod(k-2,n)+1));
    pitN=(hinN*x1')*(cos(ang)*x1+sin(ang)*x2)+(hinN*x3')*x3;
    list_P{i,6}(mod(k-1,n)+1,:)=v(list_P{i,1},:)-pitN/pitN(:,3)*v(list_P{i,1},3);
    vec2=cos(ang/2)*x2-sin(ang/2)*x1;
    hinN=-cross(biseN(mod(k-1,n)+1,:),vec2);
    list_P{i,5}(mod(k-1,n)+1,:)=v(list_P{i,1},:)-hinN/hinN(:,3)*v(list_P{i,1},3);
end
k=k+1;
    vec1=v(list_P{i,2}(mod(k-2,n)+1),:)-v(list_P{i,1}(1),:);        x3=vec1./sqrt(vec1*vec1');
    x2=cross(vec1,hinN);  x2=x2./sqrt(x2*x2');
    x1=cross(x2,vec1);  x1=x1./sqrt(x1*x1');
    ang=ang_GP(l(list_RA,1)==list_P{i,1} & l(list_RA,2)==list_P{i,2}(mod(k-2,n)+1));
    pitN=(hinN*x1')*(cos(ang)*x1+sin(ang)*x2)+(hinN*x3')*x3;
    list_P{i,6}(mod(k-1,n)+1,:)=v(list_P{i,1},:)-pitN/pitN(:,3)*v(list_P{i,1},3);
    
plot(list_P{i,5}([1:end 1],1),list_P{i,5}([1:end 1],2),'r+-')
plot(list_P{i,6}([1:end 1],1),list_P{i,6}([1:end 1],2),'r+')
n=length(list_P{i,5}(:,1));
dxf_a=DXF_lines(fid,list_P{i,5}([1:end 1],:),[(1:n)', [(2:n)';1]],dxf_a);
dxf_a=DXF_lines(fid,[list_P{i,5}(:,:);list_P{i,6}(:,:)],[(1:n)', n+(1:n)'],dxf_a);
end
DXF_end(fid);
if 0 
%%  9th generation rotating connector
for i=[16]
j3=list_P(list_N_in_P{i}(:,1),2);
seed=[];
a=1;
for k=1:length(j3)
    if ~isempty(list_P{list_N_in_P{i}(k,1),5})
        seed(a)=k;
        a=a+1;
    end
end

if length(seed)==2
    if seed(2)-seed(1)~=1
        seed=seed(end:-1:1);
    end
    n=length(list_P{list_N_in_P{i}(seed(1),1),5});
    ref1=list_P{list_N_in_P{i}(seed(1),1),5}(mod(list_N_in_P{i}(seed(1),2),n)+1,:);
    ref2=list_P{list_N_in_P{i}(seed(2),1),5}(list_N_in_P{i}(seed(2),2),:);
    %%
    % get the first hinge
    vec1=v(list_N{i,2}(seed(1)),:)-v(list_N{i,1}(1),:);
    x1=cross(vec1,ref1-v(list_N{i,1}(1),:));    x1=x1./sqrt(x1*x1');
    vec2=v(list_N{i,2}(seed(2)),:)-v(list_N{i,1}(1),:);
    x2=cross(vec2,ref2-v(list_N{i,1}(1),:));    x2=x2./sqrt(x2*x2');
    hinN=-cross(x1,x2);
    k=seed(2);
    list_N{i,5}(k,:)=v(list_N{i,1},:)-hinN/hinN(:,3)*v(list_N{i,1},3);
    % get the corners
    cornerN=v(list_N{i,3},:)-v(list_N{i,1}(1),:);   cornerN=cornerN./(sqrt(diag(cornerN*cornerN'))*ones(1,3));
    % find the pivot
    hinN=hinN./sqrt(hinN*hinN');
    pivotN=list_N{i,4}-v(list_N{i,1}(1),:);     	pivotN=pivotN./sqrt(pivotN*pivotN');
    vec1=cross(hinN,pivotN);                        vec1=vec1./sqrt(vec1*vec1');
    vec2=cross(vec1,hinN);
    j1=cornerN(k,:)*[hinN',vec1',vec2'];
    vec1=cross(hinN,cornerN(k,:));                  vec1=vec1./sqrt(vec1*vec1');
    vec2=cross(vec1,hinN);
    pivotN=[j1(1) -j1(2) j1(3)]*[hinN;vec1;vec2];
    list_N{i,4}=v(list_N{i,1},:)-pivotN/pivotN(:,3)*v(list_N{i,1},3);
    biseN=cornerN-pivotN;                          	biseN=biseN./(sqrt(diag(biseN*biseN'))*ones(1,3));
else
    ref=list_P{list_N_in_P{i}(seed,1),5}(list_N_in_P{i}(seed,2),:);
    % generating the normal vector of bisector
    pivotN=list_N{i,4}-v(list_N{i,1}(1),:);     	pivotN=pivotN./sqrt(pivotN*pivotN');
    cornerN=v(list_N{i,3},:)-v(list_N{i,1}(1),:);   cornerN=cornerN./(sqrt(diag(cornerN*cornerN'))*ones(1,3));
    biseN=cornerN-pivotN;                          	biseN=biseN./(sqrt(diag(biseN*biseN'))*ones(1,3));
    % generating the rotation polygon
    k=seed;
    vec1=v(list_N{i,2}(k),:)-v(list_N{i,1}(1),:);
    vec2=ref-v(list_N{i,1}(1),:);
    x2=cross(vec1,vec2); x2=x2./sqrt(x2*x2');
    hinN=-cross(biseN(k,:),x2);
    list_N{i,5}(k,:)=v(list_N{i,1},:)-hinN/hinN(:,3)*v(list_N{i,1},3);
end
n=length(list_N{i,2});
for k=(k+1):(k+n)
    %%
    vec1=v(list_N{i,2}(mod(k-2,n)+1),:)-v(list_N{i,1}(1),:);    x3=vec1./sqrt(vec1*vec1');
    x2=cross(vec1,hinN);  x2=x2./sqrt(x2*x2');
    x1=cross(x2,vec1);  x1=x1./sqrt(x1*x1');
    ang=ang_GP(l(list_RA,2)==list_N{i,1} & l(list_RA,1)==list_N{i,2}(mod(k-2,n)+1));
    pitN=(hinN*x1')*(cos(ang)*x1-sin(ang)*x2)+(hinN*x3')*x3;
    list_N{i,6}(mod(k-1,n)+1,:)=v(list_N{i,1},:)-pitN/pitN(:,3)*v(list_N{i,1},3);
    vec2=cos(ang/2)*x2+sin(ang/2)*x1;
    hinN=-cross(biseN(mod(k-1,n)+1,:),vec2);
    list_N{i,5}(mod(k-1,n)+1,:)=v(list_N{i,1},:)-hinN/hinN(:,3)*v(list_N{i,1},3);
end
k=k+1;
    vec1=v(list_N{i,2}(mod(k-2,n)+1),:)-v(list_N{i,1}(1),:);    x3=vec1./sqrt(vec1*vec1');
    x2=cross(vec1,hinN);  x2=x2./sqrt(x2*x2');
    x1=cross(x2,vec1);  x1=x1./sqrt(x1*x1');
    ang=ang_GP(l(list_RA,2)==list_N{i,1} & l(list_RA,1)==list_N{i,2}(mod(k-2,n)+1));
    pitN=(hinN*x1')*(cos(ang)*x1-sin(ang)*x2)+(hinN*x3')*x3;
    list_N{i,6}(mod(k-1,n)+1,:)=v(list_N{i,1},:)-pitN/pitN(:,3)*v(list_N{i,1},3);
    
plot(list_N{i,5}([1:end 1],1),list_N{i,5}([1:end 1],2),'r+-')
plot(list_N{i,6}([1:end 1],1),list_N{i,6}([1:end 1],2),'r+')
n=length(list_N{i,5}(:,1));
dxf_a=DXF_lines(fid,list_N{i,5}(:,:),[(1:n)', [(2:n)';1]],dxf_a);
dxf_a=DXF_lines(fid,[list_N{i,5}(:,:);list_N{i,6}(:,:)],[(1:n)', n+(1:n)'],dxf_a);
end
%%
end
axis equal
axis tight