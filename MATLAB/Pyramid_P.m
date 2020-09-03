function [Pyramid,Expanded_P] = Pyramid_P(i,v,l,ang_GP,list_RA,list_P,list_N,list_P_in_N)
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
    
Pyramid=list_P{i,5};
Expanded_P=list_P{i,6};
end

