%linear programming
clear;
close all;
%load sc50b;
%the problem setup is
% A*x <= b
% Aeq*x = beq
% x >=lb
% argmin f'x
% find x
Aeq = [1 0 -1 0 0 0 0 0 ;
       0 1 0 -1 -1 1 0 0 ;
       0 0 1 1 0 -1 -1 0;
       0 0 0 0 1 0 0 -1];
Aeq = sparse(Aeq);
beq = zeros(4,1);

i=[1,1,2,3,3,4,4,5];
j=[2,3,4,4,5,3,6,6];
v=[1,2,3,4,5,6,7,8];
Imx = sparse(i,j,v,6,6);
ub = [4 2 3 1 3 1 2 4]';
lb = zeros(8,1);
f = -[0 0 0 0 0 0 1 1]';

[x,fval,exitflag,output] = linprog(f,[],[],Aeq,beq,lb,ub);

cutset = mincut_dfs(1,x,ub,Imx,[]);

%% try dfs


%% study image
im1 = imread('H:\Code\Matlab\metric trainning\data\woman_g_sweater_gray_1.tiff');
im2 = imread('H:\Code\Matlab\metric trainning\data\woman_g_sweater_gray_2.tiff');
%om1 = im1(:,41:48);
%om2 = im2(:,1:8);
om1 = double(im1(end-3:end,end-3:end));
om2 = double(im2(end-3:end,1:4));
%cy = diff(:,1:end-1)+diff(:,2:end);
%cx = diff(1:end-1,:)+diff(2:end,:);
diff = abs(om1-om2);
diff = diff';
diff = diff(:);
[H,W] = size(om1);
dim = 2*(W-1)*H+2*W*(H-1)+2*H;% 2*numel(cy)+2*numel(cx)+size(om1,1)+size(om2,1);%number of edges, remeber between pixels there are two edges (undirected)
%the order of edges are column first
%Aeq = zeros(numel(om1),dim);
skip = size(om1,1); 
step = size(om1,2);

%% build connection matrix and capacity matrix
x = zeros(dim,1);%x coordinate
y = zeros(dim,1);%y coordinate
v = zeros(dim,1);%value at (x,y)
cost =zeros(dim,1);
cost = cost(:);
ind = [1:dim]';
count=skip+1;
%W=4;
%H=4;
for i=1:size(om1,1)
    for j=1:size(om1,2)
        %east
        pos = (i-1)*W+j;
        if (j~=1)
            wpos = (i-1)*W+j-1;
        else
            wpos = -1;
        end
        if (j~= W)
            epos = (i-1)*W+j+1;
        else
            epos = -1;
        end
        if (i~=1)
            npos = (i-2)*W+j;
        else 
            npos=-1;
        end
        if (i~=H)
            spos = (i)*W+j;
        else
            spos = -1;
        end
        if epos~=-1
            x(count)=pos;
            y(count)=epos;
            v(count)=1;
            cost(count) = diff(pos)+diff(epos);
            %ind(count)=count;
            count=count+1;
        end
        if wpos~=-1
            x(count)=pos;
            y(count)=wpos;
            v(count)=1;
            cost(count) = diff(pos)+diff(wpos);
            %ind(count)=count;
            count=count+1;
        end
        if npos~=-1
            x(count)=pos;
            y(count)=npos;
            v(count)=1;
            cost(count) = diff(pos)+diff(npos);
            %ind(count)=count;
            count=count+1;
        end
        if spos~=-1
            x(count)=pos;
            y(count)=spos;
            v(count)=1;
            cost(count) = diff(pos)+diff(spos);
            %ind(count)=count;
            count=count+1;
        end
    end
end
for i=1:skip
    x(i) = 0;
    y(i) = (i-1)*W+1;
    v(i) = 1;
    cost(i) = 10000;
    x(end-skip+i)=i*W;
    y(end-skip+i)=H*W+1;
    v(end-skip+i)=1;
    cost(end-skip+i)= 10000;
    
end
x=x+1;
y=y+1;
Cmx = sparse(x,y,cost,H*W+2,H*W+2);
figure
imagesc((Cmx));
Imx = sparse(x,y,ind,H*W+2,H*W+2);
%% build constraint mx A from Cmx
%A= zeros(W*H,dim);
x = zeros(dim,1);%x is id of node
y = zeros(dim,1);%y is the id of edge, which is saved in Imx, for corresponding node (i,j)
v = zeros(dim,1);%value is the 1 for incoming, -1 for outgoing
count = 1;
for i=2:size(Cmx,1)-1
    %incoming 
    %[a,b,c] = find(Imx(:,i));
    ind = Imx(find(Imx(:,i)),i);
    for j=1:length(ind)
        y(count) = ind(j);
        x(count) = i;
        v(count) = 1;
        count=count+1;
    end
    %outgoing
    ind = Imx(i,find(Imx(i,:)));
    for j=1:length(ind)
        y(count) = ind(j);
        x(count) = i;
        v(count) = -1;
        count=count+1;
    end
end
x=x-1;%eleminate the all zero line (first line for S)
A = sparse(x,y,v,H*W,dim);
figure
imagesc(A);
%% solve this LP
f = zeros(dim,1);
f(1:skip)=-1;
b = zeros(H*W,1);
lb = zeros(dim,1);
ub = cost;
[x,fval,exitflag,output] = linprog(f,[],[],A,b,lb,ub);
%% use DFS find min-cut
cutset = mincut_dfs(1,x,cost,Imx,[]);

%% now test on bigger image
close all;
clear all;
im1 = imread('H:\Code\Matlab\metric trainning\data\woman_g_sweater_gray_1.tiff');
im2 = imread('H:\Code\Matlab\metric trainning\data\woman_g_sweater_dark_2.tiff');
om1 = double(im1(:,41:48));
om2 = double(im2(:,1:8));
[H,W] = size(om1);
diff = abs(om1-om2);
%change diff to guide seam
%diff(:,4)=20;
%diff(:,5)=20;
x = 0:47;
theta = 0;%pi/12;
cx = (H-1)/2;
cy = (W-1)/2;
%b = -cx*tan(theta);
y = (x-cx)*((H-1)*tan(theta)/2)/(H-1);
y = y+cy;
y = round(y)+1;
x = x+1;
%y = round((x-1)*(W-1)/(H-1))+1;
for i=1:48
    if (y(i)==0 || y(i)>W)
        continue;
    end
    diff(x(i),y(i))=diff(x(i),y(i))/10;
end
figure 
imagesc(diff);
%
diff = diff';
diff = diff(:);

dim = 2*(W-1)*H+2*W*(H-1)+2*H;% 2*numel(cy)+2*numel(cx)+size(om1,1)+size(om2,1);%number of edges, remeber between pixels there are two edges (undirected)
%the order of edges are column first
%Aeq = zeros(numel(om1),dim);
skip = size(om1,1); 
step = size(om1,2);

% build connection matrix and capacity matrix
x = zeros(dim,1);%x coordinate
y = zeros(dim,1);%y coordinate
v = zeros(dim,1);%value at (x,y)
cost =zeros(dim,1);
cost = cost(:);
ind = [1:dim]';
count=skip+1;
%W=4;
%H=4;
for i=1:size(om1,1)
    for j=1:size(om1,2)
        %east
        pos = (i-1)*W+j;
        if (j~=1)
            wpos = (i-1)*W+j-1;
        else
            wpos = -1;
        end
        if (j~= W)
            epos = (i-1)*W+j+1;
        else
            epos = -1;
        end
        if (i~=1)
            npos = (i-2)*W+j;
        else 
            npos=-1;
        end
        if (i~=H)
            spos = (i)*W+j;
        else
            spos = -1;
        end
        if epos~=-1
            x(count)=pos;
            y(count)=epos;
            v(count)=1;
            cost(count) = diff(pos)+diff(epos);
            %ind(count)=count;
            count=count+1;
        end
        if wpos~=-1
            x(count)=pos;
            y(count)=wpos;
            v(count)=1;
            cost(count) = diff(pos)+diff(wpos);
            %ind(count)=count;
            count=count+1;
        end
        if npos~=-1
            x(count)=pos;
            y(count)=npos;
            v(count)=1;
            cost(count) = diff(pos)+diff(npos);
            %ind(count)=count;
            count=count+1;
        end
        if spos~=-1
            x(count)=pos;
            y(count)=spos;
            v(count)=1;
            cost(count) = diff(pos)+diff(spos);
            %ind(count)=count;
            count=count+1;
        end
    end
end
for i=1:skip
    x(i) = 0;
    y(i) = (i-1)*W+1;
    v(i) = 1;
    cost(i) = 10000;
    x(end-skip+i)=i*W;
    y(end-skip+i)=H*W+1;
    v(end-skip+i)=1;
    cost(end-skip+i)= 10000;
    
end
x=x+1;
y=y+1;
Cmx = sparse(x,y,cost,H*W+2,H*W+2);
figure
imagesc((Cmx));
Imx = sparse(x,y,ind,H*W+2,H*W+2);
% build constraint mx A from Cmx
%A= zeros(W*H,dim);
x = zeros(dim,1);%x is id of node
y = zeros(dim,1);%y is the id of edge, which is saved in Imx, for corresponding node (i,j)
v = zeros(dim,1);%value is the 1 for incoming, -1 for outgoing
count = 1;
for i=2:size(Cmx,1)-1
    %incoming 
    %[a,b,c] = find(Imx(:,i));
    ind = Imx(find(Imx(:,i)),i);
    for j=1:length(ind)
        y(count) = ind(j);
        x(count) = i;
        v(count) = 1;
        count=count+1;
    end
    %outgoing
    ind = Imx(i,find(Imx(i,:)));
    for j=1:length(ind)
        y(count) = ind(j);
        x(count) = i;
        v(count) = -1;
        count=count+1;
    end
end
x=x-1;%eleminate the all zero line (first line for S)
A = sparse(x,y,v,H*W,dim);
figure
imagesc(A);
% solve this LP
f = zeros(dim,1);
f(end-skip+1:end)=-1;
b = zeros(H*W,1);
lb = zeros(dim,1);
ub = cost;
[x,fval,exitflag,output] = linprog(f,[],[],A,b,lb,ub);
% use DFS find min-cut
cutset = mincut_dfs(1,x,cost,Imx,[]);

% now build the image
%om1 = double(im1(:,41:48));
%om2 = double(im2(:,1:8));
cutset(1)=[];%delete the S
cutset=cutset-1;%renumbering
mask = zeros(H,W);
x = floor((cutset-1)/W)+1;
y = mod(cutset-1,W)+1;
for i=1:length(x)
    mask(x(i),y(i))=1;
end
imagesc(mask)

overlap = om1.*mask + om2.*(1-mask);
lapf=[-1,1];
seam = conv2(mask,lapf,'same');
seam3 = zeros(size(seam,1),size(seam,2),3);
seam3(:,:,1)=seam;
%seam3(:,:,2)=seam;
%seam3(:,:,3)=seam;
rst = zeros(size(im1,1),size(im1,2)*2-8);
rst(:,1:40)=im1(:,1:40);
rst(:,end-40+1:end)=im2(:,end-40+1:end);
rst(:,41:48)= overlap;
rgbrst = gray2rgb(rst);
rgbrst_seam = rgbrst;
rgbrst_seam(:,41:48,:)=rgbrst_seam(:,41:48,:)+seam3*255;
figure
subplot(2,1,1)
imshow(uint8(rgbrst));
subplot(2,1,2)
imshow(uint8(rgbrst_seam));