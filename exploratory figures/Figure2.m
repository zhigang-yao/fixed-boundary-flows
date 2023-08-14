%% Code for the heatmap in Figure 2
%
% xiayq0121@zufe.edu.cn
% refered to Z. Yao, Y. Xia and Z. Fan, Random Fixed Boundary Flows, Journal of the American Statistical Association

% Generate data
sigma = 0.05;
r = 1/pi; theta = 0:pi/500:pi;
x0 = [r*cos(theta);r*sin(theta)];

rng(1);
noise = randn(size(x0))*sigma;

x = x0+noise;

xbar = mean(x,2);
 

% Calculate the first eigenvalues
n = 100; 
m = 200;
s = linspace(-0.5,0.5,m);
t = linspace(-0.1,0.5,n);
lambda = zeros(m,n);

for i = 1:m
    for j = 1:n
        H = bsxfun(@minus,x,[s(i);t(j)]);
        Sigma = (H*H')/size(x,2);
        e = eig(Sigma);
        lambda(i,j) = max(e);
    end
end

[S,T] = meshgrid(s,t);


%% Code for Figure2(a)
figure(1)
surf(S',T',lambda,'FaceColor','interp','LineStyle','None')
colorbar;
hold on;

scatter3(x(1,:),x(2,:),ones(1,size(x,2))*0.5,10, 'k', 'filled');
scatter3(xbar(1),xbar(2),0.5,20, 'r');

H = bsxfun(@minus,x,xbar);
Sigma = (H*H')/size(x,2);
[V,D] = eig(Sigma);
d = diag(D);
[~, idx] = max(d);
v = V(:,idx);
quiver3(xbar(1),xbar(2),0.5,-v(1),-v(2),0,0.2,'r','Linewidth',2)
view([0,0,1])

xlim([-0.5,0.5])
ylim([-0.1,0.5])
axis off;

%% Code for Figure2(b)
figure(2)

surf(S',T',lambda,'FaceColor','interp','LineStyle','None')
colorbar;
hold on;

scatter3(x(1,:),x(2,:),ones(1,size(x,2))*0.5,10, 'k', 'filled');
scatter3(x0(1,1),x0(2,1),0.5,20, 'r');
scatter3(x0(1,end),x0(2,end),0.5,20, 'r');

H = bsxfun(@minus,x,x0(:,end));
Sigma = (H*H')/size(x,2);
[V,D] = eig(Sigma);
d = diag(D);
[~, idx] = max(d);
v = V(:,idx);
quiver3(x0(1,end),x0(2,end),0.5,-v(1),-v(2),0,0.2,'r','LineWidth',2)

view([0,0,1])

xlim([-0.5,0.5])
ylim([-0.1,0.5])
axis off;
