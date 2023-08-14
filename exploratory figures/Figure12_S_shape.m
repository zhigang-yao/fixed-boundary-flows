%% Code for visualizing the S-shape data in Figure 12
%
% xiayq0121@zufe.edu.cn
% refered to Z. Yao, Y. Xia and Z. Fan, Random Fixed Boundary Flows, Journal of the American Statistical Association

%% Load S-shape data
sigma = 0.1;
theta = pi:-pi/400:-pi;
x0 = [theta;-sin(theta)];
noise = randn(size(x0))*sigma;
x = x0+noise;

xbar = mean(x,2);

%% calculate v1 and v2
H = bsxfun(@minus,x,xbar);
Sigma = (H*H')/size(x,2);
[V,D] = eig(Sigma);
d = diag(D);
[~,idx] = sort(d,'descend');
v1 = V(:,idx(1));
v2 = V(:,idx(2));
if v1(1) < 0; v1 = -v1; end
if v2(2) < 0; v2 = -v2; end

%%
%rotate the coordinate system, taking v1 as the horizontal direction 
% and v2 as the vertical direction
x = [v1'*x; v2'*x];
xbar = mean(x,2);


%% Code for Figure 12(c)
n = 10; 
m = 10;
s = linspace(min(x(1,:)),max(x(1,:)),m);
t = linspace(min(x(2,:)),max(x(2,:)),n);
W = zeros(m, n, 2);

for i = 1:m
    for j = 1:n
        H = bsxfun(@minus,x,[s(i);t(j)]);
        Sigma = (H*H')/size(x,2);
        [V,D] = eig(Sigma);
        d = diag(D);
        [~, idx] = max(d);
        v = V(:,idx);
        v = v*sign(v(1));
        W(i,j,:) = v;
    end
end

[S,T] = meshgrid(s,t);

figure(1)

quiver(S', T', W(:,:,1), W(:,:,2)); hold on;
scatter(x(1,:),x(2,:), 6, 'k', 'filled');
scatter(xbar(1),xbar(2),20,'k','LineWidth',1);

axis off;

%% Code for Figure 12(d)
n = 100; 
m = 200;
s = linspace(min(x(1,:)),max(x(1,:)),m);
t = linspace(min(x(2,:)),max(x(2,:)),n);
lambda = zeros(m, n);

for i = 1:m
    for j = 1:n
        H = bsxfun(@minus,x,[s(i);t(j)]);
        Sigma = (H*H')/size(x,2);
        [V,D] = eig(Sigma);
        d = diag(D);
        [~, idx] = max(d);
        v = V(:,idx);
        lambda(i,j) = abs(v(2));
    end
end

[S,T] = meshgrid(s,t);

figure(2)

surf(S',T',lambda,'FaceColor','interp','LineStyle','None'); 
hold on;
colorbar;
scatter3(x(1,:),x(2,:),ones(1,size(x,2))*0.5,8, 'k', 'filled');
scatter3(xbar(1),xbar(2),0.5,20, 'k', 'LineWidth',2);
view([0,0,1])

xlim([min(x(1,:)),max(x(1,:))])
ylim([min(x(2,:)),max(x(2,:))])

axis off;
