im = imread('images/tornetrask2.jpg');

im = double(im)/255;

imagesc(im)

%% 

red = im(:,:,1);
green = im(:,:,2);
blue = im(:,:,3);
figure
imagesc(red)
caxis([0 1])
colormap(gray);

%%
imagesc(im)
caxis([0 1])
%%
[m n k] = size(im);
leftim = im(:, 1:n/2,[3,2,1]);
rightim = im(:, (n/2+1):end, [1, 2, 3]);
%%
imagesc([leftim rightim])
caxis([0 1])
%%
img = [];
im=im/max(im(:));

img(:,:,1)=(im(:,:,1)<=0.25);
img(:,:,2)=(im(:,:,1)>0.25)&(im(:,:,1)<=0.5);
img(:,:,3)=(im(:,:,1)>0.5)&(im(:,:,1)<=0.75);
img(:,:,4)=(im(:,:,1)>0.75);

% Uses colors rgb + cyan
im2 = rgbimage(img);
imagesc(im2);

% Uses colors RGB, all blue is zero so they become black
figure
im3 = rgbimage(img(:,:,1:3))
imagesc(im3)

%%
x = imread('images/tornetrask2.jpg');
x = double(x)/255;
y=x./repmat(sum(x,3),[1,1,3]);
imagesc(y)

%%


%%%%%%%%%% Temperature data %%%%%%%%%%%%%%%


load T_lund.mat
plot(T_lund(:,1), T_lund(:,2))
datetick

%%
t = T_lund(:,1);  Y = T_lund(:,2);  n = length(Y);
X = [ones(n,1) sin(2*pi*t/365) cos(2*pi*t/365)];
beta = regress(Y, X);
eta = Y-X*beta;
plot(t, Y, 'b', t, X*beta, 'r');
datetick