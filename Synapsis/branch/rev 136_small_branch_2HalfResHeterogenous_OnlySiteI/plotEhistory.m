close all;clc;
E=load('fE.txt')';
%K=load('fK.txt')';
%HS=load('fHS.txt')';
%f=load('ff.txt')';
%x=load('fx.txt')';
%%
close all;clc;
figure;
hold on;

subplot(2,1,1);
%plot(x(1:end),1:length(x(1:end)));
%plot(x(1:end));

subplot(2,1,2);
hold on;
plot(0:0.1:55-0.01,E(:,floor(size(E,2)*0.3)),'r-+');
plot(0:0.1:55-0.01,E(:,floor(size(E,2)*0.7)),'b-+');
plot(0:0.1:55-0.01,E(:,end),'r-.+');
xlim([0 10]);