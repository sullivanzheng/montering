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
plot(0:0.1:55-0.01,E(:,floor(size(E,2)*0.5)),'b-+');
plot(0:0.1:55-0.01,E(:,end),'r-.+');
xlim([0 10]);
figure;
bar(0:0.1:55-0.01,E(:,floor(size(E,2)*0.5)));
xlim([0 10]);

