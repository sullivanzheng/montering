function test()
clear;close all;clc;
x = importdata('_200seg7-10,103-106_log.txt');
size(x)
subplot(4,1,1)
hist(x(:,1),1000)
xlim([0,30])
subplot(4,1,2)
plot(x(:,1))
ylim([0,80])
subplot(4,1,3)
hold on
plot(x(:,2)+x(:,3),'r')
ylim([0,50])
subplot(4,1,4)
scatter(1:length(x(:,1)),x(:,1)<0.5 & x(:,2)<20 & x(:,3)<20,0.1);
ylim([0,2])
sum(x(:,1)<0.5 & x(:,2)<20 & x(:,3)<20)/length(x(:,1))
sum(x(:,1)>30)/length(x(:,1))
% A=8;r01=5;
% Er=@(r)-A*exp(-r*r/(r01*r01)/2);
% Er(5)-Er(0.05)
% Hr=@(r)log(4*pi*r*r);
% Hr(5)-Hr(0.05)
% x=1:10
% iff(x<5,ones(10,1),0)
% PI=pi;
% r0=5;r01=5;a0=2.0/180.*PI;R0=20.0/180.*PI;
% A=8;B=17-A;C=5; 
% re=5;ree=0.05;
% 
% Ar=(-(A*exp(-re*re/(r01*r01)/2))+A*exp(-ree*ree/(r01*r01)/2))/ ...
%     (-log(4*PI*ree*ree)+log(4*PI*re*re))
% Ar=1;
% 
% E1=@(r)-A.*(iff(r>re,exp(-r.*r./(r01.*r01)/2),exp(-re.*re./(r01.*r01)/2)));
% E11=@(r)-A.*exp(-r.*r./(r01.*r01)/2);
% E2=@(r)+Ar.*(iff(r>re, ...
% 	log(4*PI.*re.*re), ...
% 	iff(r>=ree,log(4*PI.*r.*r), ...
% 	log(4*PI.*ree.*ree))));
% 
% x=0:0.1:80;
% hold on;
% plot(x,E2(x),'b');
% plot(x,E11(x),'c');
% plot(x,E1(x),'r');
% plot(x,E1(x)+E2(x),'k');
% figure;
% hold on;
% plot(x,4*pi.*x.*x.*exp(-E11(x)),'r');
% plot(x,4*pi.*x.*x.*exp(-E1(x)-E2(x)));
end

function iff=iff(in,a,b)
iff=zeros(size(in));
for i=1:length(in)
    if in(i)
        if length(a)==1
            iff(i)=a;
        else
            iff(i)=a(i);
        end
    else
        if length(b)==1
            iff(i)=b;
        else
            iff(i)=b(i);
        end
    end
end
end

