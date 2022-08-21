clc
clear all
hold off
%me6151 project plotting
T = load('F:\Studies\Mtech\Dept of Aerospace Engg\Sem2\ME6151\Project\readmesh\phi.txt');
x = T(:,1);
y = T(:,2);
z = T(:,3);
line = load('F:\Studies\Mtech\Dept of Aerospace Engg\Sem2\ME6151\Project\readmesh\line.txt');
figure(1)
plot(line(:,1),line(:,3),'.')
figure(2)
scatter(x,y,[],z,'.')
colorbar
colormap jet
figure(3)
scatter3(x,y,z,[],z,'.')
colorbar
colormap jet

        