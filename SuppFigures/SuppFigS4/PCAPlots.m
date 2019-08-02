
clear
close all

%% AR, ESR1, and KIT

pca = csvread('totPCA.csv');
pca_obs = csvread('totPCAobs.csv');

figure
scatter3(pca(:,1),pca(:,2),pca(:,3),10,[0.5 0.5 0.5],'MarkerEdgeAlpha',0.25,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerFaceAlpha',0.25)
hold on
scatter3(pca_obs(:,1),pca_obs(:,2),-pca_obs(:,3),500,'r.')
scatter3(0,0,0,'k+')
xlabel('PC 1','FontSize',22,'FontWeight','bold','FontName','Arial')
ylabel('PC 2','FontSize',22,'FontWeight','bold','FontName','Arial')
zlabel('PC 3','FontSize',22,'FontWeight','bold','FontName','Arial')
view(-13.5,15)
set(get(gca,'ylabel'), 'Rotation',-50)
set(get(gca,'xlabel'), 'Rotation',0)
set(gcf,'Position',[50 50 650 600])

xticks([-0.3 0 0.3 0.6])
yticks([-0.3 0 0.3])
zticks([-0.3 0 0.3])

title('Principle Components Analysis','FontSize',28,'FontWeight','bold','FontName','Arial')

% %% ABL
% 
% pca = csvread('ABLPCA.csv');
% pca_obs = csvread('ABLPCAobs.csv');
% 
% figure
% scatter3(pca(:,1),pca(:,2),pca(:,3),[],[0.5 0.5 0.5],'MarkerEdgeAlpha',0.25,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerFaceAlpha',0.25)
% hold on
% scatter3(pca_obs(:,1),pca_obs(:,2),pca_obs(:,3),1000,'r.')
% scatter3(0,0,0,'+')
% xlabel('Principle Component 1')
% ylabel('Principle Component 2')
% zlabel('Principle Component 3')
% 
% %% EGFR
% 
% pca = csvread('EGFRPCA.csv');
% pca_obs = csvread('EGFRPCAobs.csv');
% 
% figure
% scatter3(pca(:,1),pca(:,2),pca(:,3),[],[0.5 0.5 0.5],'MarkerEdgeAlpha',0.25,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerFaceAlpha',0.25)
% hold on
% scatter3(pca_obs(:,1),pca_obs(:,2),pca_obs(:,3),1000,'r.')
% xlabel('Principle Component 1')
% ylabel('Principle Component 2')
% zlabel('Principle Component 3')