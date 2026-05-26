function [file_deconv] = f_deconvolution(file_profile)
file_deconv = [];
vars = load(file_profile);
if ~isfield(vars, 'data_S') || isempty(vars.data_S)
   return
end
data_S = vars.data_S;
Depth = vars.Depth;
S = vars.S;
noise = vars.noise;
data_So = data_S;
Depth_o = Depth;
ind_o = find(Depth_o >= -3);
Depth_o = Depth_o(ind_o);
So = S(ind_o,:);
Sc = So;
semilogy(Depth_o,So)
ind_water = find(Depth_o >= 0);
Depth_w =  Depth_o(ind_water);
So1 = Sc(ind_water,:);
dots = length(ind_water);
num_profile = size(So1,2);
load('SystemResponse_try1.mat');
ind_land = find(Depth == 0);
Depth_land = Depth(ind_land : ind_land+dots-1);
% S_land = SR_sumNorm(ind_land : ind_land+dots-1);
S_land = SR_maxNorm(ind_land : ind_land+dots-1)';
hold on
semilogy(Depth_land,S_land)
xlim([-15 20])
%% ?1?7?1?7?0?7?0?5?1?7?1?7?1?7?1?7?1?7?1?9?1?7?1?7?1?7?1?7?0?4?1?7?1?7?1?7 F * Sc(correct) = So(observed)
Sc1 = zeros(size(So1));
for i = 1:num_profile
    F = zeros(dots, dots);
    for j = 1:dots
        F(j,1:j) = S_land(j:-1:1);
    end
    Sc1(:,i) = F\So1(:,i);
    % ?1?7?0?4?1?7?1?7?1?7?1?7?1?7?0?2?1?7?1?7?1?7?1?7?0?2?0?3?1?7?1?7?1?7?1?7<0?1?7?1?7?1?7?1?7?1?7

end
Sc(ind_water,:) = Sc1;
%% 折射率纠正
 Depth_oair = Depth_o;
% ind_depth = find(Depth_o > 0);
% theta1 = 0; 
% n2 = 1.32; 
% n1 = 1;   
% theta2 = asin(n1*sin(theta1)/n2);
% dtheta = theta1 - theta2;
% S = Depth_o(ind_depth) / cos(theta1);
% R = n1*S/n2;
% P = (R.^2 + S.^2 - 2.*R.*S.*cos(dtheta)).^(1/2);
% alpha = asin(R*sin(dtheta)./P);
% beta = (pi/2 - theta1) - alpha;
% dZ = P.*sin(beta);
% Depth_o(ind_depth) = Depth_o(ind_depth) - dZ;
%% 平滑
% start = min(find(Depth_o>=1.5)); % 1.5 m
% for i = 1:num_profile
%     Sc(start:end,i) = smoothdata(Sc(start:end,i),'movmean',5,'includenan');
% end

f = figure('Visible','on');
set(f,'units','normalized','position',[0.05 0.2 0.5 0.6]);  
set (gca,'position',[0.15  0.15  0.7  0.7] );   
for i = 1:1:num_profile
    semilogy(Depth_oair, So(:,i)/max(So(:,i)),'LineWidth',1,'Color', '#0072BD','LineWidth',0.5);
    hold on
    semilogy(Depth_oair, Sc(:,i)/max(Sc(:,i)),'LineWidth',1,'Color', '#D95319','LineWidth',0.5);
    legend('So','Sc','FontSize',10,'Location','northeast');
    xlabel('Depth','FontSize',10);
    ylabel('Normalized Signal','FontSize',10);
end

semilogy(Depth_land, S_land/max(S_land),'Color', '#7E2F8E','LineWidth',1);
ylim([10^-6 10^0]);

Sc(Sc<=0)=nan;
a = dir(file_profile);
file_deconv = [a.folder,'\','deconv_',a.name];
save(file_deconv,'Depth_o','So','Sc','data_So','noise'); %theta1?1?7?1?7?1?7?1?7?1?7


end



