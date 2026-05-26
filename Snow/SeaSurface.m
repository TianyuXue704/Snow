function [x, ind_peak] = SeaSurface(ind_peak0)
    %%
    %ind_peak = filloutliers(ind_peak0,'spline','quartiles','ThresholdFactor',1.5); %发现异常值并插补
    x0 = 1:length(ind_peak0);
    %剔除ind_peak0为1的数据点
    %原因： 部分bin没有找到光子,因此检测到的ind_peak0 = 1
    ind_one  = find(ind_peak0 == 1);
    ind_peak1 = ind_peak0;
    ind_peak1(ind_one) = [];
    x1 = x0;
    x1(ind_one) = [];
    f = figure('Visible','on');
    sgtitle('筛选海面峰值中正确的值','FontSize',16);
    set(f,'units','normalized','position',[0.1 0.1 0.8 0.8]);  %设置figure的窗口大小
    set (gca,'position',[0.1  0.1  0.85  0.8] );   %坐标轴在figure中的左边界，下边界，宽度，高度
    subplot(2,2,1);
    scatter(x0,ind_peak0,5,[0,0.447,0.741],'filled');hold on;
    scatter(x1,ind_peak1,10,[0.85,0.33,0.10],'LineWidth',1);
    set(gca,'FontSize',12);  %坐标刻度大小
    title('步骤一：剔除没有光子的bin（ind_peak0=1)','Interpreter', 'none','FontSize',16);
    xlabel('bin序号（时间方向）','FontSize',14);
    ylabel('峰值点序号（高度方向)','FontSize',14);
    
    %%
    %差分N遍,留下diff<阈值的海面数据点
    Threshold_diff = 50;  %30
    x_old = x1;
    x_new = [];
    ind_peak_old = ind_peak1;
    diffTime = 0;
    while(length(x_new) ~= length(x_old) & diffTime <2)
        if(diffTime >0)
            x_old = x_new;
            ind_peak_old = ind_peak_new;
            clearvars ind_peak_new x_new
        end
        diff = zeros(1,length(x_old)-1);
        for i = 1:length(diff)
            diff(i) = ind_peak_old(i+1) - ind_peak_old(i);
        end
        ind_remain = find(abs(diff) <= Threshold_diff);
        ind_peak_new = ind_peak_old([1 ind_remain+1]); %这里先把第1个点考虑进去
        x_new = x_old([1 ind_remain+1]);
        diffTime = diffTime + 1;
        clearvars diff ind_remain ind_peak_old
    end
    
    subplot(2,2,2);
    scatter(x0,ind_peak0,5,[0,0.447,0.741],'filled');hold on;
    scatter(x_new,ind_peak_new,10,[0.85,0.33,0.10],'LineWidth',1);
    set(gca,'FontSize',12);  %坐标刻度大小
    title('步骤二：差分：剔除突变的峰值点','FontSize',16);
    xlabel('bin序号（时间方向）','FontSize',14);
    ylabel('峰值点序号（高度方向)','FontSize',14);
    %%
    if(length(ind_peak_new)<=1)
        x = [];
        ind_peak = [];
        return;
    end
    %在差分基础上根据大的diff之间的数据点个数筛选,进行N次
    for i = 1:length(x_new)-1
        diff(i) = ind_peak_new(i+1) - ind_peak_new(i);
    end
    % Threshold_num = 15000;
    % ind = find(abs(diff) > Threshold_diff);
    % if(length(ind) > 0)  %假如ind = []，说明数据都很连续，不需要进行这步筛除工作
    %     num_btw_diff = zeros(length(ind)+1,1); %N个点将数据分为N+1段
    %     num_btw_diff(1) = ind(1);
    %     num_btw_diff(end) = length(x_new)-ind(end);
    %     for i = 2:length(num_btw_diff)-1
    %         num_btw_diff(i) = ind(i) - ind(i-1);
    %     end
    %     ind_seg = find(abs(num_btw_diff) > Threshold_num);  %useful data segment有用的数据段
    %     ind_data = [];
    %     for i = 1:length(ind_seg)
    %         if(ind_seg(i) > 1 & ind_seg(i) < length(num_btw_diff))
    %             ind_data = [ind_data   ind(ind_seg(i)-1)+1 : ind(ind_seg(i))-1];
    %         else if(ind_seg(i) == 1)
    %                 ind_data = [ind_data   1 : ind(1)];
    %             else
    %                 ind_data = [ind_data  ind(ind_seg(i)-1)+1 : length(x_new)];
    %             end
    %         end
    %     end
        x_result = x_new;
        ind_peak_result = ind_peak_new;
    % else
    %     x_result = x_new;
    %     ind_peak_result = ind_peak_new;
    % end
    subplot(2,2,3);
    scatter(x0,ind_peak0,5,[0,0.447,0.741],'filled');hold on;
    scatter(x_result,ind_peak_result,10,[0.85,0.33,0.10],'LineWidth',1);
    set(gca,'FontSize',12);  %坐标刻度大小
    title('步骤三：差分：保留连续一段波动小的峰值点','FontSize',16);
    xlabel('bin序号（时间方向）','FontSize',14);
    ylabel('峰值点序号（高度方向)','FontSize',14);

    x = x0;
    ind_peak = nan(size(ind_peak0));
    ind_peak(x_result) = ind_peak_result;    

 
end

