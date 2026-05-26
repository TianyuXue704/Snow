function alpha_Perturbation = f_solve_IS_HK(S, Depth)
%准备基本参数
SignalT = S';
[Num_of_Seg,Segment]=size(SignalT);
dz = Depth(2) - Depth(1);


dis = double(Depth);
depth = double(Depth);
%%
%基本信号信息
surface = find(dis == 0);

Signal_cor = abs(SignalT);
logSignal_cor = log(abs(Signal_cor));

%确定扰动法和Fernald法反演参考点,同时也是确定噪声范围
for i = 1:Num_of_Seg
   refPointArray = find(isnan(SignalT(i,80:Segment))) + 80;% 可能廓线前面有几个点也是nan,需要排除（这里排除约10m前）
%      refPointArray = find(SignalT(i,80:Segment)<noise(i)) + 80;
    if(isempty(refPointArray))
        refPoint(i) = min(90, length(Depth)); % 10m
        %refPoint(i) = Segment - 50;
        
    else
        refPoint(i) = min(refPointArray(1)-20, Segment);%Fernald的参考点，取信号在水面下第一次下降到0以下的深度，往水面前取一定的距离
    end
end


%% 扰动法
%
beta_Perturbation = zeros(Num_of_Seg,Segment);
surRefer = 20;

%对包含衰减信息的信号进行拟合求斜率
for i = 1:1:Num_of_Seg
    Signal_cor_fit = Signal_cor(i,surface+surRefer:refPoint(i));
    logSignal_cor_fit = logSignal_cor(i,surface+surRefer:refPoint(i));
    dis_fit = dis(surface+surRefer:refPoint(i));
    
    Y = logSignal_cor_fit';
    X = [ones(length(Y),1) dis_fit']; 
    coefficients = (X'*X)\(X'*Y);
    alpha_Perturbation(i) = -coefficients(2)/2;
    
    
end
beta_Perturbation(find(beta_Perturbation == 0)) = nan;
%经验公式
% bbp_Perturbation = 6.43*(beta_Perturbation-2.4*0.0001);
%% 

end

