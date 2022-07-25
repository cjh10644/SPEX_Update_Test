clc
clear
close all
list = dir('../LPnetlib/*.txt');

% for i = 1:length(list) 
% great cases [65]
% slowly worsen [3 67 68]
% sudden jump worse[90 51  52]
% sudden jump but ok[1 61 17]
LP_index =[52  67  65  25  1   3   17 51  61  68  90];%interesting cases
max_iter =[100 60 100 100 100  45 100 100 100 100 100];
ylb_shift=[ 4  3  4   4   4    4    4   4   4   4   4];
% LP_index =[52 65  61  68];%used cases
% max_iter =[100 100 100 100  ];
% ylb_shift=[ 4   4   4   4   ];
for index =1:length(LP_index)
    i = LP_index(index);
    list(i).name
    fRead = fopen(strcat('../LPnetlib/',list(i).name), 'r');
    A = fscanf(fRead, '%f %f %f %d %d %d %d %f %d %d %d %d %f %d %d %d %d',[17, Inf]);
    fclose(fRead);
    
    num_iter = max_iter(index);
    x_axis = 1:size(A,2);
    figure(i)
    set(gcf,'Position',[-90 3 560 560])
    ax1=subplot(4,1,1);
    title (list(i).name,'interpreter','none');
    hold(ax1,'on')
    plot(x_axis,A(4,:)+A(6,:),'r--*');
    plot(x_axis(2:end),A(9,2:end)+A(11,2:end),':.','color',[0 0.5 0]);
    plot(x_axis,A(14,:)+A(16,:),'b-.o');
    if i == 52
        xline(14,'-','color',[194 197 204]/255);
        xline(82,'-','color',[194 197 204]/255);
    end
    xlim([0,num_iter]);
    xticks(0:5:num_iter);
    xticklabels({''})
    max_val=max([max(A(4,:)+A(6,:)),...
                 max(A(9,2:end)+A(11,2:end)),...
                 max(A(14,:)+A(16,:))]);
    min_val=min([min(A(4,:)+A(6,:)),...
                 min(A(9,2:end)+A(11,2:end)),...
                 min(A(14,:)+A(16,:))]);
    expon = floor(log10(max_val))-1;
    ub = ceil(max_val/10^expon);
    lb = floor(min_val/10^expon);
    if (mod(ub+lb,2) ~= 0)
        if (ub*10^expon-max_val > min_val-lb*10^expon)
            lb = lb-1;
        else
            ub = ub+1;
        end
    end
    ylim([lb,ub]*10^expon);
    yticks([lb,(lb+ub)/2,ub]*10^expon);
    yticklabels([lb,(lb+ub)/2,ub]);
    ylb1=ylabel("nnz(\times10^"+expon+")");
    legend('DLU','lb','LUU','Orientation','horizontal','Location','best');
    box on;
    
    %total number of bit size

    ax2=subplot(4,1,2);
    hold(ax2,'on')
    plot(x_axis,A(5,:)+A(7,:),'r--*');
    plot(x_axis(2:end),A(10,2:end)+A(12,2:end),':.','color',[0 0.5 0]);
    plot(x_axis,A(15,:)+A(17,:),'b-.o');
    if i == 52
        xline(14,'-','color',[194 197 204]/255);
        xline(82,'-','color',[194 197 204]/255);
    end
    xlim([0,num_iter]);
    xticks(0:5:num_iter);
    xticklabels({''});
    max_val=max([max(A(5,:)+A(7,:)),...
                 max(A(10,2:end)+A(12,2:end)),...
                 max(A(15,:)+A(17,:))]);
    min_val=min([min(A(5,:)+A(7,:)),...
                 min(A(10,2:end)+A(12,2:end)),...
                 min(A(15,:)+A(17,:))]);
    expon = floor(log10(max_val))-1;
    ub = ceil(max_val/10^expon);
    lb = floor(min_val/10^expon);
    if (mod(ub+lb,2) ~= 0)
        if (ub*10^expon-max_val > min_val-lb*10^expon)
            lb = lb-1;
        else
            ub = ub+1;
        end
    end
    ylim([lb,ub]*10^expon);
    yticks([lb,(lb+ub)/2,ub]*10^expon);
    yticklabels([lb,(lb+ub)/2,ub]);
    ylb2=ylabel("bits(\times10^"+expon+")");
%     legend('DLU','lb','LUU','Orientation','horizontal');
    box on;
    
    % factorization time

    ax3=subplot(4,1,3);
    hold(ax3,'on')
    plot(x_axis,A(3,:),'r--*');
    plot(x_axis(2:end),A(8,2:end),':.','color',[0 0.5 0]);
    plot(x_axis,A(13,:),'b-.o');
    if i == 52
        xline(14,'-','color',[194 197 204]/255);
        xline(82,'-','color',[194 197 204]/255);
    end
    xlim([0,num_iter]);
    xticks(0:5:num_iter);
    xticklabels({''})
    ax3.YAxis.Exponent = 0;
    if i == 61
        ylim([0,130]);
        yticks([0,65,130]);
    elseif i == 52 
        ylim([0,100]);
    elseif i ==65
        ylim([0,700]);
        yticks([0,350,700]);
    elseif i == 68
        ylim([0,0.03]);
        yticks([0,0.01,0.02,0.03]);
        yticklabels({'0','0.01','0.02','0.03'});
    end
    ylb3=ylabel("time (sec)");
    %yticks(0:5:20);
    %ylim([0,20]);
%     legend('t_{DLU}','t_{lb}','t_{LUU}','Orientation','horizontal');
    box on;
    
    %sovling and searching time

    ax4=subplot(4,1,4);
    hold(ax4,'on')
     plot(x_axis-1,A(1,:),':','color',[0.6350 0.0780 0.1840]);
    plot(x_axis-1,A(2,:),'--','color',[0.3010 0.7450 0.9330]);
%     plot(x_axis-1,A(2,:),'--','color',[0.8500 0.3250 0.0980]);
    xlim([0,num_iter]);
    if i == 52
        xline(14,'-',{'iter=14'},'color',[194 197 204]/255);
        xline(82,'-',{'iter=82'},'color',[194 197 204]/255);
        ylim([0,200]);
        yticks([0 100 200])
        yticklabels({'0','100','200'})
    elseif i == 65
         xticks(0:5:100);
        ylim([3,5]);
        yticks([3 4 5]);
%         yticklabels({'0','100','200'})
    elseif i == 61
         xticks(0:5:100);
        ylim([9,25]);
        yticks([9 17 25]);
%         yticklabels({'0','100','200'})
    else
%         ylim([0,200]);
%         yticks([0 100 200])
%         yticklabels({'0','100','200'})
    end
     xticks(0:5:100);
     xticklabels({'0','','','','20','','','','40','','','','60','','','','80','','','','100'});
    xlb=xlabel("iteration index");
%     xlb.Position(2)= xlb.Position(2)-8;
    ylb4=ylabel("time (sec)");
    legend('solve','search','Orientation','horizontal','Location','best');
    box on;
    
    ylb_ypos=min([ylb1.Position(1),ylb2.Position(1),ylb3.Position(1),ylb4.Position(1)]);
    ylb1.Position(1) = ylb_ypos-ylb_shift(index);
    ylb2.Position(1) = ylb_ypos-ylb_shift(index);
    ylb3.Position(1) = ylb_ypos-ylb_shift(index);
    ylb4.Position(1) = ylb_ypos-ylb_shift(index);
    subplot_x = 0.17;
    ax1.Position(1) =subplot_x;
    ax2.Position(1) =subplot_x;
    ax3.Position(1) =subplot_x;
    ax4.Position(1) =subplot_x;
    ax1.Position(2) =ax1.Position(2)-0.01;
    ax2.Position(2) =ax2.Position(2)-0.01;
    ax3.Position(2) =ax3.Position(2)-0.01;
    ax4.Position(2) =ax4.Position(2)-0.01;
    subplot_height = 0.19;
    ax1.Position(4) =subplot_height;
    ax2.Position(4) =subplot_height;
    ax3.Position(4) =subplot_height;
    ax4.Position(4) =subplot_height;
    set(ax1,'FontSize',14)
    set(ax2,'FontSize',14)
    set(ax3,'FontSize',14)
    set(ax4,'FontSize',14)
end