clc
clear
close all
list = dir('../LPnetlib/*.txt');

% for i = 1:length(list)
for i =[52 67 65 25 1 3 17 51 61 68 90]
    close all
    list(i).name
    fRead = fopen(strcat('../LPnetlib/',list(i).name), 'r');
    A = fscanf(fRead, '%f %f %f %d %d %d %d %f %d %d %d %d %f %d %d %d %d',[17, Inf]);
    fclose(fRead);
    
    % ------------------------------------------------------------------
    % using stackedplot but work not well
    %-------------------------------------------------------------------
%     varNames = ["nnzLU_DLU","nnzLU_lb","nnzLU_LUU",...
%                 "digitLU_DLU","digitLU_lb","digitLU_LUU",...
%                 "t_DLU","t_lb","t_LUU",...
%                 "t_solve","t_search"];
%      varTypes = ["double","double","double",...
%                  "double","double","double",...
%                  "double","double","double",...
%                  "double","double",];
%     tb = table('Size',[length(A(1,:)),11],'VariableTypes',varTypes,'VariableNames',varNames);
%     tb(:,1)    =table(A(4,:)'+A(6,:)');%nnz(L+U) for DLU
%     tb(2:end,2)=table(A(9,2:end)'+A(11,2:end)');%nnz(L+U) for lb
%     tb(:,3)    =table(A(14,:)'+A(16,:)');%nnz(L+U) for LUU
%     tb(:,4)    =table(A(5,:)'+A(7,:)');%digit(L+U) for DLU
%     tb(2:end,5)=table(A(10,2:end)'+A(12,2:end)');%digit(L+U) for lb
%     tb(:,6)    =table(A(15,:)'+A(17,:)');%digit(L+U) for LUU
%     tb(:,7)    =table(A(3,:)');%time for DLU
%     tb(2:end,8)=table(A(8,2:end)');%time for lb
%     tb(:,9)    =table(A(13,:)');%time for LUU
%     tb(:,10)   =table(A(1,:)');% time for solving 3 linear equations using LUU
%     tb(:,11)   =table(A(2,:)');%time for searching
%     grouped_vars = {["nnzLU_DLU","nnzLU_lb","nnzLU_LUU"],...
%                 ["digitLU_DLU","digitLU_lb","digitLU_LUU"],...
%                 ["t_DLU","t_lb","t_LUU"],...
%                 ["t_solve","t_search"]};
%     h=stackedplot(tb,grouped_vars);
%     h.DisplayLabels = {'nnz','digit','time','time'};
%     h.XLabel = 'iteration index';
% %     h.AxesProperties(4).YLimits = [3, 10];
%     for ax_index = 1:3
%         %set the line style and marker for first three plots as ['r--*'; 'g:.';'b-.o']
%         ax = flip(findobj(h.NodeChildren, 'Type','Axes')); 
%         lines = flip(ax(ax_index).Children); 
%         lines(1).LineStyle = '--';
%         lines(1).Marker = '*';
%         lines(1).Color =[1 0 0];
%         lines(1).MarkerFaceColor =[1 0 0];
%         lines(2).LineStyle = ':';
%         lines(2).Marker = '.';
%         lines(1).Color =[0 0.5 0];
%         lines(1).MarkerFaceColor =[0 0.5 0];
%         lines(3).LineStyle = '-.';
%         lines(3).Marker = 'o';
%         lines(1).Color =[0 0 1];
%         lines(1).MarkerFaceColor =[0 0 1];
% %         h.LineProperties(ax_index).Color =  [;0 0.5 0; 0 0 1];
%     end
%     h.LineProperties(4).Color = [1 0 1;0.8500 0.3250 0.0980];
%     h.AxesProperties(2).LegendVisible = 'off';
%     h.AxesProperties(3).LegendVisible = 'off';
%     h.AxesProperties(4).LegendLabels= {'solve',  'search'};
%     h.AxesProperties(1).LegendLabels= {'DLU',  'lb',  'LUU'};
%     h.AxesProperties(1).LegendLocation = 'north';
%     legs = getfield(struct(h),'LegendHandle');
%     legs(1).Orientation = 'Horizontal';
%     legs(4).Orientation = 'Horizontal';
%     legs(1).Location = 'northoutside';
%     legs(1).Units = 'Normalize';
%     legs(1).Position = [.25, sum(h.Position([2,4])), .5, .05];
%     legs(1).Position = [.40, sum(h.Position([2,4])), .2, .05];
%     

    % ------------------------------------------------------------------
    % subplotplus
    %-------------------------------------------------------------------
    x_axis = 1:size(A,2);
    cell41={{['-g']};{['-g']};{['-g']};{['-g']}};
    C = cell41;
    [h,labelfontsize] = subplotplus(C);
    %total number of nnz of L+U
    subidx = 1;
    set(gcf,'CurrentAxes',h(subidx));
    hold on;
    plot(x_axis,A(4,:)+A(6,:),'r--*');
    plot(x_axis(2:end),A(9,2:end)+A(11,2:end),':.','color',[0 0.5 0]);
    plot(x_axis,A(14,:)+A(16,:),'b-.o');
    ylim([0,10^5]);
    xticks(0:10:100);
    xticklabels({''})
    ylb1=ylabel("nnz");
    legend('DLU','lb','LUU','Orientation','horizontal','Location','best');
    
    %total number of bit size
    subidx = subidx+1;
    set(gcf,'CurrentAxes',h(subidx));
    hold on
    plot(x_axis,A(5,:)+A(7,:),'r--*');
    plot(x_axis(2:end),A(10,2:end)+A(12,2:end),':.','color',[0 0.5 0]);
    plot(x_axis,A(15,:)+A(17,:),'b-.o');
    ylim([0,3*10^9]);
    xticks(0:10:100);
    xticklabels({''})
    ylb2=ylabel("digit");
%     legend('DLU','lb','LUU','Orientation','horizontal');
    
    % factorization time
    subidx = subidx+1;
    set(gcf,'CurrentAxes',h(subidx));
    hold on
    plot(x_axis,A(3,:),'r--*');
    plot(x_axis,A(8,:),':.','color',[0 0.5 0]);
    plot(x_axis,A(13,:),'b-.o');
    ylim([0,120]);
    xlim([0,100]);
    xticks(0:10:100);
    xticklabels({''})
    ylb3=ylabel("time (sec)");
    %yticks(0:5:20);
    %ylim([0,20]);
%     legend('t_{DLU}','t_{lb}','t_{LUU}','Orientation','horizontal');
    
    %sovling and searching time
    subidx = subidx+1;
    set(gcf,'CurrentAxes',h(subidx));
%     title('solving/searching time');
    hold on
     plot(x_axis,A(1,:),'-','color',[0.6350 0.0780 0.1840]);
    plot(x_axis,A(2,:),'--','color',[0.3010 0.7450 0.9330]);
%     plot(x_axis,A(2,:),'--','color',[0.8500 0.3250 0.0980]);
    xlim([0,100]);
    ylim([0,220]);
yticks([0 100 200])
yticklabels({'0','100','200'})
    xlabel("iteration index")
    ylb4=ylabel("time (sec)");
    legend('solve','search','Orientation','horizontal','Location','best');
    
    ylb_ypos=min([ylb1.Position(1),ylb2.Position(1),ylb3.Position(1),ylb4.Position(1)]);
    ylb1.Position(1) = ylb_ypos;
    ylb2.Position(1) = ylb_ypos;
    ylb3.Position(1) = ylb_ypos;
    ylb4.Position(1) = ylb_ypos;

    % ------------------------------------------------------------------
    %subplot
    % ------------------------------------------------------------------
    figure(20)
    set(gcf,'Position',[-600 3 560 560])
    ax1=subplot(4,1,1);
    hold(ax1,'on')
    plot(x_axis,A(4,:)+A(6,:),'r--*');
    plot(x_axis(2:end),A(9,2:end)+A(11,2:end),':.','color',[0 0.5 0]);
    plot(x_axis,A(14,:)+A(16,:),'b-.o');
    xline(14,'-','color',[194 197 204]/255);
    xline(82,'-','color',[194 197 204]/255);
    xticks(0:5:100);
    xticklabels({''})
    max_val=max([max(A(4,:)+A(6,:)),...
                 max(A(9,2:end)+A(11,2:end)),...
                 max(A(14,:)+A(16,:))]);
    min_val=min([min(A(4,:)+A(6,:)),...
                 min(A(9,2:end)+A(11,2:end)),...
                 min(A(14,:)+A(16,:))]);
    expon = floor(log10(max_val));
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
    xline(14,'-','color',[194 197 204]/255);
    xline(82,'-','color',[194 197 204]/255);
    xticks(0:5:100);
    xticklabels({''});
    max_val=max([max(A(5,:)+A(7,:)),...
                 max(A(10,2:end)+A(12,2:end)),...
                 max(A(15,:)+A(17,:))]);
    min_val=min([min(A(5,:)+A(7,:)),...
                 min(A(10,2:end)+A(12,2:end)),...
                 min(A(15,:)+A(17,:))]);
    expon = floor(log10(max_val));
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
    ylb2=ylabel("digit(\times10^"+expon+")");
%     legend('DLU','lb','LUU','Orientation','horizontal');
    box on;
    
    % factorization time

    ax3=subplot(4,1,3);
    hold(ax3,'on')
    plot(x_axis,A(3,:),'r--*');
    plot(x_axis,A(8,:),':.','color',[0 0.5 0]);
    plot(x_axis,A(13,:),'b-.o');
    xline(14,'-','color',[194 197 204]/255);
    xline(82,'-','color',[194 197 204]/255);
    ylim([0,100]);
    xlim([0,100]);
    xticks(0:5:100);
    xticklabels({''})
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
    xline(14,'-','color',[194 197 204]/255);
    xline(82,'-','color',[194 197 204]/255);
%     plot(x_axis-1,A(2,:),'--','color',[0.8500 0.3250 0.0980]);
    xlim([0,100]);
    xticks(0:5:100);
    ylim([0,200]);
    yticks([0 100 200])
    yticklabels({'0','100','200'})
    xticks([0:5:10,14,15:5:80,82,85:5:100]);
    xticklabels({'0','','','14','','20','','','','40','','','','60','','','','','82','','','','100'});
    xlabel("iteration index")
    ylb4=ylabel("time (sec)");
    legend('solve','search','Orientation','horizontal','Location','best');
    box on;
    
    ylb_ypos=min([ylb1.Position(1),ylb2.Position(1),ylb3.Position(1),ylb4.Position(1)]);
    ylb1.Position(1) = ylb_ypos-4;
    ylb2.Position(1) = ylb_ypos-4;
    ylb3.Position(1) = ylb_ypos-4;
    ylb4.Position(1) = ylb_ypos-4;
    subplot_height = 0.19;
    ax1.Position(4) =subplot_height;
    ax2.Position(4) =subplot_height;
    ax3.Position(4) =subplot_height;
    ax4.Position(4) =subplot_height;
    set(ax1,'FontSize',14)
    set(ax2,'FontSize',14)
    set(ax3,'FontSize',14)
    set(ax4,'FontSize',14)

    % ------------------------------------------------------------------
    % seperate plots
    %-------------------------------------------------------------------
    
%     x_axis = 1:size(A,2);
%     figure(1);
% %     title('number of nonzero');
%     hold on
%     plot(x_axis,A(4,:),'r--*');
%     plot(x_axis(2:end),A(9,2:end),':.','color',[0 0.5 0]);
%     plot(x_axis,A(14,:),'b-.o');
%     xlabel("iteration index", 'FontSize', 20)
%     ylabel("number of nonzero", 'FontSize', 20);
%     legend('L_{DLU}','L_{lb}','L_{LUU}','Orientation','horizontal');
%     box on;
%     %figure_handle = gca;
%     %figure_handle.YAxis.Exponent = 3;
    
%     figure(2);
% %     title('number of digit');
%     hold on
%     plot(x_axis,A(5,:),'r--*');
%     plot(x_axis(2:end),A(10,2:end),':.','color',[0 0.5 0]);
%     plot(x_axis,A(15,:),'b-.o');
%     xlabel("iteration index")
%     ylabel("number of digit");
%     legend('L_{DLU}','L_{lb}','L_{LUU}','Orientation','horizontal');
%     box on;
    
%     figure(3);
% %     title('number of nonzero');
%     hold on
%     plot(x_axis,A(6,:),'r--*');
%     plot(x_axis(2:end),A(11,2:end),':.','color',[0 0.5 0]);
%     plot(x_axis,A(16,:),'b-.o');
%     xlabel("iteration index")
%     ylabel("number of nonzero");
%     legend('U_{DLU}','U_{lb}','U_{LUU}','Orientation','horizontal');
%     box on;
%     figure_handle = gca;
%     figure_handle.YAxis.Exponent = 4;
    
%     figure(4);
% %     title('number of digit');
%     hold on
%     plot(x_axis,A(7,:),'r--*');
%     plot(x_axis(2:end),A(12,2:end),':.','color',[0 0.5 0]);
%     plot(x_axis,A(17,:),'b-.o');
%     xlabel("iteration index")
%     ylabel("number of digit");
%     legend('U_{DLU}','U_{lb}','U_{LUU}','Orientation','horizontal');
%     box on;
    
    
    figure(7);
%     title('total number of nonzero');
    hold on
    plot(x_axis,A(4,:)+A(6,:),'r--*');
    plot(x_axis(2:end),A(9,2:end)+A(11,2:end),':.','color',[0 0.5 0]);
    plot(x_axis,A(14,:)+A(16,:),'b-.o');
    xline(14,'k:',{'iter=14'});
    xlabel("iteration index")
    ylabel("total number of nonzero");
    legend('DLU','lb','LUU','Orientation','horizontal');
    box on;
    %figure_handle = gca;
    %figure_handle.YAxis.Exponent = 3;
    
    figure(8);
%     title('total number of digit');
    hold on
    plot(x_axis,A(5,:)+A(7,:),'r--*');
    plot(x_axis(2:end),A(10,2:end)+A(12,2:end),':.','color',[0 0.5 0]);
    plot(x_axis,A(15,:)+A(17,:),'b-.o');
    xlabel("iteration index")
    ylabel("total number of digit");
    legend('DLU','lb','LUU','Orientation','horizontal');
    box on;
    
%     figure(i+2);
%     title('factorization/updating time over updating time');
%     hold on
%     semilogy(x_axis,A(3,:)./A(13,:),'r--*');
%     semilogy(x_axis,A(8,:)./A(13,:),':.','color',[0 0.5 0]);
%     semilogy(x_axis,A(13,:)./A(13,:),'b-.o');
%     legend('t_{DLU}','t_{lb}','t_{LUU}');
    
    figure(5);
%     title('factorization time');
    hold on
    semilogy(x_axis,A(3,:),'r--*');
    semilogy(x_axis,A(8,:),':.','color',[0 0.5 0]);
    semilogy(x_axis,A(13,:),'b-.o');
    xline(14,'k:',{'iter=14'});
    xlabel("iteration index")
    ylabel("time (sec)");
    %yticks(0:5:20);
    %ylim([0,20]);
    legend('t_{DLU}','t_{lb}','t_{LUU}','Orientation','horizontal');
    box on;
    
    figure(6);
%     title('solving/searching time');
    hold on
     plot(x_axis-1,A(1,:));
    plot(x_axis-1,A(2,:),'color',[0.8500 0.3250 0.0980]);
    xline(14,'k:',{'iter=14'});
    xlabel("iteration index")
    ylabel("time (sec)");
    %yticks(0:5:35);
    box on;
%     legend('solving 3 linear equations','search entering and existing column');
end