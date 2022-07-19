list = dir('../LPnetlib/*.txt');
clc
close all

%for i = 1:length(list)
%for i = 1:104
for i =[5 6 17 25 35 50 51 63 65]
    list(i).name
    fRead = fopen(strcat('../LPnetlib/',list(i).name), 'r');
    A = fscanf(fRead, '%f %f %f %d %d %d %d %f %d %d %d %d %f %d %d %d %d',[17, Inf]);
    close all
    fclose(fRead);
    
    x_axis = 1:size(A,2);
    
    figure(1);
%     title('number of nonzero');
    hold on
    plot(x_axis,A(4,:),'r');
    plot(x_axis(2:end),A(9,2:end),'b');
    plot(x_axis,A(14,:),'g');
    xlabel("iteration index", 'FontSize', 20)
    ylabel("number of nonzero", 'FontSize', 20);
    legend('L_{DLU}','L_{lb}','L_{LUU}','Orientation','horizontal');
    box on;
    figure_handle = gca;
    figure_handle.YAxis.Exponent = 3;
    
    figure(2);
%     title('number of digit');
    hold on
    plot(x_axis,A(5,:),'r');
    plot(x_axis(2:end),A(10,2:end),'b');
    plot(x_axis,A(15,:),'g');
    xlabel("iteration index")
    ylabel("number of digit");
    legend('L_{DLU}','L_{lb}','L_{LUU}','Orientation','horizontal');
    box on;
    
    figure(3);
%     title('number of nonzero');
    hold on
    plot(x_axis,A(6,:),'r');
    plot(x_axis(2:end),A(11,2:end),'b');
    plot(x_axis,A(16,:),'g');
    xlabel("iteration index")
    ylabel("number of nonzero");
    legend('U_{DLU}','U_{lb}','U_{LUU}','Orientation','horizontal');
    box on;
    figure_handle = gca;
    figure_handle.YAxis.Exponent = 4;
    
    figure(4);
%     title('number of digit');
    hold on
    plot(x_axis,A(7,:),'r');
    plot(x_axis(2:end),A(12,2:end),'b');
    plot(x_axis,A(17,:),'g');
    xlabel("iteration index")
    ylabel("number of digit");
    legend('U_{DLU}','U_{lb}','U_{LUU}','Orientation','horizontal');
    box on;
    
%     figure(i+2);
%     title('factorization/updating time over updating time');
%     hold on
%     semilogy(x_axis,A(3,:)./A(13,:),'r');
%     semilogy(x_axis,A(8,:)./A(13,:),'b');
%     semilogy(x_axis,A(13,:)./A(13,:),'g');
%     legend('t_{DLU}','t_{lb}','t_{LUU}');
    
    figure(5);
%     title('factorization time');
    hold on
    semilogy(x_axis,A(3,:),'r');
    semilogy(x_axis,A(8,:),'b');
    semilogy(x_axis,A(13,:),'g');
    xlabel("iteration index")
    ylabel("time (sec)");
    %yticks(0:5:20);
    %ylim([0,20]);
    legend('t_{DLU}','t_{lb}','t_{LUU}','Orientation','horizontal');
    box on;
    
    figure(6);
%     title('solving/searching time');
    hold on
%     plot(x_axis,A(1,:));
    plot(x_axis,A(2,:),'k');
    xlabel("iteration index")
    ylabel("time (sec)");
    %yticks(0:5:35);
    box on;
%     legend('solving 3 linear equations','search entering and existing column');
end