close all
clc
clear all
list = dir('../LPnetlib/*.txt');

max_ratio = 0;
min_ratio = 1;
max_lp_ratio = 10;
min_lp_ratio = 10;
k = 0;
t1_all = [];
t3_all = [];
LU1 = [];
LU3 = [];
plotColors = jet(length(list));
% plotColors = hsv(length(list));

figure(50);
%title('factorization/updating time over updating time');
loglog([10^-5,10^5],[10^-5,10^5],'LineStyle',':','Color','k');
line([10^-5,10^5],[10^-6,10^4],'LineStyle','-.','Color','k');
line([10^-5,10^5],[10^-7,10^3],'LineStyle','--','Color','k');
line([10^-5,10^5],[10^-8,10^2],'LineStyle','-','Color','k');
xlim([10^-5, 10^5]);
xticks([10^-5, 10^-3, 10^-1,10^0,10^1,10^3,10^5]);
ylim([10^-7, 10^5]);
yticks([10^-7, 10^-5, 10^-3, 10^-1,10^0,10^1,10^3,10^5]);
xlabel("t_{DLU} (sec)");
ylabel("t_{LUU} (sec)");
legend('1x','1/10x','1/100x','1/1000x','Location','northwest','NumColumns',2);
hold on

figure(51);
loglog([10^-4,10^8],[10^-4,10^8],'LineStyle',':','Color','k');
line([10^-4,10^8],[10^-5,10^7],'LineStyle','-.','Color','k');
line([10^-4,10^8],[10^-6,10^6],'LineStyle','--','Color','k');
xlabel("T_{DLU} (sec)");
ylabel("T_{LUU} (sec)");
legend('1x','1/10x','1/100x','Location','northwest','NumColumns',3);
ylim([10^-6, 10^8]);
yticks([10^-6, 10^-4, 10^-2,10^0,10^2,10^4,10^6, 10^8]);
xlim([10^-4, 10^8]);
xticks([10^-4, 10^-2, 10^0, 10^2,10^4,10^6,10^8]);
hold on

    
for i = 1:length(list)
%     list(i).name
    fRead = fopen(strcat('../LPnetlib/',list(i).name), 'r');
    A = fscanf(fRead, '%f %f %f %d %d %d %d %f %d %d %d %d %f %d %d %d %d',[17, Inf]);
    fclose(fRead);
    
    if sum(A(13,:))~=0
        k=k+1;
        t1(k) = sum(A(3,:));
        t1_all = cat(2, t1_all, A(3,:));
        LU1 = cat(2, LU1, A(4,:)+A(6,:));
        %t2(k) = sum(A(8,2:end));
        t3(k) = sum(A(13,:));
        t3_all = cat(2, t3_all, A(13,:));
        LU3 = cat(2, LU3, A(14,:)+A(16,:));
    else
        i
    end
    if max_lp_ratio < t1(k)/t3(k)
        max_lp_ratio = t1(k)/t3(k);
        max_lp_index = i;
    end
    if min_lp_ratio > t1(k)/t3(k)
        min_lp_ratio = t1(k)/t3(k);
        min_lp_index = i;
    end

    for j = 1:size(A,2)
        if A(13,j)~=0 && max_ratio < A(3,j)/A(13,j)
            max_ratio = A(3,j)/A(13,j);
            max_index = i;
        end
        if A(13,j)~=0 && min_ratio > A(3,j)/A(13,j)
            min_ratio = A(3,j)/A(13,j);
            min_index = i;
        end
    end
    %max_ratio
    
    t1_case = A(3,:);
    t3_case = A(13,:);
    nt1_case = t1_case(t1_case.*t3_case ~= 0);
    nt3_case = t3_case(t1_case.*t3_case ~= 0);
    thisColor = plotColors(length(list)-i+1, :);
    figure(50);
    %title('factorization/updating time over updating time');
    h=scatter(nt1_case,nt3_case, 6,'MarkerEdgeColor',thisColor,...
              'MarkerFaceColor','k',...
              'LineWidth',0.5);
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    
    
    figure(51);
    %title('factorization/updating time over updating time');
    %h=loglog(t1(k),t3(k),'.', 'Color', thisColor,'MarkerSize',20);
    h=scatter(t1(k),t3(k), 18,'MarkerEdgeColor',thisColor,...
              'MarkerFaceColor','k',...
              'LineWidth',1.5);
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';
end
fprintf("worst ratio for one iteration\n");
min_ratio
list(min_index).name
fprintf("best ratio for one iteration\n");
max_ratio
list(max_index).name
fprintf("worst ratio for one lp\n");
min_lp_ratio
list(min_lp_index).name
fprintf("best ratio for one lp\n");
max_lp_ratio
list(max_lp_index).name

figure(101)
loglog(LU1./LU3,t1_all./t3_all,'b.');

% figure(1);
% semilogy(1:k,t1./t3,'r-o');
% %title('factorization/updating time over updating time');
% xticks(0: 10: 80);
% xlabel("case index");
% ylabel("time ratio");
% hold on
% semilogy(1:k,t2./t3,'b-o');
% semilogy(1:k,t3./t3,'g-o');
% legend('t_{DLU}','t_{lb}','t_{LUU}');

figure(3);
%title('factorization/updating time over updating time');
loglog(t1,t1);
ylim([10^-6, 10^8]);
yticks([10^-6, 10^-4, 10^-2,10^0,10^2,10^4,10^6, 10^8]);
xlim([10^-4, 10^8]);
xticks([10^-4, 10^-2, 10^0, 10^2,10^4,10^6,10^8]);
hold on
loglog(t1,t1/10);
loglog(t1,t1/100);
loglog(t1,t3,'b.');
xlabel("T_{DLU} (sec)");
ylabel("T_{LUU} (sec)");
legend('1x','1/10x','1/100x');

%===================================================================
% ratio = t1./t3;
% [B, I]=sort(ratio);
% numbin=20;
% maxratio = max(ratio);
% if round(maxratio,-ceil(log10(maxratio))+1) < maxratio
%     maxratio = round(maxratio,-ceil(log10(maxratio))+1)+10^(ceil(log10(maxratio))-1);
% else
%     maxratio = round(maxratio,-ceil(log10(maxratio))+1);
% end
% minratio = 0;
% binsize = (maxratio-minratio)/numbin;
% bins = minratio:binsize:maxratio;
% [t1_avg, t3_avg,bin_center,bin_count] = compute_avr_per_bin(t1, t3, bins);
% 
% figure(4);
% xlabel('time ratio for update/direct');
% yyaxis left
% h1=semilogy(ratio(I),t3(I), 'b-o');
% hold on;
% h2=semilogy(ratio(I),t1(I), 'r-*');
% ylabel('run time');
% 
% yyaxis right
% histogram(ratio,bins)
% ylabel('number of LP problems')
% legend('LU update','direct LU');
% 
% hold off
% 
% figure(5);
% xlabel('time ratio for update/direct');
% yyaxis left
% semilogy(bin_center,t3_avg, 'b-o');
% hold on;
% semilogy(bin_center,t1_avg, 'r-*');
% ylabel('run time');
% 
% yyaxis right
% histogram(ratio,bins)
% ylabel('number of LP problems')
% legend('LU update','direct LU');
% 
% hold off

%----------------------------------------------------------

% numbin=200;
% binsize = (maxratio-minratio)/numbin;
% bins = minratio:binsize:maxratio;
% [t1_avg, t3_avg,bin_center,bin_count] = compute_avr_per_bin(t1, t3, bins);
% 
% figure(6);
% xlabel('time ratio for update/direct');
% yyaxis left
% semilogy(ratio(I),t3(I), 'b-o');
% hold on;
% semilogy(ratio(I),t1(I), 'r-*');
% ylabel('run time');
% 
% yyaxis right
% histogram(ratio,bins)
% ylabel('number of LP problems')
% legend('LU update','direct LU');
% 
% hold off
% 
% figure(7);
% xlabel('time ratio for update/direct');
% yyaxis left
% semilogy(bin_center,t3_avg, 'b-o');
% hold on;
% semilogy(bin_center,t1_avg, 'r-*');
% ylabel('run time');
% 
% yyaxis right
% histogram(ratio,bins)
% ylabel('number of LP problems')
% 
% legend('LU update','direct LU');
% hold off

%===================================================================
% nt3=t3_all(t3_all~=0);
% nt1=t1_all(t3_all~=0);
% ratio = nt1./nt3;
% [B, I]=sort(ratio);
% numbin=20;
% maxratio = max(ratio);
% if round(maxratio,-ceil(log10(maxratio))+1) < maxratio
%     maxratio = round(maxratio,-ceil(log10(maxratio))+1)+10^(ceil(log10(maxratio))-1);
% else
%     maxratio = round(maxratio,-ceil(log10(maxratio))+1);
% end
% minratio = 0;
% binsize = (maxratio-minratio)/numbin;
% bins = minratio:binsize:maxratio;
% [t1_avg, t3_avg,bin_center,bin_count] = compute_avr_per_bin(nt1, nt3, bins);
% 
% figure(8);
% xlabel('time ratio for update/direct');
% yyaxis left
% semilogy(ratio(I),nt3(I), 'b-o');
% hold on;
% semilogy(ratio(I),nt1(I), 'r-*');
% ylabel('run time');
% 
% yyaxis right
% histogram(ratio,bins)
% ylabel('number of updates/factorizations')
% legend('LU update','direct LU');
% 
% hold off
% 
% figure(9);
% xlabel('time ratio for update/direct');
% yyaxis left
% semilogy(bin_center,t3_avg, 'b-o');
% hold on;
% semilogy(bin_center,t1_avg, 'r-*');
% ylabel('run time');
% 
% yyaxis right
% histogram(ratio,bins)
% ylabel('number of updates/factorizations')
% 
% legend('LU update','direct LU');
% hold off

%----------------------------------------------------------

% numbin=200;
% binsize = (maxratio-minratio)/numbin;
% bins = minratio:binsize:maxratio;
% [t1_avg, t3_avg,bin_center,bin_count] = compute_avr_per_bin(nt1, nt3, bins);
% 
% figure(10);
% xlabel('time ratio for update/direct');
% yyaxis left
% semilogy(ratio(I),nt3(I), 'b-o');
% hold on;
% semilogy(ratio(I),nt1(I), 'r-*');
% ylabel('run time');
% 
% yyaxis right
% histogram(ratio,bins)
% ylabel('number of updates/factorizations')
% 
% legend('LU update','direct LU');
% hold off
% 
% figure(11);
% xlabel('time ratio for update/direct');
% yyaxis left
% semilogy(bin_center,t3_avg, 'b-o');
% hold on;
% semilogy(bin_center,t1_avg, 'r-*');
% ylabel('run time');
% 
% yyaxis right
% histogram(ratio,bins)
% ylabel('number of updates/factorizations')
% 
% legend('LU update','direct LU');
% hold off

%********************************log bin***********************************

ratio = t1./t3;
[B, I]=sort(ratio);
%bins = [1,2,5,10,20,50,100,200,500];
bins=2.^(0:8);
[t1_avg, t3_avg,bin_center,bin_count_lp] = compute_avr_per_bin(t1, t3, bins);

% figure(14);
% xlabel('time ratio for update/direct');
% yyaxis left
% h1=semilogy(ratio(I),t3(I), 'b-o');
% hold on;
% h2=semilogy(ratio(I),t1(I), 'r-*');
% ylabel('run time');
% 
% yyaxis right
% histogram(ratio,bins)
% ylabel('number of LP problems')
% legend('LU update','direct LU');
% set(gca, "XScale", "log")
% 
% hold off

figure(15);
xticks(10.^(0: 3) );
xlabel('time ratio for T_{DLU}/T_{LUU}');
yyaxis left
h=histogram(ratio,bins);
h.Annotation.LegendInformation.IconDisplayStyle = 'off';
%semilogy(bin_center,bin_count, 'y-+');
ylabel('number of LP problems')

yyaxis right
semilogy(bin_center,t1_avg, ':*','LineWidth',1.5,'MarkerSize',8);
hold on;
semilogy(bin_center,t3_avg, '--o','LineWidth',1.5);
ylim([10^-2, 10^6]);
ylabel('average time (sec)');
legend('T_{DLU}','T_{LUU}');
set(gca, "XScale", "log")

hold off
fprintf("number of lp for each bin\n");
bin_count_lp
bin_count_lp/length(list)

%===================================================================
nt3=t3_all(t3_all.*t1_all~=0);
nt1=t1_all(t3_all.*t1_all~=0);
ratio = nt1./nt3;
[B, I]=sort(ratio);
%bins = [0.05,0.1,0.2,0.5,1,2,5,10,20,50,100,200,500,1000,2000,5000,10000,20000];
bins = 2.^(-2:20);
[t1_avg, t3_avg,bin_center,bin_count] = compute_avr_per_bin(nt1, nt3, bins);

% figure(18);
% xlabel('time ratio for update/direct');
% yyaxis left
% semilogy(ratio(I),nt3(I), 'b-o');
% ylim([10^-2, 10^6]);
% hold on;
% semilogy(ratio(I),nt1(I), 'r-*');
% ylabel('run time');
% 
% yyaxis right
% histogram(ratio,bins)
% ylabel('number of updates/factorizations')
% legend('LU update','direct LU');
% set(gca, "XScale", "log")
% 
% hold off

figure(19);
xlim([10^-1, 10^7]);
xticks([10^-1,10^0,10,10^3,10^5,10^7] );
xlabel('time ratio for t_{DLU}/t_{LUU}');
yyaxis left
h=histogram(ratio,bins);
h.Annotation.LegendInformation.IconDisplayStyle = 'off';
ax = gca;
ax.YAxis(1).Exponent = 2;
ylabel('number of column replacements')
ylim([0, 1800]);
yticks([0, 300, 600, 900,1200,1500,1800]);

yyaxis right
semilogy(bin_center,t1_avg, ':*','LineWidth',1.5,'MarkerSize',8);
hold on;
semilogy(bin_center,t3_avg, '--o','LineWidth',1.5);
ylabel('average time (sec)');
ylim([10^-2, 10^6]);
yticks([10^-2, 10^0, 10^2, 10^4,10^6]);
legend('t_{DLU}','t_{LUU}');
set(gca, "XScale", "log")
hold off
fprintf("number of iteration for each bin\n");
bin_count
bin_count/length(nt1)


figure(2);
%title('factorization/updating time over updating time');
h=loglog(nt1,nt3,'b.');
h.Annotation.LegendInformation.IconDisplayStyle = 'off';
xlim([10^-5, 10^5]);
xticks([10^-5, 10^-3, 10^-1,10^0,10^1,10^3,10^5]);
ylim([10^-7, 10^5]);
yticks([10^-7, 10^-5, 10^-3, 10^-1,10^0,10^1,10^3,10^5]);
hold on
loglog(nt1,nt1);
loglog(nt1,nt1/10);
loglog(nt1,nt1/100);
xlabel("t_{DLU} (sec)");
ylabel("t_{LUU} (sec)");
legend('1x','1/10x','1/100x');

fprintf("fix font size = 14 points\n");
fprintf("check expand axis to fill figure\n");
