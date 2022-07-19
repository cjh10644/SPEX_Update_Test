close all
clc
clear all
list = dir('../LPnetlib/*.txt');

max_ratio = 0;
k = 0;
t1_all = [];
t3_all = [];
for i = 1:length(list)
    list(i).name
    fRead = fopen(strcat('../LPnetlib/',list(i).name), 'r');
    A = fscanf(fRead, '%f %f %f %d %d %d %d %f %d %d %d %d %f %d %d %d %d',[17, Inf]);
    fclose(fRead);
    
    if sum(A(13,2:end))~=0
        k=k+1;
        t1(k) = sum(A(3,2:end));
        t1_all = cat(2, t1_all, A(3,2:end));
        t2(k) = sum(A(8,2:end));
        t3(k) = sum(A(13,2:end));
        t3_all = cat(2, t3_all, A(13,2:end));
    else
        i
    end
    for j = 1:size(A,2)
        if A(13,j)~=0 && max_ratio < A(3,j)/A(13,j)
            max_ratio = A(3,j)/A(13,j)
            A(3,j)
            A(13,j)
        end
    end
    %max_ratio
end
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
hold on
loglog(t1,t1/10);
loglog(t1,t1/100);
loglog(t1,t3,'bo');
xlabel("DLU (sec)");
ylabel("LUU (sec)");
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
[t1_avg, t3_avg,bin_center,bin_count] = compute_avr_per_bin(t1, t3, bins);

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
xlabel('time ratio for t_{DLU}/t_{LUU}');
yyaxis left
semilogy(bin_center,t3_avg, 'b-o');
ylim([10^-2, 10^6]);
hold on;
semilogy(bin_center,t1_avg, 'r-*');
ylabel('average time (sec)');

yyaxis right
histogram(ratio,bins)
%semilogy(bin_center,bin_count, 'y-+');
ylabel('number of LP problems')
legend('LUU','DLU');
set(gca, "XScale", "log")

hold off
bin_count

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
xlim([10^-2, 10^5]);
xticks(10.^(-2: 5) );
xlabel('time ratio for t_{DLU}/t_{LUU}');
yyaxis left
semilogy(bin_center,t3_avg, 'b-o');
hold on;
semilogy(bin_center,t1_avg, 'r-*');
ylabel('average time (sec)');

yyaxis right
histogram(ratio,bins)
ylabel('number of column replacements')
legend('LUU','DLU');
set(gca, "XScale", "log")
hold off


figure(2);
%title('factorization/updating time over updating time');
h=loglog(nt1,nt3,'b.');
h.Annotation.LegendInformation.IconDisplayStyle = 'off';
xlim([10^-2, 10^4]);
ylim([10^-4, 10^4]);
hold on
loglog(nt1,nt1);
loglog(nt1,nt1/10);
loglog(nt1,nt1/100);
xlabel("DLU (sec)");
ylabel("LUU (sec)");
legend('1x','1/10x','1/100x');

%fix font size = 14 points
%check expand axis to fill figure
