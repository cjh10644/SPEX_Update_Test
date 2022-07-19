close all
clc
clear all
list = dir('LPnetlib_CholUpdate/*.txt');

max_ratio = 0;
min_ratio = Inf;
k = length(list);
nodata = 0;
too_simple = 0;
too_simple1 = 0;
too_simple2 = 0;
t_fact_all = [];
t_update_all = [];
for i = 1:length(list)
    i
    list(i).name
    fRead = fopen(strcat('LPnetlib_CholUpdate/',list(i).name), 'r');
    A = fscanf(fRead, '%f %d %d %f %d %d %f %d %d %f %d %d %f %d %d %f %d %d %f %d %d',[21, 1]);
    fclose(fRead);
    
    t_fact(i) = 0;
    t_update(i) = 0;
    if (length(A) ~= 0)
        t_fact(i) = A(1)*2+A(7)+A(16);
        t_update(i) = A(4)+A(10)+A(13)+A(19);
        t_fact_all = cat(2,t_fact_all,     [A(1), A(7), A(1), A(16)]);
        t_update_all = cat(2,t_update_all, [A(4), A(10),A(13),A(19)]);
    else
        nodata = nodata+1
    end
%     figure(2);
%     loglog(t_fact,t_update,'b*');
%     hold on;
%     if i>1
%         loglog(t_fact(1:i-1),t_update(1:i-1),'bo');
%     end
%     loglog(t_fact(i),t_update(i),'r*');
%     if i == 69 || i ==70
%         figure(i);
%         loglog(t_fact,t_update,'bo');
%         hold on;
%         loglog(t_fact(i),t_update(i),'ro');
%     end

    if t_update(i) == 0
        too_simple1 = too_simple1+1
    end
    if t_fact(i) == 0
        too_simple2 = too_simple2+1
    end
    if t_update(i) == 0 || t_fact(i) == 0
        too_simple = too_simple+1;
    end
    i-too_simple
    
    if t_update(i)~=0 && t_fact(i)/t_update(i) > max_ratio
        max_ratio =t_fact(i)/t_update(i)
    end
    if t_fact(i)~=0 && t_fact(i)/t_update(i) < min_ratio
        min_ratio =t_fact(i)/t_update(i)
    end
    
end
figure(1);
%title('factorization/updating time over updating time');
xticks(0: 10: 80);
xlabel("case index");
ylabel("time ratio");
hold on
semilogy(1:k,t_fact./t_update,'b-o');
semilogy(1:k,t_update./t_update,'g-o');
legend('t_{DLU}','t_{LUU}');


figure(3);
%title('factorization/updating time over updating time');
loglog(t_fact,t_fact);
hold on
loglog(t_fact,t_fact/10);
loglog(t_fact,t_fact/100);
loglog(t_fact,t_update,'bo');
xlabel("DC (sec)");
ylabel("CU (sec)");
legend('1x','1/10x','1/100x');






%===================================================================
t_fact=t_fact(t_update~=0);
t_update=t_update(t_update~=0);
ratio = t_fact./t_update;
[B, I]=sort(ratio);
numbin=20;
maxratio = max(ratio);
if round(maxratio,-ceil(log10(maxratio))+1) < maxratio
    maxratio = round(maxratio,-ceil(log10(maxratio))+1)+10^(ceil(log10(maxratio))-1);
else
    maxratio = round(maxratio,-ceil(log10(maxratio))+1);
end
minratio = 0;
binsize = (maxratio-minratio)/numbin;
bins = minratio:binsize:maxratio;
[t1_avg, t3_avg,bin_center,~] = compute_avr_per_bin(t_fact, t_update, bins);

figure(4);
xlabel('time ratio for update/direct');
yyaxis left
h1=semilogy(ratio(I),t_update(I), 'b-o');
hold on;
h2=semilogy(ratio(I),t_fact(I), 'r-*');
ylabel('run time');

yyaxis right
histogram(ratio,bins)
ylabel('number of LP problems')
legend('Cholesky update','direct Cholesky');

hold off

figure(5);
xlabel('time ratio for update/direct');
yyaxis left
semilogy(bin_center,t3_avg, 'b-o');
hold on;
semilogy(bin_center,t1_avg, 'r-*');
ylabel('run time');

yyaxis right
histogram(ratio,bins)
ylabel('number of LP problems')
legend('Cholesky update','direct Cholesky');

hold off

%----------------------------------------------------------

numbin=200;
binsize = (maxratio-minratio)/numbin;
bins = minratio:binsize:maxratio;
[t1_avg, t3_avg,bin_center,~] = compute_avr_per_bin(t_fact, t_update, bins);

figure(6);
xlabel('time ratio for update/direct');
yyaxis left
semilogy(ratio(I),t_update(I), 'b-o');
hold on;
semilogy(ratio(I),t_fact(I), 'r-*');
ylabel('run time');

yyaxis right
histogram(ratio,bins)
ylabel('number of LP problems')
legend('Cholesky update','direct Cholesky');

hold off

figure(7);
xlabel('time ratio for update/direct');
yyaxis left
semilogy(bin_center,t3_avg, 'b-o');
hold on;
semilogy(bin_center,t1_avg, 'r-*');
ylabel('run time');

yyaxis right
histogram(ratio,bins)
ylabel('number of LP problems')

legend('Cholesky update','direct Cholesky');
hold off

%===================================================================
nt3=t_update_all(t_update_all~=0);
nt1=t_fact_all(t_update_all~=0);
ratio = nt1./nt3;
[B, I]=sort(ratio);
numbin=20;
maxratio = max(ratio);
if round(maxratio,-ceil(log10(maxratio))+1) < maxratio
    maxratio = round(maxratio,-ceil(log10(maxratio))+1)+10^(ceil(log10(maxratio))-1);
else
    maxratio = round(maxratio,-ceil(log10(maxratio))+1);
end
minratio = 0;
binsize = (maxratio-minratio)/numbin;
bins = minratio:binsize:maxratio;
[t1_avg, t3_avg,bin_center,~] = compute_avr_per_bin(nt1, nt3, bins);

figure(8);
xlabel('time ratio for update/direct');
yyaxis left
semilogy(ratio(I),nt3(I), 'b-o');
hold on;
semilogy(ratio(I),nt1(I), 'r-*');
ylabel('run time');

yyaxis right
histogram(ratio,bins)
ylabel('number of updates/factorizations')
legend('Cholesky update','direct Cholesky');

hold off

figure(9);
xlabel('time ratio for update/direct');
yyaxis left
semilogy(bin_center,t3_avg, 'b-o');
hold on;
semilogy(bin_center,t1_avg, 'r-*');
ylabel('run time');

yyaxis right
histogram(ratio,bins)
ylabel('number of updates/factorizations')

legend('Cholesky update','direct Cholesky');
hold off

%----------------------------------------------------------

numbin=200;
binsize = (maxratio-minratio)/numbin;
bins = minratio:binsize:maxratio;
[t1_avg, t3_avg,bin_center,~] = compute_avr_per_bin(nt1, nt3, bins);

figure(10);
xlabel('time ratio for update/direct');
yyaxis left
semilogy(ratio(I),nt3(I), 'b-o');
hold on;
semilogy(ratio(I),nt1(I), 'r-*');
ylabel('run time');

yyaxis right
histogram(ratio,bins)
ylabel('number of updates/factorizations')

legend('Cholesky update','direct Cholesky');
hold off

figure(11);
xlabel('time ratio for update/direct');
yyaxis left
semilogy(bin_center,t3_avg, 'b-o');
hold on;
semilogy(bin_center,t1_avg, 'r-*');
ylabel('run time');

yyaxis right
histogram(ratio,bins)
ylabel('number of updates/factorizations')

legend('Cholesky update','direct Cholesky');
hold off

%**************************log bin**********************************

t1=t_fact(t_update.*t_fact~=0);
t3=t_update(t_update.*t_fact~=0);
ratio = t1./t3;
[B, I]=sort(ratio);
%bins = [0,3,10,30,100,300,1000,3000];
%bins = [1,2,5,10,20,50,100,200,500,1000,2000];
bins=2.^(0:11);
[t1_avg, t3_avg,bin_center,~] = compute_avr_per_bin(t1, t3, bins);

figure(14);
xlabel('time ratio for update/direct');
yyaxis left
h1=semilogy(ratio(I),t3(I), 'b-o');
hold on;
h2=semilogy(ratio(I),t1(I), 'r-*');
ylabel('run time');

yyaxis right
histogram(ratio,bins)
ylabel('number of LP problems')
legend('Cholesky update','direct Cholesky');
set(gca, "XScale", "log");
hold off

figure(15);
xticks(10.^(0: 4) );
xlabel('time ratio for t_{DC}/t_{CU}');
yyaxis left
semilogy(bin_center,t3_avg, 'b-o');
ylim([10^-2, 10^6]);
yticks(10.^(-2: 2: 6));
hold on;
semilogy(bin_center,t1_avg, 'r-*');
ylabel('average time (sec)');

yyaxis right
histogram(ratio,bins)
ylabel('number of LP problems')
legend('CU','DC');
set(gca, "XScale", "log")
hold off
%===================================================================
nt3=t_update_all(t_update_all.*t_fact_all~=0);
nt1=t_fact_all(t_update_all.*t_fact_all~=0);
ratio = nt1./nt3;
[B, I]=sort(ratio);
%bins = [0.5,1,2,5,10,20,50,100,200,500,1000,2000,5000,10000,20000];
bins = 2.^(-1:15);
[t1_avg, t3_avg,bin_center,~] = compute_avr_per_bin(nt1, nt3, bins);

figure(18);
xlabel('time ratio for t_{DC}/t_{CU}');
yyaxis left
semilogy(ratio(I),nt3(I), 'b-o');
hold on;
semilogy(ratio(I),nt1(I), 'r-*');
ylabel('run time');

yyaxis right
histogram(ratio,bins)
ylabel('number of DC/CU')
legend('t_{CU}','t_{DC}');
set(gca, "XScale", "log")

hold off

figure(19);
xticks(10.^(-1: 5) );
xlabel('time ratio for t_{DC}/t_{CU}');
yyaxis left
semilogy(bin_center,t3_avg, 'b-o');
ylim([10^-2, 10^6]);
yticks(10.^(-2: 2: 6));
hold on;
semilogy(bin_center,t1_avg, 'r-*');
ylabel('average time (sec)');

yyaxis right
histogram(ratio,bins)
ylabel('number of rank-1 modifications')
legend('CU','DC');
set(gca, "XScale", "log")
hold off


figure(2);
%title('factorization/updating time over updating time');
loglog(nt1,nt1);
hold on
loglog(nt1,nt1/10);
loglog(nt1,nt1/100);
loglog(nt1,nt3,'bo');
xlabel("DC (sec)");
ylabel("CU (sec)");
legend('1x','1/10x','1/100x');

%fix font size = 14 points
%check expand axis to fill figure
