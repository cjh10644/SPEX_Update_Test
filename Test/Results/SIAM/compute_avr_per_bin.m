function [t1_avg, t3_avg,bin_center,bin_count] = compute_avr_per_bin(t1, t3, bins)
    ratio = t1./t3;
    [B, I]=sort(ratio);
    numbin = length(bins);
    min_ratio = min(ratio);
    max_ratio = max(ratio);
    if max_ratio > bins(numbin) || min_ratio < bins(1)
        error('incorrect bins, required range [%d %d]', min_ratio,max_ratio);
    end
    
    curbin = 1;
    arr_size = 0;
    t3_avg_per_bin = 0;
    t1_avg_per_bin = 0;
    bin_count = zeros(1,numbin-1);
    for i = 1:length(ratio)
        if B(i) < bins(curbin+1)
            bin_count(curbin) = bin_count(curbin)+1;
            t3_avg_per_bin = t3_avg_per_bin+t3(I(i));
            t1_avg_per_bin = t1_avg_per_bin+t1(I(i));
        else
            if bin_count(curbin) > 0
                arr_size = arr_size+1;
                t3_avg(arr_size) = t3_avg_per_bin/bin_count(curbin);
                t1_avg(arr_size) = t1_avg_per_bin/bin_count(curbin);
                bin_center(arr_size) = (bins(curbin)+bins(curbin+1))/2;
            end
            
            while B(i) >= bins(curbin+1)
                curbin = curbin+1;
            end
            t3_avg_per_bin = t3(I(i));
            t1_avg_per_bin = t1(I(i));
            bin_count(curbin) = bin_count(curbin)+1;
        end
    end
    if bin_count(curbin) > 0
        arr_size = arr_size+1;
        t3_avg(arr_size) = t3_avg_per_bin/bin_count(curbin);
        t1_avg(arr_size) = t1_avg_per_bin/bin_count(curbin);
        bin_center(arr_size) = (bins(curbin)+bins(curbin+1))/2;
    end
    bin_count = bin_count(bin_count~=0);
end