function [t1_avg, t3_avg,bin_center] = compute_avr_per_bin(t1, t3, bins)
    ratio = t1./t3;
    [B, I]=sort(ratio);
    numbin = length(bins);
    min_ratio = min(ratio);
    max_ratio = max(ratio);
    if max_ratio > bins(numbin) || min_ratio < bins(1)
        error('incorrect bins');
    end
    
    curbin = 1;
    arr_size = 0;
    t3_avg_per_bin = 0;
    t1_avg_per_bin = 0;
    case_per_bin = 0;
    %total_count = 0;
    for i = 1:length(ratio)
        if B(i) < bins(curbin+1)
            case_per_bin = case_per_bin+1;
            %total_count = total_count+1;
            t3_avg_per_bin = t3_avg_per_bin+t3(I(i));
            t1_avg_per_bin = t1_avg_per_bin+t1(I(i));
        else
            if case_per_bin > 0
                arr_size = arr_size+1;
                t3_avg(arr_size) = t3_avg_per_bin/case_per_bin;
                t1_avg(arr_size) = t1_avg_per_bin/case_per_bin;
                bin_center(arr_size) = (bins(curbin)+bins(curbin+1))/2;
            end
            case_per_bin = 1;
            %total_count = total_count+1;
            t3_avg_per_bin = t3(I(i));
            t1_avg_per_bin = t1(I(i));
            while B(i) >= bins(curbin+1)
                curbin = curbin+1;
            end
        end
    end
    if case_per_bin > 0
        arr_size = arr_size+1;
        t3_avg(arr_size) = t3_avg_per_bin/case_per_bin;
        t1_avg(arr_size) = t1_avg_per_bin/case_per_bin;
        bin_center(arr_size) = (bins(curbin)+bins(curbin+1))/2;
    end
    %[total_count length(ratio)]
end