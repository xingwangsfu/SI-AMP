function [mean_d, std_error] = errorbar_plot(errorbar_in)
[m,n] = size(errorbar_in);
for i=1:m
    tmp = errorbar_in(i,:);
    mean_d(i) = mean(tmp);
    std_error(i) = 1.96*std(tmp)/sqrt(n);
    upper_d(i) = mean_d(i) + std_error(i)*1.96;
    lower_d(i) = mean_d(i) - std_error(i)*1.96;
end
