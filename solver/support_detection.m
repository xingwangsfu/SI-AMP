function [support_idx] = support_detection(theta_interp, percentage)

theta_energy = theta_interp.^2;

threshold_energy = sum(theta_energy)*(1-percentage);

tmp_energy = 0;

[theta_energy_sorted, energy_idx] = sort(theta_energy, 'descend');
flag = 0;
i = 1;
while flag == 0
    if tmp_energy >= threshold_energy
        boundary_energy = theta_energy_sorted(i);
        flag = 1;
    end
    tmp_energy = tmp_energy + theta_energy_sorted(i);
    i = i+1;
end

support_idx = energy_idx(1:i);


    

