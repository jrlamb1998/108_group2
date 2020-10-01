close all
clear
addpath(genpath(pwd))

figure()
my_ax = gca();
hold on

[A0, L0, e_strain, e_strain1, e_stress] = load_data('data3.csv');
[t_strain,t_stress] = get_true(e_strain, e_stress);
ss_plot(my_ax, e_strain, e_stress)
ss_plot(my_ax, t_strain, t_stress)
E = get_E(e_strain,e_stress);
uts = get_uts(e_stress);
[y_strain, y_stress] = get_yield(e_strain, e_stress, E);
[f_strain, f_stress] = get_fracture(e_strain, e_stress);
strain_necking = get_strain_necking(e_strain, e_stress);
energy_fracture = get_energy(t_strain, t_stress);
reduction = get_reduction(t_strain(end));


%% Supporting Functions
function [A0, L0, e_strain, e_strain1, e_stress] = load_data(filename)
    data = readtable(filename, 'HeaderLines',0);
    A0 = table2array(data(7,2))*10e-6;
    L0 = table2array(data(4,2))*10e-3;
    e_strain = 100 * 10e-3 * table2array(data(11:end,2))./L0;
    e_stress = table2array(data(11:end,3))./A0;
    e_strain1 = 0;
    if size(data,2) > 4
        e_strain1 = table2array(data(11:end,5));
    end
end

function [t_strain,t_stress] = get_true(e_strain, e_stress)
    t_stress = e_stress.*(1 + .01*e_strain);
    t_strain = 100 * log(1 + .01*e_strain);
end

function ss_plot(my_ax, strain, stress)
    plot(my_ax, strain, stress)
end

function E = get_E(strain, stress)
    E = (stress(2)-stress(1))/(strain(2)-strain(1));
end

function uts = get_uts(e_stress)
    uts = max(e_stress);
end

function [f_strain, f_stress] = get_fracture(strain,stress)
    f_strain = strain(end-10);
    f_stress = stress(end-10);
end

function [y_strain, y_stress] = get_yield(strain,stress, E)
    for i = 1:numel(strain)
        linear_stress = 0.01*strain(i)*E - 0.02*E;
        if linear_stress > stress(i)
            break
        end
    end
    y_strain = strain(i);
    y_stress = stress(i);
end

function strain_necking = get_strain_necking(strain,stress)
    [~,max_i] = max(stress);
    strain_necking = strain(max_i);
end

function energy_fracture = get_energy(strain, stress)
    energy_fracture = sum(stress(2:end) .* diff(strain));
end

function reduction = get_reduction(ft_strain)
    ratio = exp(0.01*ft_strain);
    reduction = 100*(1-(1/ratio));
end