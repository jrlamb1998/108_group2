close all
clear
addpath(genpath(pwd))

figure()
ax = gca();

[A0, L0, e_strain, e_strain1, e_stress] = load_data('Specimen_RawData_1.csv');
[t_strain,t_stress] = get_true(ax, e_strain, e_stress);
E = get_E(e_strain,e_stress);
uts = get_uts(e_stress);
[y_strain, y_stress] = get_yield(e_strain,e_stress);
strain_necking = necking_strain(e_strain,e_stress);
energy_necking = get_energy(e_strain,e_stress);
reduction = get_reductoin(e_strain,e_stress);




e2t(ax, e_strain, e_stress);
if ~(e_strain1 == 0)
    e2t(ax, e_strain1, e_stress);
end




clear
%% Supporting Functions
function [t_strain,t_stress] = get_true(ax, e_strain, e_stress)
    plot(ax, e_strain,e_stress, 'k')
    hold on
    t_stress = e_stress.*(1 + .01*e_strain);
    t_strain = 100 * log(1 + .01*e_strain);
    plot(ax, t_strain, t_stress);
end

function [y_strain,y_stress] = get_yield(strain,stress)

function [A0,L0,e_strain,e_strain1, e_stress] = load_data(filename)
    data = readtable(filename);
    A0 = table2array(data(7,1));
    L0 = table2array(data(4,1));
    e_strain = table2array(data(11:end,2))./L0;
    e_stress = table2array(data(11:end,3))./A0;
    e_strain1 = 0;
    if size(data,2) > 4
        e_strain1 = table2array(data(11:end,5));
    end
end
