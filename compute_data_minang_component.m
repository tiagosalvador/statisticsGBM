
close all
clear
save('./StudyOrientation/data_minang')
load('ori_components')

for choice = 1:4
    ori = eval(['ori_component_',num2str(choice)]);
    N = length(ori);
    minang = [];
    for k = 1:N
        minang(k) = misorientation3d_simple(ori{k});
    end
    save_variable('data_minang',minang,['minang_',num2str(choice)])
end