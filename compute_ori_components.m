clear
save('ori_components')
dims = [256 256 256];
seeds = 1:10;
for choice = 1:4
    fprintf('choice = %d\n',choice);
    ori_component = zeros(0);
    for l = seeds
        fprintf('seed =  %d\n',l)
        load(['simulation3d-',num2str(choice),'-',num2str(dims(1)),'x',num2str(dims(2)),'x',num2str(dims(3)),'-',num2str(l),'.mat'],'grains_t0','ori_t0')
        grains = grains_t0;
        ori = ori_t0;
        
        label = zeros(dims);
        N = size(grains,1); % number of grain types
        start_progress(' - Computing Grain''s componenents')
        for k=1:N
            ind = grains{k,1}; % Pixels within a nhd. of grain.
            val = grains{k,2}; % Lev. set. vals. at those pixels.
            posind = ind(val>0); % Pixels in the interior of grain.
            component = decompose_component3d(dims,posind);
            for i = 1:length(component)
                ori_component{end+1} = ori{k};
            end
            display_progress(k,N,1);
        end
    end
    save_variable('ori_components',ori_component,['ori_component_',num2str(choice)])
end
