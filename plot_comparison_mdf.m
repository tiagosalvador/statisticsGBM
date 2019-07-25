function plot_comparison_mdf(data,NumBins,title_name,saveto,axis_bounds,legend_text)
edgeBins = linspace(0,2*pi,NumBins+1);
BinWidth = edgeBins(2)-edgeBins(1);
edgeCenters = (edgeBins(1:end-1)+BinWidth/2)/(2*pi)*360;
BinWidth = BinWidth/(2*pi)*360;
for j = 1:length(data{1})
    clf
    hold on
    for i = 1:length(data)
        working_data = data{i}{j};
        aux = ceil(working_data(:,1)*NumBins/(2*pi));
        Values = zeros(1,NumBins);
        for k = 1:NumBins
            Values(k) = sum(working_data(aux==k,2));
        end
        Values = Values/sum(Values*BinWidth);
        plot(edgeCenters,Values,'-o')
    end
    if length(legend_text) > length(data)
        plot(edgeCenters,mackenzie_function_degrees(edgeCenters))
    end
    hold off
    legend(legend_text,'Location','best')
    axis(axis_bounds)
    title(title_name)
    print(gcf, '-depsc2', [saveto,num2str(j-1)]);
end