function [xdata,ydata,edata] = plot_errorbars_comparison(data,title_name,saveto,axis_bounds,legend_text,scale,bin_width)

Values = cell(1,length(data));
BinWidth = cell(1,length(data));
edgeBins = cell(1,length(data));
edgeCenters = cell(1,length(data));
error = cell(1,length(data));
for j = 1:length(data{1})
    for i = 1:length(data)
        working_data = data{i}{j};
        if scale
            aux = working_data/mean(working_data);
        else
            aux = working_data;
        end
        if nargin < 7
            h = histogram(aux,'BinMethod','fd','Normalization','pdf');
        else
            h = histogram(aux,'BinWidth',bin_width(j),'Normalization','pdf');
        end
        Values{i} = h.Values;
        BinWidth{i} = h.BinWidth;
        edgeBins{i} = h.BinEdges;
        edgeCenters{i} = edgeBins{i}(1:end-1)+BinWidth{i}/2;
        error{i} = sqrt(h.BinCounts)/(length(working_data)*BinWidth{i});
    end
    clf
    hold on
    for i = 1:length(data)
        xdata{i,j} = edgeCenters{i};
        ydata{i,j} = Values{i};
        edata{i,j} = error{i};
        switch i
            case 1
                marker = '-';
            case 2
                marker = '-';
            case 3
                marker = '-';
        end
        errorbar(edgeCenters{i},Values{i},error{i},marker,'LineWidth',1.5)
    end
    
    hold off
    legend(legend_text,'Location','best')
    axis(axis_bounds)
    title(title_name)
    print(gcf, '-depsc2', [saveto,num2str(j-1)]);
end