clear

set(0,'defaultaxeslinewidth',1)
set(0,'defaultlinelinewidth',2)
set(0,'defaultpatchlinewidth',1)

load('data_minang')
bin_width = 1;
folder_name = './Histograms3Dcomparison/';
dims = [256 256 256];
for choice = 2:4
    file_name = ['simulation3dcomparison-',num2str(choice),'-',num2str(dims(1)),'x',num2str(dims(2)),'x',num2str(dims(3))];
    working_data = eval(['minang_',num2str(choice)])/(2*pi)*360;
    h = histogram(working_data,'BinWidth',bin_width,'Normalization','pdf');
    Values = h.Values;
    BinWidth = h.BinWidth;
    edgeBins = h.BinEdges;
    edgeCenters = edgeBins(1:end-1)+BinWidth/2;
    error = sqrt(h.BinCounts)/(length(working_data)*BinWidth);
    clf
    hold on
    errorbar(edgeCenters,Values,error,'-','LineWidth',1.5)
    xMackenzie = 0:0.1:63;
    plot(xMackenzie,mackenzie_function_degrees(xMackenzie))
    axis([0 63 0 0.05])
    title('Orientation Angle Distribution')
    hold off
    saveto = [folder_name,file_name,'-ODF_t0'];
    print(gcf, '-depsc2', saveto);
end
