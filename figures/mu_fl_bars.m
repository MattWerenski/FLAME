sizes = ["low", "med", "high"];
titles = ["11-30", "31-100", "101-300"];
mins = [0.3, 0.3, 0.3]; % yeast [0.3, 0.3, 0.3] - human [0.15, 0.2, 0.3]
maxs = [0.8, 0.75, 0.95]; % yeast [0.8, 0.75, 0.95] - human [0.5, 0.6, 0.65]
onts = ["bp", "mf", "cc"];
ont_names = ["BP","MF","CC"];
for k = 1:3
for j = 1:3
    subplot(3,3,(k-1)*3+j);
    % Example data as before
    cat1 = flame_stats.(onts(k)).(sizes(j));
    cat2 = mash_stats.(onts(k)).(sizes(j));

    model_series = [
        mean(cat1.acc) mean(cat2.acc); 
        mean(cat1.f1) mean(cat2.f1); 
        mean(cat1.auprc) mean(cat2.auprc);
    ];
    model_error = [
        std(cat1.acc) std(cat2.acc); 
        std(cat1.f1) std(cat2.f1); 
        std(cat1.auprc) std(cat2.auprc);
    ];
    b = bar(model_series, 'grouped');
    hold on
    % Calculate the number of bars in each group
    nbars = size(model_series, 2);
    % Get the x coordinate of the bars
    x = [];
    for i = 1:nbars
        x = [x ; b(i).XEndPoints];
    end
    % Plot the errorbars
    errorbar(x',model_series,model_error,'k','linestyle','none');
    set(gca, 'XTickLabel', {'Acc' 'F1' 'AUPR'});
    title(ont_names(k)+" "+titles(j));
    if j==3 && k == 3
        legend({"FLAME", "Mashup"}, 'Location', 'southeast');
    end
    ylim([mins(k), maxs(k)]);
    axis square
    hold off
end
end
suptitle('Yeast Results');