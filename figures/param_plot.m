titles = ["Accuracy", "F1", "AUPR" ];
fields = ["acc", "f1", "aupr"];
cls = [2,4,8,16,32,52];

cl_avg = @(pen, f) [mean(pen.cl2.(f)), mean(pen.cl4.(f)), mean(pen.cl8.(f)), mean(pen.cl16.(f)), mean(pen.cl32.(f)), mean(pen.cl52.(f))];

mash.acc = mean(base.acc);
mash.f1 = mean(base.f1);
mash.aupr = mean(base.aupr);

for i = 1:3
    field = fields(i);
    subplot(1,3,i);
    avg32 = cl_avg(pen32, field);
    avg64 = cl_avg(pen64, field);
    avg128 = cl_avg(pen128, field);
    avg256 = cl_avg(pen256, field);
    avg512 = cl_avg(pen512, field);
    plot(cls, mash.(field) * ones(1,length(cls)), cls, avg32, cls, avg64, cls, avg128, cls, avg256, cls, avg512);
    title(titles(i));
    all_vals = [avg32 avg64 avg128 avg256 avg512, mash.(field)];
    ylim([min(all_vals(:)) - 0.025, max(all_vals(:)) + 0.025]);
    if i == 3
       legend({'Mashup','Flame 32', 'Flame 64', 'Flame 128', 'Flame 256', 'Flame 512'}, "Location", "southeast"); 
    end
end
sgtitle("Effects of Parameter Variation");

