
for i=1:15
    load(sprintf('f%02d.mat', i));
    if exist('R100')
        filename = sprintf('./cdatafiles/F%d-R100.txt', i);
        dlmwrite(filename, R100, 'delimiter', ',', 'precision', 200);
    end

    if exist('R25')
        filename = sprintf('./cdatafiles/F%d-R25.txt', i);
        dlmwrite(filename, R25, 'delimiter', ',', 'precision', 200);
    end

    if exist('R50')
        filename = sprintf('./cdatafiles/F%d-R50.txt', i);
        dlmwrite(filename, R50, 'delimiter', ',', 'precision', 200);
    end

    if exist('p')
        filename = sprintf('./cdatafiles/F%d-p.txt', i);
        dlmwrite(filename, p, 'delimiter', ',', 'precision', 200);
    end

    if exist('s')
        filename = sprintf('./cdatafiles/F%d-s.txt', i);
        dlmwrite(filename, s, 'delimiter', ',', 'precision', 200);
    end

    if exist('w')
        filename = sprintf('./cdatafiles/F%d-w.txt', i);
        dlmwrite(filename, w, 'delimiter', ',', 'precision', 200);
    end

    if exist('xopt')
        filename = sprintf('./cdatafiles/F%d-xopt.txt', i);
        dlmwrite(filename, xopt, 'delimiter', ',', 'precision', 200);
    end
end
