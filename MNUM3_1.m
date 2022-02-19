% Część I

x_probe = 2:0.1:12;
y_probe = give_y(x_probe);

% Zakresy przeszukiwania
ranges = [3 8.3;8.3 11];
title = ["zakresie 3-8.3", "zakresie 8.3-11"]

%==========================================================================
zero_points_1 = [];
% Liczenie dla bisekcji
for i=1:length(ranges)
    [x_output, y_output] = bisection(ranges(i, 1),ranges(i, 2), 100);
    zero_points_1(i) = x_output(end);
end
draw_function(x_probe, y_probe, zero_points_1, ranges, "bisekcji");

zero_points_2 = [];
% Liczenie dla siecznych
for i=1:length(ranges)
    [x_output, y_output] = secant(ranges(i, 1),ranges(i, 2), 100);
    zero_points_2(i) = x_output(end);
end
draw_function(x_probe, y_probe, zero_points_2, ranges, "siecznych");

% Tworzenie tabel
for i=1:length(ranges)
    [bisect_x, bisect_y] = bisection(ranges(i, 1),ranges(i, 2), 100);
    [secant_x, secant_y] = secant(ranges(i, 1),ranges(i, 2), 100);
    draw_table(bisect_x, bisect_y, secant_x(3:end), secant_y(3:end), title(i));
end
%==========================================================================


% Metoda podająca pojedyńczą tabele
function out = draw_table(bisect_x, bisect_y, secant_x, secant_y, title)
    unify = max(length(bisect_x), length(secant_x));

    for i=length(bisect_x)+1:unify
        bisect_x(i) = 0;
        bisect_y(i) = 0;
    end
    
    for i=length(secant_x)+1:unify
        secant_x(i) = 0;
        secant_y(i) = 0;
    end
    tries = 1:1:unify;
    varNames = {'Iteration', 'Bisection x', 'Bisection y', 'Secant x', 'Secant y'};
    T = table(tries', bisect_x', bisect_y', secant_x', secant_y', 'VariableNames', varNames);
    title = "Wyniki metody bisekcji i metody siecznych w " + title
    disp(T)
    % Zapisywanie tabel do pliku .csv
    writetable(T, title + ".csv")
end

% Funkcja, której zera szukamy.
function out  = task_1(x)
    out = 2.3*sin(x) + 4*log(x+2) - 11;
end

% Metoda bisekcji
function [c_n, fc_n]= bisection(a, b, n)
    precision = 2^-8;
    a_n = a;
    b_n = b;
    c_n = [];
    fc_n = [];
    for j = 1:n
        c_n(j) = (a_n + b_n)/2;
        fc_n(j) = task_1(c_n(j));
        if abs(fc_n(j)) <= precision
            return 
        end
        if task_1(a_n)*fc_n(j) < 0
            b_n = c_n(j);
        else 
            a_n = c_n(j);
        end
    end
end

% Metoda siecznych
function [x_n, y_n] = secant(a, b, n)
    precision = 2^-8;
    x_n = [a, b];
    y_n = [task_1(a), task_1(b)];
    
    for j = 2:n
        x_n(j+1) = x_n(j) - y_n(j)*(x_n(j) - x_n(j-1))/(y_n(j)-y_n(j-1));
        y_n(j+1) = task_1(x_n(j+1));

        if abs(y_n(j+1)) <= precision
            return
        end
    end
end

% Funkcja zwracająca wektor y dla wektora x
function out = give_y(x)
    out = zeros(1,length(x));
    for j =1:length(x)
        out(j) = task_1(x(j));
    end
end

function draw_function(func_x,func_y, zero_points, ranges, method)
    figure;
    hold on;
    plot(func_x,func_y);
    title("Metoda  " + method, 'fontsize', 24)
    xlabel("X", 'fontsize', 18)
    ylabel("Y", 'fontsize', 18)
    lineSpecs = ["--r"; "--b"];

    for i = 1:length(zero_points)
        plot(zero_points(i), task_1(zero_points(i)), '*');
        xline(ranges(i,1), lineSpecs(i), i);
        xline(ranges(i,2), lineSpecs(i), i);
    end
    grid on;

end