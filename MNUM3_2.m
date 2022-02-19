% Część II

p = [-1, 2.5, 2.5, 1, 0.5];
x_probe = -2:0.1:4;
y_probe = give_y(x_probe);


ranges = [-400 0 300];
x_result = [];

%==========================================================================
% Odnajdowanie wszystkich pierwiatsków rzeczywistych

x_result(1) =  MM1(ranges, p, 10000);
p = deflate_simple(p, x_result(1));
x_result(2) =  MM1(ranges, p, 10000);
p = deflate_simple(p, x_result(2));
[x_result(3), x_result(4)] = solv_quadra(p);
disp(x_result')

draw_function(x_result)
%==========================================================================


% Deflacja - Schemat Hornera(prosty)
function out = deflate_simple(p, x0)
    len = length(p);
    for i = 2: len
        p(1, i) = p(1, i) + p(1, i - 1) * x0;
    end
    out = p(1, 1 : (len - 1));
end

% Rozwiązywanie równania kwadratowego
function [x1, x2] = solv_quadra(p)
    a = p(1, 1);
    b = p(1, 2);
    c = p(1, 3);

    d = sqrt(b^2 - 4 * a * c);

    x1 = (-b - d)/(2*a);
    x2 = (-b + d)/(2*a);
end

% Metoda mullera MM1
function x_new = MM1(x, p, N_iter)
    precision =  2^-8;
    x_3 = x;
    z = zeros(1, 2);
    i = 0;
    while i <= N_iter
        
        f_3 = [polyval(p, x_3(1)), polyval(p, x_3(2)), polyval(p, x_3(3))];
        % Rozwiązywanie układu
        z(1) = x_3(1) - x_3(3);
        z(2) = x_3(2) - x_3(3);
    
        a = [z(1)^2 z(1); z(2)^2 z(2)] ;
        b = [(f_3(1) - f_3(3)); (f_3(2) - f_3(3))];
        
        % gdzie a_b(1) to a, a_b(2) to b
        a_b = linsolve(a, b);
    
        c = f_3(3);
        root = sqrt(a_b(2)^2 - 4*a_b(1)*c);

        % Wybieranie rozwiązania
        if abs(a_b(2) + root) >= abs(a_b(2) - root)
            x_new =  x_3(3) + -2*c/(a_b(2) + root);
        else 
            x_new =  x_3(3) + -2*c/(a_b(2) - root);
        end
        
        if abs(polyval(p, x_new))< precision
            return
        end
        % Podmiana z najdalej znajdującym się punktem od 'x_new'
        [~, index] = max([abs(x_3(1)-x_new), abs(x_3(2)-x_new), abs(x_3(3)-x_new)]);
        if index == 1
            x_3(1) = x_new;
        elseif index == 2
            x_3(2) = x_new;
        else 
            x_3(3) = x_new;
        end
        i = i+1;
    end
    
end

% Obliczanie y dla podanych wartości funkcji x
function out = f(x)
    out = -1*x^4 + 2.5*x^3 + 2.5*x^2 + 1*x + 0.5;
end

function out = give_y(x)
    out = zeros(1,length(x));
    for j =1:length(x)
        out(j) = f(x(j));
    end
end

function draw_function(test_x)
    figure;
    hold on;
    plot(test_x, '*');
    xline(0);
    yline(0);
    title("Pierwiastki wielomianu f(x)", 'fontsize', 24 )
    xlabel("Część rzeczywista",'fontsize', 18)
    ylabel("Część urojona", 'fontsize', 18)
    grid on;
end