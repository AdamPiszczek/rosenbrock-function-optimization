clear;
clc;
close all;
format long

global coordX coordY iter
iter = 0;

% deklaracja wykresów
figure(1)
figure(2)

% ==================== Wybranie punktu startowego ====================
startingpoint = 4; % Wybierz punkt startowy z zakresu 1-4 zmieniając zmienną startingpoint
X0 = [2 0.5; 1 -1.5; -1 -1.5; -1 0.5]; % punkty startowe
x0 = [X0(startingpoint,1) X0(startingpoint,2)]; % wartość wektora punktu startowego

figure(1) %  wykres 2D konturowy funkcji z naniesionymi trajektoriam    i
hold on
axis tight
[X, Y] = meshgrid(-4:0.001:4,-4:0.001:4);
Z = plotRosenbrock(X, Y); % obliczenie wartości funkcji Rosenbrock'a
[M, c] = contourf(X,Y,log(Z),'ShowText','on');
c.Fill = 0; % wyłączenie kolorów
c.LineWidth = 0.33; % zmniejszenie rozmiaru linii dla lepszej widoczności

% Jeżeli nie ma folderu na wykresy, stwórz go
if ~exist("./wykresy", 'dir')
       mkdir("./wykresy")
end

% ==================== Wywołanie poszczególnych metod ====================
% wybierz metode wpisując w miejsce zmiennej choice numer od 1 do 4, gdzie:
% 1 - oznacza metodę Quasi-Newton
% 2 - oznacza metodę Regionu Zaufania
% 3 - oznacza metodę Regionu Zaufania z podanym hesjanem
% 4 - oznacza metodę Neldera-Meada

choice = 4;

switch choice
    case 1
        [x, value] = optimQuasiNewton(x0, startingpoint);
    case 2
        [x, value] = optimTrustRegion(x0, startingpoint);
    case 3
        [x, value] = optimTrustRegionHessian(x0, startingpoint);
    case 4
        [x, value] = optimSimplex(x0, startingpoint);
end


% !!! aby funkcja zadziałała konieczna jest instalcja: Optimization Toolbox
% źródło: https://www.google.com/search?q=Optimization+Toolbox&sourceid=chrome&ie=UTF-8
function [x, value] = optimQuasiNewton(x0, startingpoint) % funkcja odpowiedzialna za optymalizacje rozważanej funkcji metodą Quasi-Newton
    
    global coordX coordY iter
    % ustawienie odpowiednich ustawień dla solvera
    options = optimoptions(@fminunc,'OutputFcn',@outputQuasiNewton,'Display', ...
            'iter-detailed','Algorithm','quasi-newton','MaxIterations', 1000, 'StepTolerance', ...
            1e-12, 'FunctionTolerance', 1e-12);
    [x, value] = fminunc(@rosenbrock, x0, options); % wywołanie solvera

    function stop = outputQuasiNewton(x, optimValues, state)
        stop = false; % ustawienie zmiennej stop odpowiedzialnej za działanie funkcji
        if isequal(state,'init') % przygotowanie wykresów do zaznaczania kolejnych iteracji
            figure(1)
            hold on
            title('Metoda Quasi-Newton')
            xlim([-3 3])
            ylim([-3 3])
            xlabel("wartość x")
            ylabel("wartość y")

            figure(2)
            hold on
            title('Metoda Quasi-Newton (wartości funkcji celu)');
            xlabel("ilość iteracji [n]")
            ylabel("wartość funkcji [log(y)]")
            set(gca, 'YScale', 'log') % wyświetlenie osi Y w skali logarytmicznej

        elseif isequal(state,'iter') % aktualizacja wykresów w kolejnych iteracjach
            iter = iter + 1;
            figure(1)
            coordX(iter) = x(1);
            coordY(iter) = x(2);
            figure(2)
            plot(optimValues.iteration, optimValues.fval, 'bx')

        elseif isequal(state,'done') % obsługa końcowej iteracji funkcji
            figure(1)
            plot(coordX(1,1:iter), coordY(1,1:iter), 'ro-', 'MarkerSize', 4)
            xlim([min(coordX)-0.05 max(coordX)+0.05]);
            ylim([min(coordY)-0.05 max(coordY)+0.05]);
            path = char("./wykresy/QuasiNewton1_" + startingpoint + ".png");
            saveas(gcf,path);
            fprintf('Wykres nr. 1 został zapisany\n')
            hold off
            figure(2)
            path = char("./wykresy/QuasiNewton2_" + startingpoint + ".png");
            saveas(gcf,path);
            fprintf('Wykres nr. 2 został zapisany\n')
            hold off 
    
        end
    end
end

% !!! aby funkcja zadziałała konieczna jest instalacja: Optimization Toolbox
% źródło: https://www.google.com/search?q=Optimization+Toolbox&sourceid=chrome&ie=UTF-8
function [x, value] = optimTrustRegion(x0, startingpoint) % funkcja odpowiedzialna za optymalizacje rozważanej funkcji metodą Regionu Zaufania
    
    global coordX coordY iter
    % ustawienie odpowiednich ustawień dla solvera
    options = optimoptions(@fminunc,'OutputFcn',@outputTrustRegion,'Display', ...
            'iter-detailed','Algorithm','trust-region', 'SpecifyObjectiveGradient',true, ...
            'MaxIterations', 1000, 'StepTolerance', 1e-12, 'FunctionTolerance', 1e-12);
    [x, value] = fminunc(@rosenbrockwithgrad, x0, options); % wywołanie solvera
    
    function stop = outputTrustRegion(x, optimValues, state)
        stop = false; % ustawienie zmiennej stop odpowiedzialnej za działanie funkcji
        if isequal(state,'init') % przygotowanie wykresów do zaznaczania kolejnych iteracji
            figure(1)
            hold on
            title('Metoda Regionu Zaufania bez podanego hesjanu')
            xlim([-3 3])
            ylim([-3 3])
            xlabel("wartość x")
            ylabel("wartość y")

            figure(2)
            hold on
            title('Metoda Regionu Zaufania bez podanego hesjanu (wartości funkcji celu)');
            xlabel("ilość iteracji [n]")
            ylabel("wartość funkcji [log(y)]")
            set(gca, 'YScale', 'log') % wyświetlenie osi Y w formie logarytmicznej

        elseif isequal(state,'iter') % aktualizacja wykresów w kolejnych iteracjach
            iter = iter + 1;
            figure(1)
            coordX(iter) = x(1);
            coordY(iter) = x(2);
            figure(2)
            plot(optimValues.iteration, optimValues.fval, 'bx')

        elseif isequal(state,'done') % obsługa końcowej iteracji funkcji
            fprintf('done \n')
            figure(1)
            plot(coordX(1,1:iter), coordY(1,1:iter), 'ro-', 'MarkerSize', 4)
            xlim([min(coordX)-0.05 max(coordX)+0.05]);
            ylim([min(coordY)-0.05 max(coordY)+0.05]);
            path = "./wykresy/TrustRegion1_" + startingpoint + ".png";
            saveas(gcf,path);
            fprintf('Wykres nr. 1 został zapisany\n')
            hold off
            figure(2)
            path = "./wykresy/TrustRegion2_" + startingpoint + ".png";
            saveas(gcf,path);
            fprintf('Wykres nr. 2 został zapisany\n')
            hold off 
            
        end
    end
end

% !!! aby funkcja zadziałała konieczna jest instalacja: Optimization Toolbox
% źródło: https://www.google.com/search?q=Optimization+Toolbox&sourceid=chrome&ie=UTF-8
function [x, value] = optimTrustRegionHessian(x0, startingpoint) % funkcja odpowiedzialna za optymalizacje rozważanej funkcji metodą Regionu Zaufania
    
    global coordX coordY iter
    % ustawienie odpowiednich ustawień dla solvera
    options = optimoptions(@fminunc,'OutputFcn',@outputTrustRegionHessian,'Display', ...
            'iter-detailed','Algorithm','trust-region', 'SpecifyObjectiveGradient',true, ...
            'MaxIterations', 1000, 'StepTolerance', 1e-12, 'FunctionTolerance', 1e-12, 'HessianFcn', 'objective');
    [x, value] = fminunc(@rosenbrockwithhes, x0, options); % wywołanie solvera
    
    function stop = outputTrustRegionHessian(x, optimValues, state)
        stop = false; % ustawienie zmiennej stop odpowiedzialnej za działanie funkcji
        if isequal(state,'init') % przygotowanie wykresów do zaznaczania kolejnych iteracji
            figure(1)
            hold on
            title('Metoda Regionu Zaufania z podanym hesjanem')
            xlim([-3 3])
            ylim([-3 3])
            xlabel("wartość x")
            ylabel("wartość y")

            figure(2)
            hold on
            title('Metoda Regionu Zaufania z podanym hesjanem (wartości funkcji celu)');
            xlabel("ilość iteracji [n]")
            ylabel("wartość funkcji [log(y)]")
            set(gca, 'YScale', 'log') % wyświetlenie osi Y w formie logarytmicznej

        elseif isequal(state,'iter') % aktualizacja wykresów w kolejnych iteracjach
            iter = iter + 1;
            figure(1)
            coordX(iter) = x(1);
            coordY(iter) = x(2);
            figure(2)
            plot(optimValues.iteration, optimValues.fval, 'bx')

        elseif isequal(state,'done') % obsługa końcowej iteracji funkcji
            fprintf('done \n')
            figure(1)
            plot(coordX(1,1:iter), coordY(1,1:iter), 'ro-', 'MarkerSize', 4)
            xlim([min(coordX)-0.05 max(coordX)+0.05]);
            ylim([min(coordY)-0.05 max(coordY)+0.05]);
            path = "./wykresy/TrustRegionHessian1_" + startingpoint + ".png";
            saveas(gcf,path);
            fprintf('Wykres nr. 1 został zapisany\n')
            hold off
            figure(2)
            path = "./wykresy/TrustRegionHessian2_" + startingpoint + ".png";
            saveas(gcf,path);
            fprintf('Wykres nr. 2 został zapisany\n')
            hold off 
        end
    end
end


function [x, value] = optimSimplex(x0, startingpoint) % funkcja odpowiedzialna za optymalizacje rozważanej funkcji metodą Neldera-Meada
    
    global coordX coordY iter
    
    options = optimset('OutputFcn', @outputSimplex,...
            'MaxFunEvals', 1000, 'TolX', 1e-12, 'TolFun', 1e-12);
    [x, value, exitflag, output] = fminsearch(@rosenbrock, x0, options); % wywołanie solvera
    
    function stop = outputSimplex(x, optimValues, state)
        stop = false; % ustawienie zmiennej stop odpowiedzialnej za działanie funkcji
        if isequal(state,'init') % przygotowanie wykresów do zaznaczania kolejnych iteracji
            figure(1)
            hold on
            title('Metoda Neldera-Meada')
            xlim([-3 3])
            ylim([-3 3])
            xlabel("wartość x")
            ylabel("wartość y")

            figure(2)
            hold on
            title('Metoda Neldera-Meada (wartości funkcji celu)');
            xlabel("ilość iteracji [n]")
            ylabel("wartość funkcji [log(y)]")
            set(gca, 'YScale', 'log') % wyświetlenie osi Y w formie logarytmicznej

        elseif isequal(state,'iter') % aktualizacja wykresów w kolejnych iteracjach
            iter = iter + 1;
            figure(1)
            coordX(iter) = x(1);
            coordY(iter) = x(2);
            figure(2)
            plot(optimValues.iteration, optimValues.fval, 'bx')

        elseif isequal(state,'done') % obsługa końcowej iteracji funkcji
            fprintf('done \n')
            figure(1)
            plot(coordX(1,1:iter), coordY(1,1:iter), 'ro-', 'MarkerSize', 4)
            xlim([min(coordX)-0.05 max(coordX)+0.05]);
            ylim([min(coordY)-0.05 max(coordY)+0.05]);
            path = "./wykresy/Nelder-Mead1_" + startingpoint + ".png";
            saveas(gcf,path);
            fprintf('Wykres nr. 1 został zapisany\n')
            hold off
            figure(2)
            path = "./wykresy/Nelder-Mead2_" + startingpoint + ".png";
            saveas(gcf,path);
            fprintf('Wykres nr. 2 został zapisany\n')
            hold off 
            
        end
    end
    output
end


% Wartości stałych a = 0 i b = -0.5, dla funkcji Rosenbrock'a („bananowej"):
% f(x) = (1-x+a)^2 + 100[y-b-(x-a)^2]^2

function f = plotRosenbrock(x, y) % funkcja potrzebna do stworzenia wykresu funkcji Rosenbrock'a
    f = (1-x).^2 + 100 * (y + 0.5 - x.^2).^2;
end

function f = rosenbrock(x) % właściwa funkcja użyta do optymalizacji funkcji Rosenbrock'a
    f = (1-x(1)).^2 + 100 * (x(2) + 0.5 - x(1).^2).^2;
end

function [f,g] = rosenbrockwithgrad(x) % właściwa funkcja użyta do optymalizacji funkcji Rosenbrock'a ze zdefiniowanym gradientem
    f = (1-x(1)).^2 + 100 * (x(2) + 0.5 - x(1).^2).^2;

    if nargout > 1 % Jeżeli jest więcej niż jeden argument, który zwraca funkcja
        g = [400*((-x(2)-0.495)*x(1)+x(1).^3-0.005); % obliczony gradient, według wzoru: https://wikimedia.org/api/rest_v1/media/math/render/svg/d632a346cd0677aef80d9fa32f476a5b5bf4dc58
        200*(x(2)-x(1).^2+0.5)];
    end
end

function [f, g, H] = rosenbrockwithhes(x) % funkcja użyta do optymalizacji funkcji Rosenbrock'a ze zdefiniowanym gradientem i podanym hesjanem
    f = (1 - x(1)).^2 + 100 * (x(2) + 0.5 - x(1).^2).^2;

    if nargout > 1

        g = [400 * ((-x(2) - 0.495) * x(1) + x(1).^3 - 0.005);
        200 * (x(2) - x(1).^2 + 0.5)];

        if nargout > 2

            H = [-400*(-1 * x(1).^2 + x(2) + 0.5) + 800 * x(1).^2 + 2, -400 * x(1);
            -400 * x(1), 200];
        end 
    end
end