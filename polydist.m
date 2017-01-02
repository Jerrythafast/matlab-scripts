function distances = polydist(p, points, draw)
%POLYDIST Distance of points to a polynomial.
%   POLYDIST(P, POINTS) returns the Euclidean distances of the
%   points in matrix POINTS to the polynomial whose coefficients
%   are the elements of vector P.
%
%   POINTS must be an n-by-2 matrix in which each row corresponds
%   to the [x, y] coordinates of a point in a two-dimensional space.
%
%   The return value is a column vector of length n containing the
%   Euclidean distance of each input point to the polynomial.
%
%   POLYDIST(P, POINTS, 1) additionally plots the points,
%   the polynomial, and the shortest distance lines.
%
%   See also POLYFIT, POLYVAL.
%
%   Author: Jerry Hoogenboom
%   Release: 1.0
%   Release date: 2016-12-28

if nargin < 3;
    draw = 0;
end;
if ~isvector(p);
    error('P must be a vector of polynomial coefficients.');
end;
if ~ismatrix(points) || size(points, 2) ~= 2;
    error('Points must be a 2-D array/matrix with two columns.');
end;

if draw == 1;
    % Create figure and draw the points.
    figure; hold on; grid on;
    plot(points(:, 1), points(:, 2), '*');
end;

distances = zeros([size(points, 1), 1]);
for i = 1:size(points, 1);
    % Squared Euclidean distance of point (x0, y0) to the curve is:
    % distpoly = (x - x0)^2 + (y - y0)^2
    x_part = [1, -2 * points(i, 1), points(i, 1)^2];  % (x - x0)^2
    y_part = p;  % y
    y_part(end) = y_part(end) - points(i, 2);  % y - y0
    y_part = conv2(y_part, y_part);  % (y - y0)^2

    % Make sure the x and y parts can be added together (are same length).
    if length(y_part) < 3;
        y_part = padarray(y_part, [0, 3-length(y_part)], 'pre');
    elseif length(y_part) > 3;
        x_part = padarray(x_part, [0, length(y_part)-3], 'pre');
    end;
    distpoly = x_part + y_part;

    % Closest point on the polynomial is a real (non-complex) root of the
    % derivative of the distance polynomial.
    r = roots(polyder(distpoly));
    r(imag(r) ~= 0) = [];  % Gets rid of complex roots.

    % One of the roots is a minimizer of the distance.  Since distpoly
    % gives squared Euclidean distances, we take the square root of it.
    distances(i) = sqrt(min(polyval(distpoly, r)));

    if draw == 1;
        % Draw shortest distance lines.
        x = r(find(distances(i) == sqrt(polyval(distpoly, r)), 1));
        plot([points(i, 1), x], [points(i, 2), polyval(p, x)], '-r');
        text(mean([points(i, 1), x]), ...
            mean([points(i, 2), polyval(p, x)]), ...
            num2str(distances(i)), ...
            'FontSize', 6);
    end;
end;

if draw == 1;
    % Set both axes to the same limits and draw the polynomial.
    a = axis();
    axis(repmat([min(a(1), a(3)), max(a(2), a(4))], 1, 2));
    xgrid = linspace(min(a(1), a(3)), max(a(2), a(4)));
    plot(xgrid, polyval(p, xgrid));
end;

end