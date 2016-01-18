from graphics import *
# https://en.wikipedia.org/wiki/Hilbert_curve
#
# Based on Python 2.7
# Requires graphics.py and package python-tk


def point_to_d(point, n):
    s = n/2
    d = 0
    x = point.x
    y = point.y
    while s > 0:
        rx = (x & s) > 0
        ry = (y & s) > 0
        d += s * s * ((3 * rx) ^ ry)
        x, y = rotate(x, y, rx, ry, s)
        s /= 2
    return d


def d_to_point(d, n):
    s = 1
    temp = d
    x = 0
    y = 0
    while s < n:
        rx = 1 & (temp/2)
        ry = 1 & (temp ^ rx)
        x, y = rotate(x, y, rx, ry, s)
        x += s * rx
        y += s * ry
        temp /= 4
        s *= 2
    return Point(x, y)


def rotate(x, y, rx, ry, n):
    if ry == 0:
        if rx == 1:
            x = n-1 - x
            y = n-1 - y
        return y, x
    else:
        return x, y


def transform(point, coeff=5, xo=5, yo=50):
    return Point(point.x * coeff + xo, point.y * coeff + yo)


def main():
    win = GraphWin('Hilbert Curve', 400, 400)

    n = 35
    length = 4095

    for d in xrange(0, length):  # xrange(0, pow(2, n-1)-1):  #
        line = Line(transform(d_to_point(d, n)), transform(d_to_point(d+1, n)))
        line.setWidth(1)
        # luminance = d*256/length
        # line.setFill(color_rgb(luminance, luminance, luminance))
        line.draw(win)

    message = Text(Point(win.getWidth()/2, 20), 'Click anywhere to quit.')
    message.draw(win)
    win.getMouse()
    win.close()

main()
