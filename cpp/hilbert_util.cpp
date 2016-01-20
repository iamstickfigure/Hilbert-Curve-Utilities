#include<vector>
#include<math.h>
//#include<CImg.h>
#include "CImg.h"
using namespace cimg_library;

/*
 * https://en.wikipedia.org/wiki/Hilbert_curve
 *
 */

//rotate/flip a quadrant appropriately
void rot(int n, int *x, int *y, int rx, int ry) {
    if (ry == 0) {
        if (rx == 1) {
            *x = n-1 - *x;
            *y = n-1 - *y;
        }

        //Swap x and y
        int t  = *x;
        *x = *y;
        *y = t;
    }
}

//convert (x,y) to d
int xy2d (int n, int x, int y) {
    int rx, ry, s, d=0;
    for (s=n/2; s>0; s/=2) {
        rx = (x & s) > 0;
        ry = (y & s) > 0;
        d += s * s * ((3 * rx) ^ ry);
        rot(s, &x, &y, rx, ry);
    }
    return d;
}

//convert d to (x,y)
void d2xy(int n, int d, int *x, int *y) {
    int rx, ry, s, t=d;
    *x = *y = 0;
    for (s=1; s<n; s*=2) {
        rx = 1 & (t/2);
        ry = 1 & (t ^ rx);
        rot(s, x, y, rx, ry);
        *x += s * rx;
        *y += s * ry;
        t /= 4;
    }
}

std::vector<std::pair<int, int>> generate() {
    std::vector<std::pair<int, int>> hilbert;
    int nlog = 8;
    int n = (int)pow(2, nlog);
    int length = n*n;
    int x, y;

    for(int d = 0; d < length; d++) {
        d2xy(n, d, &x, &y);
        hilbert.push_back(std::make_pair(x, y));
    }
    return hilbert;
}

// http://stackoverflow.com/questions/2374959/algorithm-to-convert-any-positive-integer-to-an-rgb-value
// Based off of the graph provided.
void colorbar(int d, int max, unsigned char (&color)[3]) {
    int interval = max/8;
    if(d < interval) {
        color[0] = 0;
        color[1] = 0;
        color[2] = 127 + 127*d/interval;
    }
    else if(d < interval*3) {
        color[0] = 0;
        color[1] = 127*(d-interval)/interval;
        color[2] = 255;
    }
    else if(d < interval*5) {
        color[0] = 127*(d-3*interval)/interval;
        color[1] = 255;
        color[2] = 255 - 127*(d-3*interval)/interval;
    }
    else if(d < interval*7) {
        color[0] = 255;
        color[1] = 255 - 127*(d-5*interval)/interval;
        color[2] = 0;
    }
    else {
        color[0] = 255 - 127*(d-7*interval)/interval;
        color[1] = 0;
        color[2] = 0;
    }
}

void draw(std::vector<std::pair<int, int>>& curve) {
    int scale = 1;
    int left=10,top=10;
    CImg<unsigned char> image(400,400,1,3,0);
    for(int d = 0; d < curve.size()-1; d++) {
        unsigned char color[3] = {0, 0, 0};
        colorbar(d, curve.size(), color);
        image.draw_line(left + curve[d].first * scale, top + curve[d].second * scale, left + curve[d+1].first * scale, top + curve[d+1].second * scale, color);
    }
    CImgDisplay draw_disp(image);

    while (!draw_disp.is_closed()) {
        draw_disp.wait();
    }
}

int main()
{
    std::vector<std::pair<int, int>> curve = generate();
    draw(curve);

    return 0;
}
