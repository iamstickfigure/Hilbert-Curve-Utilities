#include<graphics.h>
#include<vector>

/*
 * https://en.wikipedia.org/wiki/Hilbert_curve
 *
 * The only documentation of libgraph that I could find:
 * http://www.programmingsimplified.com/c/graphics.h
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
    int n = 35;
    int length = 4095;
    int x, y;

    for(int d = 0; d < length; d++) {
        d2xy(n, d, &x, &y);
        hilbert.push_back(std::make_pair(x, y));
    }
    return hilbert;
}

void draw(std::vector<std::pair<int, int>>& curve) {
    int n = 35;
    int length = 4095;
    int scale = 5;
    int left=100,top=100;

    for(int d = 0; d < curve.size()-1; d++) {
//        setcolor(d*16 / length);
        setcolor(d % 16 + 1);
        line(left + curve[d].first * scale, top + curve[d].second * scale, left + curve[d+1].first * scale, top + curve[d+1].second * scale);
    }
}

int main()
{
    int gd = DETECT,gm,left=100,top=100;
    initgraph(&gd,&gm,NULL);

    std::vector<std::pair<int, int>> curve = generate();
    draw(curve);

    delay(50000);
    closegraph();
    return 0;
}
