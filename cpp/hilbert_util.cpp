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

int main()
{
    int gd = DETECT,gm,left=100,top=100;
    initgraph(&gd,&gm,NULL);

    int n = 35;
    int length = 4095;
    int scale = 1;
    int x1, y1, x2, y2;

    for(int d = 0; d < length; d++) {
//        setcolor(d*16 / length);
        setcolor(d % 16);
        d2xy(n, d, &x1, &y1);
        d2xy(n, d+1, &x2, &y2);
        line(left + x1 * scale, top + y1 * scale, left + x2 * scale, top + y2 * scale);
    }

    delay(5000);
    closegraph();
    return 0;
}
