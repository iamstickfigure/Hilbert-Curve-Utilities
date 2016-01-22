#include<vector>
#include<cmath>
#include<limits>
#include<math.h>
//#include<CImg.h>
#include "CImg.h"
#include "portaudio.h"
#include "sndfile.hh"
#include "sndfile.h"
using namespace cimg_library;

#define NUM_SECONDS   (50)
//#define SAMPLE_RATE   (8000)
#define SAMPLE_RATE   (44100)
#define FRAMES_PER_BUFFER  (64)

#ifndef M_PI
#define M_PI  (3.14159265)
#endif

//#define TABLE_SIZE   (20000)
#define TABLE_SIZE   (200000)

typedef struct
{
    float data[TABLE_SIZE];
    int left_phase;
    int right_phase;
    char message[20];
} soundData;

// paCallback and StreamFinished partially from PortAudio examples (paex_sine.cpp)
/* This routine will be called by the PortAudio engine when audio is needed.
** It may called at interrupt level on some machines so don't do anything
** that could mess up the system like calling malloc() or free().
*/
static int paCallback( const void *inputBuffer, void *outputBuffer,
                           unsigned long framesPerBuffer,
                           const PaStreamCallbackTimeInfo* timeInfo,
                           PaStreamCallbackFlags statusFlags,
                           void *userData )
{
    soundData *sound = (soundData*)userData;
    float *out = (float*)outputBuffer;  // This output buffer is where to send the sound data.
    unsigned long i;

    (void) timeInfo; /* Prevent unused variable warnings. */
    (void) statusFlags;
    (void) inputBuffer;

    for( i=0; i<framesPerBuffer; i++ )
    {
        *out++ = sound->data[sound->left_phase];  /* left */
        *out++ = sound->data[sound->right_phase];  /* right */
        sound->left_phase += 1;
        if( sound->left_phase >= TABLE_SIZE ) sound->left_phase -= TABLE_SIZE;
        sound->right_phase += 1;
        if( sound->right_phase >= TABLE_SIZE ) sound->right_phase -= TABLE_SIZE;
    }

    return paContinue;
}

/*
 * This routine is called by portaudio when playback is done.
 */
static void StreamFinished( void* userData )
{
    soundData *sound = (soundData *) userData;
    printf( "Stream Completed: %s\n", sound->message );
}

int playSound(soundData& sound) {

    PaStreamParameters outputParameters;
    PaStream *stream;
    PaError err;
    err = Pa_Initialize();
    if( err != paNoError ) goto error;

    outputParameters.device = Pa_GetDefaultOutputDevice(); /* default output device */
    if (outputParameters.device == paNoDevice) {
        fprintf(stderr,"Error: No default output device.\n");
        goto error;
    }
    outputParameters.channelCount = 2;       /* stereo output */
    outputParameters.sampleFormat = paFloat32; /* 32 bit floating point output */
    outputParameters.suggestedLatency = Pa_GetDeviceInfo( outputParameters.device )->defaultLowOutputLatency;
    outputParameters.hostApiSpecificStreamInfo = NULL;

    err = Pa_OpenStream(
            &stream,
            NULL, /* no input */
            &outputParameters,
            SAMPLE_RATE,
            FRAMES_PER_BUFFER,
            paClipOff,      /* we won't output out of range samples so don't bother clipping them */
            paCallback,
            &sound );
    if( err != paNoError ) goto error;

    sprintf( sound.message, "No Message" );
    err = Pa_SetStreamFinishedCallback( stream, &StreamFinished );
    if( err != paNoError ) goto error;

    err = Pa_StartStream( stream );
    if( err != paNoError ) goto error;

    printf("Play for %d seconds.\n", NUM_SECONDS );
    Pa_Sleep( NUM_SECONDS * 1000 );

    err = Pa_StopStream( stream );
    if( err != paNoError ) goto error;

    err = Pa_CloseStream( stream );
    if( err != paNoError ) goto error;

    Pa_Terminate();
    printf("Test finished.\n");

    return err;
    error:
    Pa_Terminate();
    fprintf( stderr, "An error occured while using the portaudio stream\n" );
    fprintf( stderr, "Error number: %d\n", err );
    fprintf( stderr, "Error message: %s\n", Pa_GetErrorText( err ) );
    return err;
}

/*
 * https://en.wikipedia.org/wiki/Hilbert_curve
 *
 * Compile with:
 * g++11 hilbert_util.cpp -o hilbert_util -O2 -L/usr/X11R6/lib -lm -lpthread -lX11
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

typedef std::vector<std::pair<int, int>> Path;

Path generate() {
    Path hilbert;
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

void soundMap(int d, int max, unsigned char r, unsigned char g, unsigned char b, float (&data)[TABLE_SIZE]) {
    // SAMPLE_RATE samples per second
    // TABLE_SIZE total samples
    // Map the f values between 20 Hz to 20 kHz
    // sin(x*f*pi/(SAMPLE_RATE))
    float f = 20 + d*20000/max;
    float a = sqrt(pow(r,2) + pow(g,2) + pow(b,2))/255.0;
    for(int i = 0; i < TABLE_SIZE; ++i) {
        data[i] += (float) sin(i * f * M_PI / (SAMPLE_RATE)) * a;
        if(data[i] == std::numeric_limits<float>::infinity())
            printf("data[%d] overflowed\n", i);
    }
}

void generateSound(CImg<unsigned char> &image, Path curve) {
    soundData sound = { {}, 0, 0, {} };
//    printf("sound vals are: %f, %f, %f\n", sound.data[0], sound.data[5], sound.data[10]);
    for(int d = 0; d < curve.size(); ++d) {
        if(d % 10000 == 0)
            printf("%d out of %d\n", d, (int)curve.size());
        soundMap(d, curve.size(), image(curve[d].first, curve[d].second, 0), image(curve[d].first, curve[d].second, 1),
                 image(curve[d].first, curve[d].second, 2), sound.data);
    }
    SndfileHandle file("out.wav", SFM_WRITE, SF_FORMAT_WAV | SF_FORMAT_FLOAT, 2, SAMPLE_RATE);
    file.write(sound.data, TABLE_SIZE) ;
//    playSound(sound);
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

void draw(Path& curve) {
    int scale = 1;
    int left=0,top=0;
//    CImg<unsigned char> image("color_bars_1121.jpg");//400,400,1,3,0);
    CImg<unsigned char> image("square_leiss256.jpg");//400,400,1,3,0);

//    generateSound(image, curve);

    for(int d = 0; d < curve.size()-1; d++) {
        unsigned char color[3] = {0, 0, 0};
        colorbar(d, curve.size(), color);
        image.draw_line(left + curve[d].first * scale, top + curve[d].second * scale, left + curve[d+1].first * scale,
                        top + curve[d+1].second * scale, color, 0.5);
    }
    CImgDisplay draw_disp(image);

    while (!draw_disp.is_closed()) {
//        if (draw_disp.mouse_x()>=0 && draw_disp.mouse_y()>=0) { // Mouse pointer is over the image
//
//            const int x = draw_disp.mouse_x(), y = draw_disp.mouse_y();
//            printf("(%d, %d) - %d, %d, %d\n", x, y, image(x, y, 0), image(x, y, 1), image(x, y, 2));
//        }
        draw_disp.wait();
    }
}

int main() {
    Path curve = generate();
    draw(curve);

    return 0;
}
