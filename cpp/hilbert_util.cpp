#include<vector>
#include<cmath>
#include<limits>
#include<math.h>
#include<thread>
#include "CImg.h"
#include "portaudio.h"
#include "sndfile.hh"
#include "sndfile.h"
using namespace cimg_library;

#define LOG_N       (9)

#define NUM_SECONDS   (50)
//#define SAMPLE_RATE   (8000)
#define SAMPLE_RATE   (44100)
//#define SAMPLE_RATE   (88200)

#define FRAMES_PER_BUFFER  (64)

#define NUM_THREADS     (10)

#define LOW_FREQUENCY   (500.0)
#define HIGH_FREQUENCY  (1760.0)
// 20.0 to 20000.0   or  110.0 to 1760.0

#define GRAPHICS_ON   // Comment out to turn off graphics
#define SOUND_OUTPUT_ON   // Comment out to turn off all sound output
//#define VERBOSE

#ifndef M_PI
#define M_PI  (3.14159265)
#endif

//#define TABLE_SIZE   (20000)
#define TABLE_SIZE   (200000)

#ifdef SOUND_OUTPUT_ON
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

#endif

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
    int n = (int)pow(2, LOG_N);
    int length = n*n;
    int x, y;

    for(int d = 0; d < length; d++) {
        d2xy(n, d, &x, &y);
        hilbert.push_back(std::make_pair(x, y));
    }
    return hilbert;
}

typedef struct {
    float h,s,v;
} hsv_t;

hsv_t rgb_to_hsv(int R, int G, int B) {
    using std::max;
    using std::min;
    float h, s, v, r = R/255.0, g = G/255.0, b = B/255.0;
    float cmax = max(max(r, g), b);
    float cmin = min(min(r, g), b);
    float cdelta = cmax - cmin;

    if(R == G && G == B)
        h = 0;
    else if(r > b) {
        if(r > g) // cmax == r
            h = 60.0*fmod((g-b)/cdelta, 6);
        else // cmax == g
            h = 60.0*((b-r)/cdelta + 2);
    }
    else // cmax == b
        h = 60.0*((r-g)/cdelta + 4);

    if(R == 0 && G == 0 && B == 0) // All of them are zero
        s = 0;
    else
        s = cdelta/cmax;

    v = cmax;

    return hsv_t{h, s, v};
}

float frequencyMap(int d, int max) {
    return LOW_FREQUENCY + d*HIGH_FREQUENCY/max;
}

float frequencySaturationMap(float s) {
    return LOW_FREQUENCY + s*HIGH_FREQUENCY;
}

#ifdef SOUND_OUTPUT_ON  // If sound output is off, don't compile this


//float soundMap(int d, int max, unsigned char r, unsigned char g, unsigned char b, float (&data)[TABLE_SIZE]) {
//    // SAMPLE_RATE samples per second
//    // TABLE_SIZE total samples
//    // Map the f values between 20 Hz to 20 kHz (Maybe 110 to 1760 instead?)
//    // sin(x*f*pi/(SAMPLE_RATE))
//    float max_amp = 0;
//    float f = frequencyMap(d, max));
//    float a = sqrt(pow(r,2) + pow(g,2) + pow(b,2))/255.0;
//    for(int i = 0; i < TABLE_SIZE; ++i) {
//        data[i] += (float) sin((i-10000) * f * M_PI / (SAMPLE_RATE)) * a;
//        if(data[i] > max_amp)
//            max_amp = data[i];
//        if(data[i] == std::numeric_limits<float>::infinity())
//            printf("data[%d] overflowed\n", i);
//    }
//
//    #ifdef VERBOSE
//        printf("Frequency: %f Hz    Amplitude: %f\n", f, a);
//    #endif
//
//    return max_amp; // Returns the largest amplitude present
//}

//void soundPartition(int t, Path &curve, CImg<unsigned char> &image, soundData &sound, float &max_amp, int &d) {
//    int interval = curve.size() / NUM_THREADS;
//    for(d = interval*t; d < interval*(t+1); ++d) {
//        float max_amp_local = soundMap(d, curve.size(), image(curve[d].first, curve[d].second, 0),
//                                       image(curve[d].first, curve[d].second, 1),
//                                       image(curve[d].first, curve[d].second, 2), sound.data);
////        float max_amp_local = soundMap(d, curve.size(), 255, 255, 255, sound.data); // White test (Full on all frequencies)
////        float max_amp_local = (d==(int)curve.size()/100)?soundMap(d, curve.size(), 255, 255, 255, sound.data):0;
//        if(max_amp_local > max_amp)
//            max_amp = max_amp_local;
//    }
//}

float pulse(int location, float f, float a, float width, float (&data)[TABLE_SIZE]) {
    // high frequency sine (f) * (low frequency sine + 1)
    float max_amp = 0;
    int offset = 0;
    for (int i = (int)((location > width/2)?(location-width/2):0); i < (int)width/2+location && i < TABLE_SIZE; ++i) {
        data[i] += (float) sin((i-offset) * f * M_PI / (SAMPLE_RATE)) * a * (cos(2*M_PI*i/width)+1)/2;
        if(data[i] > max_amp)
            max_amp = data[i];
        if(data[i] == std::numeric_limits<float>::infinity())
            printf("data[%d] overflowed\n", i);
    }
    return max_amp;
}

float pulsePartition(int t, Path &curve, CImg<unsigned char> &image, soundData &sound, float &max_amp, int &d) {
    int interval = curve.size() / NUM_THREADS;
    float f, a, local_max = 0;
    int location, r, g, b;
    hsv_t hsv;

    for(d = interval*t; d < interval*(t+1); ++d) {
        r = image(curve[d].first, curve[d].second, 0);
        g = image(curve[d].first, curve[d].second, 1);
        b = image(curve[d].first, curve[d].second, 2);
        hsv = rgb_to_hsv(r, g, b);
        f = frequencySaturationMap(hsv.s);
        a = sqrt(pow(r,2) + pow(g,2) + pow(b,2));
//        f = frequencySaturationMap(1.0); // Max saturation test
//        a = sqrt(pow(255,2) + pow(255,2) + pow(255,2));  // White test
        location = (int)(d*(float)TABLE_SIZE/(float)curve.size());
        local_max = pulse(location, f, a, 1000, sound.data);
        if(local_max > max_amp)
            max_amp = local_max;
    }
}

void generateSound(CImg<unsigned char> &image, Path &curve, int (&progress)[NUM_THREADS]) {
    soundData sound = { {}, 0, 0, {} };
    float max_amp = 0;
    std::thread threads[NUM_THREADS];
//    for(int d = 0; d < curve.size(); ++d) {
//        if(d % 10000 == 0)
//            printf("%d out of %d\n", d, (int)curve.size());
//        float max_amp_local = soundMap(d, curve.size(), image(curve[d].first, curve[d].second, 0),
//                                       image(curve[d].first, curve[d].second, 1),
//                                       image(curve[d].first, curve[d].second, 2), sound.data);
////        float max_amp_local = soundMap(d, curve.size(), 255, 255, 255, sound.data); // White test (Full on all frequencies)
////        float max_amp_local = (d==(int)curve.size()/100)?soundMap(d, curve.size(), 255, 255, 255, sound.data):0;
//        if(max_amp_local > max_amp)
//            max_amp = max_amp_local;
//    }
    for(int t = 0; t < NUM_THREADS; ++t) {
//        threads[t] = std::thread(soundPartition, t, std::ref(curve), std::ref(image),
//                                 std::ref(sound), std::ref(max_amp), std::ref(progress[t]));
        threads[t] = std::thread(pulsePartition, t, std::ref(curve), std::ref(image),
                                 std::ref(sound), std::ref(max_amp), std::ref(progress[t]));
    }
    for(int t = 0; t < NUM_THREADS; ++t) {
        threads[t].join();
    }
    for(int i = 0; i < TABLE_SIZE; ++i) {
        sound.data[i] /= max_amp;
    }
    SndfileHandle file("out.wav", SFM_WRITE, SF_FORMAT_WAV | SF_FORMAT_FLOAT, 2, SAMPLE_RATE);
    file.write(sound.data, TABLE_SIZE) ;
    playSound(sound);
}

#endif

// http://stackoverflow.com/questions/2374959/algorithm-to-convert-any-positive-integer-to-an-rgb-value
// Based off of the graph provided.

#ifdef GRAPHICS_ON // If graphics are turned off, don't compile this.

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

void draw(Path& curve, CImg<unsigned char> &image, int (&progress)[NUM_THREADS]) {
    int scale = 1;
    int left=0,top=0;
    CImg<unsigned char> freq(400, 400, 1, 3, 0);

//    for(int d = 0; d < curve.size()-1; d++) {
//        unsigned char color[3] = {0, 0, 0};
//        colorbar(d, curve.size(), color);
//        image.draw_line(left + curve[d].first * scale, top + curve[d].second * scale, left + curve[d+1].first * scale,
//                        top + curve[d+1].second * scale, color, 0.5);
//    }
    CImgDisplay draw_disp(image);
    CImgDisplay freq_disp(freq);
    const unsigned char color[] = { 255,255,255 };

    while (!draw_disp.is_closed()) {
        freq.fill(0);
        char* message = new char[50];
        if (draw_disp.mouse_x()>=0 && draw_disp.mouse_y()>=0) { // Mouse pointer is over the image

            const int x = draw_disp.mouse_x(), y = draw_disp.mouse_y();
            int d = xy2d((int)pow(2, LOG_N), x, y);
            hsv_t hsv = rgb_to_hsv(image(x, y, 0), image(x, y, 1), image(x, y, 2));
            sprintf(message, "(%d, %d) d=%d hsv=%4.2f,%4.2f,%4.2f f=%4.2fHz\n", x, y, d, hsv.h, hsv.s, hsv.v, frequencyMap(d, curve.size()));
            freq.draw_text(10, 10, message, color, 0, 0.8f, 24);
        }
        int interval = curve.size() / NUM_THREADS;
        for(int t = 0; t < NUM_THREADS; ++t) {
            sprintf(message, "Thread %d: %4.2f%%\n", t, 100.0*(progress[t]-interval*t)/(float)interval);
            freq.draw_text(10, 34 + 24*t, message, color, 0, 0.8f, 24);
        }
        freq.display(freq_disp);
    }
}
#endif

int main() {
    Path curve = generate();
//    CImg<unsigned char> image("color_bars_512.jpg");
    CImg<unsigned char> image("square_leiss.jpg");

    int progress[NUM_THREADS];

    #ifdef SOUND_OUTPUT_ON
        #ifdef GRAPHICS_ON
            std::thread graphics(draw, std::ref(curve), std::ref(image), std::ref(progress));
            std::thread sounds(generateSound, std::ref(image), std::ref(curve), std::ref(progress));

            sounds.join();
            graphics.join();
        #else
            generateSound(image, curve, progress);
        #endif
    #else
    #ifdef GRAPHICS_ON
        draw(curve, image);
    #endif
    #endif

    return 0;
}
