#include <math.h>
#include "hmatrix.h"


#define true 1
#define false 0


void hw_u32_array_to_float(unsigned int input[MSGLEN], float *result);
float hw_u32_to_float(unsigned int val);


void calculateBits(float p0[MSGLEN], float p1[MSGLEN], unsigned int *bits, short *c){
    int change = 0;
    for (int i = 0; i < MSGLEN; ++i) {
        if (p0[i] > p1[i]){
            if (bits[i] == 1) {
                change += 1;
            }
            bits[i] = 0;
        } else {
            if (bits[i] == 0) {
                change += 1;
            }
            bits[i] = 1;
        }
    }
    *c = change;
}

void sumProductDecoding(unsigned int max_iter, unsigned int end, unsigned int intmsg[MSGLEN], unsigned int output[MSGLEN], unsigned int *iterations) {
    #pragma HLS INTERFACE s_axilite port=return bundle=CRTLS
    #pragma HLS INTERFACE s_axilite port=max_iter bundle=CRTLS
    #pragma HLS INTERFACE s_axilite port=end bundle=CRTLS
    #pragma HLS INTERFACE s_axilite port=intmsg bundle=CRTLS
    #pragma HLS INTERFACE s_axilite port=output bundle=CRTLS
    #pragma HLS INTERFACE s_axilite port=iterations bundle=CRTLS

    float p0[MSGLEN];
    float p1[MSGLEN];
    float message[MSGLEN];
    float P0[MSGLEN][HROWS];
    float P1[MSGLEN][HROWS];
    short i;
    short i2;
    short j;
    short j2;
    short change = 0;
    // init as max+1 so you know how it finished.
    *iterations = max_iter+1;
    int same = 0;

    // convert message:
    hw_u32_array_to_float(intmsg, message);


    // init arrays
    for (i = 0; i < MSGLEN; ++i) {
        p1[i] = 1 / (1 + exp((2 * message[i]) / SIGMA));
        p0[i] = 1 - p1[i];
    }

    // get original message
    calculateBits(p0, p1, output, &change);
    // init P matrix's
    for (i = 0; i < HROWS; ++i) {
        for ( j = 0; j < MSGLEN; ++j) {

            P0[j][i] = p0[j];
            P1[j][i] = p1[j];
        }
    }

    // actual algorithm
    for (int iter = 0; iter < max_iter; ++iter) {
        // DeltaP
    	float deltaP[MSGLEN][HROWS];
        for ( i = 0; i < HROWS; ++i) {
            for ( j = 0; j < MSGLEN; ++j) {
                deltaP[j][i] = P0[j][i] - P1[j][i];
            }
        }

        // Q:s
        float deltaQ[MSGLEN][HROWS] = {0};
        float Q0[MSGLEN][HROWS] = {0};
        float Q1[MSGLEN][HROWS] = {0};
        for (i = 0; i < HROWS; ++i) {
            for (j = 0; j < MSGLEN; ++j) {
                deltaQ[j][i] = 1;
                for (i2 = 0; i2 < HROWS; ++i2) {
                    for (j2 = 0; j2 < MSGLEN; ++j2) {
                        if (j != j2 && h_matrix[j][i] == h_matrix[j2][i2]) {
                            deltaQ[j][i] *= deltaP[j2][i2];
                        }
                    }
                }

            }
        }

        for (i = 0; i < HROWS; ++i) {
            for (j = 0; j < MSGLEN; ++j) {
                Q0[j][i] = (float)(0.5 * (1.0+deltaQ[j][i]));
                Q1[j][i] = (float)(0.5 * (1.0-deltaQ[j][i]));
            }
        }

        for (j = 0; j < MSGLEN; ++j) {
            for (i = 0; i < HROWS; ++i) {
                P0[j][i] = p0[j];
                P1[j][i] = p1[j];
                for (i2 = 0; i2 < HROWS; ++i2) {
                    if (i != i2) {
                        P0[j][i] *= Q0[j][i2];
                        P1[j][i] *= Q1[j][i2];
                    }
                }
            }
        }

        //scale
        for ( j = 0; j < MSGLEN; ++j) {
            for ( i = 0; i < HROWS; ++i) {
            	float sum = (P0[j][i] + P1[j][i]);
                if (sum != 0) {
                    P0[j][i] /= sum;
                    P1[j][i] /= sum;
                }
            }
        }
        //new P's
        for (j = 0; j < MSGLEN; ++j) {
            for (i = 0; i < HROWS; ++i) {
                for ( i2 = 0; i2 < HROWS; ++i2) {
                    p0[j] *= Q0[j][i2];
                    p1[j] *= Q1[j][i2];
                }

            }
        }

        //scale
        for ( j = 0; j < MSGLEN; ++j) {
        	float sum = (p0[j] + p1[j]);
            if (sum != 0) {
                p0[j] /= sum;
                p1[j] /= sum;
            }
        }

        // Check ending conditions
        calculateBits(p0, p1, output, &change);

        if (change == 0){
            same += 1;
            if (same == end){
                *iterations = iter+1;
                break;
            }
        } else {
            same = 0;
        }
    }
}

void hw_u32_array_to_float(unsigned int input[MSGLEN], float *result) {
    for (int i = 0; i < MSGLEN; ++i){
        result[i] = hw_u32_to_float(input[i]);
    }
}

float hw_u32_to_float(unsigned int val) {
    // thanks, Mehdi
    union {
        float val_float;
        unsigned char bytes[4];
    } data;
    data.bytes[3] = (val >> (8 * 3)) & 0xff;
    data.bytes[2] = (val >> (8 * 2)) & 0xff;
    data.bytes[1] = (val >> (8 * 1)) & 0xff;
    data.bytes[0] = (val >> (8 * 0)) & 0xff;
    return data.val_float;
}
