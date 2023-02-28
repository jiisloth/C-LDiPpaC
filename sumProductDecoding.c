#include <math.h>

#define SIGMA 1

#define MSGLEN 50
#define HROWS 20

#define true 1
#define false 0



int calculateBits(float p0[], float p1[], int *bits){
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
    return change;
}


void sumProductDecoding(int max_iter, int end, int h_matrix[HROWS][MSGLEN], float message[MSGLEN], int output[MSGLEN], int *iterations) {
	#pragma HLS INTERFACE s_axilite port=return bundle=CRTLS
	#pragma HLS INTERFACE s_axilite port=output bundle=CRTLS
    #pragma HLS INTERFACE s_axilite port=message bundle=CRTLS
    #pragma HLS INTERFACE s_axilite port=h_matrix bundle=CRTLS
    #pragma HLS INTERFACE s_axilite port=iterations bundle=CRTLS

	float p0[MSGLEN];
	float p1[MSGLEN];
	float P0[HROWS][MSGLEN];
	float P1[HROWS][MSGLEN];
    int i;
    int i2;
    int j;
    int j2;

    //returns -1 if max_iter is reached
    *iterations = -1;
    int same = 0;

    // init arrays
    for (i = 0; i < MSGLEN; ++i) {
        p1[i] = 1 / (1 + exp((2 * message[i]) / SIGMA));
        p0[i] = 1 - p1[i];
    }
    // get original message
    calculateBits(p0, p1, output);
    // init P matrix's
    for (i = 0; i < HROWS; ++i) {
        for ( j = 0; j < MSGLEN; ++j) {
            P0[i][j] = h_matrix[i][j] * p0[j];
            P1[i][j] = h_matrix[i][j] * p1[j];
        }
    }
    // actual algorithm
    for (int iter = 0; iter < max_iter; ++iter) {
        // DeltaP
        float deltaP[HROWS][MSGLEN];
        for ( i = 0; i < HROWS; ++i) {
            for ( j = 0; j < MSGLEN; ++j) {
                deltaP[i][j] = P0[i][j] - P1[i][j];
            }
        }
        // Q:s
        float deltaQ[HROWS][MSGLEN];
        float Q0[HROWS][MSGLEN];
        float Q1[HROWS][MSGLEN];
        for (i = 0; i < HROWS; ++i) {
            for ( j = 0; j < MSGLEN; ++j) {
            	float product = 1;
                if (deltaP[i][j] != 0){
                    for ( j2 = 0; j2 < MSGLEN; ++j2) {
                        if (j != j2 && deltaP[i][j2] != 0){
                            product = product * deltaP[i][j2];
                        }
                        deltaQ[i][j] = product;
                        Q0[i][j] = (0.5 * (1+product));
                        Q1[i][j] = (0.5 * (1-product));
                    }
                } else {
                    deltaQ[i][j] = 0;
                    Q0[i][j] = 0;
                    Q1[i][j] = 0;
                }
            }
        }
        //new P's
        float np0[MSGLEN] = {0};
        float np1[MSGLEN] = {0};
        for ( i = 0; i < HROWS; ++i) {
            for ( j = 0; j < MSGLEN; ++j) {
            	float product = 1;
            	float product2 = 1;
                if (Q0[i][j] != 0) {
                    for ( i2 = 0; i2 < HROWS; ++i2) {
                        if (Q0[i2][j] != 0) {
                            product2 = product2 * Q0[i2][j];
                            if (i != i2) {
                                product = product * Q0[i2][j];
                            }
                        }
                    }
                    np0[j] = p0[j] * product2;
                    P0[i][j] = p0[j] * product;
                } else {
                    P0[i][j] = 0;
                }
                product = 1;
                product2 = 1;
                if (Q1[i][j] != 0) {
                    for ( i2 = 0; i2 < HROWS; ++i2) {
                        if (Q1[i2][j] != 0) {
                            product2 = product2 * Q1[i2][j];
                            if (i != i2) {
                                product = product * Q1[i2][j];
                            }
                        }
                    }
                    np1[j] = p1[j] * product2;
                    P1[i][j] = p1[j] * product;
                } else {
                    P1[i][j] = 0;
                }
            }
        }
        for ( i = 0; i < MSGLEN; ++i) {
            p0[i] = np0[i];
            p1[i] = np1[i];
        }
        //scale
        for ( j = 0; j < MSGLEN; ++j) {
            for ( i = 0; i < HROWS; ++i) {
                if (P0[i][j] + P1[i][j] != 0) {
                	float factor = 1 / (P0[i][j] + P1[i][j]);
                    P0[i][j] = P0[i][j] * factor;
                    P1[i][j] = P1[i][j] * factor;
                }
            }
            if (p0[j] + p1[j] != 0) {
            	float factor = 1 / (p0[j] + p1[j]);
                p0[j] = p0[j]*factor;
                p1[j] = p1[j]*factor;
            }

        }
        if (calculateBits(p0, p1, output) == 0){
            same += 1;
            if (same == end){
                *iterations = iter+1;
            }
        } else {
            same = 0;
        }
    }
}
