#include <math.h>
#include "hmatrix.h"
#include <stdio.h>


#define HARDMAXITER 16

void hw_u32_array_to_float(unsigned int input[BLOCKSIZE], float *result);
void hw_u32_to_float(unsigned int val, float *res);


void calculateBits(float prob0[BLOCKSIZE], float prob1[BLOCKSIZE], int *bits, short *c){
    int change = 0;
    for (int i = 0; i < BLOCKSIZE; ++i) {
        if (prob0[i] > prob1[i]){
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

void sumProductDecoding(int max_iter, int endsame, float message[BLOCKSIZE], int decodedmsg[BLOCKSIZE], int *iterations) {

#pragma HLS INTERFACE s_axilite port=iterations bundle=CTRLS
#pragma HLS INTERFACE s_axilite port=decodedmsg bundle=CTRLS
#pragma HLS INTERFACE s_axilite port=message bundle=CTRLS
#pragma HLS INTERFACE s_axilite port=endsame bundle=CTRLS
#pragma HLS INTERFACE s_axilite port=max_iter bundle=CTRLS
#pragma HLS INTERFACE s_axilite port=return bundle=CTRLS

    float prob0[BLOCKSIZE];
    float prob1[BLOCKSIZE];

    //float messagef[BLOCKSIZE];
    int messaged[BLOCKSIZE];

    float deltaP[BLOCKSIZE][PARITYCHECKS];

    float Q0[BLOCKSIZE][PARITYCHECKS];
    float Q1[BLOCKSIZE][PARITYCHECKS];

    short i;
    short i2;
    short j;
    short j2;
    short change = 0;


    int maxi = max_iter;
    int itercount = 0;

    // convert message:
    //hw_u32_array_to_float(message, messagef);


    // init arrays
    for (j = 0; j < BLOCKSIZE; ++j) {
        messaged[j] = 0;
        float prob = (float)(1 / (1 + expf((2 * message[j]) / (float)SIGMA)));
        prob1[j] = prob;
        prob0[j] = (float)(1 - prob);
        for (i = 0; i < PARITYCHECKS; ++i) {
            // DeltaP - First time out of the main loop..
            deltaP[j][i] = 1 - prob - prob;
        }
    }

    // get original message
    calculateBits(prob0, prob1, messaged, &change);

    // init P matrix's
    /*
    for ( j = 0; j < BLOCKSIZE; ++j) {
        for (i = 0; i < PARITYCHECKS; ++i) {
            // DeltaP - First time out of the main loop..
            deltaP[j][i] = (float)(prob0[j] - prob1[j]);
        }
    }*/
    int same = 0;
    int endonsames = endsame;

    // actual algorithm
    for (int iter = 0; iter < HARDMAXITER; ++iter) {
    	if (iter < maxi){
			// Q:s
			for (j = 0; j < BLOCKSIZE; ++j){
				for (i = 0; i < PARITYCHECKS; ++i) {
					float  deltaQ = 1;
					int row = h_matrix[j][i];
					for (j2 = 0; j2 < BLOCKSIZE; ++j2) {
						for (i2 = 0; i2 < PARITYCHECKS; ++i2) {
							if (h_matrix[j2][i2] == row && j != j2) {
								deltaQ *= deltaP[j2][i2];
							}
						}
					}
					Q0[j][i] = (float)(1.0+deltaQ);
					Q1[j][i] = (float)(1.0-deltaQ);
				}
			}
			// getting probabilities
			for (j = 0; j < BLOCKSIZE; ++j) {
				for (i = 0; i < PARITYCHECKS; ++i) {
					float prob0matrix = prob0[j];
					float prob1matrix = prob1[j];
					for (i2 = 0; i2 < PARITYCHECKS; ++i2) {
						if (i != i2) {
							prob0matrix = (float)(prob0matrix * Q0[j][i2]);
							prob1matrix = (float)(prob1matrix * Q1[j][i2]);
						}
					}
					//scale
					float sum = (float)(prob0matrix + prob1matrix);
					if (sum != 0) {
						prob0matrix = (float)(prob0matrix / sum);
						prob1matrix = (float)(prob1matrix / sum);
					}
					//delta P
					deltaP[j][i] = prob0matrix - prob1matrix;
				}
			}
			//new P's
			for (j = 0; j < BLOCKSIZE; ++j) {
				float p0 = prob0[j];
				float p1 = prob1[j];
				for (i = 0; i < PARITYCHECKS; ++i) {
					for ( i2 = 0; i2 < PARITYCHECKS; ++i2) {
						p0 = (float)(p0 * Q0[j][i2]);
						p1 = (float)(p1 * Q1[j][i2]);
					}

				}
				//scale
				float sum = (float)(p0 + p1);
				if (sum != 0) {
					p0 = (float)(p0 / sum);
					p1 = (float)(p1 / sum);
				}
				prob0[j] = p0;
				prob1[j] = p1;

			}
			// Check ending conditions
			calculateBits(prob0, prob1, messaged, &change);

			if (change == 0){
				if (same >= endonsames-1){
					itercount = iter+1;
					break;
				}
				same = same + 1;
			} else {
				same = 0;
			}
    	}
    }
    // set endvalues.
    *iterations = itercount;
    for (j = 0; j < BLOCKSIZE; ++j) {
        decodedmsg[j] = messaged[j];
    }
}


void hw_u32_array_to_float(unsigned int input[BLOCKSIZE], float *result) {
    for (int i = 0; i < BLOCKSIZE; ++i){
        float res = 0;
        hw_u32_to_float(input[i], &res);
        result[i] = res;
    }
}


void hw_u32_to_float(unsigned int val, float *res) {
    union {
        float val_float;
        unsigned char bytes[4];
    } data;
    data.bytes[3] = (val >> (8 * 3)) & 0xff;
    data.bytes[2] = (val >> (8 * 2)) & 0xff;
    data.bytes[1] = (val >> (8 * 1)) & 0xff;
    data.bytes[0] = (val >> (8 * 0)) & 0xff;
    *res = data.val_float;
}
