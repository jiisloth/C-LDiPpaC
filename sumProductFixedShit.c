#include <math.h>
#include "hmatrix.h"
#include <stdio.h>

#define SHIFT 20
#define DIVSHIFT 10 // shift/2
#define ONE 1048576 //1 << SHIFT

void hw_u32_array_to_float(unsigned int input[BLOCKSIZE], float *result);
void hw_u32_to_float(unsigned int val, float *res);


void calculateBits(int prob0[BLOCKSIZE], int prob1[BLOCKSIZE], unsigned int *bits, short *c){
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

void sumProductDecoding(unsigned int max_iter, unsigned int endsame, int message[BLOCKSIZE], unsigned int decodedmsg[BLOCKSIZE], unsigned int *iterations) {

    #pragma HLS INTERFACE s_axilite port=return bundle=CTRLS
    #pragma HLS INTERFACE s_axilite port=max_iter bundle=CTRLS
	#pragma HLS INTERFACE s_axilite port=endsame bundle=CTRLS
    #pragma HLS INTERFACE s_axilite port=message bundle=CTRLS
    #pragma HLS INTERFACE s_axilite port=decodedmsg bundle=CTRLS
    #pragma HLS INTERFACE s_axilite port=iterations bundle=CTRLS

	int prob0[BLOCKSIZE];
	int prob1[BLOCKSIZE];

	float messagef[BLOCKSIZE];
	unsigned int messaged[BLOCKSIZE];

    int deltaP[BLOCKSIZE][PARITYCHECKS];


    int Q0[BLOCKSIZE][PARITYCHECKS];
    int Q1[BLOCKSIZE][PARITYCHECKS];

    short i;
    short i2;
    short j;
    short j2;
    short change = 0;

    // init as max+1 so you know how it finished.
    short itercount = max_iter+1;

    // convert message:
    //hw_u32_array_to_float(message, messagef);




    // init arrays
    for (i = 0; i < BLOCKSIZE; ++i) {
    	messagef[i] = (float)message[i]/65536.0;
    	messaged[i] = 0;
    	float prob = (1 / (1  + expf((2 * messagef[i] / SIGMA))));
    	prob1[i] = (int) (prob * ONE);
        prob0[i] = (int)(ONE - prob1[i]);
    }
    // get original message
    calculateBits(prob0, prob1, messaged, &change);

    for (i = 0; i < BLOCKSIZE; ++i) {
        //printf("%f %d %d %d\n", messagef[i], messaged[i], prob1[i], prob0[i]);
    }

    printf("DeltaP:\n");
    // init P matrix's
    for ( j = 0; j < BLOCKSIZE; ++j) {
    	for (i = 0; i < PARITYCHECKS; ++i) {
            // DeltaP - First time out of the main loop..
            deltaP[j][i] = prob0[j] - prob1[j];
            //printf("%d ", deltaP[j][i]);
        }
    }
    printf("\n");
    unsigned int same = 0;

    // actual algorithm
    for (unsigned int iter = 0; iter < max_iter; ++iter) {
        // Q:s
        printf("DeltaQ:\n");
    	for (j = 0; j < BLOCKSIZE; ++j){
    		for (i = 0; i < PARITYCHECKS; ++i) {
    			int deltaQ = ONE;
    			int row = h_matrix[j][i];
				for (j2 = 0; j2 < BLOCKSIZE; ++j2) {
					if (j2 != j){
						for (i2 = 0; i2 < PARITYCHECKS; ++i2) {
							if (h_matrix[j2][i2] == row) {
								deltaQ = ((deltaQ >> 2) * (deltaP[j2][i2] >> 2)) >> (SHIFT-4);
							}
						}
					}
				}
                printf("%d ", deltaQ);
				Q0[j][i] = (1 << SHIFT) + deltaQ;
				Q1[j][i] = (1 << SHIFT) - deltaQ;
			}
		}
        printf("\n");

        printf("pP:\n");
		// getting probabilities
        for (j = 0; j < BLOCKSIZE; ++j) {
            for (i = 0; i < PARITYCHECKS; ++i) {
                int prob0matrix = prob0[j];
                int prob1matrix = prob1[j];
                for (i2 = 0; i2 < PARITYCHECKS; ++i2) {
                    if (i != i2) {
                        prob0matrix = (prob0matrix * Q0[j][i2]) >> SHIFT;
                        prob1matrix = (prob1matrix * Q1[j][i2]) >> SHIFT;
                    }
                }
                printf("%d %d\n", prob0matrix, prob1matrix);
                //scale
                if (prob0matrix + prob1matrix != 0){
                	int sum = (int)ONE/(prob0matrix + prob1matrix);
					if (sum != 0) {
						prob0matrix = (int)(prob0matrix * sum);
						prob1matrix = (int)(prob1matrix * sum);
					}
                }
                printf("  %d %d\n", prob0matrix, prob1matrix);
                //delta P
                deltaP[j][i] = prob0matrix - prob1matrix;
            }
        }
        //new P's
        for (j = 0; j < BLOCKSIZE; ++j) {
            for (i = 0; i < PARITYCHECKS; ++i) {
                for ( i2 = 0; i2 < PARITYCHECKS; ++i2) {
                    prob0[j] = prob0[j] * Q0[j][i2] >> SHIFT;
                    prob1[j] = prob1[j] * Q1[j][i2] >> SHIFT;
                }

            }
            //scale
            int sum = prob0[j] + prob1[j] >>DIVSHIFT;
            if (sum != 0) {
                prob0[j] = (int)(prob0[j] / sum) >> DIVSHIFT;
                prob1[j] = (int)(prob1[j] / sum) >> DIVSHIFT;
            }
        }
        // Check ending conditions
        calculateBits(prob0, prob1, messaged, &change);

        if (change == 0){
            if (same >= endsame-1){
                itercount = iter+1;
                break;
            }
            same = same + 1;
        } else {
            same = 0;
        }
    }
    // set endvalues.
    for (i = 0; i < BLOCKSIZE; ++i) {
    	decodedmsg[i] = messaged[i];
    }
    *iterations = itercount;
    return;
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
