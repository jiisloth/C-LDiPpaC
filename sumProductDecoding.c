#include <math.h>
#include "hmatrix.h"
#include <stdio.h>


#define HARDMAXITER 16


void calculateBits(float prob0[BLOCKSIZE], float prob1[BLOCKSIZE], volatile int *bits, short *c){
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

    float prob0[BLOCKSIZE] = {0};
    float prob1[BLOCKSIZE] = {0};

    //float messagef[BLOCKSIZE];
    int messaged[BLOCKSIZE] = {0};

    float deltaP[BLOCKSIZE][PARITYCHECKS] = {0};
#pragma HLS ARRAY_PARTITION dim=2 type=complete variable=deltaP

    float Q0[BLOCKSIZE][PARITYCHECKS] = {0};
    float Q1[BLOCKSIZE][PARITYCHECKS] = {0};

    short i = 0;
    short i2 = 0;
    short j = 0;
    short j2 = 0;
    //short change = 0;


    int maxi = max_iter;
    int itercount = 0;



    // init arrays
    init:for (j = 0; j < BLOCKSIZE; ++j) {
        float prob = (float)(1 / (1 + expf((2 * message[j]) / (float)SIGMA)));
        prob1[j] = prob;
        prob0[j] = (float)(1 - prob);

        //init msg.
        if (prob < 0.5){
            messaged[j] = 0;
        } else {
            messaged[j] = 1;
        }

        init_deltaP:for (i = 0; i < PARITYCHECKS; ++i) {
            // DeltaP - First time out of the main loop..
            deltaP[j][i] = 1 - prob - prob;
        }
    }

    // get original message
    //calculateBits(prob0, prob1, messaged, &change);



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
    main_loop:for (int iter = 0; iter < HARDMAXITER; ++iter) {

    	if (iter < maxi){
			// Q:s
			Qs:for (j = 0; j < BLOCKSIZE; ++j){
				for (i = 0; i < PARITYCHECKS; ++i) {
					float  deltaQ = 1;
#pragma HLS PIPELINE II=7
					int row = h_matrix[j][i];
					get_deltaPs:for (j2 = 0; j2 < BLOCKSIZE; ++j2) {
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
			probabilities:for (j = 0; j < BLOCKSIZE; ++j) {
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
						float factor = (float)1 / sum;
						prob0matrix = (float)(prob0matrix * factor);
						prob1matrix = (float)(prob1matrix * factor);
					}
					//delta P
					deltaP[j][i] = prob0matrix - prob1matrix;
				}
			}
			//new P's
			new_ps:for (j = 0; j < BLOCKSIZE; ++j) {
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
			//calculateBits(prob0, prob1, messaged, &change);

		    int change = 0;
		    check_result:for (int i = 0; i < BLOCKSIZE; ++i) {
		        if (prob0[i] > prob1[i]){
		            if (messaged[i] == 1) {
		                change += 1;
		            }
		            messaged[i] = 0;
		        } else {
		            if (messaged[i] == 0) {
		                change += 1;
		            }
		            messaged[i] = 1;
		        }
		    }


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
    endloop:for (j = 0; j < BLOCKSIZE; ++j) {
        decodedmsg[j] = messaged[j];
    }
}

