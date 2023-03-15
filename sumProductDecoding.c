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

#pragma HLS INTERFACE s_axilite port=iterations //bundle=CTRLS
#pragma HLS INTERFACE s_axilite port=decodedmsg //bundle=CTRLS
#pragma HLS INTERFACE s_axilite port=message //bundle=CTRLS
#pragma HLS INTERFACE s_axilite port=endsame //bundle=CTRLS
#pragma HLS INTERFACE s_axilite port=max_iter //bundle=CTRLS
#pragma HLS INTERFACE s_axilite port=return //bundle=CTRLS

    float prob0[BLOCKSIZE];
    float prob1[BLOCKSIZE];

    //float messagef[BLOCKSIZE];
    volatile int messaged[BLOCKSIZE];

    float deltaP[BLOCKSIZE][PARITYCHECKS];
#pragma HLS ARRAY_PARTITION dim=2 type=complete variable=deltaP


    int i = 0;
    int i2 = 0;
    int j = 0;
    int j2 = 0;
    //short change = 0;


    int maxi = max_iter;
    int itercount = 0;



    // init arrays
    init:
	for (j = 0; j < BLOCKSIZE; ++j) {
        float prob = (float)(1 / (1 + expf((2 * message[j]) / (float)SIGMA)));
        prob1[j] = prob;
        prob0[j] = (float)(1 - prob);

        //init msg.
        if (prob < 0.5){
            messaged[j] = 0;
        } else {
            messaged[j] = 1;
        }

        init_deltaP:
		for (i = 0; i < PARITYCHECKS; ++i) {
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
    main_loop:
	for (int iter = 0; iter < HARDMAXITER; ++iter) {
    	if (iter < maxi){
    	    float Q0[BLOCKSIZE][PARITYCHECKS];
    	    float Q1[BLOCKSIZE][PARITYCHECKS];
			// Q:s
    	    float tQ[CNODES];
#pragma HLS ARRAY_PARTITION type=complete variable=tQ
    	    initQ:
    	    for (int c = 0; c < CNODES; ++c){
    	    	tQ[c] = 1;
    	    }
			Cj:
    	    for (j = 0; j < BLOCKSIZE; ++j){
				Ci:
				for (i = 0; i < PARITYCHECKS; ++i) {
#pragma HLS PIPELINE II=5
					tQ[h_matrix[j][i]] *= deltaP[j][i];
				}
			}

			Qj:
			for (j = 0; j < BLOCKSIZE; ++j){
				Qi:
				for (i = 0; i < PARITYCHECKS; ++i) {
					float  deltaQ = (float)(tQ[h_matrix[j][i]]/deltaP[j][i]);
					Q0[j][i] = (float)(1.0+deltaQ);
					Q1[j][i] = (float)(1.0-deltaQ);
				}
			}
			// getting probabilities
			dPj:
			for (j = 0; j < BLOCKSIZE; ++j) {
				dPi:
				for (i = 0; i < PARITYCHECKS; ++i) {
					float prob0matrix = prob0[j];
					float prob1matrix = prob1[j];
					pmi:
					for (i2 = 0; i2 < PARITYCHECKS; ++i2) {
						if (i != i2) {
							prob0matrix = (float)(prob0matrix * Q0[j][i2]);
							prob1matrix = (float)(prob1matrix * Q1[j][i2]);
						}
					}
					//scale
					float sum = (float)(prob0matrix + prob1matrix);
					if (sum != 0) {
						float factor = (float)(1 / sum);
						prob0matrix = (float)(prob0matrix * factor);
						prob1matrix = (float)(prob1matrix * factor);
					}
					//delta P
					deltaP[j][i] = prob0matrix - prob1matrix;
				}
			}
			//new P's
			pj:
			for (j = 0; j < BLOCKSIZE; ++j) {
				float p0 = prob0[j];
				float p1 = prob1[j];
				pi:
				for (i = 0; i < PARITYCHECKS; ++i) {
					pif:
					for ( i2 = 0; i2 < PARITYCHECKS; ++i2) {
						p0 = (float)(p0 * Q0[j][i2]);
						p1 = (float)(p1 * Q1[j][i2]);
					}

				}
				//scale
				float sum = (float)(p0 + p1);
				if (sum != 0) {
					float factor = (float)(1 / sum);
					p0 = (float)(p0 * factor);
					p1 = (float)(p1 * factor);
				}
				prob0[j] = p0;
				prob1[j] = p1;

			}
			// Check ending conditions
			//calculateBits(prob0, prob1, messaged, &change);

		    int change = 0;
		    check_result:
			for (i = 0; i < BLOCKSIZE; ++i) {
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
	if (itercount == 0 && maxi > HARDMAXITER){
		// Tells that Hard max was reached.
		itercount == -HARDMAXITER;
	}
    *iterations = itercount;
    endloop:
	for (j = 0; j < BLOCKSIZE; ++j) {
        decodedmsg[j] = messaged[j];
    }
}

