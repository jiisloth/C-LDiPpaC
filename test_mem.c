#include <stdio.h>
#include "testcases.h"
#include "Testbench.h"

#define printfv(v, ...) do{if(verbose_lvl>=(v)){printf(__VA_ARGS__);}}while(0);

#define MAXITER 10
#define ENDSAME 3
#define SUMPRODUCT 0
#define HARDDECISION 1


#define TESTCASES 4 // 0 = all
#define CASEOFFSET 0

int verbose_lvl = 2; //-1: Nothing, 0: Final results, 1: Individual test results, 2: Individual test outputs
int runsingle = 0; //0 for all, others for case #

int algorithm = SUMPRODUCT; // Remember to comment/uncomment from run_case

int main() {
    if (algorithm == SUMPRODUCT) {
        printfv(0, "Running Sum Product Decoding..\n");
    } else if (algorithm == HARDDECISION) {
        printfv(0, "Running Hard Decision Decoding..\n");
    } else {
        printf("Wrong algorithm configuration!\n");
        return 1;
    }
    int cases_to_test = TESTCASES;
    if (cases_to_test > TESTCASECOUNT || cases_to_test < 1) {
        cases_to_test = TESTCASECOUNT;
        if (TESTCASES > 0) {
            printf("No enough test cases.. Running all %d cases\n", TESTCASECOUNT);
        }
    }

    int errors = 0;
    if (runsingle <= 0) {
        for (int t = 0; t < cases_to_test; t++) {
            int tc = t + CASEOFFSET;
            tc = tc % TESTCASECOUNT;

            printfv(1, "\nTest #%d, Test case: %d, ", t, tc + 1);

            errors += run_case(*cases[tc]);

        }
        if (errors == cases_to_test) {
            printfv(0, "\nAll %d tests FAILED!\n", cases_to_test);
        } else if (errors == 0) {
            printfv(0, "\nAll %d tests PASSED!\n", cases_to_test);
        } else {
            printfv(0, "\nTests PASSED: %3d/%d\n", cases_to_test - errors, cases_to_test);
            printfv(0, "Tests FAILED: %3d/%d\n", errors, cases_to_test);
        }
    } else if (runsingle <= TESTCASECOUNT) {
        printfv(1, "\nTest case: %d, ", runsingle);
        errors += run_case(*cases[runsingle - 1]);
        if (errors == 0) {
            printfv(0, "\nTest passed!\n");
        } else {
            printfv(0, "\nTest failed!\n");
        }
    }
    return 0;
}

int run_case(Testcase custard) {
    // custard = case but cant use case as variable...

    //get binary message and initial errors..
    int errorsinmessage = 0;

    int i;
    for (i = 0; i < BLOCKSIZE; i++) {
        if (custard.receivedhard[i] != custard.GoldRef[i]) {
            errorsinmessage += 1;
        }
    }
    printfv(1, "Flipped bits in the message: %d.\n", errorsinmessage);

    int output[BLOCKSIZE] = {0};
    int itercount = 0;
    int *iterations = &itercount;

    int endsame = ENDSAME;


    printfv(1, "Running...\n");
    if (algorithm == SUMPRODUCT) {
        sumProductDecoding(MAXITER, endsame, custard.receivedfloat, output, iterations);
    } else if (algorithm == HARDDECISION) {
        //hardDecisionDecoding(MAXITER, custard.receivedhard, output, iterations);
    }

    if (itercount == 0) {
        printfv(1, "- Max iterations reached. (%d)\n", MAXITER);
    } else {
        printfv(1, "- Finished on iteration # %d.\n", itercount);
    }
    int crcresult = -1;
    if (crc_bits[0] == 1) {
        crcresult = test_crc(output);
        if (crcresult == 1) {
            printfv(1, "- CRC failed\n");
        } else {
            printfv(1, "- CRC succesfull\n");
        }
    }

    int errors = 0;
    int errorline[BLOCKSIZE] = {0};

    printfv(1, "- Comparing to GoldRef: ");

    printfv(2, "\n");

    int goldrefresult = 0;
    for (i = 0; i < BLOCKSIZE; i++) {
        if (output[i] != custard.GoldRef[i]) {
            errors += 1;
            errorline[i] = 1;
            goldrefresult = 1;
        }
        printfv(2, "%d", output[i]);
    }
    if (errors > 0) {
        if (verbose_lvl >= 2) {
            printfv(2, "\n");
            for (i = 0; i < BLOCKSIZE; i++) {
                if (errorline[i] == 1) {
                    printfv(2, "*");
                } else {
                    printfv(2, " ");
                }
            }
            printfv(2, "\n");
            printfv(2, "  - Decoding failed with %d errors!\n", errors);
        } else {
            printfv(1, "Fail! (%d errors.)\n", errors);
        }
    } else {
        if (verbose_lvl >= 2) {
            printfv(2, "  - Decoding succesfull!\n");
        } else {
            printfv(1, "OK!\n");
        }
    }
    if (crcresult == -1) {
        return goldrefresult;
    } else {
        if (crcresult != goldrefresult) {
            if (crcresult == 1) {
                printfv(1, "  - CRC gave false negative!\n");
            } else {
                printfv(1, "  - CRC gave false positive!\n");
            }
        }
        return crcresult;
    }
}

int test_crc(int output[BLOCKSIZE]) {
    int crc_padded_len = MSGLEN + CRCLEN - 1;
    int test[BLOCKSIZE] = {0};
    int c = 0;
    int i;

    for (i = 0; i < crc_padded_len; ++i) {
        test[i] = output[i];
    }

    while (c < MSGLEN) {
        for (i = c; i < MSGLEN; ++i) {
            if (test[i] == 1) {
                c = i;
                break;
            }
            c = -1;
        }
        if (c == -1) {
            break;
        }
        for (i = 0; i < CRCLEN; ++i) {
            test[c + i] ^= crc_bits[i];
        }
    }
    for (i = 0; i < crc_padded_len; ++i) {
        if (test[i] == 1) {
            return 1;
        }
    }
    return 0;
}


#include <math.h>
#include "hmatrix.h"


#define HARDMAXITER 16


void calculateBits(float prob0[BLOCKSIZE], float prob1[BLOCKSIZE], volatile int *bits, short *c) {
    int change = 0;
    for (int i = 0; i < BLOCKSIZE; ++i) {
        if (prob0[i] > prob1[i]) {
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
    init:
    for (j = 0; j < BLOCKSIZE; ++j) {
        float prob = (float) (1 / (1 + expf((2 * message[j]) / (float) SIGMA)));
        prob1[j] = prob;
        prob0[j] = (float) (1 - prob);

        //init msg.
        if (prob < 0.5) {
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

        if (iter < maxi) {
            // Q:s
            Qs:
            for (j = 0; j < BLOCKSIZE; ++j) {
                for (i = 0; i < PARITYCHECKS; ++i) {
                    float deltaQ = 1;
#pragma HLS PIPELINE II=7
                    int row = h_matrix[j][i];
                    get_deltaPs:
                    for (j2 = 0; j2 < BLOCKSIZE; ++j2) {
                        for (i2 = 0; i2 < PARITYCHECKS; ++i2) {
                            if (h_matrix[j2][i2] == row && j != j2) {
                                deltaQ *= deltaP[j2][i2];
                            }
                        }
                    }
                    Q0[j][i] = (float) (1.0 + deltaQ);
                    Q1[j][i] = (float) (1.0 - deltaQ);
                }
            }
            // getting probabilities
            probabilities:
            for (j = 0; j < BLOCKSIZE; ++j) {
                for (i = 0; i < PARITYCHECKS; ++i) {
                    float prob0matrix = prob0[j];
                    float prob1matrix = prob1[j];
                    for (i2 = 0; i2 < PARITYCHECKS; ++i2) {
                        if (i != i2) {
                            prob0matrix = (float) (prob0matrix * Q0[j][i2]);
                            prob1matrix = (float) (prob1matrix * Q1[j][i2]);
                        }
                    }
                    //scale
                    float sum = (float) (prob0matrix + prob1matrix);
                    if (sum != 0) {
                        float factor = (float) 1 / sum;
                        prob0matrix = (float) (prob0matrix * factor);
                        prob1matrix = (float) (prob1matrix * factor);
                    }
                    //delta P
                    deltaP[j][i] = prob0matrix - prob1matrix;
                }
            }
            //new P's
            new_ps:
            for (j = 0; j < BLOCKSIZE; ++j) {
                float p0 = prob0[j];
                float p1 = prob1[j];
                for (i = 0; i < PARITYCHECKS; ++i) {
                    for (i2 = 0; i2 < PARITYCHECKS; ++i2) {
                        p0 = (float) (p0 * Q0[j][i2]);
                        p1 = (float) (p1 * Q1[j][i2]);
                    }

                }
                //scale
                float sum = (float) (p0 + p1);
                if (sum != 0) {
                    p0 = (float) (p0 / sum);
                    p1 = (float) (p1 / sum);
                }
                prob0[j] = p0;
                prob1[j] = p1;

            }
            // Check ending conditions
            //calculateBits(prob0, prob1, messaged, &change);

            int change = 0;
            check_result:
            for (i = 0; i < BLOCKSIZE; ++i) {
                if (prob0[i] > prob1[i]) {
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


            if (change == 0) {
                if (same >= endonsames - 1) {
                    itercount = iter + 1;
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
    endloop:
    for (j = 0; j < BLOCKSIZE; ++j) {
        decodedmsg[j] = messaged[j];
    }
}

