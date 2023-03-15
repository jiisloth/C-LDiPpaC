#include <stdio.h>
#include "testcases.h"
#include "Testbench.h"
#define printfv(v,...) do{if(verbose_lvl>=(v)){printf(__VA_ARGS__);}}while(0);


#include <math.h>
#include "hmatrix.h"


void OLDsumProductDecoding(unsigned int max_iter, unsigned int end, unsigned int intmsg[BLOCKSIZE], unsigned int output[BLOCKSIZE], unsigned int *iterations);
void oldhardDecisionDecoding(unsigned int max_iter, unsigned int message[BLOCKSIZE], unsigned int output[BLOCKSIZE], unsigned int *iterations);
void hw_u32_array_to_float(unsigned int input[BLOCKSIZE], float *result);
float hw_u32_to_float(unsigned int val);


int compare(Testcase custard);

#define MAXITER 1
#define ENDSAME 2
#define SUMPRODUCT 0
#define HARDECISION 1

short int h_matrix_all[0][0];

int verbose_lvl = 0; //-1: Nothing, 0: Final results, 1: Individual test results, 2: Individual test outputs
int runsingle = 0; //-1 for all

int algorithm = SUMPRODUCT;

int main(){

    if (algorithm == SUMPRODUCT) {
        printfv(0,"Running Sum Product Decoding..\n");
    } else if (algorithm == HARDECISION) {
        printfv(0,"Running Hard Decision Decoding..\n");
    }
    int errors = 0;
    if (runsingle < 0) {
        for (int t = 0; t < TESTCASECOUNT; t++){

            printfv(1,"\nTestcase: %d, ", t);

            errors += run_case(*cases[t]);

        }
        if (errors == TESTCASECOUNT){
            printfv(0,"\nAll %d tests FAILED!\n", TESTCASECOUNT);
        } else if (errors == 0){
            printfv(0,"\nAll %d tests PASSED!\n", TESTCASECOUNT);
        } else {
            printfv(0,"\nTests PASSED: %3d/%d\n", TESTCASECOUNT-errors, TESTCASECOUNT);
            printfv(0,"\nTests FAILED: %3d/%d\n", errors, TESTCASECOUNT);
        }
    } else if (runsingle < TESTCASECOUNT) {
        errors += run_case(*cases[runsingle]);
        if (errors == 0) {
            printfv(0,"\nTest passed!\n");
        } else {
            printfv(0,"\nTest failed!\n");
        }
    }
    return 0;
}

int run_case(Testcase custard){
    // custard = case but cant use case as variable...

    //get binary message and initial errors..
    int errorsinmessage = 0;

    int i;
    for(i = 0; i < BLOCKSIZE; i++) {
        if (custard.receivedhard[i] != custard.GoldRef[i]) {
            errorsinmessage += 1;
        }
    }
    printfv(1,"Error bits: %d.\n", errorsinmessage);

    unsigned int output[BLOCKSIZE];
    unsigned int iterations = -1;



    printfv(1,"Running...\n");
    if (algorithm == SUMPRODUCT){
        unsigned int message[BLOCKSIZE];
        float_array_to_u32(custard.receivedfloat, message);
        sumProductDecoding(MAXITER, ENDSAME,  message, output, &iterations);
    } else if (algorithm == HARDECISION) {
        hardDecisionDecoding(MAXITER, custard.receivedhard, output, &iterations);
    }
    if (iterations == MAXITER+1){
        printfv(1,"- Max iterations reached. (%d)\n", MAXITER);
    } else {
        printfv(1,"- Finished on iteration # %d.\n", iterations);
    }
    int crcresult = -1;
    if (crc_bits[0] == 1) {
        crcresult = test_crc(output);
        if (crcresult == 1) {
            printfv(1,"CRC failed\n");
        } else {
            printfv(1,"CRC succesfull\n");
        }
    }

    int errors = 0;
    int errorline[BLOCKSIZE] = {0};

    printfv(1,"- Comparing to GoldRef: ");

    printfv(2,"\n");

    int goldrefresult = 0;
    for(i = 0; i < BLOCKSIZE; i++){
        if(output[i] != custard.GoldRef[i]){
            errors += 1;
            errorline[i] = 1;
            goldrefresult = 1;
        }
        printfv(2,"%d", output[i]);
    }
    if (errors > 0){
        if (verbose_lvl >= 2){
            printfv(2,"\n");
            for (i = 0; i < BLOCKSIZE; i++) {
                if (errorline[i] == 1) {
                    printfv(2,"*");
                } else {
                    printfv(2," ");
                }
            }
            printfv(2,"\n");
            printfv(2,"  - Decoding failed with %d errors!\n", errors);
        } else {
            printfv(1,"Fail! (%d errors.)\n", errors);
        }
    } else {
        if (verbose_lvl >= 2){
            printfv(2,"  - Decoding succesfull!\n");
        }  else {
            printfv(1,"OK!\n");
        }
    }
    if (crcresult == -1) {
        return goldrefresult;
    } else {
        if (crcresult != goldrefresult){
            if (crcresult == 1){
                printfv(1,"  - CRC gave false negative!\n");
            } else {
                printfv(1,"  - CRC gave false positive!\n");
            }
        }
        return crcresult;
    }
}


/*
int test_crcOLD(unsigned int output[MSGLEN]){
    int crc_len = sizeof(crc_bits)/ sizeof(int);
    int crc_padded_len = OGMESLEN + crc_len - 1;
    unsigned int test[MSGLEN] = {0};
    int c = 0;
    int i = 0;
    int j = 0;

    for (i = 0; i < crc_padded_len; ++i){
        test[i] = output[i];
    }

    while (c < OGMESLEN){
        for(i = c; i < OGMESLEN; ++i){
            if (test[i] == 1){
                c = i;
                break;
            }
            c = -1;
        }
        if (c == -1){
            break;
        }
        for(i = 0; i < crc_len; ++i){
            test[c+i] ^= crc_bits[i];
        }
    }
    for (i = 0; i < crc_padded_len; ++i){
        if (test[i] == 1){
            return 1;
        }
    }
    return 0;
}

*/
int test_crc(unsigned int output[BLOCKSIZE]){
    int crc_padded_len = MSGLEN + CRCLEN - 1;
    unsigned int test[BLOCKSIZE] = {0};
    int c = 0;
    int i;

    for (i = 0; i < crc_padded_len; ++i){
        test[i] = output[i];
    }

    while (c < MSGLEN){
        for(i = c; i < MSGLEN; ++i){
            if (test[i] == 1){
                c = i;
                break;
            }
            c = -1;
        }
        if (c == -1){
            break;
        }
        for(i = 0; i < CRCLEN; ++i){
            test[c+i] ^= crc_bits[i];
        }
    }
    for (i = 0; i < crc_padded_len; ++i){
        if (test[i] == 1){
            return 1;
        }
    }
    return 0;
}

void float_array_to_u32(float input[BLOCKSIZE], unsigned int *result) {
    for (int i = 0; i < BLOCKSIZE; ++i){
        result[i] = float_to_u32(input[i]);
    }
}

void u32_array_to_float(unsigned int input[BLOCKSIZE], float *result) {
    for (int i = 0; i < BLOCKSIZE; ++i){
        result[i] = u32_to_float(input[i]);
    }
}

unsigned int float_to_u32(float val){
    // thanks, Mehdi
    unsigned int result;
    union float_bytes {
        float v;
        unsigned char bytes[4];
    }data;
    data.v = val;
    result = (data.bytes[3] << 24) + (data.bytes[2] << 16) + (data.bytes[1] << 8) + (data.bytes[0]);
    return result;
}

float u32_to_float(unsigned int val) {
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





void calculateBits(float p0[BLOCKSIZE], float p1[BLOCKSIZE], unsigned int *bits, short *c){
    int change = 0;
    for (int i = 0; i < BLOCKSIZE; ++i) {
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


void sumProductDecoding(unsigned int max_iter, unsigned int endsame, unsigned int message[BLOCKSIZE], unsigned int decodedmsg[BLOCKSIZE], unsigned int *iterations) {

#pragma HLS INTERFACE s_axilite port=return bundle=CTRLS
#pragma HLS INTERFACE s_axilite port=max_iter bundle=CTRLS
#pragma HLS INTERFACE s_axilite port=endsame bundle=CTRLS
#pragma HLS INTERFACE s_axilite port=message bundle=CTRLS
#pragma HLS INTERFACE s_axilite port=decodedmsg bundle=CTRLS
#pragma HLS INTERFACE s_axilite port=iterations bundle=CTRLS

    float prob0[BLOCKSIZE];
    float prob1[BLOCKSIZE];

    float messagef[BLOCKSIZE];
    unsigned int messaged[BLOCKSIZE];

    float deltaP[BLOCKSIZE][PARITYCHECKS];


    float Q0[BLOCKSIZE][PARITYCHECKS];
    float Q1[BLOCKSIZE][PARITYCHECKS];

    short i;
    short i2;
    short j;
    short j2;
    short change = 0;

    // init as max+1 so you know how it finished.
    short itercount = max_iter+1;

    // convert message:
    hw_u32_array_to_float(message, messagef);


    // init arrays
    for (i = 0; i < BLOCKSIZE; ++i) {
        messaged[i] = 0;
        prob1[i] = (float)(1 / (1 + expf((2 * messagef[i]) / (float)SIGMA)));
        prob0[i] = (float)(1 - prob1[i]);
        printf("%f ", prob1[i]);
    }
    printf("\n");

    // get original message
    calculateBits(prob0, prob1, messaged, &change);

    // init P matrix's
    for ( j = 0; j < BLOCKSIZE; ++j) {
        for (i = 0; i < PARITYCHECKS; ++i) {
            // DeltaP - First time out of the main loop..
            deltaP[j][i] = (float)(prob0[j] - prob1[j]);
            printf("%f ", deltaP[j][i]);
        }
    }
    printf("\n");
    unsigned int same = 0;

    // actual algorithm
    for (unsigned int iter = 0; iter < max_iter; ++iter) {
        // Q:s
        for (j = 0; j < BLOCKSIZE; ++j){
            for (i = 0; i < PARITYCHECKS; ++i) {
                float  deltaQ = 1;
                int row = h_matrix[j][i];
                for (j2 = 0; j2 < BLOCKSIZE; ++j2) {
                    if (j2 != j){
                        for (i2 = 0; i2 < PARITYCHECKS; ++i2) {
                            if (h_matrix[j2][i2] == row) {
                                deltaQ *= deltaP[j2][i2];
                            }
                        }
                    }
                }
                printf("%f ", deltaQ);
                Q0[j][i] = (float)(1.0+deltaQ);
                Q1[j][i] = (float)(1.0-deltaQ);
            }
        }
        printf("\n");
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
                deltaP[j][i] = (float)(prob0matrix - prob1matrix);
            }
        }
        //new P's
        for (j = 0; j < BLOCKSIZE; ++j) {
            for (i = 0; i < PARITYCHECKS; ++i) {
                for ( i2 = 0; i2 < PARITYCHECKS; ++i2) {
                    prob0[j] = (float)(prob0[j] * Q0[j][i2]);
                    prob1[j] = (float)(prob1[j] * Q1[j][i2]);
                }

            }
            //scale
            float sum = (float)(prob0[j] + prob1[j]);
            if (sum != 0) {
                prob0[j] = (float)(prob0[j] / sum);
                prob1[j] = (float)(prob1[j] / sum);
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


void AnothersumProductDecoding(unsigned int max_iter, unsigned int end, unsigned int intmsg[BLOCKSIZE], unsigned int output[BLOCKSIZE], unsigned int *iterations) {

#pragma HLS INTERFACE s_axilite port=return bundle=CRTLS
#pragma HLS INTERFACE s_axilite port=max_iter bundle=CRTLS
#pragma HLS INTERFACE s_axilite port=end bundle=CRTLS
#pragma HLS INTERFACE s_axilite port=intmsg bundle=CRTLS
#pragma HLS INTERFACE s_axilite port=output bundle=CRTLS
#pragma HLS INTERFACE s_axilite port=iterations bundle=CRTLS

    float p0[BLOCKSIZE];
    float p1[BLOCKSIZE];
    float message[BLOCKSIZE];

    float P0[BLOCKSIZE][PARITYCHECKS];
    float P1[BLOCKSIZE][PARITYCHECKS];

    float deltaP[BLOCKSIZE][PARITYCHECKS];

    float deltaQ[BLOCKSIZE][PARITYCHECKS];

    float Q0[BLOCKSIZE][PARITYCHECKS];
    float Q1[BLOCKSIZE][PARITYCHECKS];

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
    for (i = 0; i < BLOCKSIZE; ++i) {
        p1[i] = 1 / (1 + expf((2 * message[i]) / (float)SIGMA));
        p0[i] = 1 - p1[i];
    }

    // get original message
    calculateBits(p0, p1, output, &change);

    // init P matrix's
    for (i = 0; i < PARITYCHECKS; ++i) {
        for ( j = 0; j < BLOCKSIZE; ++j) {

            P0[j][i] = p0[j];
            P1[j][i] = p1[j];
        }
    }
    // actual algorithm
    for (int iter = 0; iter < max_iter; ++iter) {
        // DeltaP
        for ( i = 0; i < PARITYCHECKS; ++i) {
            for ( j = 0; j < BLOCKSIZE; ++j) {
                deltaP[j][i] = P0[j][i] - P1[j][i];
                printf("%f ", deltaP[j][i]);
            }
        }
        printf("\n");

        // Q:s
        for (i = 0; i < PARITYCHECKS; ++i) {
            for (j = 0; j < BLOCKSIZE; ++j) {
                deltaQ[j][i] = 1;
                for (i2 = 0; i2 < PARITYCHECKS; ++i2) {
                    for (j2 = 0; j2 < BLOCKSIZE; ++j2) {
                        if (j != j2 && h_matrix[j][i] == h_matrix[j2][i2]) {
                            deltaQ[j][i] *= deltaP[j2][i2];
                        }
                    }
                }

            }
        }

        for (i = 0; i < PARITYCHECKS; ++i) {
            for (j = 0; j < BLOCKSIZE; ++j) {
                Q0[j][i] = (float)(0.5 * (1.0+deltaQ[j][i]));
                Q1[j][i] = (float)(0.5 * (1.0-deltaQ[j][i]));
            }
        }

        for (j = 0; j < BLOCKSIZE; ++j) {
            for (i = 0; i < PARITYCHECKS; ++i) {
                P0[j][i] = p0[j];
                P1[j][i] = p1[j];
                for (i2 = 0; i2 < PARITYCHECKS; ++i2) {
                    if (i != i2) {
                        P0[j][i] *= Q0[j][i2];
                        P1[j][i] *= Q1[j][i2];
                    }
                }
            }
        }
        break; ///TEST

        //scale
        for ( j = 0; j < BLOCKSIZE; ++j) {
            for ( i = 0; i < PARITYCHECKS; ++i) {
                float sum = (P0[j][i] + P1[j][i]);
                if (sum != 0) {
                    P0[j][i] /= sum;
                    P1[j][i] /= sum;
                }
            }
        }
        //new P's
        for (j = 0; j < BLOCKSIZE; ++j) {
            for (i = 0; i < PARITYCHECKS; ++i) {
                for ( i2 = 0; i2 < PARITYCHECKS; ++i2) {
                    p0[j] *= Q0[j][i2];
                    p1[j] *= Q1[j][i2];
                }

            }
        }

        //scale
        for ( j = 0; j < BLOCKSIZE; ++j) {
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
        break;
    }
}

void hw_u32_array_to_float(unsigned int input[BLOCKSIZE], float *result) {
    for (int i = 0; i < BLOCKSIZE; ++i){
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




void OLDsumProductDecoding(unsigned int max_iter, unsigned int end, unsigned int intmsg[BLOCKSIZE], unsigned int output[BLOCKSIZE], unsigned int *iterations) {
#pragma HLS INTERFACE s_axilite port=return bundle=CRTLS
#pragma HLS INTERFACE s_axilite port=max_iter bundle=CRTLS
#pragma HLS INTERFACE s_axilite port=end bundle=CRTLS
#pragma HLS INTERFACE s_axilite port=intmsg bundle=CRTLS
#pragma HLS INTERFACE s_axilite port=output bundle=CRTLS
#pragma HLS INTERFACE s_axilite port=iterations bundle=CRTLS

    float p0[BLOCKSIZE];
    float p1[BLOCKSIZE];
    float message[BLOCKSIZE];
    float P0[CNODES][BLOCKSIZE];
    float P1[CNODES][BLOCKSIZE];
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

    printf("\n>A: ");
    for(i = 0; i < 5; i++){
        printf("%.4f ", message[i]);
    }


    for (i = 0; i < BLOCKSIZE; ++i) {
        p1[i] = (float)(1 / (1 + exp((2 * message[i]) / SIGMA)));
        p0[i] = 1 - p1[i];
    }
    printf("\n>B: ");
    for(i = 0; i < 6; i++){
        printf("%.4f ", p1[i]);
    }

    // get original message
    calculateBits(p0, p1, output, &change);
    // init P matrix's
    for (i = 0; i < CNODES; ++i) {
        for ( j = 0; j < BLOCKSIZE; ++j) {
            P0[i][j] = (float)h_matrix_all[i][j] * p0[j];
            P1[i][j] = (float)h_matrix_all[i][j] * p1[j];
        }
    }
    printf("\n>C: ");
    for(i = 0; i < CNODES; i++){
        if (P1[i][0] != 0){
            printf("%.4f ", P1[i][0]);
        }
    }
    // actual algorithm
    for (int iter = 0; iter < max_iter; ++iter) {
        printf("\n> %d: ", iter);
        // DeltaP
        float deltaP[CNODES][BLOCKSIZE];
        for ( i = 0; i < CNODES; ++i) {
            for ( j = 0; j < BLOCKSIZE; ++j) {
                deltaP[i][j] = P0[i][j] - P1[i][j];
            }
        }
        printf("\n>  A: ");
        for(i = 0; i < CNODES; i++){
            if (deltaP[i][0] != 0){
                printf("%.4f ", deltaP[i][0]);
            }
        }

        // Q:s
        float deltaQ[CNODES][BLOCKSIZE] = {0};
        float Q0[CNODES][BLOCKSIZE] = {0};
        float Q1[CNODES][BLOCKSIZE] = {0};
        for (i = 0; i < CNODES; ++i) {
            for ( j = 0; j < BLOCKSIZE; ++j) {
                float product = 1;
                if (deltaP[i][j] != 0){
                    for ( j2 = 0; j2 < BLOCKSIZE; ++j2) {
                        if (j != j2 && deltaP[i][j2] != 0){
                            product = product * deltaP[i][j2];
                        }
                        deltaQ[i][j] = product;
                        Q0[i][j] = (float)(0.5 * (1+product));
                        Q1[i][j] = (float)(0.5 * (1-product));
                    }
                } else {
                    deltaQ[i][j] = 0;
                    Q0[i][j] = 0;
                    Q1[i][j] = 0;
                }
            }
        }
        printf("\n>  B: ");
        for(i = 0; i < CNODES; i++){
            if (deltaQ[i][0] != 0){
                printf("%.4f ", deltaQ[i][0]);
            }
        }
        printf("\n>  C: ");
        for(i = 0; i < CNODES; i++){
            if (Q1[i][0] != 0){
                printf("%.4f ", Q1[i][0]);
            }
        }
        //new P's
        float np0[BLOCKSIZE] = {0};
        float np1[BLOCKSIZE] = {0};
        for ( i = 0; i < CNODES; ++i) {
            for ( j = 0; j < BLOCKSIZE; ++j) {
                float product = 1;
                float product2 = 1;
                if (Q0[i][j] != 0) {
                    for ( i2 = 0; i2 < CNODES; ++i2) {
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
                    for ( i2 = 0; i2 < CNODES; ++i2) {
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
        for ( i = 0; i < BLOCKSIZE; ++i) {
            p0[i] = np0[i];
            p1[i] = np1[i];
        }

        printf("\n>  D: ");
        for(i = 0; i < 6; i++){
            printf("%.4f ", p1[i]);
        }

        //scale
        for ( j = 0; j < BLOCKSIZE; ++j) {
            for ( i = 0; i < CNODES; ++i) {
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
        printf("\n>  E: ");
        for(i = 0; i < 6; i++){
            printf("%.4f ", p1[i]);
        }
        // Check ending conditions
        calculateBits(p0, p1, output, &change);

        printf("\n>  F: ");
        for(i = 0; i < 80; i++){
            printf("%d", output[i]);
        }
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


void PrintssumProductDecoding(unsigned int max_iter, unsigned int end, unsigned int intmsg[BLOCKSIZE], unsigned int output[BLOCKSIZE], unsigned int *iterations) {
#pragma HLS INTERFACE s_axilite port=return bundle=CRTLS
#pragma HLS INTERFACE s_axilite port=max_iter bundle=CRTLS
#pragma HLS INTERFACE s_axilite port=end bundle=CRTLS
#pragma HLS INTERFACE s_axilite port=intmsg bundle=CRTLS
#pragma HLS INTERFACE s_axilite port=output bundle=CRTLS
#pragma HLS INTERFACE s_axilite port=iterations bundle=CRTLS

    float p0[BLOCKSIZE];
    float p1[BLOCKSIZE];
    float message[BLOCKSIZE];
    float P0[BLOCKSIZE][PARITYCHECKS];
    float P1[BLOCKSIZE][PARITYCHECKS];
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

    printf("\n>A: ");
    for(i = 0; i < 5; i++){
        printf("%.4f ", message[i]);
    }

    // init arrays
    for (i = 0; i < BLOCKSIZE; ++i) {
        p1[i] = 1 / (1 + exp((2 * message[i]) / SIGMA));
        p0[i] = 1 - p1[i];
    }
    printf("\n>B: ");
    for(i = 0; i < 6; i++){
        printf("%.4f ", p1[i]);
    }
    // get original message
    calculateBits(p0, p1, output, &change);
    // init P matrix's
    for (i = 0; i < PARITYCHECKS; ++i) {
        for ( j = 0; j < BLOCKSIZE; ++j) {

            P0[j][i] = p0[j];
            P1[j][i] = p1[j];
        }
    }
    printf("\n>C: ");
    for(i = 0; i < PARITYCHECKS; i++){
        printf("%.4f ", P1[0][i]);
    }
    // actual algorithm
    for (int iter = 0; iter < max_iter; ++iter) {
        printf("\n> %d: ", iter);
        // DeltaP
        float deltaP[BLOCKSIZE][PARITYCHECKS];
        for ( i = 0; i < PARITYCHECKS; ++i) {
            for ( j = 0; j < BLOCKSIZE; ++j) {
                deltaP[j][i] = P0[j][i] - P1[j][i];
            }
        }
        printf("\n>  A: ");
        for(i = 0; i < PARITYCHECKS; i++){
            printf("%.4f ", deltaP[0][i]);
        }

        // Q:s
        float deltaQ[BLOCKSIZE][PARITYCHECKS] = {0};
        float Q0[BLOCKSIZE][PARITYCHECKS] = {0};
        float Q1[BLOCKSIZE][PARITYCHECKS] = {0};
        for (i = 0; i < PARITYCHECKS; ++i) {
            for (j = 0; j < BLOCKSIZE; ++j) {
                deltaQ[j][i] = 1;
                for (i2 = 0; i2 < PARITYCHECKS; ++i2) {
                    for (j2 = 0; j2 < BLOCKSIZE; ++j2) {
                        if (j != j2 && h_matrix[j][i] == h_matrix[j2][i2]) {
                            deltaQ[j][i] *= deltaP[j2][i2];
                        }
                    }
                }

            }
        }
        printf("\n>  B: ");
        for(i = 0; i < PARITYCHECKS; i++){
            printf("%.4f ", deltaQ[0][i]);
        }
        for (i = 0; i < PARITYCHECKS; ++i) {
            for (j = 0; j < BLOCKSIZE; ++j) {
                Q0[j][i] = (float)(0.5 * (1.0+deltaQ[j][i]));
                Q1[j][i] = (float)(0.5 * (1.0-deltaQ[j][i]));
            }
        }
        printf("\n>  C: ");
        for(i = 0; i < PARITYCHECKS; i++){
            printf("%.4f ", Q1[0][i]);
        }


        for (j = 0; j < BLOCKSIZE; ++j) {
            for (i = 0; i < PARITYCHECKS; ++i) {
                P0[j][i] = p0[j];
                P1[j][i] = p1[j];
                for (i2 = 0; i2 < PARITYCHECKS; ++i2) {
                    if (i != i2) {
                        P0[j][i] *= Q0[j][i2];
                        P1[j][i] *= Q1[j][i2];
                    }
                }
            }
        }

        //scale
        for ( j = 0; j < BLOCKSIZE; ++j) {
            for ( i = 0; i < PARITYCHECKS; ++i) {
                float sum = (P0[j][i] + P1[j][i]);
                if (sum != 0) {
                    P0[j][i] /= sum;
                    P1[j][i] /= sum;
                }
            }
        }
        //new P's
        for (j = 0; j < BLOCKSIZE; ++j) {
            for (i = 0; i < PARITYCHECKS; ++i) {
                for ( i2 = 0; i2 < PARITYCHECKS; ++i2) {
                    p0[j] *= Q0[j][i2];
                    p1[j] *= Q1[j][i2];
                }

            }
        }
        printf("\n>  D: ");
        for(i = 0; i < 6; i++){
            printf("%.4f ", p1[i]);
        }

        //scale
        for ( j = 0; j < BLOCKSIZE; ++j) {
            float sum = (p0[j] + p1[j]);
            if (sum != 0) {
                p0[j] /= sum;
                p1[j] /= sum;
            }
        }

        printf("\n>  E: ");
        for(i = 0; i < 6; i++){
            printf("%.4f ", p1[i]);
        }

        // Check ending conditions
        calculateBits(p0, p1, output, &change);

        printf("\n>  F: ");
        for(i = 0; i < 80; i++){
            printf("%d", output[i]);
        }

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




#define true 1
#define false 0

void hardDecisionDecoding(unsigned int max_iter, unsigned int message[BLOCKSIZE], unsigned int output[BLOCKSIZE], unsigned int *iterations){

    #pragma HLS INTERFACE s_axilite port=return bundle=CRTLS
    #pragma HLS INTERFACE s_axilite port=max_iter bundle=CRTLS
    #pragma HLS INTERFACE s_axilite port=message bundle=CRTLS
    #pragma HLS INTERFACE s_axilite port=output bundle=CRTLS
    #pragma HLS INTERFACE s_axilite port=iterations bundle=CRTLS


    *iterations = max_iter+1;

    short int c;
    short int i;
    short int j;

    unsigned int v_nodes[BLOCKSIZE];
    for (j = 0; j < BLOCKSIZE; ++j) {
        output[j] = message[j];
        v_nodes[j] = message[j];// padding to get away from negatives and keep uint
    }

    for (unsigned int iter = 0; iter < max_iter; ++iter) {
        short int satisfied = true;
        for (c = 0; c < CNODES; ++c) {
            unsigned int c_node = 0;
            for (j = 0; j < BLOCKSIZE; ++j) {
                for (i = 0; i < PARITYCHECKS; ++i) {
                    if (h_matrix[j][i] == c) {
                        c_node ^= output[j];
                    }
                }
            }
            if (c_node == 1){
                satisfied = false;
            }
            for (j = 0; j < BLOCKSIZE; ++j) {
                for (i = 0; i < PARITYCHECKS; ++i) {
                    if (h_matrix[j][i] == c) {
                        v_nodes[j] += output[j]^c_node;
                    }
                }
            }
        }

        if (satisfied){
            *iterations = iter +1;
            break;
        }

        for (j = 0; j < BLOCKSIZE; ++j) {
            if (v_nodes[j] > (PARITYCHECKS+1)/2){
                output[j] = 1;
            } else if (v_nodes[j] < (PARITYCHECKS+1)/2){
                output[j] = 0;
            }
            v_nodes[j] = output[j];
        }
    }
}





void oldhardDecisionDecoding(unsigned int max_iter, unsigned int message[BLOCKSIZE], unsigned int output[BLOCKSIZE], unsigned int *iterations){

#pragma HLS INTERFACE s_axilite port=return bundle=CRTLS
#pragma HLS INTERFACE s_axilite port=max_iter bundle=CRTLS
#pragma HLS INTERFACE s_axilite port=message bundle=CRTLS
#pragma HLS INTERFACE s_axilite port=output bundle=CRTLS
#pragma HLS INTERFACE s_axilite port=iterations bundle=CRTLS

    //returns -1 if max_iter is reached
    *iterations = max_iter+1;

    unsigned int i = 0;
    unsigned int j = 0;
    unsigned int v_nodes[BLOCKSIZE];

    for (i = 0; i < BLOCKSIZE; ++i) {
        output[i] = message[i];
        v_nodes[i] = message[i];
    }

    for (unsigned int iter = 0; iter < max_iter; ++iter) {
        short int satisfied = true;
        unsigned int c_nodes[CNODES] = {0};
        for (i = 0; i < CNODES; ++i) {
            for (j = 0; j < BLOCKSIZE; ++j) {
                if (h_matrix_all[i][j] == 1){
                    c_nodes[i] ^= output[j];
                }
            }
            for (j = 0; j < BLOCKSIZE; ++j) {
                if (h_matrix_all[i][j] == 1) {
                    v_nodes[j] += output[j] ^ c_nodes[i];
                }
            }
            if (c_nodes[i] == 1){
                satisfied = false;
            }
        }
        if (satisfied == 1){
            *iterations = iter +1;
            break;
        }
        for (j = 0; j < BLOCKSIZE; ++j) {
            unsigned int div = 1;
            for (i = 0; i < CNODES; ++i) {
                div += h_matrix_all[i][j];
            }
            float v = (float)v_nodes[j] /(float) div;
            if (v > 0.5){
                v_nodes[j] = 1;
            } else {
                v_nodes[j] = 0;
            }
        }
        for (i = 0; i < BLOCKSIZE; ++i) {
            output[i] = v_nodes[i];
        }
    }
}
