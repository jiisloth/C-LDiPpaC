#include <stdio.h>
#include "testcases.h"
#include "Testbench.h"

#define printfv(v,...) do{if(verbose_lvl>=(v)){printf(__VA_ARGS__);}}while(0);

#define MAXITER 10
#define ENDSAME 3
#define SUMPRODUCT 0
#define HARDDECISION 1


#define TESTCASES 4 // 0 = all
#define CASEOFFSET 0

int verbose_lvl = 2; //-1: Nothing, 0: Final results, 1: Individual test results, 2: Individual test outputs
int runsingle = 0; //0 for all, others for case #

int algorithm = SUMPRODUCT; // Remember to comment/uncomment from run_case

int main(){
    if (algorithm == SUMPRODUCT) {
        printfv(0,"Running Sum Product Decoding..\n");
    } else if (algorithm == HARDDECISION) {
        printfv(0,"Running Hard Decision Decoding..\n");
    } else {
        printf("Wrong algorithm configuration!\n");
        return 1;
    }
    int cases_to_test = TESTCASES;
    if (cases_to_test > TESTCASECOUNT || cases_to_test < 1){
    	cases_to_test = TESTCASECOUNT;
    	if (TESTCASES > 0){
    		printf("No enough test cases.. Running all %d cases\n", TESTCASECOUNT);
    	}
    }

    int errors = 0;
    if (runsingle <= 0) {
        for (int t = 0; t < cases_to_test; t++){
        	int tc = t + CASEOFFSET;
        	tc = tc % TESTCASECOUNT;

            printfv(1,"\nTest #%d, Test case: %d, ", t, tc+1);

            errors += run_case(*cases[tc]);

        }
        if (errors == cases_to_test){
            printfv(0,"\nAll %d tests FAILED!\n", cases_to_test);
        } else if (errors == 0){
            printfv(0,"\nAll %d tests PASSED!\n", cases_to_test);
        } else {
            printfv(0,"\nTests PASSED: %3d/%d\n", cases_to_test-errors, cases_to_test);
            printfv(0,"Tests FAILED: %3d/%d\n", errors, cases_to_test);
        }
    } else if (runsingle <= TESTCASECOUNT) {
        printfv(1,"\nTest case: %d, ", runsingle);
        errors += run_case(*cases[runsingle-1]);
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
    printfv(1,"Flipped bits in the message: %d.\n", errorsinmessage);

    int output[BLOCKSIZE] = {0};
    int itercount = 0;
    int *iterations = &itercount;

    int endsame = ENDSAME;



    printfv(1,"Running...\n");
    if (algorithm == SUMPRODUCT){
    	//unsigned int message[BLOCKSIZE];
        //float_array_to_u32(custard.receivedfloat, message);
        sumProductDecoding(MAXITER, endsame, custard.receivedfloat, output, iterations);
    } else if (algorithm == HARDDECISION) {
        //hardDecisionDecoding(MAXITER, custard.receivedhard, output, iterations);
    }

    if (itercount < 0){
        printfv(1,"- Hardware maximum iterations reached. (%d)\n", -itercount);
    }
    else if (itercount == 0){
        printfv(1,"- Max iterations reached. (%d)\n", MAXITER);
    } else {
        printfv(1,"- Finished on iteration # %d.\n", itercount);
    }
    int crcresult = -1;
    if (crc_bits[0] == 1) {
        crcresult = test_crc(output);
        if (crcresult == 1) {
            printfv(1,"- CRC failed\n");
        } else {
            printfv(1,"- CRC succesfull\n");
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
int test_crc(int output[BLOCKSIZE]){
    int crc_padded_len = MSGLEN + CRCLEN - 1;
    int test[BLOCKSIZE] = {0};
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

void float_array_to_fixed(float input[BLOCKSIZE], int *result) {
    for (int i = 0; i < BLOCKSIZE; ++i){
        result[i] = (input[i] * (1 << 16));
    }
}


