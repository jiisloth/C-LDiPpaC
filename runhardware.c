#include <stdio.h>
#include "platform.h"
#include "xil_printf.h"
#include <xparameters.h>
#include <xsumproductdecoding.h>
#include <xharddecisiondecoding.h>
#include "xtime_l.h"

#include <stdio.h>
#include "testcases.h"
#include "Testbench.h"


#define printfv(v,...) do{if(verbose_lvl>=(v)){printf(__VA_ARGS__);}}while(0);

#define MAXITER 10
#define ENDSAME 3

#define BOTH -1
#define SUMPRODUCT 0
#define HARDDECISION 1

void run_sumProductDecoding(unsigned int max_iter, unsigned int endsame, unsigned int message[BLOCKSIZE], unsigned int decodedmsg[BLOCKSIZE], unsigned int *iterations);
void run_hardDecisionDecoding(unsigned int max_iter, unsigned int message[BLOCKSIZE], unsigned int output[BLOCKSIZE], unsigned int *iterations);


int runs[2] = {SUMPRODUCT,HARDDECISION};
int verbose_lvl = 1; //-2: Nothing, -1: Combined Final results, 0: Algorithm Final results, 1: Individual test results, 2: Individual test outputs
int runsingle = -1; //-1 for all

int algorithm = SUMPRODUCT; // Remember to comment/uncomment from run_case

int main(){
    init_platform();
    int error = hardware_test();
    cleanup_platform();
    if (error != 0){
        printfv(-1, "\nError: %d ", error)
    }
    if (error == 1){
        printfv(-1, "Failed all tests.\n");
    } else if (error == 2){
        printfv(-1, "Configuration failure.\n");
    }
    return error;
}

XSumproductdecoding doSumproductdecoding;
XSumproductdecoding_Config *doSumproductdecoding_cfg;

XHarddecisiondecoding doHarddecisiondecoding;
XHarddecisiondecoding_Config *doHarddecisiondecoding_cfg;


int hardware_test(){
    Result *results[TESTCASECOUNT*2];

    int status;

    XTime gbl_time_before_test;
    XTime gbl_time_after_test;

    int resline = 0;

    int runcount = algorithm+1;
    int a = algorithm;
    if (algorithm == BOTH) {
        runcount = 2;
        a = 0;
    }

    for (int alg = a; alg < runcount; alg++) {
        int errors = 0;
        if (alg == SUMPRODUCT) {
            doSumproductdecoding_cfg = XSumproductdecoding_LookupConfig(XPAR_SUMPRODUCTDECODING_0_DEVICE_ID);
            if (!doSumproductdecoding_cfg) {
                printf("Error loading conf for spd\n");
                return 2;
            }
            status = XSumproductdecoding_CfgInitialize(&doSumproductdecoding, doSumproductdecoding_cfg);
            if (status != XST_SUCCESS){
                printf("Error initializing spd\n");
                return 2;
            }
            printfv(0, "Running Sum Product Decoding..\n");
        } else if (alg == HARDDECISION) {
            doHarddecisiondecoding_cfg = XHarddecisiondecoding_LookupConfig(XPAR_HARDDECISIONDECODING_0_DEVICE_ID);
            if (!doHarddecisiondecoding_cfg) {
                printf("Error loading conf for hdd\n");
                return 2;
            }
            status = XHarddecisiondecoding_CfgInitialize(&doHarddecisiondecoding, doHarddecisiondecoding_cfg);
            if (status != XST_SUCCESS){
                printf("Error initializing hdd\n");
                return 2;
            }
            printfv(0, "Running Hard Decision Decoding..\n");
        } else {
            printf("Wrong algorithm configuration!\n");
            return 2;
        }
        if (runsingle < 0) {
            for (int t = 0; t < TESTCASECOUNT; t++) {
                Result res = {-1,0,0,0,0,0,0};
                printfv(1, "\nTestcase: %d, ", t);
                int ret = hw_run_case(*cases[t], alg, t, res);
                if (ret != 0){
                    return ret;
                }
                *results[resline] = res;
                if (res.crc_error > 0 || res.gf_error > 0){
                    errors += 1;
                }
            }
            if (errors == TESTCASECOUNT) {
                printfv(0, "\nAll %d tests FAILED!\n", TESTCASECOUNT);
            } else if (errors == 0) {
                printfv(0, "\nAll %d tests PASSED!\n", TESTCASECOUNT);
            } else {
                printfv(0, "\nTests PASSED: %3d/%d\n", TESTCASECOUNT - errors, TESTCASECOUNT);
                printfv(0, "Tests FAILED: %3d/%d\n", errors, TESTCASECOUNT);
            }
        } else if (runsingle < TESTCASECOUNT) {
            Result res = {-1,0,0,0,0,0,0};

            int ret = hw_run_case(*cases[runsingle], alg, runsingle, res);
            if (ret != 0){
                return ret;
            }
            *results[resline] = res;
            if (res.crc_error > 0 || res.gf_error > 0){
                errors += 1;
            }
            if (errors == 0) {
                printfv(0, "\nTest passed!\n");
            } else {
                printfv(0, "\nTest failed!\n");
            }
        }
    }
    if (verbose_lvl > -2){
        print_results(*results);
    }
    return 0;
}

void print_results(Result results[TESTCASECOUNT*2]){
    Result spd_total = {0,0,0,0,0,0,0};
    Result hdd_total = {1,0,0,0,0,0,0};

    printf("\n######################################\n\n");
    printf("ALGORITHM TESTS ERRORS CRC_FAILS ITERATIONS ITERMAXED DURATION\n");
    for(int i = 0; i < BLOCKSIZE; i++) {
        if (results[i].alg == SUMPRODUCT){
            printf("SPD       %5d  %5d     %5d      %5d     %5d   %3.2f\n", results[i].tcase, results[i].gf_error, results[i].crc_error, results[i].iter, results[i].max_iter, results[i].time);
            spd_total.tcase += 1;
            if (results[i].gf_error > 0){
                spd_total.gf_error += 1;
                if (results[i].crc_error == 0){
                    spd_total.crc_error += 1;
                }
            } else if (results[i].crc_error == 1){
                spd_total.gf_error += 1;
                spd_total.crc_error += 1;
            }
            spd_total.iter += results[i].iter;
            spd_total.max_iter += results[i].max_iter;
            spd_total.time += results[i].time;
        } else if (results[i].alg == HARDDECISION){
            printf("HDD       %5d  %5d     %5d      %5d     %5d   %3.2f\n", results[i].tcase, results[i].gf_error, results[i].crc_error, results[i].iter, results[i].max_iter, results[i].time);
            hdd_total.tcase += 1;
            if (results[i].gf_error > 0){
                hdd_total.gf_error += 1;
                if (results[i].crc_error == 0){
                    hdd_total.crc_error += 1;
                }
            } else if (results[i].crc_error == 1){
                hdd_total.gf_error += 1;
                hdd_total.crc_error += 1;
            }
            hdd_total.iter += results[i].iter;
            hdd_total.max_iter += results[i].max_iter;
            hdd_total.time += results[i].time;
        } else {
            break;
        }
    }
    printf("\n\n######################################\n\n");
    printf("ALGORITHM TESTS ERRORS CRC_FAILS ITERATIONS ITERMAXED DURATION\n");
    printf("SPD       %5d  %5d     %5d      %5d     %5d   %3.2f\n", spd_total.tcase, spd_total.gf_error, spd_total.crc_error, spd_total.iter, spd_total.max_iter, spd_total.time);
    printf("HDD       %5d  %5d     %5d      %5d     %5d   %3.2f\n", hdd_total.tcase, hdd_total.gf_error, hdd_total.crc_error, hdd_total.iter, hdd_total.max_iter, hdd_total.time);
    printf("TOTAL     %5d  %5d     %5d      %5d     %5d   %3.2f\n", hdd_total.tcase + spd_total.tcase, hdd_total.gf_error + spd_total.gf_error, hdd_total.crc_error + spd_total.crc_error, hdd_total.iter + spd_total.iter, hdd_total.max_iter + spd_total.max_iter, hdd_total.time + spd_total.time);
    int totc = hdd_total.tcase + spd_total.tcase;
    printf("TOTAL%%           %2.2f     %2.2f      %2.2f     %2.2f   %3.2f\n", (float)(hdd_total.gf_error + spd_total.gf_error)/totc, (float)(hdd_total.crc_error + spd_total.crc_error)/totc, (float)(hdd_total.iter + spd_total.iter)/totc, (float)(hdd_total.max_iter + spd_total.max_iter)/totc, (float)(hdd_total.time + spd_total.time)/totc);

}

int hw_run_case(Testcase custard, int alg, int caseid, Result resultline){
    // custard = case but cant use case as variable...

    XTime gbl_time_before_test;
    XTime gbl_time_after_test;

    //get initial errors..
    int errorsinmessage = 0;

    int i;
    for(i = 0; i < BLOCKSIZE; i++) {
        if (custard.receivedhard[i] != custard.GoldRef[i]) {
            errorsinmessage += 1;
        }
    }
    printfv(1,"Error bits: %d.\n", errorsinmessage);

    unsigned int output[BLOCKSIZE] = {0};
    unsigned int iterations = -1;

    unsigned int max_iter = MAXITER;
    unsigned int endsame = ENDSAME;


    printfv(1,"Running...\n");
    if (algorithm == SUMPRODUCT){
        unsigned int message[BLOCKSIZE];
        float_array_to_u32(custard.receivedfloat, message);

        XTime_GetTime(&gbl_time_before_test);
        run_sumProductDecoding(max_iter, endsame, message, output, &iterations);
        XTime_GetTime(&gbl_time_after_test);
    } else if (algorithm == HARDDECISION) {
        XTime_GetTime(&gbl_time_before_test);
        run_hardDecisionDecoding(max_iter, custard.receivedhard, output, &iterations);
        XTime_GetTime(&gbl_time_after_test);
    }
    int real_iter = iterations;
    int max_i_f = 0;
    if (iterations == max_iter+1){
        real_iter = max_iter;
        max_i_f = 1;
        printfv(1,"- Max iterations reached. (%d)\n", max_iter);
    } else {
        printfv(1,"- Finished on iteration # %d.\n", iterations);
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

    int biterrors = 0;
    int errorline[BLOCKSIZE] = {0};

    printfv(1,"- Comparing to GoldRef: ");

    printfv(2,"\n");

    int goldrefresult = 0;
    for(i = 0; i < BLOCKSIZE; i++){
        if(output[i] != custard.GoldRef[i]){
            biterrors += 1;
            errorline[i] = 1;
            goldrefresult = 1;
        }
        printfv(2,"%d", output[i]);
    }
    if (biterrors > 0){
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
            printfv(2,"  - Decoding failed with %d errors!\n", biterrors);
        } else {
            printfv(1,"Fail! (%d errors.)\n", biterrors);
        }
    } else {
        if (verbose_lvl >= 2){
            printfv(2,"  - Decoding succesfull!\n");
        }  else {
            printfv(1,"OK!\n");
        }
    }
    if (crcresult != -1)  {
        if (crcresult != goldrefresult){
            if (crcresult == 1){
                printfv(1,"  - CRC gave false negative!\n");
            } else {
                printfv(1,"  - CRC gave false positive!\n");
            }
        }
    }
    float time = (float) (gbl_time_after_test - gbl_time_before_test)/(COUNTS_PER_SECOND);
    resultline = (Result){alg,caseid,biterrors,crcresult,real_iter,max_i_f,time};

    return 0;
}
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



int run_sumProductDecoding(unsigned int max_iter, unsigned int endsame, unsigned int message[BLOCKSIZE], unsigned int decodedmsg[BLOCKSIZE], unsigned int *iterations) {

    XSumproductdecoding_Set_max_iter(&doSumproductdecoding, max_iter);
    XSumproductdecoding_Set_end_r(&doSumproductdecoding, end);
    XSumproductdecoding_Write_h_matrix_Words(&doSumproductdecoding, 0, h_matrix, HROWS*MSGLEN);
    XSumproductdecoding_Write_message_Words(&doSumproductdecoding, 0, message, MSGLEN);
    XSumproductdecoding_Write_output_r_Words(&doSumproductdecoding, 0, output, MSGLEN);


    XSumproductdecoding_Start(&doSumproductdecoding);
    while(!XSumproductdecoding_IsDone(&doSumproductdecoding));


    XSumproductdecoding_Read_output_r_Words(&doSumproductdecoding, 0, output, MSGLEN);
    *iterations = XSumproductdecoding_Get_iterations(&doSumproductdecoding);
    return 0;
}

int run_hardDecisionDecoding(unsigned int max_iter, unsigned int message[BLOCKSIZE], unsigned int output[BLOCKSIZE], unsigned int *iterations){
    XHarddecisiondecoding_Set_max_iter(&doHarddecisiondecoding, max_iter);
    XHarddecisiondecoding_Write_h_matrix_Words(&doHarddecisiondecoding, 0, h_matrix, HROWS*MSGLEN);
    XHarddecisiondecoding_Write_message_Words(&doHarddecisiondecoding, 0, message, MSGLEN);
    XHarddecisiondecoding_Write_output_r_Words(&doHarddecisiondecoding, 0, output, MSGLEN);

    XHarddecisiondecoding_Start(&doHarddecisiondecoding);
    while(!XHarddecisiondecoding_IsDone(&doHarddecisiondecoding));

    XHarddecisiondecoding_Read_output_r_Words(&doHarddecisiondecoding, 0, output, MSGLEN);
    *iterations = XHarddecisiondecoding_Get_iterations(&doHarddecisiondecoding);

    return 0;
}

