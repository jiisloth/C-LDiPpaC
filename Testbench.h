void sumProductDecoding(int max_iter, int endsame, float message[BLOCKSIZE], int decodedmsg[BLOCKSIZE], int *iterations);
void hardDecisionDecoding(int max_iter, int message[BLOCKSIZE], int output[BLOCKSIZE], int *iterations);

int run_case(Testcase custard);
int hardware_test();
int test_crc(int output[BLOCKSIZE]);

unsigned int float_to_u32(float val);
float u32_to_float(unsigned int val);
void float_array_to_fixed(float input[BLOCKSIZE], int *result);
void float_array_to_u32(float input[BLOCKSIZE], unsigned int *result);
void u32_array_to_float(unsigned int input[BLOCKSIZE], float *result);

typedef struct Result {
    int alg;
    int tcase;
    int gf_error;
    int crc_error;
    int iter;
    int max_iter;
    float time;
} Result;

int hw_run_case(Testcase custard, int alg, int caseid, Result resultline);
void print_results(Result *results);
