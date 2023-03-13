#define BLOCKSIZE 1024
#define PARITYCHECKS 5
#define CNODES 320
#define SIGMA 1.199526

const short int h_matrix[BLOCKSIZE][PARITYCHECKS] = {{19, 121, 168, 236, 264}, {19, 65, 166, 213, 269}, {19, 113, 181, 212, 291}, {19, 125, 149, 237, 311}, {16, 125, 177, 221, 297}, {20, 94, 187, 243, 316}, {20, 88, 157, 199, 284}, {20, 96, 165, 234, 314}, {20, 118, 139, 198, 309}, {20, 75, 165, 203, 281}, {20, 96, 184, 229, 258}, {20, 114, 183, 251, 274}, {20, 69, 129, 246, 290}, {20, 68, 142, 251, 299}, {20, 111, 161, 223, 317}, {20, 115, 159, 218, 276}, {20, 93, 152, 208, 281}, {20, 73, 162, 231, 283}, {20, 126, 153, 224, 306}, {16, 84, 145, 254, 285}, {16, 87, 150, 232, 292}, {21, 94, 158, 193, 311}, {21, 91, 153, 217, 276}, {21, 88, 171, 192, 272}, {21, 98, 190, 216, 257}, {21, 105, 179, 251, 265}, {21, 96, 148, 251, 294}, {21, 108, 154, 247, 261}, {21, 72, 137, 194, 311}, {21, 108, 149, 216, 295}, {21, 65, 191, 200, 284}, {21, 103, 128, 234, 268}, {17, 103, 180, 232, 276}, {21, 121, 137, 223, 310}, {21, 92, 159, 255, 309}, {21, 110, 151, 227, 297}, {17, 92, 184, 217, 260}, {22, 89, 139, 254, 303}, {22, 103, 160, 214, 279}, {22, 73, 131, 209, 304}, {22, 69, 157, 251, 309}, {22, 111, 150, 229, 269}, {22, 114, 168, 195, 276}, {22, 126, 156, 211, 262}, {22, 73, 131, 239, 309}, {22, 109, 137, 204, 298}, {22, 84, 174, 216, 273}, {22, 77, 159, 215, 290}, {22, 83, 170, 245, 289}, {22, 77, 173, 208, 302}, {22, 92, 141, 210, 286}, {22, 94, 150, 234, 289}, {17, 104, 177, 241, 294}, {23, 102, 178, 245, 315}, {23, 79, 128, 249, 271}, {23, 69, 173, 227, 257}, {23, 127, 172, 199, 300}, {23, 117, 157, 211, 317}, {23, 92, 140, 242, 304}, {23, 125, 151, 246, 260}, {23, 76, 131, 225, 304}, {23, 121, 143, 232, 308}, {23, 89, 145, 194, 260}, {23, 66, 129, 204, 280}, {23, 94, 141, 238, 272}, {23, 120, 176, 250, 290}, {23, 90, 139, 218, 312}, {23, 88, 167, 247, 308}, {17, 81, 151, 228, 281}, {24, 122, 146, 252, 298}, {24, 77, 154, 194, 269}, {24, 104, 147, 219, 299}, {24, 107, 186, 220, 299}, {24, 98, 188, 239, 274}, {24, 119, 148, 226, 297}, {24, 76, 160, 202, 286}, {24, 101, 185, 222, 304}, {24, 111, 161, 243, 298}, {24, 66, 170, 245, 297}, {24, 86, 156, 224, 293}, {24, 80, 190, 216, 278}, {24, 89, 171, 205, 272}, {24, 88, 139, 218, 277}, {24, 114, 176, 210, 294}, {17, 79, 188, 215, 293}, {25, 68, 183, 227, 278}, {25, 99, 165, 194, 264}, {25, 102, 190, 250, 280}, {25, 124, 160, 208, 262}, {25, 80, 177, 217, 308}, {25, 73, 163, 236, 308}, {25, 68, 160, 216, 256}, {25, 81, 155, 246, 296}, {25, 101, 164, 225, 311}, {25, 110, 184, 251, 274}, {25, 117, 180, 226, 281}, {25, 101, 147, 253, 272}, {25, 106, 140, 216, 304}, {25, 111, 164, 224, 306}, {25, 72, 155, 199, 276}, {17, 79, 178, 231, 295}, {26, 78, 166, 254, 313}, {26, 66, 156, 246, 261}, {26, 86, 168, 204, 308}, {26, 89, 132, 213, 268}, {26, 82, 171, 245, 277}, {26, 117, 150, 237, 275}, {26, 70, 157, 227, 286}, {26, 83, 130, 202, 256}, {26, 68, 162, 252, 318}, {26, 122, 172, 199, 318}, {26, 126, 144, 200, 300}, {26, 89, 149, 234, 304}, {26, 94, 140, 193, 314}, {26, 86, 141, 202, 305}, {26, 95, 162, 242, 295}, {17, 123, 131, 231, 266}, {27, 76, 166, 207, 282}, {27, 88, 147, 251, 305}, {27, 114, 162, 221, 258}, {27, 89, 171, 226, 312}, {27, 85, 134, 210, 275}, {27, 118, 138, 243, 283}, {27, 111, 180, 255, 306}, {27, 88, 158, 244, 270}, {27, 88, 163, 247, 260}, {27, 73, 130, 230, 307}, {27, 94, 149, 230, 270}, {27, 109, 187, 236, 260}, {27, 125, 130, 246, 266}, {27, 79, 131, 198, 310}, {27, 111, 157, 231, 306}, {17, 89, 177, 240, 317}, {28, 83, 185, 194, 295}, {28, 85, 186, 249, 281}, {28, 86, 188, 195, 291}, {28, 126, 179, 202, 312}, {28, 74, 169, 225, 317}, {28, 101, 132, 247, 307}, {28, 121, 128, 193, 312}, {28, 100, 160, 254, 276}, {28, 113, 179, 248, 302}, {28, 84, 167, 230, 269}, {28, 95, 178, 216, 309}, {28, 71, 133, 245, 259}, {28, 116, 160, 192, 290}, {28, 121, 156, 255, 263}, {28, 98, 149, 235, 266}, {17, 77, 157, 244, 267}, {29, 96, 189, 229, 312}, {29, 102, 170, 245, 313}, {29, 90, 171, 249, 310}, {29, 119, 163, 197, 262}, {29, 102, 147, 233, 297}, {29, 99, 144, 247, 259}, {29, 64, 142, 215, 313}, {29, 114, 130, 242, 288}, {29, 120, 136, 194, 271}, {29, 91, 134, 192, 262}, {29, 74, 156, 245, 262}, {29, 65, 154, 219, 287}, {29, 127, 183, 254, 267}, {29, 98, 180, 232, 279}, {29, 121, 190, 228, 267}, {17, 112, 183, 242, 306}, {30, 113, 151, 221, 280}, {30, 72, 132, 234, 311}, {30, 74, 137, 225, 316}, {30, 72, 184, 228, 268}, {30, 101, 151, 226, 292}, {30, 109, 152, 218, 284}, {30, 87, 144, 218, 278}, {30, 124, 143, 243, 277}, {30, 90, 144, 236, 280}, {30, 94, 188, 232, 312}, {30, 121, 147, 242, 263}, {30, 106, 180, 228, 267}, {30, 117, 136, 192, 304}, {30, 110, 151, 229, 311}, {30, 64, 148, 209, 281}, {17, 83, 138, 221, 281}, {31, 113, 148, 216, 277}, {31, 70, 152, 203, 291}, {31, 98, 169, 224, 272}, {31, 66, 154, 238, 287}, {31, 73, 153, 202, 296}, {31, 100, 173, 236, 260}, {31, 79, 176, 220, 293}, {31, 105, 139, 207, 279}, {31, 81, 148, 243, 314}, {31, 127, 130, 219, 301}, {31, 79, 178, 206, 282}, {31, 83, 135, 192, 272}, {31, 72, 185, 203, 317}, {31, 100, 163, 213, 260}, {31, 127, 140, 231, 284}, {17, 112, 184, 254, 269}, {32, 110, 148, 233, 288}, {32, 85, 144, 206, 283}, {32, 108, 153, 199, 309}, {32, 110, 175, 211, 274}, {32, 120, 188, 204, 279}, {32, 123, 189, 241, 302}, {32, 114, 143, 248, 290}, {32, 81, 178, 211, 295}, {32, 70, 131, 196, 303}, {32, 116, 166, 217, 307}, {32, 110, 146, 214, 276}, {32, 78, 140, 244, 273}, {32, 65, 140, 252, 297}, {32, 98, 184, 208, 319}, {32, 75, 188, 224, 305}, {17, 96, 138, 214, 282}, {33, 119, 163, 211, 295}, {33, 122, 187, 205, 317}, {33, 83, 152, 198, 318}, {33, 109, 137, 215, 301}, {33, 118, 140, 201, 298}, {33, 88, 175, 224, 277}, {33, 83, 168, 197, 259}, {33, 101, 150, 219, 293}, {33, 65, 139, 209, 262}, {33, 92, 141, 229, 270}, {33, 82, 165, 198, 293}, {33, 71, 179, 223, 292}, {33, 81, 165, 202, 283}, {33, 122, 141, 202, 306}, {33, 119, 155, 232, 275}, {17, 78, 191, 248, 294}, {34, 69, 153, 222, 280}, {34, 115, 173, 242, 260}, {34, 90, 172, 226, 317}, {34, 85, 141, 199, 305}, {34, 89, 172, 225, 278}, {34, 76, 162, 205, 307}, {34, 95, 159, 243, 261}, {34, 76, 151, 206, 266}, {34, 91, 179, 199, 307}, {34, 75, 179, 208, 318}, {34, 103, 142, 227, 309}, {34, 127, 161, 232, 263}, {34, 64, 132, 248, 302}, {34, 101, 132, 196, 259}, {34, 116, 185, 238, 312}, {17, 118, 163, 221, 262}, {35, 85, 153, 215, 264}, {35, 75, 184, 213, 270}, {35, 72, 176, 193, 307}, {35, 105, 190, 212, 294}, {35, 101, 185, 207, 298}, {35, 118, 137, 193, 269}, {35, 105, 136, 231, 263}, {35, 84, 131, 244, 305}, {35, 120, 158, 215, 266}, {35, 104, 160, 198, 312}, {35, 81, 133, 196, 307}, {35, 104, 182, 243, 285}, {35, 115, 183, 207, 307}, {35, 103, 168, 249, 316}, {35, 109, 141, 228, 307}, {16, 122, 178, 247, 280}, {36, 79, 164, 215, 282}, {36, 123, 176, 237, 277}, {36, 88, 177, 209, 260}, {36, 115, 130, 250, 304}, {36, 70, 188, 251, 268}, {36, 64, 175, 235, 271}, {36, 127, 140, 251, 272}, {36, 124, 135, 225, 313}, {36, 76, 157, 206, 300}, {36, 71, 190, 254, 308}, {36, 80, 139, 204, 299}, {36, 103, 152, 243, 309}, {36, 122, 148, 216, 286}, {36, 121, 146, 221, 280}, {36, 98, 133, 230, 273}, {18, 124, 149, 238, 315}, {37, 68, 178, 243, 303}, {37, 79, 144, 214, 271}, {37, 107, 177, 238, 305}, {37, 97, 153, 194, 277}, {37, 117, 187, 241, 301}, {37, 116, 135, 243, 283}, {37, 69, 146, 224, 304}, {37, 118, 168, 241, 301}, {37, 74, 185, 208, 266}, {37, 124, 150, 216, 271}, {37, 86, 187, 233, 318}, {37, 97, 129, 223, 270}, {37, 94, 135, 200, 257}, {37, 87, 183, 244, 263}, {37, 91, 188, 217, 285}, {18, 84, 142, 217, 285}, {38, 127, 144, 209, 317}, {38, 122, 169, 197, 273}, {38, 110, 162, 244, 285}, {38, 69, 137, 241, 295}, {38, 66, 189, 203, 263}, {38, 86, 188, 250, 313}, {38, 97, 159, 198, 286}, {38, 123, 129, 232, 269}, {38, 75, 144, 211, 275}, {38, 118, 181, 194, 299}, {38, 95, 177, 241, 287}, {38, 100, 155, 237, 300}, {38, 93, 152, 229, 289}, {38, 80, 151, 217, 289}, {38, 65, 145, 244, 301}, {18, 105, 147, 233, 265}, {39, 122, 145, 205, 312}, {39, 94, 149, 254, 310}, {39, 106, 158, 211, 303}, {39, 78, 135, 206, 308}, {39, 84, 177, 222, 313}, {39, 125, 177, 208, 315}, {39, 119, 141, 226, 284}, {39, 109, 171, 255, 277}, {39, 102, 176, 211, 299}, {39, 113, 134, 230, 311}, {39, 120, 161, 248, 306}, {39, 116, 180, 193, 316}, {39, 110, 172, 235, 313}, {39, 84, 133, 226, 302}, {39, 100, 133, 193, 315}, {18, 124, 162, 202, 318}, {40, 101, 180, 211, 301}, {40, 90, 151, 207, 300}, {40, 93, 159, 235, 306}, {40, 82, 189, 206, 279}, {40, 103, 140, 209, 259}, {40, 98, 161, 218, 276}, {40, 96, 129, 222, 300}, {40, 73, 155, 207, 309}, {40, 81, 133, 203, 263}, {40, 124, 160, 205, 277}, {40, 100, 145, 224, 275}, {40, 84, 142, 198, 258}, {40, 82, 164, 249, 291}, {40, 95, 144, 205, 302}, {40, 93, 186, 207, 287}, {18, 99, 188, 239, 312}, {41, 74, 181, 214, 296}, {41, 112, 190, 201, 257}, {41, 76, 134, 235, 257}, {41, 119, 187, 240, 257}, {41, 93, 152, 194, 307}, {41, 71, 173, 220, 279}, {41, 95, 159, 215, 257}, {41, 69, 135, 201, 281}, {41, 107, 191, 247, 294}, {41, 64, 177, 196, 275}, {41, 83, 154, 215, 297}, {41, 64, 166, 230, 306}, {41, 67, 157, 209, 289}, {41, 64, 175, 236, 267}, {41, 93, 167, 215, 286}, {18, 99, 169, 229, 263}, {42, 94, 175, 247, 289}, {42, 116, 139, 224, 290}, {42, 116, 134, 225, 262}, {42, 109, 149, 233, 265}, {42, 114, 164, 247, 274}, {42, 65, 172, 235, 286}, {42, 98, 133, 243, 288}, {42, 93, 141, 244, 264}, {42, 106, 190, 205, 292}, {42, 69, 147, 240, 274}, {42, 113, 186, 195, 307}, {42, 126, 129, 248, 261}, {42, 112, 158, 208, 280}, {42, 79, 146, 198, 310}, {42, 110, 150, 195, 292}, {18, 124, 151, 239, 308}, {43, 83, 150, 195, 308}, {43, 120, 191, 239, 273}, {43, 116, 135, 239, 273}, {43, 96, 164, 235, 273}, {43, 70, 130, 253, 295}, {43, 71, 158, 218, 306}, {43, 124, 128, 206, 272}, {43, 66, 156, 232, 295}, {43, 84, 161, 217, 310}, {43, 113, 152, 205, 316}, {43, 104, 153, 217, 307}, {43, 107, 145, 243, 305}, {43, 119, 134, 210, 310}, {43, 82, 159, 253, 311}, {43, 68, 155, 233, 304}, {18, 65, 181, 208, 256}, {44, 105, 182, 232, 288}, {44, 127, 190, 246, 315}, {44, 123, 136, 212, 265}, {44, 91, 142, 230, 265}, {44, 112, 147, 197, 296}, {44, 105, 155, 238, 292}, {44, 94, 135, 233, 313}, {44, 92, 155, 255, 290}, {44, 104, 134, 228, 303}, {44, 127, 190, 234, 303}, {44, 90, 146, 193, 262}, {44, 112, 187, 219, 279}, {44, 115, 167, 249, 268}, {44, 117, 189, 202, 307}, {44, 107, 162, 217, 274}, {18, 106, 185, 245, 313}, {45, 71, 180, 240, 312}, {45, 78, 186, 251, 293}, {45, 103, 151, 228, 283}, {45, 123, 137, 246, 273}, {45, 108, 131, 213, 284}, {45, 80, 180, 227, 290}, {45, 99, 140, 199, 261}, {45, 120, 155, 220, 303}, {45, 114, 129, 255, 288}, {45, 79, 149, 249, 285}, {45, 85, 153, 230, 315}, {45, 85, 146, 230, 316}, {45, 87, 185, 232, 262}, {45, 120, 152, 195, 267}, {45, 107, 166, 227, 307}, {18, 117, 128, 200, 272}, {46, 121, 154, 235, 319}, {46, 116, 172, 250, 277}, {46, 91, 129, 215, 299}, {46, 115, 145, 219, 264}, {46, 121, 171, 198, 285}, {46, 118, 172, 232, 272}, {46, 82, 157, 237, 278}, {46, 71, 153, 235, 286}, {46, 81, 135, 214, 292}, {46, 110, 189, 247, 267}, {46, 89, 186, 243, 289}, {46, 91, 173, 197, 284}, {46, 112, 178, 243, 256}, {46, 114, 135, 224, 257}, {46, 121, 159, 255, 317}, {18, 115, 150, 235, 272}, {47, 125, 166, 193, 281}, {47, 64, 147, 226, 288}, {47, 116, 165, 195, 291}, {47, 127, 172, 241, 303}, {47, 76, 188, 192, 279}, {47, 116, 172, 192, 283}, {47, 79, 175, 197, 260}, {47, 126, 154, 223, 315}, {47, 82, 148, 237, 294}, {47, 97, 175, 230, 293}, {47, 94, 140, 210, 267}, {47, 119, 163, 216, 274}, {47, 68, 174, 228, 318}, {47, 86, 156, 231, 303}, {47, 85, 133, 219, 313}, {18, 75, 182, 248, 281}, {48, 70, 177, 249, 305}, {48, 103, 169, 241, 276}, {48, 87, 179, 233, 314}, {48, 125, 143, 202, 259}, {48, 104, 135, 243, 264}, {48, 127, 183, 229, 315}, {48, 84, 164, 202, 269}, {48, 112, 131, 221, 275}, {48, 65, 141, 200, 261}, {48, 98, 176, 224, 265}, {48, 66, 143, 239, 286}, {48, 108, 175, 237, 259}, {48, 113, 144, 208, 284}, {48, 111, 158, 250, 276}, {48, 116, 129, 215, 310}, {18, 100, 180, 198, 298}, {49, 82, 149, 222, 285}, {49, 98, 165, 210, 314}, {49, 108, 144, 194, 309}, {49, 80, 188, 200, 271}, {49, 77, 128, 210, 274}, {49, 106, 142, 231, 317}, {49, 90, 145, 197, 287}, {49, 91, 171, 193, 271}, {49, 109, 130, 236, 265}, {49, 106, 141, 253, 295}, {49, 105, 174, 213, 314}, {49, 107, 183, 196, 261}, {49, 115, 135, 247, 315}, {49, 125, 153, 246, 265}, {49, 75, 188, 215, 284}, {18, 113, 138, 217, 304}, {50, 80, 181, 228, 283}, {50, 80, 145, 222, 288}, {50, 90, 180, 212, 311}, {50, 77, 187, 213, 271}, {50, 76, 168, 218, 276}, {50, 109, 179, 249, 285}, {50, 110, 186, 234, 274}, {50, 124, 146, 242, 265}, {50, 101, 182, 240, 295}, {50, 90, 182, 219, 301}, {50, 89, 190, 207, 270}, {50, 115, 132, 242, 270}, {50, 110, 184, 253, 287}, {50, 85, 171, 225, 302}, {50, 105, 189, 198, 287}, {18, 65, 183, 197, 314}, {51, 106, 132, 196, 289}, {51, 126, 138, 237, 271}, {51, 93, 138, 244, 300}, {51, 105, 149, 224, 268}, {51, 87, 155, 208, 279}, {51, 97, 167, 212, 265}, {51, 73, 173, 223, 301}, {51, 113, 191, 201, 315}, {51, 112, 172, 201, 317}, {51, 68, 174, 200, 258}, {51, 101, 164, 230, 315}, {51, 104, 139, 222, 267}, {51, 114, 179, 195, 256}, {51, 73, 148, 217, 280}, {51, 100, 181, 197, 300}, {16, 102, 165, 252, 270}, {52, 112, 169, 218, 297}, {52, 74, 146, 245, 296}, {52, 72, 132, 193, 260}, {52, 110, 182, 242, 289}, {52, 81, 189, 216, 297}, {52, 120, 132, 227, 257}, {52, 92, 145, 192, 300}, {52, 126, 191, 196, 281}, {52, 122, 170, 210, 295}, {52, 104, 176, 195, 316}, {52, 105, 168, 219, 282}, {52, 97, 159, 207, 256}, {52, 117, 182, 206, 284}, {52, 70, 135, 246, 305}, {52, 64, 147, 223, 261}, {19, 100, 165, 240, 271}, {53, 126, 169, 235, 264}, {53, 88, 156, 201, 275}, {53, 126, 185, 240, 288}, {53, 81, 169, 203, 263}, {53, 67, 159, 214, 311}, {53, 99, 139, 207, 297}, {53, 77, 136, 214, 293}, {53, 123, 131, 244, 274}, {53, 67, 179, 193, 291}, {53, 117, 138, 242, 292}, {53, 127, 151, 225, 267}, {53, 71, 156, 252, 264}, {53, 94, 134, 204, 292}, {53, 88, 133, 220, 269}, {53, 116, 159, 213, 262}, {19, 108, 137, 209, 276}, {54, 83, 130, 218, 312}, {54, 67, 191, 206, 310}, {54, 91, 173, 209, 296}, {54, 77, 161, 231, 263}, {54, 95, 180, 226, 268}, {54, 90, 186, 225, 256}, {54, 65, 170, 240, 261}, {54, 67, 176, 241, 318}, {54, 105, 170, 220, 295}, {54, 71, 136, 231, 277}, {54, 75, 155, 207, 318}, {54, 64, 167, 255, 270}, {54, 74, 183, 238, 282}, {54, 93, 137, 214, 308}, {54, 95, 191, 252, 286}, {19, 75, 179, 225, 296}, {55, 82, 165, 245, 287}, {55, 64, 186, 227, 309}, {55, 72, 181, 236, 271}, {55, 108, 133, 198, 306}, {55, 106, 139, 198, 256}, {55, 67, 177, 201, 308}, {55, 92, 181, 193, 285}, {55, 70, 189, 213, 258}, {55, 108, 166, 200, 274}, {55, 109, 148, 196, 317}, {55, 74, 150, 220, 313}, {55, 97, 174, 210, 269}, {55, 93, 167, 201, 313}, {55, 96, 140, 220, 307}, {55, 91, 160, 246, 258}, {19, 87, 164, 216, 314}, {56, 70, 145, 249, 267}, {56, 107, 161, 241, 305}, {56, 66, 136, 211, 284}, {56, 118, 173, 241, 280}, {56, 101, 158, 242, 306}, {56, 78, 143, 210, 270}, {56, 72, 185, 193, 261}, {56, 70, 137, 195, 296}, {56, 78, 170, 206, 276}, {56, 88, 153, 225, 262}, {56, 82, 190, 204, 297}, {56, 66, 156, 224, 265}, {56, 67, 164, 212, 293}, {56, 113, 172, 208, 263}, {56, 80, 134, 252, 262}, {19, 111, 158, 228, 277}, {57, 70, 165, 237, 258}, {57, 104, 182, 231, 319}, {57, 119, 138, 202, 317}, {57, 78, 179, 251, 316}, {57, 77, 183, 223, 259}, {57, 99, 173, 237, 279}, {57, 92, 181, 224, 280}, {57, 76, 138, 207, 264}, {57, 86, 166, 230, 291}, {57, 109, 148, 209, 268}, {57, 80, 133, 225, 271}, {57, 83, 149, 237, 308}, {57, 77, 155, 233, 272}, {57, 119, 146, 252, 314}, {57, 66, 137, 196, 319}, {19, 113, 152, 209, 289}, {58, 87, 134, 250, 296}, {58, 83, 133, 238, 280}, {58, 68, 129, 204, 258}, {58, 107, 187, 213, 261}, {58, 74, 170, 203, 282}, {58, 97, 151, 230, 314}, {58, 125, 191, 239, 290}, {58, 89, 174, 226, 275}, {58, 71, 166, 207, 294}, {58, 119, 162, 250, 279}, {58, 83, 160, 195, 305}, {58, 73, 184, 193, 304}, {58, 67, 156, 248, 294}, {58, 78, 189, 197, 275}, {58, 117, 128, 217, 306}, {19, 127, 164, 233, 263}, {59, 118, 183, 234, 316}, {59, 126, 149, 214, 282}, {59, 103, 174, 236, 277}, {59, 97, 181, 194, 272}, {59, 76, 168, 201, 272}, {59, 94, 187, 214, 267}, {59, 90, 130, 214, 276}, {59, 123, 181, 222, 274}, {59, 119, 163, 248, 288}, {59, 68, 153, 211, 264}, {59, 118, 150, 233, 293}, {59, 74, 189, 196, 288}, {59, 102, 188, 239, 314}, {59, 102, 175, 206, 300}, {59, 82, 142, 203, 288}, {19, 71, 175, 204, 264}, {60, 69, 131, 245, 261}, {60, 82, 184, 212, 291}, {60, 121, 169, 254, 280}, {60, 102, 172, 221, 288}, {60, 80, 147, 211, 314}, {60, 107, 144, 249, 281}, {60, 90, 148, 192, 266}, {60, 77, 169, 240, 305}, {60, 84, 170, 210, 299}, {60, 80, 169, 255, 256}, {60, 90, 163, 225, 315}, {60, 100, 156, 232, 259}, {60, 67, 153, 231, 278}, {60, 85, 173, 197, 269}, {60, 126, 147, 238, 278}, {19, 100, 179, 218, 312}, {61, 77, 147, 233, 312}, {61, 123, 181, 197, 282}, {61, 80, 133, 198, 275}, {61, 116, 178, 244, 270}, {61, 92, 130, 228, 298}, {61, 93, 161, 214, 319}, {61, 80, 167, 196, 285}, {61, 101, 128, 240, 303}, {61, 122, 149, 199, 315}, {61, 81, 152, 202, 319}, {61, 111, 175, 195, 258}, {61, 102, 174, 240, 306}, {61, 98, 148, 219, 271}, {61, 64, 171, 210, 302}, {61, 84, 134, 249, 279}, {19, 98, 141, 222, 283}, {62, 75, 143, 203, 318}, {62, 86, 187, 246, 281}, {62, 105, 182, 195, 290}, {62, 103, 140, 250, 294}, {62, 77, 163, 218, 298}, {62, 98, 129, 197, 261}, {62, 112, 136, 208, 288}, {62, 92, 148, 255, 256}, {62, 66, 189, 240, 266}, {62, 85, 145, 205, 296}, {62, 72, 171, 219, 278}, {62, 74, 139, 223, 260}, {62, 107, 182, 226, 318}, {62, 113, 175, 222, 283}, {62, 117, 185, 202, 308}, {19, 119, 129, 253, 256}, {63, 115, 143, 250, 301}, {63, 93, 176, 236, 258}, {63, 70, 168, 234, 260}, {63, 91, 146, 235, 256}, {63, 116, 143, 208, 295}, {63, 85, 188, 253, 268}, {63, 117, 157, 223, 285}, {63, 75, 132, 205, 285}, {63, 79, 157, 192, 319}, {63, 123, 176, 199, 287}, {63, 95, 175, 219, 293}, {63, 70, 132, 213, 291}, {63, 82, 138, 218, 305}, {63, 64, 152, 251, 284}, {63, 112, 166, 255, 274}, {0, 79, 178, 254, 287}, {0, 81, 139, 201, 282}, {0, 99, 173, 246, 300}, {0, 80, 180, 231, 278}, {0, 96, 160, 200, 286}, {0, 122, 185, 249, 269}, {0, 70, 143, 234, 286}, {0, 83, 145, 239, 271}, {0, 106, 130, 252, 266}, {0, 103, 130, 230, 316}, {0, 79, 163, 255, 259}, {0, 101, 136, 226, 257}, {0, 105, 183, 255, 278}, {0, 102, 162, 239, 315}, {0, 122, 163, 201, 289}, {0, 122, 186, 196, 293}, {1, 66, 181, 240, 295}, {1, 88, 158, 253, 269}, {1, 69, 186, 247, 294}, {1, 68, 136, 200, 298}, {1, 87, 154, 203, 282}, {1, 79, 151, 241, 317}, {1, 80, 164, 216, 294}, {1, 66, 168, 253, 274}, {1, 86, 190, 252, 311}, {1, 123, 183, 228, 298}, {1, 123, 187, 220, 257}, {1, 127, 128, 250, 272}, {1, 78, 146, 225, 298}, {1, 106, 142, 253, 291}, {1, 108, 145, 218, 263}, {1, 71, 156, 234, 299}, {2, 123, 174, 237, 318}, {2, 120, 154, 218, 261}, {2, 123, 154, 234, 303}, {2, 107, 160, 201, 266}, {2, 110, 165, 212, 301}, {2, 85, 160, 248, 303}, {2, 81, 154, 204, 259}, {2, 99, 147, 253, 279}, {2, 86, 185, 207, 279}, {2, 87, 164, 230, 289}, {2, 123, 187, 239, 310}, {2, 99, 169, 205, 310}, {2, 66, 134, 226, 266}, {2, 97, 144, 245, 277}, {2, 106, 157, 201, 270}, {2, 117, 171, 215, 294}, {3, 107, 142, 204, 310}, {3, 124, 191, 250, 302}, {3, 91, 157, 238, 261}, {3, 103, 157, 212, 305}, {3, 118, 132, 234, 273}, {3, 73, 171, 235, 276}, {3, 84, 177, 215, 263}, {3, 108, 154, 227, 258}, {3, 75, 146, 254, 293}, {3, 93, 145, 213, 287}, {3, 120, 173, 206, 289}, {3, 104, 137, 231, 283}, {3, 69, 139, 194, 264}, {3, 74, 189, 208, 260}, {3, 69, 130, 212, 292}, {3, 95, 137, 221, 319}, {4, 96, 148, 208, 275}, {4, 114, 181, 254, 319}, {4, 104, 167, 192, 259}, {4, 126, 166, 212, 290}, {4, 125, 178, 209, 268}, {4, 111, 168, 192, 310}, {4, 64, 143, 207, 292}, {4, 76, 137, 234, 293}, {4, 118, 172, 237, 263}, {4, 103, 170, 220, 267}, {4, 82, 164, 242, 300}, {4, 125, 144, 236, 267}, {4, 110, 177, 247, 292}, {4, 95, 189, 253, 292}, {4, 84, 154, 229, 311}, {4, 115, 143, 224, 303}, {5, 92, 162, 213, 311}, {5, 96, 158, 209, 278}, {5, 94, 152, 200, 290}, {5, 121, 144, 200, 284}, {5, 69, 156, 236, 267}, {5, 74, 163, 211, 318}, {5, 67, 136, 242, 262}, {5, 126, 184, 223, 319}, {5, 100, 143, 242, 264}, {5, 87, 156, 210, 317}, {5, 114, 133, 229, 304}, {5, 99, 129, 196, 299}, {5, 85, 138, 250, 319}, {5, 72, 174, 245, 301}, {5, 92, 143, 227, 278}, {5, 92, 189, 223, 265}, {6, 115, 168, 252, 259}, {6, 106, 142, 230, 312}, {6, 67, 177, 237, 287}, {6, 109, 189, 222, 309}, {6, 85, 178, 219, 318}, {6, 125, 131, 245, 292}, {6, 122, 144, 247, 282}, {6, 100, 162, 195, 296}, {6, 87, 168, 199, 273}, {6, 121, 166, 212, 268}, {6, 81, 167, 238, 265}, {6, 121, 174, 198, 313}, {6, 114, 165, 210, 265}, {6, 81, 129, 218, 270}, {6, 74, 147, 233, 310}, {6, 73, 168, 203, 289}, {7, 108, 178, 220, 275}, {7, 73, 146, 210, 258}, {7, 90, 191, 211, 273}, {7, 100, 156, 254, 277}, {7, 123, 179, 240, 301}, {7, 78, 154, 227, 257}, {7, 113, 150, 194, 283}, {7, 83, 170, 220, 304}, {7, 111, 172, 199, 271}, {7, 72, 134, 220, 297}, {7, 103, 138, 236, 258}, {7, 68, 138, 212, 285}, {7, 78, 149, 250, 319}, {7, 120, 163, 232, 283}, {7, 101, 168, 235, 316}, {7, 105, 191, 199, 294}, {8, 77, 136, 203, 288}, {8, 86, 184, 241, 316}, {8, 122, 173, 204, 319}, {8, 96, 176, 206, 278}, {8, 89, 182, 239, 270}, {8, 106, 132, 241, 311}, {8, 76, 178, 220, 284}, {8, 117, 175, 197, 266}, {8, 104, 151, 226, 282}, {8, 125, 153, 221, 302}, {8, 103, 174, 212, 306}, {8, 92, 169, 205, 280}, {8, 75, 170, 233, 291}, {8, 65, 159, 219, 316}, {8, 107, 158, 205, 289}, {8, 118, 128, 244, 278}, {9, 65, 159, 254, 263}, {9, 118, 191, 211, 258}, {9, 95, 174, 248, 274}, {9, 89, 142, 228, 278}, {9, 100, 138, 197, 256}, {9, 91, 184, 229, 288}, {9, 87, 187, 192, 285}, {9, 88, 153, 246, 318}, {9, 86, 167, 244, 259}, {9, 72, 171, 204, 290}, {9, 111, 149, 221, 297}, {9, 69, 154, 240, 279}, {9, 127, 128, 229, 260}, {9, 125, 152, 249, 266}, {9, 99, 166, 238, 269}, {9, 69, 167, 238, 299}, {10, 125, 185, 209, 268}, {10, 102, 132, 200, 269}, {10, 83, 138, 254, 264}, {10, 106, 139, 243, 272}, {10, 114, 155, 232, 309}, {10, 97, 161, 244, 265}, {10, 102, 152, 214, 265}, {10, 79, 150, 248, 299}, {10, 112, 131, 196, 286}, {10, 91, 171, 231, 297}, {10, 124, 160, 221, 268}, {10, 110, 162, 222, 308}, {10, 96, 177, 239, 285}, {10, 112, 186, 237, 262}, {10, 92, 158, 203, 302}, {10, 94, 139, 195, 281}, {11, 95, 140, 222, 256}, {11, 75, 182, 199, 256}, {11, 78, 154, 236, 290}, {11, 102, 191, 248, 281}, {11, 86, 128, 213, 279}, {11, 73, 147, 251, 296}, {11, 89, 187, 209, 286}, {11, 117, 169, 239, 257}, {11, 124, 137, 252, 299}, {11, 108, 173, 248, 296}, {11, 117, 132, 217, 317}, {11, 67, 178, 227, 310}, {11, 95, 155, 244, 266}, {11, 99, 174, 241, 292}, {11, 102, 170, 196, 302}, {11, 118, 191, 253, 309}, {12, 76, 164, 205, 291}, {12, 108, 142, 225, 312}, {12, 97, 150, 227, 266}, {12, 115, 180, 253, 282}, {12, 111, 176, 226, 278}, {12, 78, 160, 229, 311}, {12, 111, 170, 198, 314}, {12, 90, 158, 242, 319}, {12, 71, 159, 222, 287}, {12, 124, 134, 220, 275}, {12, 107, 148, 221, 302}, {12, 86, 146, 248, 286}, {12, 111, 169, 223, 281}, {12, 87, 131, 193, 314}, {12, 74, 188, 212, 291}, {12, 67, 133, 216, 291}, {13, 105, 191, 235, 264}, {13, 100, 161, 240, 299}, {13, 120, 182, 247, 293}, {13, 119, 186, 216, 276}, {13, 67, 142, 210, 313}, {13, 120, 165, 236, 258}, {13, 88, 173, 199, 298}, {13, 107, 166, 204, 319}, {13, 87, 180, 217, 300}, {13, 72, 169, 246, 284}, {13, 71, 161, 223, 288}, {13, 78, 186, 206, 311}, {13, 96, 172, 227, 294}, {13, 108, 131, 222, 314}, {13, 68, 161, 199, 284}, {13, 91, 133, 232, 258}, {14, 76, 136, 192, 266}, {14, 111, 128, 255, 309}, {14, 72, 130, 234, 308}, {14, 116, 174, 229, 305}, {14, 71, 164, 233, 296}, {14, 65, 140, 253, 297}, {14, 104, 143, 220, 316}, {14, 88, 152, 205, 309}, {14, 99, 179, 202, 283}, {14, 93, 142, 228, 302}, {14, 96, 187, 235, 298}, {14, 78, 136, 245, 280}, {14, 98, 134, 252, 305}, {14, 68, 185, 200, 259}, {14, 64, 161, 228, 286}, {14, 73, 157, 222, 268}, {15, 86, 167, 197, 294}, {15, 126, 129, 223, 302}, {15, 82, 155, 215, 270}, {15, 99, 158, 236, 275}, {15, 114, 162, 246, 298}, {15, 113, 182, 238, 291}, {15, 90, 181, 250, 283}, {15, 102, 184, 214, 270}, {15, 106, 130, 194, 313}, {15, 112, 179, 242, 273}, {15, 109, 175, 233, 313}, {15, 125, 183, 206, 304}, {15, 70, 175, 254, 264}, {15, 97, 176, 246, 290}, {15, 127, 158, 251, 303}, {15, 104, 141, 237, 282}, {16, 120, 170, 249, 295}, {16, 66, 128, 253, 301}, {16, 87, 163, 195, 267}, {16, 95, 143, 239, 271}, {16, 68, 141, 217, 297}, {16, 65, 166, 207, 268}, {16, 119, 140, 221, 261}, {16, 81, 190, 221, 298}, {16, 115, 132, 202, 273}, {16, 121, 129, 234, 273}, {16, 98, 131, 251, 303}, {17, 97, 141, 200, 273}, {18, 114, 185, 205, 302}, {19, 73, 182, 252, 262}, {20, 97, 128, 244, 293}, {20, 89, 135, 204, 296}, {21, 67, 163, 241, 299}, {21, 93, 138, 214, 307}, {22, 109, 170, 211, 300}, {23, 96, 145, 213, 283}, {24, 84, 160, 245, 299}, {25, 91, 165, 192, 260}, {26, 115, 181, 255, 280}, {27, 119, 180, 209, 287}, {28, 76, 151, 229, 301}, {29, 70, 167, 252, 287}, {30, 124, 165, 192, 257}, {31, 82, 184, 251, 315}, {32, 108, 134, 223, 300}, {33, 112, 135, 213, 300}, {34, 64, 146, 248, 301}, {35, 101, 190, 203, 291}, {36, 122, 174, 219, 289}, {37, 66, 171, 238, 259}, {38, 95, 135, 224, 268}, {39, 99, 167, 238, 275}, {40, 74, 142, 194, 308}, {41, 77, 176, 227, 256}, {42, 67, 136, 204, 290}, {43, 71, 155, 228, 296}, {44, 75, 186, 203, 287}, {45, 77, 150, 212, 314}, {46, 69, 141, 221, 316}, {47, 72, 186, 201, 277}, {48, 120, 150, 219, 257}, {49, 126, 176, 250, 310}, {50, 85, 183, 226, 318}, {51, 68, 162, 249, 269}, {52, 75, 182, 255, 316}, {53, 65, 161, 196, 257}, {54, 84, 159, 247, 315}, {55, 109, 167, 200, 306}, {56, 109, 184, 201, 303}, {57, 93, 161, 231, 273}, {58, 97, 178, 194, 298}, {59, 113, 157, 206, 292}, {60, 124, 128, 229, 319}, {61, 115, 136, 201, 282}, {62, 89, 162, 252, 304}, {63, 104, 190, 203, 301}};
