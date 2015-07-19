
#include <stdlib.h>
#include <math.h>

int getnum(int decaynum){ // translate decay hashes to the internal number scheme
    int res = 0;
    switch(decaynum){
        case -522186: res = -17; break;
        case -522167: res = -16; break;
        case -522056: res = -15; break;
        case -522054: res = -14; break;
        case -521964: res = -13; break;
        case -521875: res = -12; break;
        case -521856: res = -11; break;
        case -521855: res = -10; break;
        case -521854: res = -9 ; break;
        case -521852: res = -8 ; break;
        case -521833: res = -7 ; break;
        case -521764: res = -6 ; break;
        case -521745: res = -5 ; break;
        case -521743: res = -4 ; break;
        case -521742: res = -3 ; break;
        case -521634: res = -2 ; break;
        case -521632: res = -1 ; break;
        case  521632: res = 1 ; break;
        case  521634: res = 2 ; break;
        case  521742: res = 3 ; break;
        case  521743: res = 4 ; break;
        case  521745: res = 5 ; break;
        case  521764: res = 6 ; break;
        case  521833: res = 7 ; break;
        case  521852: res = 8 ; break;
        case  521854: res = 9 ; break;
        case  521855: res = 10; break;
        case  521856: res = 11; break;
        case  521875: res = 12; break;
        case  521964: res = 13; break;
        case  522054: res = 14; break;
        case  522056: res = 15; break;
        case  522167: res = 16; break;
        case  522186: res = 17; break;
            
        case -512157: res = -13; break;
        case -512046: res = -12; break;
        case -512044: res = -11; break;
        case -511975: res = -10; break;
        case -511846: res = -9; break;
        case -511845: res = -8; break;
        case -511844: res = -7; break;
        case -511842: res = -6; break;
        case -511735: res = -5; break;
        case -511733: res = -4; break;
        case -511624: res = -3; break;
        case -511622: res = -2; break;
        case -511532: res = -1; break;
        case 511532: res = 1; break;
        case 511622: res = 2; break;
        case 511624: res = 3; break;
        case 511733: res = 4; break;
        case 511735: res = 5; break;
        case 511842: res = 6; break;
        case 511844: res = 7; break;
        case 511845: res = 8; break;
        case 511846: res = 9; break;
        case 511975: res = 10; break;
        case 512044: res = 11; break;
        case 512046: res = 12; break;
        case 512157: res = 13; break;
            
        case 512175: res = 14; break;
        case 511753: res = 15; break;
            
        default: break;
    }
    return res;
}

double pol4(double x, const double *p){ // 4th order polynomial evaluation
    return p[0] + x*(p[1] + x*(p[2] + x*(p[3] + x*(p[4]))));
}

double tagcorr(int decaynum, double nbout){
    // input variables:
    //
    // decaynum is tag decay mode (Ekpfullrecon::decay() method output)
    // nbout is neural net output variable (Ekpfullrecon::NBout() method output)
    //
    // return value:
    //
    // Multiplicative MC event correction. Note that calibration was
    // made under assumption that Upsilon(4S) decays into B+B-(B0B0)
    // with BR = 51.6%(48.4%)
    
    static const double _pbc[][5] = {
        {0.803102, -0.116015, 0.014123, 0.006996, 0.000620},
        {0.819491, 0.029608, 0.054351, 0.010438, 0.000779},
        {0.887012, 0.310010, 0.103320, 0.007104, 0.000106},
        {0.626101, -0.143045, 0.007881, 0.005269, 0.000418},
        {0.599316, -0.178898, 0.009490, 0.007289, 0.000542},
        {0.900000, 0.000000, 0.000000, 0.000000, 0.000000},
        {0.237321, -0.495781, -0.045058, 0.004273, 0.000567},
        {0.657237, -0.373295, -0.103654, -0.012458, -0.000557},
        {0.592405, -0.261191, -0.050523, -0.005440, -0.000249},
        {1.063210, 0.445836, 0.099481, 0.001242, -0.000400},
        {0.644721, -0.120312, -0.031039, -0.004366, -0.000086},
        {1.175940, -0.187520, -0.199063, -0.033703, -0.001613},
        {0.349156, -0.768838, -0.052804, 0.009078, 0.000654},
        {0.557068, -0.238081, -0.019655, 0.003448, 0.000380},
        {0.706730, -0.029007, 0.025267, 0.006106, 0.000442},
        {0.684215, 0.037921, 0.034631, 0.004792, 0.000242},
        {0.778064, -0.184310, -0.007253, 0.002446, 0.000219}
    };
    
    static const double _pbci[] = {1.057080, 0.918807, 0.978150, 0.733643,
        0.727512, 0.952388, 0.485972, 0.671686,
        0.984174, 0.812720, 0.928933, 1.094490,
        0.914915, 0.672623, 0.563267, 0.516600, 0.609034};
    
    static const double _pb0[][5] = {
        {1.693560, 1.135230, 0.364091, 0.037053, 0.001183},
        {0.797147, -0.180640, -0.006953, 0.007928, 0.000889},
        {0.849117, -0.293022, -0.122973, -0.021244, -0.001328},
        {0.843642, 0.037873, 0.007734, -0.004009, -0.000410},
        {0.861456, -0.084080, -0.001069, 0.003512, 0.000337},
        {0.145721, -1.166530, -0.311187, -0.026577, -0.000623},
        {0.447494, -0.360701, -0.092223, -0.009766, -0.000239},
        {1.335940, 0.842397, 0.245318, 0.022089, 0.000612},
        {0.721848, -0.253327, -0.026691, 0.001022, 0.000089},
        {1.634700, 1.559520, 0.585662, 0.097255, 0.006251},
        {0.702499, -0.147479, -0.001137, 0.005302, 0.000470},
        {0.675704, -0.290131, -0.069624, -0.008028, -0.000363},
        {0.792166, 0.018745, 0.022773, 0.002513, 0.000104},
        {0.552737, -0.351686, -0.174651, -0.027102, -0.001183},
        {1.597510, 2.063240, 2.005750, 0.593154, 0.053781}
    };
    
    static const double _pb0i[] = {0.896533, 0.930940, 0.957721, 0.905519,
        1.088230, 0.952349, 0.805700, 1.018510,
        0.761165, 1.032340, 0.554638, 0.523624,
        0.402844, 0.380570, 0.739893};
    
    double kweight = 1;
    const int num = getnum(decaynum);
    if(!num) return kweight;
    
    const bool b0 = (abs(decaynum)/10000) == 51;
    
    if(!b0){ // calibration was made under assumption that Upsilon(4S)
        // decays into B+B-(B0B0) with BR = 51.6%(48.4%)
        kweight *= 0.516/0.5;
    } else {
        kweight *= 0.484/0.5;
    }
    
    const int i = abs(num) - 1;
    const double x = floor(2*log(nbout))/2;
    
    double nbout_corr;
    if(!b0){
        nbout_corr = _pbci[i]*pol4(x,_pbc[i]);
    } else {
        nbout_corr = _pb0i[i]*pol4(x,_pb0[i]);
    }
    if(nbout_corr<0) nbout_corr = 0.5; // due to polynom oscillation -- keep correction always positive
    
    kweight *= nbout_corr;
    
    return kweight;
}
