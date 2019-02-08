#include "mpi_functions.h"

const int Blur[S+2][S+2] = {
    {1,1,1},
    {1,1,1},
    {1,1,1}
};

const int Sharpen[S+2][S+2] = {
    {0,-1,0},
    {-1,5,-1},
    {0,-1,0}
};

const int Edge[S+2][S+2] = {
    {0,-1,0},
    {-1,4,-1},
    {0,-1,0}
};

const int EdgeH[S+2][S+2] = {
    {0,0,0},
    {-1,2,-1},
    {0,0,0}
};

const int EdgeV[S+2][S+2] = {
    {0,-1,0},
    {0,2,0},
    {0,-1,0}
};

const int GradientH[S+2][S+2] = {
    {-1,-1,-1},
    {0,0,0},
    {1,1,1}
};

const int GradientV[S+2][S+2] = {
    {-1,0,1},
    {-1,0,1},
    {-1,0,1}
};

const int SobelH[S+2][S+2] = {
    {1,2,1},
    {0,0,0},
    {-1,-2,-1}
};

const int SobelV[S+2][S+2] = {
    {1,0,-1},
    {2,0,-2},
    {1,0,-1}
};

const int Emboss[S+2][S+2] = {
    {-2,-1,0},
    {-1,1,1},
    {0,1,2}
};
