#include "MRCB.h"
int main()
{  
    char file[1024] = "";
    prcopt_t prcopt = prcopt_default;   /* processing options */
    solopt_t solopt = solopt_default;   /* solution options   */
    filopt_t filopt = { "" };
    int config_loaded = 0;
#ifdef _WIN32
    sprintf(file, "..\\conf\\config.dat");
#else
    sprintf(file, "../conf/config.dat");
#endif

    resetsysopts();
    if (loadopts(file, sysopts)) {
        config_loaded = 1;
    }

    if (!config_loaded) {
#ifdef _WIN32
        sprintf(file, "..\\..\\conf\\config.dat");
#else
        sprintf(file, "../../conf/config.dat");
#endif
        resetsysopts(); 
        if (loadopts(file, sysopts)) {
            config_loaded = 1;
        }
    }

    if (!config_loaded) {
        printf("Error: Cannot find config file in both paths:\n");
#ifdef _WIN32
        printf("1. ..\\conf\\config.dat\n");
        printf("2. ..\\..\\conf\\config.dat\n");
#else
        printf("1. ../conf/config.dat\n");
        printf("2. ../../conf/config.dat\n");
#endif
        printf("Recommendation: Use an absolute path in main.c to ensure reliable reading of config.dat.\n");
        return 0;
    }
    getsysopts(&prcopt, &solopt, &filopt);
    postpos(&prcopt, &solopt, filopt);
    return 0;
}