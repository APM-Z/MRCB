#include "MRCB.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

// Platform-specific headers
#ifdef _WIN32
#include <windows.h>
#include <io.h>
#else
#include <dirent.h>
#include <unistd.h>
#endif

#define MAXFILESTR 128

extern int initIfile(Ifile* ifile)
{
    int i = 0;
    for (i = 0; i < NIfile; i++) ifile->idx[i] = 0;
    // Observation files
    for (i = 0; i < MAXOBSFILE; i++) {
        if (!(ifile->obsfile[i] = (char*)malloc(MAXFILESTR))) {
            for (; i >= 0; i--) free(ifile->obsfile[i]);
            return 0;
        }
    }
    // Navigation files
    for (i = 0; i < MAXPROFILE; i++) {
        if (!(ifile->rnxfile[i] = (char*)malloc(MAXFILESTR))) {
            for (; i >= 0; i--) free(ifile->rnxfile[i]);
            return 0;
        }
    }
    return 1;
}

extern int freeIfile(Ifile* ifile)
{
    int i = 0;

    for (i = 0; i < MAXOBSFILE; i++) {
        free(ifile->obsfile[i]);
        free(ifile->rnxfile[i]);
    }
    return 1;
}

static int dotofstr(char* str, int n)
{
    int i = 0;
    int pos = 0;

    for (i = 0; i < n; i++) {
        if (str[i] == '.')
        {
            pos = i;
        }
    }

    return pos ? pos : i;
}

// Unified directory scanning function
static int scan_directory(const char* dir_path, Ifile* ifile, int file_type)
{
    char fname[1024] = { 0 };
    int dot = 0;

#ifdef _WIN32
    WIN32_FIND_DATAA findFileData;
    HANDLE hFind;
    char searchPath[MAX_PATH];

    // Build search path: dir_path\*
    snprintf(searchPath, sizeof(searchPath), "%s\\*", dir_path);
    hFind = FindFirstFileA(searchPath, &findFileData);  
    if (hFind == INVALID_HANDLE_VALUE) {
        printf("open dir failed: %s\n", dir_path);
        return 0;
    }

    do {
        // Skip "." and ".." directories
        if (strcmp(findFileData.cFileName, ".") == 0 ||
            strcmp(findFileData.cFileName, "..") == 0) {
            continue;
        }
        strcpy(fname, findFileData.cFileName);
#else
    DIR* dir = NULL;
    struct dirent* entry;

    if ((dir = opendir(dir_path)) == NULL) {
        printf("open dir failed: %s\n", dir_path);
        return -1;
    }

    while ((entry = readdir(dir))) {
        strcpy(fname, entry->d_name);
#endif

        dot = dotofstr(fname, sizeof(fname));
        if (dot >= 1024) continue;

        // Observation files (.o or .O)
        if (file_type == 0) {
            if (strstr(fname + dot, "o") != NULL || strstr(fname + dot, "O") != NULL) {
                snprintf(ifile->obsfile[ifile->idx[0]], MAXFILESTR, "%s%s", dir_path, fname);
                ifile->idx[0]++;
            }
        }
        // Navigation files (.N, .n, .P, .p)
        else if (file_type == 1) {
            if (strstr(fname + dot, "N") != NULL || strstr(fname + dot, "n") != NULL ||
                strstr(fname + dot, "P") != NULL || strstr(fname + dot, "p") != NULL) {
                snprintf(ifile->rnxfile[ifile->idx[1]], MAXFILESTR, "%s%s", dir_path, fname);
                ifile->idx[1]++;
            }
        }

#ifdef _WIN32
    } while (FindNextFileA(hFind, &findFileData)); 
    FindClose(hFind);
#else
    }
    closedir(dir);
#endif

    return 1;
}

extern int getpath(Ifile* ifile, prcopt_t* opt, filopt_t filopt)
{
    char obs_dir[120] = "", nav_dir[120] = "";
    char obs[120] = "", nav[120] = "";

    if (opt->CONTINUE == 1) {
        strcpy(obs_dir, filopt.ObsDir);
        strcpy(nav_dir, filopt.NavDir);

        // Scan observation files directory
        if (scan_directory(obs_dir, ifile, 0) < 0) {
            printf("open obs dir failed!");
            return 0;
        }

        // Scan navigation files directory
        if (scan_directory(nav_dir, ifile, 1) < 0) {
            printf("open nav dir failed!");
            return 0;
        }
    }
    else {
        strcpy(obs, filopt.obs);
        strcpy(nav, filopt.nav);
        if (strstr(obs, "o") != NULL || strstr(obs, "O") != NULL) {
            snprintf(ifile->obsfile[ifile->idx[0]], MAXFILESTR, "%s", obs);
            ifile->idx[0]++;
        }
        if (strstr(nav, "N") != NULL || strstr(nav, "n") != NULL ||
            strstr(nav, "P") != NULL || strstr(nav, "p") != NULL) {
            snprintf(ifile->rnxfile[ifile->idx[1]], MAXFILESTR, "%s", nav);
            ifile->idx[1]++;
        }
    }

    /* Check if files are found */
    if (ifile->idx[0] == 0) {
        printf("lack obs file\n");
        return 0;
    }
    else if (ifile->idx[1] == 0) {
        printf("lack nav file\n");
        return 0;
    }
    return 1;
}