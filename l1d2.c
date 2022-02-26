/*
  l1d2.c ... generates scaling data for calculating the largest
             Lyapunov exponent and the correlation dimension

             Copyright (c) 1999. Michael T. Rosenstein.

This program is free software; you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation; either version 2 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program; if not, write to the Free Software Foundation, Inc., 59 Temple
Place - Suite 330, Boston, MA 02111-1307, USA.

You may contact the author by e-mail (mtr@cs.umass.edu) or postal mail
(c/o Department of Computer Science, University of Massachusetts, 140 Governors
Drive, Amherst, MA 01003-4610 USA).  For updates to this software, please visit
PhysioNet (http://www.physionet.org/).

	     reference: M.T. Rosenstein, J.J. Collins, C.J. De Luca,
                        A practical method for calculating largest
                        Lyapunov exponents from small data sets,
                        Physica D 65:117-134, 1993.

             email contact: mtr@cs.umass.edu	     
*/

#include <ctype.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define IOTA			10e-15
#define MAX_LN_R		12
#define MIN_LN_R		-12
#define N_LN_R			600

typedef struct
{
    char  fileName[16];
    long  startIndex, stopIndex;
    int   seriesN, m, J, W, divergeT;
} test;
	
double	**AllocateDMatrix(long nRows, long nCols);
void	ComputeSlopes(void);
void	FreeDMatrix(double **mat, long nRows);
void	GenerateTemplateFile(void);
long	GetData(char *fileName, int tsN, long start, long stop);
void	LineFit(double *data, double n, double *m, double *b, double *rr);
void	PercentDone(int percentDone);
void	ProcessTest(int testN);
void	ReadTestFile(char *fileName);
void	SaveL1Results(char *fileRoot);
void	SaveD2Results(char *fileRoot);

int     gNTests, gMaxDivergeT;
double  *gData, *gNDivergence;
double  **gDivergence, **gCSum;
test    *gTest;

int    main(void)
{
    int   i;
    char  str[256];

    printf("\n");
    printf("*** L1D2: generates scaling information for L1 and D2 ***\n");
    printf("          v1.0 Copyright (c) 1999 M.T. Rosenstein\n");
    printf("                                  (mtr@cs.umass.edu)\n\n");
    printf("          reference: M.T. Rosenstein, J.J. Collins, C.J. De Luca,\n");
    printf("                     A practical method for calculating largest\n");
    printf("                     Lyapunov exponents from small data sets,\n");
    printf("                     Physica D 65:117-134, 1993.\n\n");

    GenerateTemplateFile();

    printf("Enter test file name: ");
    scanf("%s", str);
    ReadTestFile(str);

    printf("\nEnter output file root (no extension): ");
    scanf("%s", str);

    /* allocate the divergence and correlation sum arrays */
    for (i = 0; i < gNTests; i++)
        if (gTest[i].divergeT > gMaxDivergeT)
            gMaxDivergeT = 1 + gTest[i].divergeT;
    gNDivergence = (double*)calloc(gMaxDivergeT, sizeof(double));
    gDivergence = AllocateDMatrix(gMaxDivergeT, 2 * gNTests);
    gCSum = AllocateDMatrix(N_LN_R, 2 * gNTests);

    for (i = 0; i < gNTests; i++)
        ProcessTest(i);

    ComputeSlopes();
    printf("\n");
    SaveL1Results(str);
    SaveD2Results(str);

    printf("\nSuccess!\n");
    return(0);
}

double** AllocateDMatrix(long nRows, long nCols)
{
    long   i;
    double** mat;

    /* allocate space for the row pointers */
    mat = (double**)calloc(nRows, sizeof(double*));
    if (mat == NULL)
    {
        printf("OUT OF MEMORY: AllocateDMatrix(%ld, %ld)\n\n", nRows, nCols);
        exit(1);
    };

    /* allocate space for each row pointer */
    for (i = 0; i < nRows; i++)
        mat[i] = (double*)calloc(nCols, sizeof(double));
    if (mat[i - 1] == NULL)
    {
        printf("OUT OF MEMORY: AllocateDMatrix(%ld, %ld)\n\n", nRows, nCols);
        exit(1);
    };

    return(mat);
}

void   ComputeSlopes(void)
{
    int    i, i2, j;
    double k, m, b, rr;
    double* data;

    data = (double*)calloc(N_LN_R > gMaxDivergeT ? N_LN_R : gMaxDivergeT,
        sizeof(double));
    if (data == NULL)
    {
        printf("OUT OF MEMORY: ComputeSlopes\n\n");
        exit(1);
    };

    for (i = 0; i < gNTests; i++)
    {
        i2 = i + gNTests;

        /*** work on correlation dimension first ***/

        k = (double)N_LN_R / (MAX_LN_R - MIN_LN_R);

        /* pack the array column into the local data array */
        for (j = 0; j < N_LN_R; j++)
            data[j] = gCSum[j][i];
        /* compute the 7-point slopes */
        for (j = 3; j < N_LN_R - 3; j++)
        {
            LineFit(data + j - 3, 7, &m, &b, &rr);
            gCSum[j][i2] = k * m;
        };
        /* handle the edges */
        LineFit(data, 5, &m, &b, &rr); gCSum[2][i2] = k * m;
        LineFit(data + N_LN_R - 5, 5, &m, &b, &rr); gCSum[N_LN_R - 3][i2] = k * m;
        LineFit(data, 3, &m, &b, &rr); gCSum[1][i2] = k * m;
        LineFit(data + N_LN_R - 3, 3, &m, &b, &rr); gCSum[N_LN_R - 2][i2] = k * m;
        gCSum[0][i2] = k * (data[1] - data[0]);
        gCSum[N_LN_R - 1][i2] = k * (data[N_LN_R - 1] - data[N_LN_R - 2]);

        /*** now work on divergence data ***/

        /* pack the array column into the local data array */
        for (j = 0; j < gMaxDivergeT; j++)
            data[j] = gDivergence[j][i];
        /* compute the 7-point slopes */
        for (j = 3; j < gMaxDivergeT - 3; j++)
        {
            LineFit(data + j - 3, 7, &m, &b, &rr);
            gDivergence[j][i2] = m;
        };
        /* handle the edges */
        LineFit(data, 5, &m, &b, &rr); gDivergence[2][i2] = m;
        LineFit(data + gMaxDivergeT - 5, 5, &m, &b, &rr);
        gDivergence[gMaxDivergeT - 3][i2] = m;
        LineFit(data, 3, &m, &b, &rr);
        gDivergence[1][i2] = m;
        LineFit(data + gMaxDivergeT - 3, 3, &m, &b, &rr);
        gDivergence[gMaxDivergeT - 2][i2] = m;
        gDivergence[0][i2] = data[1] - data[0];
        gDivergence[gMaxDivergeT - 1][i2] = data[gMaxDivergeT - 1] -
            data[gMaxDivergeT - 2];
    };
}

void   FreeDMatrix(double** mat, long nRows)
{
    long   i;

    /* free space for each row pointer */
    for (i = nRows - 1; i >= 0; i--)
        free(mat[i]);

    /* free space for the row pointers */
    free(mat);
}

void   GenerateTemplateFile(void)
{
    FILE* outFile;

    outFile = fopen("sample.l1d2", "w");

    fprintf(outFile, "* Header info starts with an asterisk.\n*\n");
    fprintf(outFile, "* Each line of the test file contains the name of a data file followed\n");
    fprintf(outFile, "* by the parameters for the test:\n");
    fprintf(outFile, "*   file_name series# startIndex stopIndex m J W divergeT\n*\n");
    fprintf(outFile, "*      file_name  = name of the data file\n");
    fprintf(outFile, "*      series#    = time series number to use for delay reconstruction\n");
    fprintf(outFile, "*      startIndex = index of first data point to read (usually 1)\n");
    fprintf(outFile, "*      stopIndex  = index of last data point to read\n");
    fprintf(outFile, "*                      (enter 0 for maximum)\n");
    fprintf(outFile, "*      m          = embedding dimension\n");
    fprintf(outFile, "*      J          = delay in samples\n");
    fprintf(outFile, "*      W          = window size for skipping temporally close nearest neighbors\n");
    fprintf(outFile, "*      divergT    = total divergence time in samples\n");
    fprintf(outFile, "*   example: lorenz.dat 1 1 0 5 7 100 300\n");

    fclose(outFile);
}

long   GetData(char* fileName, int tsN, long start, long stop)
{
    long   i, j, len;
    long   nHead, nCols, nRows, nPts;
    char   str[1024];
    double dummy;
    FILE* inFile;

    /* try to open the data file */
    inFile = fopen(fileName, "r");
    if (inFile == NULL)
    {
        printf("FILE ERROR: GetData(%s)\n\n", fileName);
        exit(1);
    };

    /* skip the header */
    nHead = 0;
    fgets(str, 1024, inFile);
    while (str[0] == '*')
    {
        nHead++;
        fgets(str, 1024, inFile);
    };

    /* figure out the number of columns */
    len = strlen(str);
    i = 0;
    nCols = 0;
    while (i < len)
    {
        if (!isspace(str[i]))
        {
            nCols++;
            while (!isspace(str[i]) && i < len)
                i++;
            while (isspace(str[i]) && i < len)
                i++;
        }
        else
            i++;
    };

    /* figure out the number of rows; assume there's at least one */
    nRows = 1;
    while (fgets(str, 1024, inFile) != NULL)
        nRows++;

    if (stop < start)
        stop = nRows;
    else if (stop > nRows)
        stop = nRows;
    if (start<1 || start>stop - 3)
        start = 1;
    nPts = stop - start + 1;
    gData = calloc(nPts, sizeof(double));

    /* now read the time series data */
    rewind(inFile);
    for (i = 0; i < nHead; i++)
        fgets(str, 1024, inFile);

    for (i = 1; i < start; i++)
        for (j = 0; j < nCols; j++)
            fscanf(inFile, "%lf", &dummy);
    for (i = 0; i < nPts; i++)
    {
        for (j = 0; j < tsN; j++)
            fscanf(inFile, "%lf", &dummy);
        gData[i] = dummy;
        for (; j < nCols; j++)
            fscanf(inFile, "%lf", &dummy);
    };
    fclose(inFile);

    return(nPts);
}

void   LineFit(double* data, double n, double* m, double* b, double* rr)
{
    int    i;
    double sx, sy, sxy, sx2, sy2;
    double x, y, k, mTemp, bTemp, rrTemp;

    sx = sy = sxy = sx2 = sy2 = 0;
    for (i = 0; i < n; i++)
    {
        x = i;
        y = data[i];
        sx += x; sy += y;
        sx2 += x * x; sy2 += y * y;
        sxy += x * y;
    };
    k = n * sx2 - sx * sx;
    mTemp = (n * sxy - sx * sy) / k;
    bTemp = (sx2 * sy - sx * sxy) / k;
    k = sy * sy / n;
    if (k == sy2)
        rrTemp = 1.0;
    else
    {
        rrTemp = (bTemp * sy + mTemp * sxy - k) / (sy2 - k);
        rrTemp = 1.0 - (1.0 - rrTemp) * (n - 1.0) / (n - 2.0);
    };
    *m = mTemp;
    *b = bTemp;
    *rr = rrTemp;
}

void   PercentDone(int percentDone)
{
    static last = 100;

    if (percentDone < last)
    {
        last = 0;
        printf("0");
        fflush(stdout);
    }
    else if (percentDone > last && percentDone % 2 == 0)
    {
        last = percentDone;
        if (percentDone % 10 == 0)
            printf("%d", percentDone / 10);
        else
            printf(".");
        fflush(stdout);
    };
}

void   ProcessTest(int testN)
{
    long   m, J, W, divergeT, neighborIndex, maxIndex;
    long   i, j, k, T, CSumIndex, percentDone;
    long   nPts, nCompletedPairs = 0, nVectors;
    char* isNeighbor;
    double distance, d;
    double k1, k2, temp, kNorm;

    printf("\nProcessing test %d of %d:  ", testN + 1, gNTests);
    fflush(stdout);

    m = gTest[testN].m;
    J = gTest[testN].J;
    W = gTest[testN].W;
    divergeT = gTest[testN].divergeT;
    nPts = GetData(gTest[testN].fileName, gTest[testN].seriesN,
        gTest[testN].startIndex, gTest[testN].stopIndex);

    k1 = (double)N_LN_R / (MAX_LN_R - MIN_LN_R);
    k1 *= 0.5; /* accounts for the SQUARED distances later on */
    k2 = N_LN_R / 2;

    nVectors = nPts - J * (m - 1);

    isNeighbor = (char*)calloc(nVectors, sizeof(char));
    if (isNeighbor == NULL)
    {
        printf("\nOUT OF MEMORY: ProcessTest\n\n");
        exit(1);
    };

    /* clear the divergence arrays */
    for (i = 0; i < gMaxDivergeT; i++)
        gNDivergence[i] = gDivergence[i][testN] = 0;

    /* loop through vectors */
    i = 0;
    while (i < nVectors)
    {
        percentDone = 100.0 * nCompletedPairs / nVectors * 2 + 0.5;
        percentDone = 100.0 * i / nVectors + 0.5;
        PercentDone(percentDone);

        if (!isNeighbor[i])
        {
            distance = 10e10;

            /* find the nearest neighbor for the vector at i */
            /* ignore temporally close neighbors using W */
            if (i > W)
                for (j = 0; j < i - W; j++)
                {
                    /* calculate distance squared */
                    d = 0;
                    for (k = 0; k < m; k++)
                    {
                        T = k * J;
                        temp = gData[i + T] - gData[j + T];
                        temp *= temp;
                        d += temp;
                    };
                    d += IOTA;

                    /* map the squared distance to array position */
                    CSumIndex = k1 * log(d) + k2;
                    if (CSumIndex < 0)
                        CSumIndex = 0;
                    if (CSumIndex >= N_LN_R)
                        CSumIndex = N_LN_R - 1;

                    /* increment the correlation sum array */
                    gCSum[CSumIndex][testN]++;

                    /* now compare to current nearest neighbor */
                    /* use IOTA above to ignore nbrs that are too close */
                    if (d < distance)
                    {
                        distance = d;
                        neighborIndex = j;
                    };
                };
            if (i < nVectors - W)
                for (j = i + W; j < nVectors; j++)
                {
                    d = 0;
                    for (k = 0; k < m; k++)
                    {
                        T = k * J;
                        temp = gData[i + T] - gData[j + T];
                        temp *= temp;
                        d += temp;
                    };
                    d += IOTA;

                    CSumIndex = k1 * log(d) + k2;
                    if (CSumIndex < 0)
                        CSumIndex = 0;
                    if (CSumIndex >= N_LN_R)
                        CSumIndex = N_LN_R - 1;

                    gCSum[CSumIndex][testN]++;

                    if (d < distance)
                    {
                        distance = d;
                        neighborIndex = j;
                    };
                };
            isNeighbor[neighborIndex] = 1;

            /* track divergence */
            for (j = 0; j <= divergeT; j++)
            {
                maxIndex = nPts - m * J - j - 1;
                if (i < maxIndex && neighborIndex < maxIndex)
                {
                    d = 0;
                    for (k = 0; k < m; k++)
                    {
                        T = k * J + j;
                        temp = gData[i + T] - gData[neighborIndex + T];
                        temp *= temp;
                        d += temp;
                    };
                    d += IOTA;
                    gNDivergence[j]++;
                    temp = 0.5 * log(d);
                    gDivergence[j][testN] += temp;
                };
            };
            nCompletedPairs++;
        };
        i++;
    };

    /* integrate the correlation sum array to get the correlation sum curve */
    for (i = 1; i < N_LN_R; i++)
        gCSum[i][testN] += gCSum[i - 1][testN];

    /* next normalize values */
    kNorm = 1.0 / gCSum[N_LN_R - 1][testN];
    for (i = 0; i < N_LN_R; i++)
        gCSum[i][testN] *= kNorm;

    /* now take natural log of the values */
    for (i = 0; i < N_LN_R; i++)
    {
        temp = gCSum[i][testN];
        if ((temp < 0.000045) || (temp > 0.990050))
            gCSum[i][testN] = 0;
        else
            gCSum[i][testN] = log(temp);
    };

    /* now take care of Lyapunovv average */
    for (i = 0; i <= divergeT; i++)
        if (gNDivergence[i] > 0)
            gDivergence[i][testN] /= gNDivergence[i];

    free(isNeighbor);
    free(gData);
}

void   ReadTestFile(char* fileName)
{
    int    i;
    int    nHead, nRows;
    char   str[1024];
    FILE* inFile;

    printf("\nReading Test File...\n");

    /* try to open the data file */
    inFile = fopen(fileName, "r");
    if (inFile == NULL)
    {
        printf("FILE ERROR: ReadTestFile(%s)\n\n", fileName);
        exit(1);
    };

    /* skip the header */
    nHead = 0;
    fgets(str, 1024, inFile);
    while (str[0] == '*')
    {
        nHead++;
        fgets(str, 1024, inFile);
    };

    /* figure out the number of rows; assume there's at least one */
    nRows = 1;
    while (fgets(str, 1024, inFile) != NULL && !isspace(str[0]))
        nRows++;
    gNTests = nRows;

    /* allocate the test array */
    gTest = (test*)calloc(gNTests, sizeof(test));
    if (gTest == NULL)
    {
        printf("OUT OF MEMORY: ReadTestFile(%d)\n\n", gNTests);
        exit(1);
    };

    printf("detected %d %s\n", gNTests, gNTests == 1 ? "test" : "tests");

    /* rewind the file and skip the header */
    rewind(inFile);
    for (i = 0; i < nHead; i++)
        fgets(str, 1024, inFile);

    for (i = 0; i < gNTests; i++)
        fscanf(inFile, "%s %d %ld %ld %d %d %d %d\n",
            gTest[i].fileName, &gTest[i].seriesN, &gTest[i].startIndex,
            &gTest[i].stopIndex, &gTest[i].m, &gTest[i].J, &gTest[i].W,
            &gTest[i].divergeT);

    fclose(inFile);
}

void   SaveD2Results(char* fileRoot)
{
    int    i, i1, i2, testN, keepGoing;
    char   str[256];
    double k1, k2;
    FILE* outFile;

    printf("\nSaving data for correlation dimension...\n");

    sprintf(str, "%s.d2", fileRoot);
    outFile = fopen(str, "w");

    k1 = (double)(MAX_LN_R - MIN_LN_R) / N_LN_R;
    k2 = MIN_LN_R;

    /* don't save rows of just zeros */
    keepGoing = 1;
    i1 = 0;
    while (keepGoing)
    {
        for (testN = 0; testN < gNTests; testN++)
            if (gCSum[i1][testN] != 0)
            {
                keepGoing = 0;
                break;
            };
        i1 += keepGoing;
    };
    i1--;
    if (i1 < 0 || i1 >= N_LN_R)
        i1 = 0;
    keepGoing = 1;
    i2 = N_LN_R - 1;
    while (keepGoing)
    {
        for (testN = 0; testN < gNTests; testN++)
            if (gCSum[i2][testN] != 0)
            {
                keepGoing = 0;
                break;
            };
        i2 -= keepGoing;
    };
    i2++;
    if (i2 < 0 || i2 >= N_LN_R)
        i2 = N_LN_R - 1;

    /* write the data */
    for (i = i1; i < i2 + 1; i++)
    {
        fprintf(outFile, "%lf\t", k1 * i + k2);
        for (testN = 0; testN < gNTests; testN++)
            fprintf(outFile, "%lf\t", gCSum[i][testN]);

        /* write slope data */
        fprintf(outFile, "%lf\t", k1 * i + k2);
        for (; testN < 2 * gNTests - 1; testN++)
            fprintf(outFile, "%lf\t", gCSum[i][testN]);
        fprintf(outFile, "%lf\n", gCSum[i][testN]);
    };

    fclose(outFile);
}

void    SaveL1Results(char* fileRoot)
{
    int   i, testN;
    char  str[256];
    FILE* outFile;

    printf("\nSaving data for largest Lyapunov exponent...\n");

    sprintf(str, "%s.l1", fileRoot);
    outFile = fopen(str, "w");

    for (i = 0; i < gMaxDivergeT; i++)
    {
        fprintf(outFile, "%d\t", i);
        for (testN = 0; testN < gNTests; testN++)
            if (i <= gTest[testN].divergeT)
                fprintf(outFile, "%lf\t", gDivergence[i][testN]);
            else
                fprintf(outFile, "\t");

        /* write slope data */
        fprintf(outFile, "%d\t", i);
        for (; testN < 2 * gNTests - 1; testN++)
            if (i <= gTest[testN - gNTests].divergeT)
                fprintf(outFile, "%lf\t", gDivergence[i][testN]);
            else
                fprintf(outFile, "\t");
        if (i <= gTest[testN - gNTests].divergeT)
            fprintf(outFile, "%lf\n", gDivergence[i][testN]);
        else
            fprintf(outFile, "\n");
    };

    fclose(outFile);
}




