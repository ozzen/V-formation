#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <float.h>
#include "pso.h"

#define PI 3.14159265358979323846

const int Num = 7;
// Number of time stpes
const int Steps = 50;
// Initial range

// Initial configs

const double Init_box = 5.0; 

const double InitVmin = 0.25;
const double InitVmax = 0.75;

// Minimum distance for collision freedom
const double Dmin = 1;

int Ph = 1;
const int MaxPh = 5;
// wing span
const double w = 1.0;
// y_opt for upwash
const double d0 = 1.0;
// x_opt is 2w-lambda 
const double lambda = 0.5 - PI / 8.0;
// angle of clear view cone
const double angle = PI / 6.0;

// Gaussian params for upwash
const double u_sigma1 = 5.0;
const double u_sigma2 = 5.0;
const double d_sigma1 = 1.0 / 0.3;
const double d_sigma2 = 1.0 / 0.7;

// bound on acceleration w.r.t velocity
const double delta = 1.0;

int currentConf;
int currentRun;
const int Clones = 5;

static inline double randd(double min, double max) {
    return ((double) rand() / (double) RAND_MAX) * (max - min) + min;
}

static inline double norm(double x, double y) {
    return sqrt(x * x + y * y);
}


static inline double mvnpdf(double x, double y, double a, double b) {
    return exp(-0.5 * (a * x * x + b * y * y));
}

static inline double dot(double x1, double y1, double x2, double y2) {
    return x1*x2 + y1*y2;
}

double v_matching(double *vx, double *vy, int num) {
    double sum = 0.0f;
    for (int i = 0; i < num; i++) {
        for (int j = i + 1; j < num; j++) {
            double diff = norm(vx[i] - vx[j], vy[i] - vy[j])
                / (norm(vx[i], vy[i]) + norm(vx[j], vy[j]));
            sum += diff * diff;
        }
    }
    return sum;
}

// minimum distance
double cpa(double x1, double y1, double x2, double y2, double vx1, double vy1, double vx2, double vy2)
{
    double dvx = vx1 - vy1, dvy = vx2 - vy2;
    double dv2 = dot(dvx,dvy,dvx,dvy);
    double cpatime = 0;
    if (dv2 > 1e-8) {
        double wx = x1 - x2, wy = y1 - y2;
        cpatime = -dot(wx, wy, dvx, dvy) / dv2;
    }
    if (cpatime < 0 || cpatime > 1)
        return INFINITY;
    else
        return norm(x1 - x2 + cpatime * (vx1 - vx2), y1 - y2 + cpatime * (vy1 - vy2));
}

//
void trim(double *x, double *y, double mag)
{
    if (norm(*x, *y) <= mag)
        return;
    double theta = atan2(*y, *x);
    *x = mag * cos(theta);
    *y = mag * sin(theta);
    printf("trim\n");
}

void init(double x[Steps][Num], double y[Steps][Num], double vx[Steps][Num],
        double vy[Steps][Num]) {
    int isCollision;
    do {
        isCollision = 0;
        for (int i = 0; i < Steps; i++) {
            for (int j = 0; j < Num; j++) {
                if (i != 0) {
                    x[i][j] = y[i][j] = vx[i][j] = vy[i][j] = 0;
                } else {
                    x[i][j] = randd(0, Init_box);
                    y[i][j] = randd(0, Init_box);
                    vx[i][j] = randd(InitVmin, InitVmax);
                    vy[i][j] = randd(InitVmin, InitVmax);
                }
            }
        }
        for (int bi = 0; bi < Num; bi++) {
            for (int bj = bi + 1; bj < Num; bj++) {
                if (norm(x[0][bi]-x[0][bj], y[0][bi]-y[0][bj]) < Dmin)
                    isCollision = 1;
            }
        }
    } while (isCollision);
}


double calculateFitness(double nextVX[Num], double nextVY[Num], double nextX[Num], double nextY[Num], double firstVX[Num], double firstVY[Num], void *params, int forTheRecord)
{
    flock_info *info = (flock_info *) params;
    double obstacle = 0.0, benefit = 0.0, ca = 0.0;
    double blocks[Num][2];
    double px = 0.0, py = 0.0, k = 0.0, A = 0.0, B = 0.0, C = 0.0, side = 0.0, h_dis = 0.0,
           v_dis = 0.0, sm = 0.0, dot_prod = 0.0, ub_j = 0.0;
    for (int i = 0; i < Num; i++) {
        memset(blocks, 0, sizeof(double) * Num * 2);
        A = nextVX[i];
        B = nextVY[i];
        C = -nextVY[i] * nextY[i] - nextVX[i] * nextX[i];
        ub_j = 0;
        for (int j = 0; j < Num; j++) {
            if (j != i) {
                if(forTheRecord == 0){
                    if (cpa(info->cx[i], info->cy[i], info->cx[j], info->cy[j], 
                                firstVX[i], firstVY[i], firstVX[j], firstVY[j]) < Dmin) {
                        ca = INFINITY;
                        break;
                    }
                }
                if (nextVX[i] == 0.0) {
                    px = nextX[j];
                    py = nextY[i];
                } else if (nextVY[i] == 0.0) {
                    px = nextX[i];
                    py = nextY[j];
                } else {
                    k = -nextVX[i] / nextVY[i];
                    px = (k * nextX[i] + nextX[j] / k + nextY[j] - nextY[i])
                        / (k + 1.0 / k);
                    py = -1.0 / k * (px - nextX[j]) + nextY[j];
                }

                side = A * nextX[j] + B * nextY[j] + C;
                h_dis = norm(px - nextX[i], py - nextY[i]);
                v_dis = fabs(side) / norm(A, B);

                if (side >= 0.0 && (h_dis < w || (h_dis - w) / v_dis < tan(angle))) {
                    blocks[j][0] = atan(v_dis / (h_dis + w));
                    blocks[j][1] = atan2(v_dis, h_dis - w);
                    if (blocks[j][0] < PI / 2.0 - angle)
                        blocks[j][0] = PI / 2.0 - angle;
                    if (blocks[j][1] > PI / 2.0 + angle)
                        blocks[j][1] = PI / 2.0 + angle;
                    obstacle += (blocks[j][1] - blocks[j][0]) / (angle*2);
                }

                sm = erf((h_dis - (w - lambda)) * sqrt(2.0) * 8.0);
                dot_prod = (nextVX[i] * nextVX[j] + nextVY[i] * nextVY[j])
                    / (norm(nextVX[i], nextVY[i]) * norm(nextVX[j], nextVY[j]));
                if (side > 0.0 && h_dis >= w - lambda)
                    ub_j += dot_prod * sm
                        * mvnpdf(h_dis - (2.0 * w - lambda), v_dis - d0,
                                u_sigma1, u_sigma2);
                else if (side >= 0.0 && h_dis < w - lambda)
                    ub_j += sm * mvnpdf(h_dis, v_dis, d_sigma1, d_sigma2);
            }
        }

        if(ub_j < 1.0){
            benefit += ub_j;
        }
        else{
            benefit += 1.0;
        }

        //benefit += MIN(ub_j,1);
    }
    return pow(v_matching(nextVX, nextVY, Num),2) + 4*pow(obstacle,2) + pow(Num - 1.0 - benefit,2) + ca;
}


double flock_fit(double *va, int dim, void *params, int *ph) {
    flock_info *info = (flock_info *) params;

    double curFitness;
    double bestFitness = DBL_MAX;
    double nextX[Num], nextY[Num], nextVX[Num], nextVY[Num], prevX[Num], prevY[Num], prevVX[Num], prevVY[Num], firstVX[Num], firstVY[Num];
    for (int j = 0; j < Ph; j++){
        for (int i = 0; i < Num; i++) {
            if(j == 0){
                nextVX[i] = info->cvx[i] + va[i] * cos(va[i + Ph * Num]);
                nextVY[i] = info->cvy[i] + va[i] * sin(va[i + Ph * Num]);
                nextX[i] = info->cx[i] + nextVX[i];
                nextY[i] = info->cy[i] + nextVY[i];
                firstVX[i] = nextVX[i];
                firstVY[i] = nextVY[i];
            }
            else{
                nextVX[i] = prevVX[i] + va[i+j*Num] * cos(va[i + (Ph + j) * Num]);
                nextVY[i] = prevVY[i] + va[i+j*Num] * sin(va[i + (Ph + j) * Num]);
                nextX[i] = prevX[i] + nextVX[i];
                nextY[i] = prevY[i] + nextVY[i];
            }
            prevVX[i] = nextVX[i];
            prevVY[i] = nextVY[i];
            prevX[i] = nextX[i];
            prevY[i] = nextY[i];
        }
        curFitness = calculateFitness(nextVX, nextVY, nextX, nextY, firstVX, firstVY, params, 0);
        if(curFitness <= bestFitness){
            *ph = j+1;
            bestFitness = curFitness;		
        }
    }

    return bestFitness;

}

void disp(int t, double x[Steps][Num], double y[Steps][Num],
        double vx[Steps][Num], double vy[Steps][Num]) {
    for (int i = 0; i < Num; i++) {
        printf("%.4f\t%.4f\t%.4f\t%.4f\n", x[t][i], y[t][i], vx[t][i],
                vy[t][i]);
    }
}

char* concat(char *s1, char *s2){
    char *result = malloc(strlen(s1) + strlen(s2)+1);
    strcpy(result, s1);
    strcat(result, s2);
    return result;
}

void save(double x[Steps][Num], double y[Steps][Num], double vx[Steps][Num],
        double vy[Steps][Num]) {
    FILE *fp;

    char numstr[21];
    sprintf(numstr, "%i", currentConf);
    char numstr2[21];
    sprintf(numstr2, "%i", currentRun);
    fp = fopen(concat(concat(concat(concat("config", numstr),"_"),numstr2), ".txt"), "w");

    for (int t = 0; t < Steps; t++) {
        for (int i = 0; i < Num; i++) {
            fprintf(fp, "%f\t", x[t][i]);
        }
        for (int i = 0; i < Num; i++) {
            fprintf(fp, "%f\t", y[t][i]);
        }
        for (int i = 0; i < Num; i++) {
            fprintf(fp, "%f\t", vx[t][i]);
        }
        for (int i = 0; i < Num; i++) {
            fprintf(fp, "%f\t", vy[t][i]);
        }
        fprintf(fp, "\n");
    }
    fclose(fp);
}

void saveFitness(int timestep, double fitness){
    //	printf("fitness at step %i: %.15f\n",timestep, fitness);
    char numstr[21];
    sprintf(numstr, "%i", currentConf);
    char* filename = concat(concat("fit", numstr), ".txt");
    FILE *fp;
    fp = fopen(filename, "a");//strcat(strcat("fit", to_string(currentConf)),".txt"), "w");
    if(timestep == 0){
        if(currentRun == 0){
            fprintf(fp,"%i;",0);
        }	
        else{
            fprintf(fp,"\n%i;",currentRun);
        }
    }
    fprintf(fp, "%f;",fitness);
    fclose(fp);
}

void flock_pso(double x[Steps][Num], double y[Steps][Num],
        double vx[Steps][Num], double vy[Steps][Num]) {
    objFcn obj_fun = flock_fit;
    flock_info init_info;
    init_info.cx = x[0];
    init_info.cy = y[0];
    init_info.cvx = vx[0];
    init_info.cvy = vy[0];
    init_info.step = 0;
    double initFitness = calculateFitness(vx[0], vy[0], x[0], y[0], vx[0], vy[0], &init_info, 1);
    double curFitness = initFitness;
    printf("curFitness: %f\n", curFitness);
    int t = 0; int v_reached = 0; int k = 0;
    for (t = 0; t < Steps - 1 && v_reached == 0;) {
        pso_result sol[Clones];
        int best_sol_idx = 0;
        double best_sol = DBL_MAX;
        int best_sol_ph = 1;
        flock_info info;
        info.cx = x[t];
        info.cy = y[t];
        info.cvx = vx[t];
        info.cvy = vy[t];
        info.step = t;
        for (k = 0; k < Clones; ++k) {
            for (Ph = 1; Ph <= MaxPh; ++Ph) {
                pso_options options;

                double nvars = 2 * Num * Ph;
                settings(&options, nvars);
                sol[k].gbest = malloc(options.dim * sizeof(double));

                options.lb = malloc(nvars * sizeof(double));
                options.ub = malloc(nvars * sizeof(double));
                for (int i = 0; i < nvars; i++) {
                    if (i < Ph * Num) {
                        options.lb[i] = 0.0;
                        options.ub[i] = 0.0;
                    } else {
                        options.lb[i] = 0.0;
                        options.ub[i] = 2.0*PI;
                    }
                }

                for (int i = 0; i < Num; i++) {
                    options.ub[i] = delta * norm(vx[t][i], vy[t][i]);
                }
                for (int p = 1; p < Ph; p++) {
                    memcpy(options.ub+p*Num, options.ub, Num * sizeof(double));
                }

                pso(obj_fun, &info, &sol[k], &options);
                if (curFitness - sol[k].error >=  curFitness/(Steps - t) || Ph == MaxPh) {
                    if (sol[k].error < best_sol) {
                        best_sol = sol[k].error;
                        best_sol_idx = k;
                        best_sol_ph = Ph;
                    }
                    free(options.lb);
                    free(options.ub);
                    break;
                } else {
                    free(options.lb);
                    free(options.ub);
                }
            }
        }
        for (int p = 1; p <= 1 && t < Steps - 1; p++, t++) {
            for (int i = 0; i < Num; i++) {
                double ax = sol[best_sol_idx].gbest[i+(p-1)*Num] * cos(sol[best_sol_idx].gbest[i + (best_sol_ph+p-1) * Num]);
                double ay = sol[best_sol_idx].gbest[i+(p-1)*Num] * sin(sol[best_sol_idx].gbest[i + (best_sol_ph+p-1) * Num]);
                vx[t + 1][i] = vx[t][i] + ax;
                vy[t + 1][i] = vy[t][i] + ay;
                x[t + 1][i] = x[t][i] + vx[t+1][i];
                y[t + 1][i] = y[t][i] + vy[t+1][i];
            }
        }
        curFitness = calculateFitness(vx[t],vy[t], x[t], y[t], vx[t], vy[t], &info, 0);
        printf("found. sol.error: %f ph:%d, step: %d, curFitness: %f\n", sol[best_sol_idx].error, best_sol_ph, t, curFitness);
        for (k = 0; k < Clones; ++k) {
            free(sol[k].gbest);
        }
        if(curFitness <= 1e-4) {
            v_reached = 1;
            break;  
        }
        if(t == 0){
            disp(t, x, y, vx, vy);
        }
    }
    printf("step: %d\n", t);
    //   for (; t < Steps - 1; t++) {
    //     for (int i = 0; i < Num; i++) {
    //       vx[t + 1][i] = vx[t][i];
    //       vy[t + 1][i] = vy[t][i];
    //       x[t + 1][i] = x[t][i] + vx[t + 1][i];
    //       y[t + 1][i] = y[t][i] + vy[t + 1][i];
    //     }
    //   }
}

void help(){
    printf("======= HELP ========\n First argument: [int] unique number of current CONFIGURATION tried\n Second argument: [int] number of RUNS performed by this current configuration\n ======================\n");
}

int main(int argc, char *argv[]) {

    if (argc > 1){
        int helpcalled = 0;
        for(int i = 1; i < argc; i++){
            if(strcmp(argv[i], "-help")==0){
                help();
                helpcalled = 1;
            }
        }
        if(helpcalled == 0){
            if(argc > 2){	
                currentConf = atoi(argv[1]);
                int runs = atoi(argv[2]);
                for(int k = 0; k < runs; k++){ //iterate through the amount of runs per configuration
                    currentRun = k;
                    srand(currentConf+time(NULL));
                    double x[Steps][Num], y[Steps][Num], vx[Steps][Num], vy[Steps][Num];
                    int loading = 0;
                    if (!loading){
                        init(x, y, vx, vy);
                        //saveConf(0,x,y,vx,vy);
                    }
                    disp(0, x, y, vx, vy);
                    srand(time(NULL));
                    clock_t begin, end;
                    double time_spent;
                    begin = clock();

                    flock_pso(x, y, vx, vy);
                    //disp(Steps - 1, x, y, vx, vy);
                    save(x, y, vx, vy);

                    //saveConf(Steps-1,x,y,vx,vy);

                    end = clock();
                    time_spent = (double) (end - begin) / CLOCKS_PER_SEC;
                    printf("duration: %f\n", time_spent);
                }
            }
            else
                help();
        }
    }
    else{
        help();
    }
}
