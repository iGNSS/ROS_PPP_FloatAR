/*------------------------------------------------------------------------------
 * ppp_ar.c : ppp ambiguity resolution
 *
 * reference :
 *    [1] H.Okumura, C-gengo niyoru saishin algorithm jiten (in Japanese),
 *        Software Technology, 1991
 *
 *          Copyright (C) 2012-2015 by T.TAKASU, All rights reserved.
 *
 * version : $Revision:$ $Date:$
 * history : 2013/03/11  1.0  new
 *           2015/05/15  1.1  refine complete algorithms
 *           2015/05/31  1.2  delete WL-ambiguity resolution by ILS
 *                            add PAR (partial ambiguity resolution)
 *           2015/11/26  1.3  support option opt->pppopt=-TRACE_AR
 *           2015/12/25  1.4  disable GPS-QZS ambiguity resolution
 *                            fix problem on zero-divide for ratio-factor
 *-----------------------------------------------------------------------------*/
#include "rtklib.h"

static const char rcsid[] = "$Id:$";

#define MIN_AMB_RES 4     /* min number of ambiguities for ILS-AR */
#define MIN_ARC_GAP 300.0 /* min arc gap (s) */
#define CONST_AMB 0.003   /* constraint to fixed ambiguity */

#define LOG_PI 1.14472988584940017 /* log(pi) */
#define SQRT2 1.41421356237309510  /* sqrt(2) */

#define MIN(x, y) ((x) < (y) ? (x) : (y))
#define MAX(x, y) ((x) > (y) ? (x) : (y))
#define SQR(x) ((x) * (x))
#define ROUND(x) (int)floor((x) + 0.5)
#define SWAP_I(x, y)  \
    do                \
    {                 \
        int _tmp = x; \
        x = y;        \
        y = _tmp;     \
    } while (0)
#define SWAP_D(x, y)     \
    do                   \
    {                    \
        double _tmp = x; \
        x = y;           \
        y = _tmp;        \
    } while (0)

/* number and index of ekf states */
#define NF(opt) ((opt)->ionoopt == IONOOPT_IFLC ? 1 : (opt)->nf)
#define NP(opt) ((opt)->dynamics ? 9 : 3)
#define NC(opt) (NSYS)
#define NT(opt) ((opt)->tropopt < TROPOPT_EST ? 0 : ((opt)->tropopt == TROPOPT_EST ? 1 : 3))
#define NI(opt) ((opt)->ionoopt == IONOOPT_EST ? MAXSAT : 0)
/*#define ND(opt)     ((opt)->nf>=3?1:0)
#define NR(opt)     (NP(opt)+NC(opt)+NT(opt)+NI(opt)+ND(opt)+1)*/
#define ND(opt) ((opt)->nf >= 2 ? 2 : 0)
#define NR(opt) (NP(opt) + NC(opt) + NT(opt) + NI(opt) + ND(opt))
#define IB(s, f, opt) (NR(opt) + MAXSAT * (f) + (s)-1)
/* wave length of LC (m) -----------------------------------------------------*/
static double lam_LC(int i, int j, int k)
{
    const double f1 = FREQL1, f2 = FREQL2, f5 = FREQL5;

    return CLIGHT / (i * f1 + j * f2 + k * f5);
}
/* wave length of LC (m) -----------------------------------------------------*/
static double lam_LC_CMP(int i, int j, int k)
{
    const double f1 = FREQ1_BDS, f2 = FREQ2_BDS, f5 = FREQ3_BDS;

    return CLIGHT / (i * f1 + j * f2 + k * f5);
}
/* complementaty error function (ref [1] p.227-229) --------------------------*/
static double q_gamma(double a, double x, double log_gamma_a);
static double p_gamma(double a, double x, double log_gamma_a)
{
    double y, w;
    int i;

    if (x == 0.0)
        return 0.0;
    if (x >= a + 1.0)
        return 1.0 - q_gamma(a, x, log_gamma_a);
    y = w = exp(a * log(x) - x - log_gamma_a) / a;
    for (i = 1; i < 100; i++)
    {
        w *= x / (a + i);
        y += w;
        if (fabs(w) < 1E-15)
            break;
    }
    return y;
}
static double q_gamma(double a, double x, double log_gamma_a)
{
    double y, w, la = 1.0, lb = x + 1.0 - a, lc;
    int i;

    if (x < a + 1.0)
        return 1.0 - p_gamma(a, x, log_gamma_a);
    w = exp(-x + a * log(x) - log_gamma_a);
    y = w / lb;
    for (i = 2; i < 100; i++)
    {
        lc = ((i - 1 - a) * (lb - la) + (i + x) * lb) / i;
        la = lb;
        lb = lc;
        w *= (i - 1 - a) / i;
        y += w / la / lb;
        if (fabs(w / la / lb) < 1E-15)
            break;
    }
    return y;
}
static double f_erfc(double x)
{
    return x >= 0.0 ? q_gamma(0.5, x * x, LOG_PI / 2.0) : 1.0 + p_gamma(0.5, x * x, LOG_PI / 2.0);
}
/* confidence function of integer ambiguity ----------------------------------*/
static double conf_func(int N, double B, double var)
{
    double x, p = 1.0, sig = sqrt(var);
    int i;

    x = fabs(B - N);
    for (i = 1; i < 6; i++)
    {
        p -= f_erfc((i - x) / (SQRT2 * sig)) - f_erfc((i + x) / (SQRT2 * sig));
    }
    return p;
}
/* generate satellite SD (single-difference) ---------------------------------*/
static int gen_sat_sd(rtk_t *rtk, const obsd_t *obs, int n, const int *exc,
                      const double *azel, int f, int *sat1, int *sat2, int *frq)
{
#if 0
    const int sys[]={SYS_GPS|SYS_QZS,SYS_GAL,SYS_CMP,0};
#else
    const int sys[] = {SYS_GPS, SYS_GAL, SYS_BDS, 0};
#endif
    double elmask, el[MAXOBS];
    int i, j, k, m, ns = 0, sat[MAXOBS];

    elmask = MAX(rtk->opt.elmaskar, rtk->opt.elmin);

    for (i = 0; sys[i]; i++)
    { /* for each system */

        /* sort by elevation angle */
        for (j = m = 0; j < n; j++)
        {
            if (exc[j] || !(satsys(obs[j].sat, NULL) & sys[i]) ||
                !rtk->ssat[obs[j].sat - 1].vsat[f] || azel[1 + j * 2] < elmask)
                continue;
            if (obs[j].L[f] == 0.0)
                continue;
            sat[m] = obs[j].sat;
            el[m++] = azel[1 + j * 2];
        }
        for (j = 0; j < m - 1; j++)
            for (k = j + 1; k < m; k++)
            {
                if (el[j] >= el[k])
                    continue;
                SWAP_I(sat[j], sat[k]);
                SWAP_D(el[j], el[k]);
            }
        /* generate SD referenced to max elevation angle */
        for (j = 1; j < m; j++)
        {
            sat1[ns] = sat[j];
            sat2[ns] = sat[0];
            frq[ns++] = f;
        }
    }
    return ns; /* # of SD */
}
/* filter EWL-ambiguity ------------------------------------------------------*/
static void filter_EWL(rtk_t *rtk, const obsd_t *obs, int n, const int *exc,
                       const nav_t *nav, const double *azel)
{
    const double *lam;
    ambc_t *amb;
    double lamE, lamN, EW, var, K;
    int i, sat, sys;

    for (i = 0; i < n; i++)
    {
        sat = obs[i].sat;
        sys = satsys(sat, NULL);
        if (sys != SYS_GPS && sys != SYS_QZS)
            continue;
        amb = rtk->ambc + sat - 1;
        lam = nav->lam[sat - 1];

        if (rtk->ssat[sat - 1].slip[1] || rtk->ssat[sat - 1].slip[2] ||
            fabs(timediff(amb->epoch[0], obs[i].time)) > MIN_ARC_GAP)
        {
            amb->n[1] = 0;
        }
        if (exc[i] || azel[1 + 2 * i] < rtk->opt.elmin || lam[1] == 0.0 || lam[2] == 0.0 ||
            obs[i].L[1] == 0.0 || obs[i].P[2] == 0.0 ||
            obs[i].L[1] == 0.0 || obs[i].P[2] == 0.0)
            continue;

        /* EMW-LC and variance */
        lamE = lam[1] * lam[2] / (lam[2] - lam[1]);
        lamN = lam[1] * lam[2] / (lam[2] + lam[1]);
        EW = lamE * (obs[i].L[1] - obs[i].L[2]) -
             lamN * (obs[i].P[1] / lam[1] + obs[i].P[2] / lam[2]);
        var = SQR(0.7 * rtk->opt.eratio[0]) *
              (SQR(rtk->opt.err[1]) + SQR(rtk->opt.err[2] / sin(azel[1 + 2 * i])));

        /* filter EWL-ambiguity */
        if (amb->n[1] <= 0)
        {
            amb->LC[1] = EW;
            amb->LCv[1] = SQR(10.0);
            amb->n[1] = 1;
            amb->epoch[1] = obs[i].time;
        }
        else if (SQR(EW - amb->LC[1]) <= amb->LCv[1] * 25.0)
        {
            K = amb->LCv[1] / (amb->LCv[1] + var);
            amb->LC[1] += K * (EW - amb->LC[1]);
            amb->LCv[1] -= K * amb->LCv[1];
            amb->n[1]++;
            amb->epoch[1] = obs[i].time;
        }
    }
}
/* filter WL-ambiguity -------------------------------------------------------*/
static void filter_WL(rtk_t *rtk, const obsd_t *obs, int n, const int *exc,
                      const nav_t *nav, const double *azel)
{
    const double *lam;
    ambc_t *amb;
    double lamW, lamN, MW, var, K, C1, C2, P[NFREQ];
    int i, l, sat;

    for (i = 0; i < n; i++)
    {
        sat = obs[i].sat;
        amb = rtk->ambc + sat - 1;
        lam = nav->lam[sat - 1];
        l = satsys(sat, NULL) == SYS_GAL ? 2 : 1; /* L1/L2 or L1/L5 */

        if (rtk->ssat[sat - 1].slip[0] || rtk->ssat[sat - 1].slip[l] ||
            fabs(timediff(amb->epoch[0], obs[i].time)) > MIN_ARC_GAP)
        {
            amb->n[0] = 0;
        }
        if (exc[i] || azel[1 + 2 * i] < rtk->opt.elmin || lam[0] == 0.0 || lam[l] == 0.0 ||
            obs[i].L[0] == 0.0 || obs[i].P[0] == 0.0 ||
            obs[i].L[l] == 0.0 || obs[i].P[l] == 0.0)
            continue;

        /* P1-P2 dcb correction (P1->Pc,P2->Pc)  */
        C1 = SQR(lam[l]) / (SQR(lam[l]) - SQR(lam[0]));
        C2 = SQR(lam[0]) / (SQR(lam[l]) - SQR(lam[0]));
        P[0] = obs[i].P[0] + C2 * nav->cbias[obs[i].sat - 1][0];
        P[l] = obs[i].P[l] + C1 * nav->cbias[obs[i].sat - 1][0];

        /* MW-LC and variance */
        lamW = lam[0] * lam[l] / (lam[l] - lam[0]);
        lamN = lam[0] * lam[l] / (lam[l] + lam[0]);
        MW = lamW * (obs[i].L[0] - obs[i].L[l]) -
             lamN * (P[0] / lam[0] + P[l] / lam[l]);
        var = SQR(0.7 * rtk->opt.eratio[0]) *
              (SQR(rtk->opt.err[1]) + SQR(rtk->opt.err[2]) / sin(azel[1 + 2 * i]));

        /* filter WL-ambiguity */
        if (amb->n[0] <= 0)
        {
            amb->LC[0] = MW;
            amb->LCv[0] = SQR(10.0);
            amb->n[0] = 1;
            amb->epoch[0] = obs[i].time;
        }
        else if (SQR(MW - amb->LC[0]) <= amb->LCv[0] * 25.0 * 10)
        {
            K = amb->LCv[0] / (amb->LCv[0] + var);
            amb->LC[0] += K * (MW - amb->LC[0]);
            amb->LCv[0] -= K * amb->LCv[0];
            amb->n[0]++;
            amb->epoch[0] = obs[i].time;
        }
    }
}
/* ambiguity resolution by WL/NL for iono-free LC ----------------------------*/
static int ppp_amb_IFLC(rtk_t *rtk, const obsd_t *obs, int n, int *exc,
                        const nav_t *nav, const double *azel, double *x,
                        double *P)
{
    const double *lam;
    ambc_t *amb1, *amb2;
    double lamN, lamW, lamE, C1, C2, Be, Bw, B1, varw, var1;
    double v[MAXOBS], R[MAXOBS * MAXOBS] = {0}, *H, *thres = rtk->opt.thresar;
    int i, j, k, l, info, ns, nw = 0, na = 0, Ne, Nw, N1, sat1[MAXOBS], sat2[MAXOBS], frq[MAXOBS];
    double tow;
    int week;
    tow = time2gpst(rtk->sol.time, &week);

    /* filter EWL-ambiguity */
    filter_EWL(rtk, obs, n, exc, nav, azel);

    /* filter WL-ambiguity */
    filter_WL(rtk, obs, n, exc, nav, azel);

    /* generate satellite SD */
    if (!(ns = gen_sat_sd(rtk, obs, n, exc, azel, 0, sat1, sat2, frq)))
        return 0;

    H = zeros(rtk->nx, ns);

    for (i = 0; i < ns; i++)
    {
        lam = nav->lam[sat1[i] - 1];
        l = satsys(sat1[i], NULL) == SYS_GAL ? 2 : 1; /* L1/L2 or L1/L5 */
        lamE = lam[1] * lam[2] / (lam[2] - lam[1]);
        lamW = lam[0] * lam[l] / (lam[l] - lam[0]);
        lamN = lam[0] * lam[l] / (lam[l] + lam[0]);
        C1 = SQR(lam[l]) / (SQR(lam[l]) - SQR(lam[0]));
        C2 = -SQR(lam[0]) / (SQR(lam[l]) - SQR(lam[0]));

        j = IB(sat1[i], 0, &rtk->opt);
        k = IB(sat2[i], 0, &rtk->opt);
        amb1 = rtk->ambc + sat1[i] - 1;
        amb2 = rtk->ambc + sat2[i] - 1;
        if (!amb1->n[0] || !amb2->n[0])
            continue;

        if (amb1->n[1] && amb2->n[1])
        {
            Be = (amb1->LC[1] - amb2->LC[1]) / lamE;
            Ne = ROUND(Be);
#if 0
            trace(2,"%s sat=%2d-%2d:Be=%13.4f Fe=%7.4f\n",time_str(obs[0].time,0),
                  sat1[i],sat2[i],Be,Be-Ne);
#endif
        }
        /* round WL- and L1-ambiguity */
        Bw = (amb1->LC[0] - amb2->LC[0]) / lamW;
        Bw += nav->wlbias[sat1[i] - 1] - nav->wlbias[sat2[i] - 1]; /* correct WL-bias sat2 is ref */
        Nw = ROUND(Bw);
        B1 = (x[j] - x[k] + C2 * lam[l] * Nw) / lamN;
        N1 = ROUND(B1);
        varw = (amb1->LCv[0] + amb2->LCv[0]) / SQR(lamW);
        var1 = (P[j + j * rtk->nx] + P[k + k * rtk->nx] - 2.0 * P[j + k * rtk->nx]) / SQR(lamN);

        if (conf_func(Nw, Bw, varw) < thres[1] || fabs(Nw - Bw) > thres[2])
            continue;
        nw++;
        /*------added by xiang for output the errors debug only those pass fix of WL------*/

        /*narrow lane using rounding*/
        if (conf_func(N1, B1, var1) < thres[1] || fabs(N1 - B1) > thres[3])
            continue;

        /* constraint to fixed LC-ambiguty */
        H[j + na * rtk->nx] = 1.0;
        H[k + na * rtk->nx] = -1.0;
        v[na++] = C1 * lam[0] * N1 + C2 * lam[l] * (N1 - Nw) - (rtk->x[j] - rtk->x[k]);

        /* update fix flags */
        rtk->ssat[sat1[i] - 1].fix[0] = rtk->ssat[sat2[i] - 1].fix[0] = 2;
    }
    rtk->sol.age = (float)nw;   /* # of WL-fixed */
    rtk->sol.ratio = (float)na; /* # of NL-fixed */

    if (na <= 0)
    {
        free(H);
        return 0;
    }
    for (i = 0; i < na; i++)
        R[i + i * na] = SQR(CONST_AMB);

    /* update states with fixed LC-ambiguity constraints */
    if ((info = filter(x, P, H, v, R, rtk->nx, na)))
    {
        trace(1, "filter error (info=%d)\n", info);
        free(H);
        return 0;
    }
    free(H);
    return 1;
}

/* fixed solution ------------------------------------------------------------*/
static int fix_sol(rtk_t *rtk, const int *sat1, const int *sat2,
                   const double *NC, int n)
{
    double *v, *H, *R;
    int i, j, k, info;

    trace(2, "solution fix!\n");
    if (n <= 0)
        return 0;

    v = zeros(n, 1);
    H = zeros(rtk->nx, n);
    R = zeros(n, n);

    /* constraints to fixed ambiguities */
    for (i = 0; i < n; i++)
    {
        j = IB(sat1[i], 0, &rtk->opt);
        k = IB(sat2[i], 0, &rtk->opt);
        v[i] = NC[i] - (rtk->x[j] - rtk->x[k]);
        H[j + i * rtk->nx] = 1.0;
        H[k + i * rtk->nx] = -1.0;
        R[i + i * n] = SQR(CONST_AMB);
    }
    /* update states with constraints */
    if ((info = filter(rtk->x, rtk->P, H, v, R, rtk->nx, n)))
    {
        trace(1, "filter error (info=%d)\n", info);
        free(v);
        free(H);
        free(R);
        return 0;
    }
    /* set solution */
    for (i = 0; i < rtk->nx; i++)
    {
        rtk->xa[i] = rtk->x[i];
        for (j = 0; j < rtk->nx; j++)
        {
            rtk->Pa[i + j * rtk->nx] = rtk->Pa[j + i * rtk->nx] = rtk->P[i + j * rtk->nx];
        }
    }
    rtk->na = n;
    /* set flags */
    for (i = 0; i < n; i++)
    {
        rtk->ambc[sat1[i] - 1].flags[sat2[i] - 1] = 1;
        rtk->ambc[sat2[i] - 1].flags[sat1[i] - 1] = 1;
    }

    free(v);
    free(H);
    free(R);
    return 1;
}
/* linear dependency check ---------------------------------------------------*/
static int is_depend(int sat1, int sat2, int *flgs, int *max_flg)
{
    int i;

    if (flgs[sat1 - 1] == 0 && flgs[sat2 - 1] == 0)
    {
        flgs[sat1 - 1] = flgs[sat2 - 1] = ++(*max_flg);
    }
    else if (flgs[sat1 - 1] == 0 && flgs[sat2 - 1] != 0)
    {
        flgs[sat1 - 1] = flgs[sat2 - 1];
    }
    else if (flgs[sat1 - 1] != 0 && flgs[sat2 - 1] == 0)
    {
        flgs[sat2 - 1] = flgs[sat1 - 1];
    }
    else if (flgs[sat1 - 1] > flgs[sat2 - 1])
    {
        for (i = 0; i < MAXSAT; i++)
            if (flgs[i] == flgs[sat2 - 1])
                flgs[i] = flgs[sat1 - 1];
    }
    else if (flgs[sat1 - 1] < flgs[sat2 - 1])
    {
        for (i = 0; i < MAXSAT; i++)
            if (flgs[i] == flgs[sat1 - 1])
                flgs[i] = flgs[sat2 - 1];
    }
    else
        return 0; /* linear depenent */
    return 1;
}
/* ppp partial ambiguity resolutions ------------------------------------------*/
static int ppp_amb_IFLAM_PAR(rtk_t *rtk, const nav_t *nav, int m, int n, int *sat1, int *sat2, int *Nw, double *z, double *N1, const double *Q)
{
    const double *lam;
    double lamN, lamW;
    double *D, *E, *B1, *Q_, s[2];
    int i, j, k, l, info, del, p;
    double var[MAXOBS];
    char sys[32];

    /* init */
    del = 1;
    p = m - del;

    /* extract variance */
    for (i = 0; i < m; i++)
    {
        var[i] = Q[i + i * m];
    }

    /* sort by variance */
    for (i = 0; i < m; i++)
        for (j = 1; j < m - i; j++)
        {
            if (var[j] >= var[j - 1])
                continue;
            SWAP_I(sat1[j], sat1[j - 1]);
            SWAP_I(sat2[j], sat2[j - 1]);
            SWAP_I(Nw[j], Nw[j - 1]);
            SWAP_D(N1[j], N1[j - 1]);
            SWAP_D(z[j], z[j - 1]);
            SWAP_D(var[j], var[j - 1]);
        }

    /* iteration */
    while (p > 4)
    {
        D = zeros(rtk->nx, m);
        E = mat(m, rtk->nx);
        B1 = zeros(p, 1);
        Q_ = mat(p, p);

        p = m - del;
        for (i = 0; i < p; i++)
        {
            if (satsys(sat1[i], NULL) == SYS_GPS)
            {
                lamN = lam_LC(1, 1, 0);
                lamW = lam_LC(1, -1, 0);
            }
            else if (satsys(sat1[i], NULL) == SYS_GAL)
            {
                lamN = lam_LC(1, 0, 1);
                lamW = lam_LC(1, 0, -1);
            }
            else if (satsys(sat1[i], NULL) == SYS_BDS)
            {
                lamN = lam_LC_CMP(1, 0, 1);
                lamW = lam_LC_CMP(1, 0, -1);
            }

            j = IB(sat1[i], 0, &rtk->opt);
            k = IB(sat2[i], 0, &rtk->opt);

            /* narrow-line ambiguity transformation matrix */
            D[j + i * rtk->nx] = 1.0 / lamN;
            D[k + i * rtk->nx] = -1.0 / lamN;
            B1[i] = z[i];
            N1[i] = ROUND(z[i]);
        }
        /* covariance of narrow-lane ambiguities */
        matmul("TN", p, rtk->nx, rtk->nx, 1.0, D, rtk->P, 0.0, E);
        matmul("NN", p, p, rtk->nx, 1.0, E, D, 0.0, Q_);

        if ((info = lambda(p, 2, B1, Q_, N1, s)) || (s[0] <= 0.0))
        {
            del++;
        }
        else
        {
            rtk->sol.ratio = (float)(MIN(s[1] / s[0], 999.9));
            if (rtk->sol.ratio < rtk->opt.thresar[0])
            {
                trace(2, "ratio test error: ratio=%f ,threshold=%f \n", rtk->sol.ratio, rtk->opt.thresar[0]);
                del++;
            }
            else
            {

                free(D);
                free(E);
                free(B1);
                free(Q_);
                return del;
            }
        }
        free(D);
        free(E);
        free(B1);
        free(Q_);
    }
    return 0;
}
/* fix narrow-lane ambiguity by ILS ------------------------------------------*/
static int fix_amb_ILS(rtk_t *rtk, const nav_t *nav, int *sat1, int *sat2, int *NW, int n)
{
    double C1, C2, *B1, *N1, *NC, *D, *E, *Q, s[2], lam_NL, lam1, lam2;
    int i, j = 0, k = 0, m = 0, info, stat, del = 0, flgs[MAXSAT] = {0}, max_flg = 0, l = 0;
    double *lam = {0}, var1;
    int PARerror;

    /*added by xiang inorder for output ambiguity 2015-11-16*/
    int week;
    double tow;

    tow = time2gpst(rtk->sol.time, &week);
    /*end of addition*/

    B1 = zeros(n, 1);
    N1 = zeros(n, 2);
    D = zeros(rtk->nx, n);
    E = mat(n, rtk->nx);
    Q = mat(n, n);
    NC = mat(n, 1);

    for (i = 0; i < n; i++)
    {
        lam = nav->lam[sat1[i] - 1];
        l = satsys(sat1[i], NULL) == SYS_GAL ? 2 : 1; /* L1/L2 or L1/L5 */
        lam_NL = lam[0] * lam[l] / (lam[l] + lam[0]);
        C1 = SQR(lam[l]) / (SQR(lam[l]) - SQR(lam[0]));
        C2 = -SQR(lam[0]) / (SQR(lam[l]) - SQR(lam[0]));
        /* check linear independency */
        if (!is_depend(sat1[i], sat2[i], flgs, &max_flg))
            continue;

        j = IB(sat1[i], 0, &rtk->opt);
        k = IB(sat2[i], 0, &rtk->opt);

        /* float narrow-lane ambiguity (cycle) */
        B1[m] = (rtk->x[j] - rtk->x[k] + C2 * lam[l] * NW[i]) / lam_NL;

        N1[m] = ROUND(B1[m]);
        var1 = (rtk->P[j + j * rtk->nx] + rtk->P[k + k * rtk->nx] - 2.0 * rtk->P[j + k * rtk->nx]) / SQR(lam_NL);

        /* validation of narrow-lane ambiguity conf_func(N1[m], B1[m], var1)<rtk->opt.thresar[1] ||*/
        if (fabs(N1[m] - B1[m]) > rtk->opt.thresar[3])
            continue; /*here maybe is the threshold for NL at 0.1, 0.25 is for WL*/

        /* narrow-lane ambiguity transformation matrix */
        D[j + m * rtk->nx] = 1.0 / lam_NL;
        D[k + m * rtk->nx] = -1.0 / lam_NL;

        sat1[m] = sat1[i];
        sat2[m] = sat2[i];
        NW[m++] = NW[i];
    }
    trace(2, " %d NL amb is close to integer\n", m);
    if (m < 3){
        return 0;
    }

    /* covariance of narrow-lane ambiguities */
    matmul("TN", m, rtk->nx, rtk->nx, 1.0, D, rtk->P, 0.0, E);
    matmul("NN", m, m, rtk->nx, 1.0, E, D, 0.0, Q);

    /* integer least square */
    if ((info = lambda(m, 2, B1, Q, N1, s)))
    {
        trace(2, "lambda error: info=%d\n", info);
        return 0;
    }
    if (s[0] <= 0.0)
        return 0;

    rtk->sol.ratio = (float)(MIN(s[1] / s[0], 999.9));

    /* varidation by ratio-test */
    if (rtk->opt.thresar[0] > 0.0 && rtk->sol.ratio < rtk->opt.thresar[0])
    {
        trace(2, "ratio= %f , do PAR\n",rtk->sol.ratio);
        del = ppp_amb_IFLAM_PAR(rtk, nav, m, m, sat1, sat2, NW, B1, N1, Q);
        if (del == 0)
        {
            free(B1);
            free(N1);
            free(D);
            free(E);
            free(Q);
            free(NC);
            trace(2, "PAR fail, del = 0\n");
            return 0;
        }
    }
    
    /*|| (m - del < 4 && rtk->opt.navsys > 1)*/
    if ((m - del < 3 && rtk->opt.navsys >= 1) || (m - del < 4 && rtk->opt.navsys > 1))
    {
        PARerror = 1;
    }
    if (((double)del / m >= 0.5))
        PARerror = 2;
    if ((PARerror == 1) || (PARerror == 2))
    {
        trace(2, "PAR not enough sat = %d\n", m-del);
        free(B1);
        free(N1);
        free(D);
        free(E);
        free(Q);
        free(NC);
        return 0;
    }
    trace(2, "varidation ok: %s n=%2d ratio=%8.3f\n", time_str(rtk->sol.time, 0), m - del,
          rtk->sol.ratio);

    /*output the fix between-satellites ambiguity(m) xiang 2015-11-14*/

    for (i = 0; i < m - del; i++)
    {
        j = IB(sat1[i], 0, &rtk->opt);
        k = IB(sat2[i], 0, &rtk->opt);
        /*fprintf(stderr, "$N1_fix %d,%.3f,%d,%3d,%3d,%.2f,%.2f,%.2f\n", week, tow, rtk->sol.stat,
                sat1[i], sat2[i], B1[i], N1[i + m], rtk->sol.ratio);*/
    }
    /*end of addition*/

    /* narrow-lane to iono-free ambiguity */
    for (i = 0; i < m - del; i++)
    {
        NC[i] = C1 * lam[0] * N1[i] + C2 * lam[l] * (N1[i] - NW[i]);
    }
    /* fixed solution */
    stat = fix_sol(rtk, sat1, sat2, NC, m - del);
    trace(2, "fix success %d %d\n", m, m - del);
    free(B1);
    free(N1);
    free(D);
    free(E);
    free(Q);
    free(NC);

    return stat;
}

static int ppp_amb_IFLC_LAMBDA(rtk_t *rtk, const obsd_t *obs, int n, int *exc,
                               const nav_t *nav, const double *azel, double *x,
                               double *P)
{
    const double *lam;
    ambc_t *amb1, *amb2;
    double lamN, lamW, lamE, C1, C2, Be, Bw, B1, varw, var1;
    double v[MAXOBS], R[MAXOBS * MAXOBS] = {0}, *H, *thres = rtk->opt.thresar;
    int i, j, k, l, info, ns, nw = 0, na = 0, Ne, Nw, N1, sat1[MAXOBS], sat2[MAXOBS], satN1[MAXOBS], satN2[MAXOBS], frq[MAXOBS];
    double tow, p;
    int week, stat = 0, *NW;
    tow = time2gpst(rtk->sol.time, &week);

    NW = imat(n * n, 1);
    /* filter EWL-ambiguity */
    filter_EWL(rtk, obs, n, exc, nav, azel);

    /* filter WL-ambiguity */
    filter_WL(rtk, obs, n, exc, nav, azel);

    /* generate satellite SD */
    if (!(ns = gen_sat_sd(rtk, obs, n, exc, azel, 0, sat1, sat2, frq)))
        return 0;

    H = zeros(rtk->nx, ns);

    for (i = 0; i < ns; i++)
    {
        lam = nav->lam[sat1[i] - 1];
        l = satsys(sat1[i], NULL) == SYS_GAL ? 2 : 1; /* L1/L2 or L1/L5 */
        lamE = lam[1] * lam[2] / (lam[2] - lam[1]);
        lamW = lam[0] * lam[l] / (lam[l] - lam[0]);
        lamN = lam[0] * lam[l] / (lam[l] + lam[0]);
        C1 = SQR(lam[l]) / (SQR(lam[l]) - SQR(lam[0]));
        C2 = -SQR(lam[0]) / (SQR(lam[l]) - SQR(lam[0]));

        j = IB(sat1[i], 0, &rtk->opt);
        k = IB(sat2[i], 0, &rtk->opt);
        amb1 = rtk->ambc + sat1[i] - 1;
        amb2 = rtk->ambc + sat2[i] - 1;
        if (!amb1->n[0] || !amb2->n[0])
            continue;

        if (amb1->n[1] && amb2->n[1])
        {
            Be = (amb1->LC[1] - amb2->LC[1]) / lamE;
            Ne = ROUND(Be);
#if 0
			trace(2, "%s sat=%2d-%2d:Be=%13.4f Fe=%7.4f\n", time_str(obs[0].time, 0),
				sat1[i], sat2[i], Be, Be - Ne);
#endif
        }
        /* round WL- and L1-ambiguity */
        Bw = (amb1->LC[0] - amb2->LC[0]) / lamW;
        Bw += nav->wlbias[sat1[i] - 1] - nav->wlbias[sat2[i] - 1]; /* correct WL-bias sat2 is ref */
        Nw = ROUND(Bw);
        B1 = (x[j] - x[k] + C2 * lam[l] * Nw) / lamN;

        varw = (amb1->LCv[0] + amb2->LCv[0]) / SQR(lamW);
        p = conf_func(Nw, Bw, varw);
        if (p < thres[1] || fabs(Nw - Bw) > thres[2])
            continue;
        NW[nw] = Nw;
        satN1[nw] = sat1[i];
        satN2[nw] = sat2[i];
        nw++;

        /* constraint to fixed LC-ambiguty */
        /* H[j + na*rtk->nx] = 1.0; */
        /* H[k + na*rtk->nx] = -1.0; */
        /* v[na++] = C1*lam[0] * N1 + C2*lam[l] * (N1 - Nw) - (rtk->x[j] - rtk->x[k]); */

        /* update fix flags */
        /*fprintf(stderr, "$NW_fix %d,%.3f,%d,%3d,%3d,%.2f,%d,%.2f,%.4f\n", week, tow, rtk->sol.stat,
                sat1[i], sat2[i], Bw, Nw, Bw - Nw, p);*/
    }
    trace(2, " %f WL fix success = %d\n", tow ,nw);
    if (nw < 4)
        return 0;

    /*fix narrowlane ambiguity using LAMBDA*/
    stat = fix_amb_ILS(rtk, nav, satN1, satN2, NW, nw);

    rtk->sol.age = (float)nw;        /* # of WL-fixed */
    rtk->sol.ratio = (float)rtk->na; /* # of NL-fixed rtk->sol.ratio;//*/

    return stat;
}
/* write debug trace for PAR -------------------------------------------------*/
static void write_trace1(rtk_t *rtk, const double *Z, const double *a,
                         const double *Q, int na, const int *sat1,
                         const int *sat2, const int *frq)
{
    const char freq[] = "125678";
    char buff[1024], s[32], *p = buff;
    int i, j;

    trace(2, "EPOCH=%s NFIX=%d\n", time_str(rtk->sol.time, 0), rtk->nfix);

    for (i = 0, p = buff; i < na; i++)
    {
        satno2id(sat1[i], s);
        p += sprintf(p, "%s ", s);
    }
    trace(2, "     %s          Z*a     STD\n", buff);
    for (i = 0, p = buff; i < na; i++)
    {
        satno2id(sat2[i], s);
        p += sprintf(p, "%s ", s);
    }
    trace(2, "     %s         (cyc)   (cyc)\n", buff);
    for (i = 0, p = buff; i < na; i++)
    {
        p += sprintf(p, "L%c  ", freq[frq[i]]);
    }
    trace(2, "      %s\n", buff);
    for (i = na - 1; i >= 0; i--)
    {
        p = buff;
        p += sprintf(p, "%3d: ", na - i);
        for (j = 0; j < na; j++)
            p += sprintf(p, "%3.0f ", Z[j + i * na]);
        p += sprintf(p, "%14.3f %7.3f", a[i], sqrt(Q[i + i * na]));
        trace(2, "%s\n", buff);
    }
    trace(2, "%3s: %7s %9s (%9s/%9s) [ N1 N2 ... NN ]\n", "FIX", "STD-POS", "RATIO",
          "S1", "S2");
}
static void write_trace2(rtk_t *rtk, const double *x, const double *P,
                         const double *a, const double *N, const double *D,
                         int na, const double *s)
{
    double *xp, *Pp, b[256], R[256 * 256] = {0}, std[3];
    char buff[1024], *p = buff;
    int i;

    xp = mat(rtk->nx, 1);
    Pp = mat(rtk->nx, rtk->nx);
    matcpy(xp, x, rtk->nx, 1);
    matcpy(Pp, P, rtk->nx, rtk->nx);
    for (i = 0; i < na; i++)
    {
        b[i] = N[i] - a[i];
        R[i + i * na] = SQR(CONST_AMB);
    }
    for (i = na - 1; i >= 0; i--)
    {
        p += sprintf(p, "%s%d", i == na - 1 ? "" : " ", (int)N[i]);
    }
    if (!filter(xp, Pp, D, b, R, rtk->nx, na))
    {
        for (i = 0; i < 3; i++)
            std[i] = sqrt(Pp[i + i * rtk->nx]);
        trace(2, "%3d: %7.3f %9.3f (%9.3f/%9.3f) [%s]\n", na, norm(std, 3),
              MIN(99999.999, s[1] / s[0]), s[0], s[1], buff);
    }
    free(xp);
    free(Pp);
}
/* decorrelate search space --------------------------------------------------*/
static int decorr_space(rtk_t *rtk, double *a, double *Q, double *D, int na,
                        const int *sat1, const int *sat2, const int *frq)
{
    double *W = mat(na, rtk->nx), *Z = eye(na);

    /* lambda reduction */
    if (lambda_reduction(na, Q, Z))
    {
        free(W);
        free(Z);
        return 0;
    }
    /* a=Z'*a, Q=Z'*Q*Z, D'=Z'*D' */
    matmul("TN", na, 1, na, 1.0, Z, a, 0.0, W);
    matcpy(a, W, na, 1);
    matmul("TN", na, na, na, 1.0, Z, Q, 0.0, W);
    matmul("NN", na, na, na, 1.0, W, Z, 0.0, Q);
    matmul("NN", rtk->nx, na, na, 1.0, D, Z, 0.0, W);
    matcpy(D, W, rtk->nx, na);

    /*if (strstr(rtk->opt.pppopt, "-TRACE_AR"))
    {
        write_trace1(rtk, Z, a, Q, na, sat1, sat2, frq);
    }*/
    write_trace1(rtk, Z, a, Q, na, sat1, sat2, frq);
    free(W);
    free(Z);
    return 1;
}
/* shrink search space -------------------------------------------------------*/
static int shrink_space(double *a, double *Q, double *H, int is, int na, int nx)
{
    int i, j, n = 0, *index = imat(na, 1);

    for (i = is; i < na; i++)
    {
        index[n++] = i;
    }
    for (i = 0; i < n; i++)
    {
        a[i] = a[index[i]];
        for (j = 0; j < n; j++)
            Q[j + i * n] = Q[index[j] + index[i] * na];
        matcpy(H + i * nx, H + index[i] * nx, nx, 1);
    }
    free(index);
    return n;
}
/* freq-ambiguity resolution by ILS ------------------------------------------*/
static int ppp_amb_ILS_FRQ(rtk_t *rtk, const obsd_t *obs, int n, const int *exc,
                           const nav_t *nav, const double *azel,
                           const double *x, const double *P, double *D,
                           double *a, double *N)
{
    const double *lam, thres_fact[] = {3.0, 2.5, 2.0, 1.5, 1.2};
    double *W, *Q, s[2] = {0}, thres = 0.0;
    int i, j, k, na = 0, sat1[MAXOBS * NFREQ], sat2[MAXOBS * NFREQ], frq[MAXOBS * NFREQ];
    int trace_AR = 0;

    if (strstr(rtk->opt.pppopt, "-TRACE_AR"))
        trace_AR = 1;

    /* generate satellite SD */
    for (i = 0; i < rtk->opt.nf; i++)
    {
        if (rtk->opt.freqopt == 1 && i == 1)
            continue;
        na += gen_sat_sd(rtk, obs, n, exc, azel, i, sat1 + na, sat2 + na, frq + na);
    }
    if (na <= 0)
        return 0;

    /* freq-ambiguity SD-matrix */
    for (i = 0; i < na; i++)
    {
        lam = nav->lam[sat1[i] - 1];
        j = IB(sat1[i], frq[i], &rtk->opt);
        k = IB(sat2[i], frq[i], &rtk->opt);
        D[j + i * rtk->nx] = 1.0 / lam[frq[i]];
        D[k + i * rtk->nx] = -1.0 / lam[frq[i]];
    }
    /* a=D'*x, Q=D'*P*D */
    W = mat(rtk->nx, na);
    Q = mat(na, na);
    matmul("TN", na, 1, rtk->nx, 1.0, D, x, 0.0, a);
    matmul("TN", na, rtk->nx, rtk->nx, 1.0, D, P, 0.0, W);
    matmul("NN", na, na, rtk->nx, 1.0, W, D, 0.0, Q);

    /* decorrelate search space */
    if (!decorr_space(rtk, a, Q, D, na, sat1, sat2, frq))
    {
        free(W);
        free(Q);
        return 0;
    }
    for (i = 0; i < rtk->opt.armaxiter && na >= MIN_AMB_RES; i++)
    {

        s[0] = s[1] = 0.0;

        /* integer least-square (a->N) */
        if (lambda_search(na, 2, a, Q, N, s) || s[0] <= 0.0)
        {
            free(W);
            free(Q);
            return 0;
        }
        if (trace_AR)
        {
            write_trace2(rtk, x, P, a, N, D, na, s);
        }
        thres = rtk->opt.thresar[0];
        if (na - MIN_AMB_RES < 5)
            thres *= thres_fact[na - MIN_AMB_RES];

        /* validation by ratio-test */
        if (s[1] / s[0] >= thres)
            break;

        /* shrink search space */
        na = shrink_space(a, Q, D, 1, na, rtk->nx);
    }
    free(W);
    free(Q);

    rtk->sol.ratio = s[0] == 0.0 ? 0.0 : (float)MIN(s[1] / s[0], 999.9);
    rtk->sol.thres = (float)thres;

    if (i >= rtk->opt.armaxiter || na < MIN_AMB_RES)
        return 0;

    /* update fix flags */
    for (i = 0; i < na; i++)
    {
        rtk->ssat[sat1[i] - 1].fix[frq[i]] = rtk->ssat[sat2[i] - 1].fix[frq[i]] = 2;
    }
    return na;
}
/* ambiguity resolution by ILS (integer-least-square) ------------------------*/
static int ppp_amb_ILS(rtk_t *rtk, const obsd_t *obs, int n, int *exc,
                       const nav_t *nav, const double *azel, double *x,
                       double *P)
{
    double *D, *R, a[MAXOBS * NFREQ], N[MAXOBS * NFREQ * 2];
    int i, na, info;

    D = zeros(rtk->nx, n * NF(&rtk->opt));

    /* freq-ambiguity resolution */
    if (!(na = ppp_amb_ILS_FRQ(rtk, obs, n, exc, nav, azel, x, P, D, a, N)))
    {
        free(D);
        return 0;
    }
    R = zeros(na, na);

    /* update states with integer ambiguity constraints */
    for (i = 0; i < na; i++)
    {
        a[i] = N[i] - a[i];
        R[i + i * na] = SQR(CONST_AMB);
    }
    if ((info = filter(x, P, D, a, R, rtk->nx, na)))
    {
        trace(1, "filter error (info=%d)\n", info);
    }
    free(D);
    free(R);
    return info ? 0 : 1;
}
/* carrier-phase bias correction by fcb --------------------------------------*/

static void corr_phase_bias_fcb(obsd_t *obs, int n, const nav_t *nav)
{
    int i, j, k;

    for (i = 0; i < nav->nf; i++)
    {
        if (timediff(nav->fcb[i].te, obs[0].time) < -1E-3)
            continue;
        if (timediff(nav->fcb[i].ts, obs[0].time) > 1E-3)
            break;
        for (j = 0; j < n; j++)
        {
            for (k = 0; k < NFREQ; k++)
            {
                if (obs[j].L[k] == 0.0)
                    continue;
                obs[j].L[k] -= nav->fcb[i].bias[obs[j].sat - 1][k];
            }
        }
        return;
    }
}
/* carrier-phase bias correction by ssr --------------------------------------*/
static void corr_phase_bias_ssr(obsd_t *obs, int n, const nav_t *nav)
{
    double lam;
    int i, j, code, pbias_num = 0, pbias_satnum = 0, pbias_num_sat = 0;
    char id[32];
    for (i = 0; i < n; i++)
    {
        pbias_num_sat = pbias_num;
        for (j = 0; j < NFREQ; j++)
        {

            if (!(code = obs[i].code[j]))
                continue;
            if ((lam = nav->lam[obs[i].sat - 1][j]) == 0.0)
                continue;

            /* correct phase bias (cyc) */
            obs[i].L[j] -= nav->ssr[obs[i].sat - 1].pbias[code - 1] / lam;
            satno2id(obs[i].sat, id);
            if (nav->ssr[obs[i].sat - 1].pbias[code - 1] != 0)
            {
                /*trace(2, "%s %d carrier phase corrected =%f %f\n", id, code, nav->ssr[obs[i].sat - 1].pbias[code - 1]/lam, obs[i].L[j]);*/
                pbias_num++;
            }
        }
        if (pbias_num > pbias_num_sat)
            pbias_satnum++;
    }
    /*trace(2, "%d sats %d freq pbias corrected\n", pbias_satnum, pbias_num);*/
}
/* ambiguity resolution in ppp -----------------------------------------------*/
extern int ppp_ar(rtk_t *rtk, const obsd_t *obs, int n, int *exc,
                  const nav_t *nav, const double *azel, double *x, double *P)
{
    obsd_t *obs_tmpt = {0};

    if (n <= 0 || rtk->opt.modear < ARMODE_CONT)
        return 0;
    /* carrier-phase bias correction */
    if (nav->nf >= 0)
    {
        corr_phase_bias_ssr(obs, n, nav);
    }

    /* ambiguity resolution by WL/NL for iono-free LC */
    if (rtk->opt.ionoopt == IONOOPT_IFLC)
    {
        return ppp_amb_IFLC_LAMBDA(rtk, obs, n, exc, nav, azel, x, P);
    }
    /* ambiguity resolution by ILS */
    else
    {
        return ppp_amb_ILS(rtk, obs, n, exc, nav, azel, x, P);
    }
}
