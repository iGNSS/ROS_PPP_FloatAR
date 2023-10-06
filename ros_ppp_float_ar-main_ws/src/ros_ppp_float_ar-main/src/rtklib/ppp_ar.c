/*------------------------------------------------------------------------------
* ppp_ar.c : ppp ambiguity resolution
*
* options : -DREV_WL_FCB reversed polarity of WL FCB
*
* reference :
*    [1] H.Okumura, C-gengo niyoru saishin algorithm jiten (in Japanese),
*        Software Technology, 1991
*
*          Copyright (C) 2012-2013 by T.TAKASU, All rights reserved.
*
* version : $Revision:$ $Date:$
* history : 2013/03/11 1.0  new
*-----------------------------------------------------------------------------*/
#include "rtklib.h"

static const char rcsid[]="$Id:$";

/* constants/macros ----------------------------------------------------------*/

#define MIN_ARC_GAP     300.0       /* min arc gap (s) */
#define CONST_AMB       0.001       /* constraint to fixed ambiguity */
#define THRES_RES       0.3         /* threashold of residuals test (m) */
#define LOG_PI          1.14472988584940017 /* log(pi) */
#define SQRT2           1.41421356237309510 /* sqrt(2) */

#define MIN(x,y)        ((x)<(y)?(x):(y))
#define SQR(x)          ((x)*(x))
#define ROUND(x)        (int)floor((x)+0.5)

#define SWAP_I(x,y)     do {int _z=x; x=y; y=_z;} while (0)
#define SWAP_D(x,y)     do {double _z=x; x=y; y=_z;} while (0)

#define NP(opt)     ((opt)->dynamics?9:3) /* number of pos solution */
#define IC(s,opt)   (NP(opt)+(s))      /* state index of clocks (s=0:gps,1:glo) */
#define IT(opt)     (IC(0,opt)+NSYS)   /* state index of tropos */
#define NR(opt)     (IT(opt)+((opt)->tropopt<TROPOPT_EST?0:((opt)->tropopt==TROPOPT_EST?1:3)))
                                       /* number of solutions */
#define IB(s,opt)   (NR(opt)+(s)-1)    /* state index of phase bias */

/* select narrow-lane FCB --------------------------------------------------------*/
//static double selfcb(gtime_t time, int sat,const nav_t *nav)
//{
//	double t,tmax,tmin;
//	int i,j=-1,sys,sel=0;
//
//	sys=satsys(sat,NULL);
//	tmax=600;// FCB sample==10 min
//	tmin=tmax+1.0;
//
//	for (i=0;i<nav->nf;i++) {
//        if ((t = fabs(timediff(nav->fcb[i].ts, time))) > tmax) {
//            continue;
//        }
//		if (t<=tmin) {
//            j=i; tmin=t;
//        } /* toe closest to time */
//	}
//	if (j<0) {
//        return 10;
//    }
//	return nav->fcb[j].bias[sat-1][0];
//}
/* wave length of LC (m) -----------------------------------------------------*/
static double lam_LC(int i, int j, int k)
{
    const double f1=FREQ1,f2=FREQ2,f5=FREQ5;
    
    return CLIGHT/(i*f1+j*f2+k*f5);
}
/* wave length of LC (m) -----------------------------------------------------*/
static double lam_LC_CMP(int i, int j, int k)
{
    const double f1=FREQ1_CMP,f2=FREQ2_CMP,f5=FREQ3_CMP;
    
    return CLIGHT/(i*f1+j*f2+k*f5);
}
/* carrier-phase LC (m) ------------------------------------------------------*/
static double L_LC(int i, int j, int k, const double *L)
{
    const double f1=FREQ1,f2=FREQ2,f5=FREQ5;
    double L1,L2,L5;
    
    if ((i&&!L[0])||(j&&!L[1])||(k&&!L[2])) {
        return 0.0;
    }
    L1=CLIGHT/f1*L[0];
    L2=CLIGHT/f2*L[1];
    L5=CLIGHT/f5*L[2];
    return (i*f1*L1+j*f2*L2+k*f5*L5)/(i*f1+j*f2+k*f5);
}
/* carrier-phase LC (m) ------------------------------------------------------*/
static double L_LC_CMP(int i, int j, int k, const double *L)
{
    const double f1= FREQ1_CMP,f2=FREQ2_CMP,f5=FREQ3_CMP;
    double L1,L2,L5;
    
    if ((i&&!L[0])||(j&&!L[1])||(k&&!L[2])) {
        return 0.0;
    }
    L1=CLIGHT/f1*L[0];
    L2=CLIGHT/f2*L[1];
    L5=CLIGHT/f5*L[2];
    return (i*f1*L1+j*f2*L2+k*f5*L5)/(i*f1+j*f2+k*f5);
}
/* pseudorange LC (m) --------------------------------------------------------*/
static double P_LC(int i, int j, int k, const double *P)
{
    const double f1=FREQ1,f2=FREQ2,f5=FREQ5;
    double P1,P2,P5;
    
    if ((i&&!P[0])||(j&&!P[1])||(k&&!P[2])) {
        return 0.0;
    }
    P1=P[0];
    P2=P[1];
    P5=P[2];
    return (i*f1*P1+j*f2*P2+k*f5*P5)/(i*f1+j*f2+k*f5);
}
/* pseudorange LC (m) --------------------------------------------------------*/
static double P_LC_CMP(int i, int j, int k, const double *P)
{
    const double f1=FREQ1_CMP,f2=FREQ2_CMP,f5=FREQ3_CMP;
    double P1,P2,P5;
    
    if ((i&&!P[0])||(j&&!P[1])||(k&&!P[2])) {
        return 0.0;
    }
    P1=P[0];
    P2=P[1];
    P5=P[2];
    return (i*f1*P1+j*f2*P2+k*f5*P5)/(i*f1+j*f2+k*f5);
}
/* noise variance of LC (m) --------------------------------------------------*/
static double var_LC(int i, int j, int k, double sig)
{
    const double f1=FREQ1,f2=FREQ2,f5=FREQ5;
    
    return (SQR(i*f1)+SQR(j*f2)+SQR(k*f5))/SQR(i*f1+j*f2+k*f5)*SQR(sig);
}
/* noise variance of LC (m) --------------------------------------------------*/
static double var_LC_CMP(int i, int j, int k, double sig)
{
    const double f1=FREQ1_CMP,f2=FREQ2_CMP,f5=FREQ3_CMP;
    
    return (SQR(i*f1)+SQR(j*f2)+SQR(k*f5))/SQR(i*f1+j*f2+k*f5)*SQR(sig);
}
/* complementaty error function (ref [1] p.227-229) --------------------------*/
static double q_gamma(double a, double x, double log_gamma_a);
static double p_gamma(double a, double x, double log_gamma_a)
{
    double y,w;
    int i;
    
    if (x==0.0) return 0.0;
    if (x>=a+1.0) return 1.0-q_gamma(a,x,log_gamma_a);
    
    y=w=exp(a*log(x)-x-log_gamma_a)/a;
    
    for (i=1;i<100;i++) {
        w*=x/(a+i);
        y+=w;
        if (fabs(w)<1E-15) break;
    }
    return y;
}
static double q_gamma(double a, double x, double log_gamma_a)
{
    double y,w,la=1.0,lb=x+1.0-a,lc;
    int i;
    
    if (x<a+1.0) return 1.0-p_gamma(a,x,log_gamma_a);
    w=exp(-x+a*log(x)-log_gamma_a);
    y=w/lb;
    for (i=2;i<100;i++) {
        lc=((i-1-a)*(lb-la)+(i+x)*lb)/i;
        la=lb; lb=lc;
        w*=(i-1-a)/i;
        y+=w/la/lb;
        if (fabs(w/la/lb)<1E-15) break;
    }
    return y;
}
static double f_erfc(double x)
{
    return x>=0.0?q_gamma(0.5,x*x,LOG_PI/2.0):1.0+p_gamma(0.5,x*x,LOG_PI/2.0);
}
/* confidence function of integer ambiguity ----------------------------------*/
static double conffunc(int N, double B, double sig)
{
    double x,p=1.0;
    int i;
    
    x=fabs(B-N);
    for (i=1;i<8;i++) {
        p-=f_erfc((i-x)/(SQRT2*sig))-f_erfc((i+x)/(SQRT2*sig));
    }
    return p;
}
/* average LC ----------------------------------------------------------------*/
static void average_LC(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav,
                       const double *azel)
{
    ambc_t *amb;
	double *lam;
    double LC1,LC2,LC3,var1,var2,var3,sig,eratio,P[NFREQ],L[NFREQ],dmp[3]={0.0},C1,C2;
    if (n > 12) {
        int zyf = 1;
    }
    int i,j,sat,sys,index,jj;
    
    for (i=0;i<n;i++) {
        sat=obs[i].sat;index=0;
		lam=nav->lam[sat-1];
		sys=satsys(sat,NULL);

        
        if (azel[1+2*i]<rtk->opt.elmin) continue;
        
        //if (satsys(sat,NULL)!=SYS_GPS) continue;// support multi-GNSS
        
		/* P1-C1,P2-C2 dcb correction (C1->P1,C2->P2) */
		for(j=0;j<NFREQ;j++){
			P[j]=obs[i].P[j];L[j]=obs[i].L[j];
			if (obs[i].code[j]==CODE_L1C) {
			P[j]+=nav->cbias[obs[i].sat-1][1];}
		}

		/* remove obs when data is lost */
		jj=(sys&(SYS_GAL|SYS_SBS|SYS_CMP))?2:1;
		if(!P[0]||!P[jj]||!L[0]||!L[jj]) index=1;

        /* BDS2 code pseudorange multipath*/
        if (sys == SYS_CMP) {
            BDmultipath(obs+i, azel[1+2*i], dmp);
            P[0] += dmp[0];
            P[jj] += dmp[jj];
        }
		/* P1-P2 dcb correction (P1->Pc,P2->Pc) */
	    C1= SQR(lam[jj])/(SQR(lam[jj])-SQR(lam[0]));
        C2=-SQR(lam[0])/(SQR(lam[jj])-SQR(lam[0]));
        P[0]-=C2*nav->cbias[obs[i].sat-1][0];
        P[jj]+=C1*nav->cbias[obs[i].sat-1][0];
        /* triple-freq carrier and code LC (m) */
        if(sys!=SYS_CMP){
            LC1=L_LC(1,-1, 0,obs[i].L)-P_LC(1,1,0,P);// MW Combination
		    LC2=L_LC(1,0, -1,obs[i].L)-P_LC(1,0,1,P);// L1/L5 for galileo
        }
        else {
            LC1=L_LC_CMP(1,-1, 0,obs[i].L)-P_LC_CMP(1,1,0,P);// MW Combination
		    LC2=L_LC_CMP(1,0, -1,obs[i].L)-P_LC_CMP(1,0,1,P);// add L1/L5 
        }

       /* LC2=L_LC(0, 1,-1,obs[i].L)-P_LC(0,1,1,P);*/
		LC3=L_LC(1,-6, 5,obs[i].L)-P_LC(1,1,0,P);

        sig=sqrt(SQR(rtk->opt.err[1])+SQR(rtk->opt.err[2]/sin(azel[1+2*i])));
        
        /* measurement noise variance (m) */
		eratio=rtk->opt.eratio[0];
        if(sys!=SYS_CMP){
            var1=var_LC(1,1,0,sig*eratio);//sig*rtk->opt.eratio[0]
		    var2=var_LC(1,0,1,sig*eratio);// add L1/L5 for galileo
           //var2=var_LC(0,1,1,sig*eratio);//sig*rtk->opt.eratio[0]
            var3=var_LC(1,1,0,sig*eratio);//sig*rtk->opt.eratio[0]
        }
        else {
            var1=var_LC_CMP(1,1,0,sig*eratio);//sig*rtk->opt.eratio[0]
		    var2=var_LC_CMP(1,0,1,sig*eratio);// add L1/L5
            var3=var_LC_CMP(1,1,0,sig*eratio);//sig*rtk->opt.eratio[0]
        }

        amb=rtk->ambc+sat-1;
		if (satsys(obs[i].sat, NULL) == SYS_GPS){
			amb->L4 = obs[i].L[0] * lam[0] - obs[i].L[1] * lam[1];
		}
		else
		{
			amb->L4 = obs[i].L[0] * lam[0] - obs[i].L[2] * lam[2];
		}

		// ��ʼ��
        if (fabs(timediff(amb->epoch[0],obs[0].time))>MIN_ARC_GAP) {
            
            amb->n[0]=amb->n[1]=amb->n[2]=0.0;
            amb->LC[0]=amb->LC[1]=amb->LC[2]=0.0;
            amb->LCv[0]=amb->LCv[1]=amb->LCv[2]=0.0;
            amb->fixcnt=0;
            for (j=0;j<MAXSAT;j++) amb->flags[j]=0;
        }

		// ����̽��
        if (rtk->ssat[sat-1].slip[0]||rtk->ssat[sat-1].slip[1]||
            rtk->ssat[sat-1].slip[2]||index) {
            
            amb->n[0]=amb->n[1]=amb->n[2]=0.0;
            amb->LC[0]=amb->LC[1]=amb->LC[2]=0.0;
            amb->LCv[0]=amb->LCv[1]=amb->LCv[2]=0.0;
            amb->fixcnt=0;
            for (j=0;j<MAXSAT;j++) amb->flags[j]=0;
			continue;
        }

	
        /* averaging */
        if (LC1) {
            amb->n[0]+=1.0;
            amb->LC [0]+=(LC1 -amb->LC [0])/amb->n[0];
            amb->LCv[0]+=(var1-amb->LCv[0])/amb->n[0];
			//trace(2,"%s  %02d  %8.3f %5d\n",time_str(rtk->sol.time,0), sat, amb->LC[0], amb->n[0]);
        }
        if (LC2) {
            amb->n[1]+=1.0;
            amb->LC [1]+=(LC2 -amb->LC [1])/amb->n[1];
            amb->LCv[1]+=(var2-amb->LCv[1])/amb->n[1];
			//trace(2, "%s  %02d  %8.3f\n", time_str(rtk->sol.time, 0), sat, amb->LC[1]);
        }
        if (LC3) {
            amb->n[2]+=1.0;
            amb->LC [2]+=(LC3 -amb->LC [2])/amb->n[2];
            amb->LCv[2]+=(var3-amb->LCv[2])/amb->n[2];
        }
        amb->epoch[0]=obs[0].time;
    }
}
/* fix wide-lane ambiguity ---------------------------------------------------*/
static int fix_amb_WL(rtk_t *rtk, const nav_t *nav, int sat1, int sat2, int *NW)
{
    ambc_t *amb1,*amb2;
    double BW=0.0,vW=0.0,lam_WL=0.0;
    char id[64];
    char idd[64];
    satno2id(sat1, id); satno2id(sat2, idd);
    
	if(satsys(sat1, NULL)==SYS_GPS){
        lam_WL=lam_LC(1,-1,0);
    }  
	else if(satsys(sat1,NULL)==SYS_GAL){
        lam_WL=lam_LC(1,0,-1);// support GAL
    }
    else if(satsys(sat1,NULL)==SYS_CMP){
        lam_WL=lam_LC_CMP(1,0,-1);// support BDS
    }
		
    amb1=rtk->ambc+sat1-1;
    amb2=rtk->ambc+sat2-1;

    /* amb1->n can be set to 20 to constrain */
    if ((!amb1->n[0]||!amb2->n[0])&&satsys(sat1, NULL)==SYS_GPS) return 0;
	if ((!amb1->n[1]||!amb2->n[1])&&satsys(sat1, NULL)==SYS_GAL) return 0;
    if ((!amb1->n[1]||!amb2->n[1])&&satsys(sat1, NULL)==SYS_CMP) return 0;


    /* wide-lane ambiguity */
#ifndef REV_WL_FCB
    if(satsys(sat1, NULL)==SYS_GPS){
        BW=(amb1->LC[0]-amb2->LC[0])/lam_WL;
        }
	else if(satsys(sat1,NULL)==SYS_GAL){
        BW=(amb1->LC[1]-amb2->LC[1])/lam_WL;
    }
    else if(satsys(sat1,NULL)==SYS_CMP){
        BW=(amb1->LC[1]-amb2->LC[1])/lam_WL;
    }
#else
    BW=(amb1->LC[0]-amb2->LC[0])/lam_WL;
#endif
    *NW=ROUND(BW);
    
    /* variance of wide-lane ambiguity */
	if(satsys(sat1, NULL)==SYS_GPS){
        vW=(amb1->LCv[0]/amb1->n[0]+amb2->LCv[0]/amb2->n[0])/SQR(lam_WL);
    }   
    else if(satsys(sat1,NULL)==SYS_GAL){
        vW=(amb1->LCv[1]/amb1->n[1]+amb2->LCv[1]/amb2->n[1])/SQR(lam_WL);
    }
    else if(satsys(sat1,NULL)==SYS_CMP){
        vW=(amb1->LCv[1]/amb1->n[1]+amb2->LCv[1]/amb2->n[1])/SQR(lam_WL);
    }
    /* validation of integer wide-lane ambigyity */
    int WLfixflag = (fabs(*NW - BW) <= rtk->opt.thresar[2] && conffunc(*NW, BW, sqrt(vW)) >= rtk->opt.thresar[1]);
    if (WLfixflag == 0) {
        trace(2, "WL fix error, fabs=%f, conffunc=%f %s %s \n)", fabs(*NW - BW), conffunc(*NW, BW, sqrt(vW)), id ,idd);
    }

    return WLfixflag;
}
/* linear dependency check ---------------------------------------------------*/
static int is_depend(int sat1, int sat2, int *flgs, int *max_flg)
{
    int i;
    
    if (flgs[sat1-1]==0&&flgs[sat2-1]==0) {
        flgs[sat1-1]=flgs[sat2-1]=++(*max_flg);
    }
    else if (flgs[sat1-1]==0&&flgs[sat2-1]!=0) {
        flgs[sat1-1]=flgs[sat2-1];
    }
    else if (flgs[sat1-1]!=0&&flgs[sat2-1]==0) {
        flgs[sat2-1]=flgs[sat1-1];
    }
    else if (flgs[sat1-1]>flgs[sat2-1]) {
        for (i=0;i<MAXSAT;i++) if (flgs[i]==flgs[sat2-1]) flgs[i]=flgs[sat1-1];
    }
    else if (flgs[sat1-1]<flgs[sat2-1]) {
        for (i=0;i<MAXSAT;i++) if (flgs[i]==flgs[sat1-1]) flgs[i]=flgs[sat2-1];
    }
    else {
        int zyf = 1;
        return 0;
    }/* linear depenent */
    return 1;
}
/* select fixed ambiguities --------------------------------------------------*/
static int sel_amb(int *sat1, int *sat2, double *N, double *var, int n)
{
    int i,j,flgs[MAXSAT]={0},max_flg=0;
    
    /* sort by variance */
    for (i=0;i<n;i++) for (j=1;j<n-i;j++) {
        if (var[j]>=var[j-1]) continue;
        SWAP_I(sat1[j],sat1[j-1]);
        SWAP_I(sat2[j],sat2[j-1]);
        SWAP_D(N[j],N[j-1]);
        SWAP_D(var[j],var[j-1]);
    }
    /* select linearly independent satellite pair */
    for (i=j=0;i<n;i++) {
        if (!is_depend(sat1[i],sat2[i],flgs,&max_flg)) continue;
        sat1[j]=sat1[i];
        sat2[j]=sat2[i];
        N[j]=N[i];
        var[j++]=var[i];
    }
    return j;
}
/* fixed solution ------------------------------------------------------------*/
static int fix_sol(rtk_t *rtk, const int *sat1, const int *sat2,
                   const double *NC, int n)
{
    double *v,*H,*R;
    int i,j,k,info,sys;
    
    if (n<=0) return 0;
    
    v=zeros(n,1); H=zeros(rtk->nx,n); R=zeros(n,n);
    
    /* constraints to fixed ambiguities */
    for (i=0;i<n;i++) {
        j=IB(sat1[i],&rtk->opt);
        k=IB(sat2[i],&rtk->opt);
        v[i]=NC[i]-(rtk->x[j]-rtk->x[k]);
        H[j+i*rtk->nx]= 1.0;
        H[k+i*rtk->nx]=-1.0;
        R[i+i*n]=SQR(CONST_AMB);
    }
    /* update states with constraints */
    if ((info=filter(rtk->x,rtk->P,H,v,R,rtk->nx,n))) {
        trace(1,"filter error (info=%d)\n",info);
        free(v); free(H); free(R);
        return 0;
    }
    /* set solution */
    for (i=0;i<rtk->na;i++) {
        rtk->xa[i]=rtk->x[i];
        for (j=0;j<rtk->na;j++) {
            rtk->Pa[i+j*rtk->na]=rtk->Pa[j+i*rtk->na]=rtk->P[i+j*rtk->nx];
        }
    }
    /* set flags */
    for (i=0;i<n;i++) {
        rtk->ambc[sat1[i]-1].flags[sat1[i]-1]=1;
        rtk->ambc[sat2[i]-1].flags[sat2[i]-1]=1;
    }

    /* record SD ambiguity number with AR */
    rtk->sol.gpsns=0;
    rtk->sol.galns=0;
    rtk->sol.bdsns=0;
    for (i=0;i<n;i++) {
        sys=satsys(sat1[i],NULL);
        if (sys==1) rtk->sol.gpsns++;
        if (sys==8) rtk->sol.galns++;
        if (sys==32) rtk->sol.bdsns++;
    }
    free(v); free(H); free(R);
    return 1;
}
/* fix narrow-lane ambiguity by rounding -------------------------------------*/
static int fix_amb_ROUND(rtk_t *rtk, int *sat1, int *sat2, const int *NW, int n)
{
    double C1,C2,B1,v1,BC,v,vc,*NC,*var,lam_NL=lam_LC(1,1,0),lam1,lam2;
    int i,j,k,m=0,N1,stat;
    
    lam1=lam_carr[0]; lam2=lam_carr[1];
    
    C1= SQR(lam2)/(SQR(lam2)-SQR(lam1));
    C2=-SQR(lam1)/(SQR(lam2)-SQR(lam1));
    
    NC=zeros(n,1); var=zeros(n,1);
    
    for (i=0;i<n;i++) {
        j=IB(sat1[i],&rtk->opt);
        k=IB(sat2[i],&rtk->opt);
        
        /* narrow-lane ambiguity */
        B1=(rtk->x[j]-rtk->x[k]+C2*lam2*NW[i])/lam_NL;
        N1=ROUND(B1);
        
        /* variance of narrow-lane ambiguity */
        var[m]=rtk->P[j+j*rtk->nx]+rtk->P[k+k*rtk->nx]-2.0*rtk->P[j+k*rtk->nx];
        v1=var[m]/SQR(lam_NL);
        
        /* validation of narrow-lane ambiguity */
        if (fabs(N1-B1)>rtk->opt.thresar[2]||
            conffunc(N1,B1,sqrt(v1))<rtk->opt.thresar[1]) {
            continue;
        }
        /* iono-free ambiguity (m) */
        BC=C1*lam1*N1+C2*lam2*(N1-NW[i]);
        
        /* check residuals */
        v=rtk->ssat[sat1[i]-1].resc[0]-rtk->ssat[sat2[i]-1].resc[0];
        vc=v+(BC-(rtk->x[j]-rtk->x[k]));
        if (fabs(vc)>THRES_RES) continue;
        
        sat1[m]=sat1[i];
        sat2[m]=sat2[i];
        NC[m++]=BC;
    }
    /* select fixed ambiguities by dependancy check */
    m=sel_amb(sat1,sat2,NC,var,m);
    
    /* fixed solution */
    stat=fix_sol(rtk,sat1,sat2,NC,m);
    
    free(NC); free(var);
    
    return stat&&m>=3;
}
/* ppp partial ambiguity resolutions ------------------------------------------*/
static int ppp_amb_IFLAM_PAR(rtk_t *rtk, const nav_t *nav, int m, int n, int *sat1, int *sat2, int *Nw, double *z, double *N1, const double *Q) {
	const double *lam;
	double lamN, lamW;
	double *D, *E, *B1, *Q_, s[2];
	int i, j, k, l, info, del, p;
	double var[MAXOBS];
	char sys[32];
    double boot=0.0;

	/* init */
	del = 1; p = m - del;

	/* extract variance */
	for (i = 0; i < m; i++) {
		var[i] = Q[i + i * m];
	}

	/* sort by variance */
	for (i = 0; i < m; i++)for (j = 1; j < m - i; j++) {
		if (var[j] >= var[j - 1])continue;
		SWAP_I(sat1[j], sat1[j - 1]);
		SWAP_I(sat2[j], sat2[j - 1]);
		SWAP_I(Nw[j], Nw[j - 1]);
		SWAP_D(N1[j], N1[j - 1]);
		SWAP_D(z[j], z[j - 1]);
		SWAP_D(var[j], var[j - 1]);
	}

	/* iteration */
	while (p > 4) {
		D = zeros(rtk->nx, m);
		E = mat(m, rtk->nx);
		B1 = zeros(p, 1);
		//N1 = zeros(p, 2); //must remove
		Q_ = mat(p, p);

		p = m - del;
		for (i = 0; i < p; i++) {
		    if(satsys(sat1[i], NULL)==SYS_GPS)
		    {lamN=lam_LC(1,1,0);lamW=lam_LC(1,-1,0);}
		    else if(satsys(sat1[i],NULL)==SYS_GAL) // support galileo
		    {lamN=lam_LC(1,0,1);lamW=lam_LC(1,0,-1);}
            else if(satsys(sat1[i],NULL)==SYS_CMP) // support BDS
		    {lamN=lam_LC_CMP(1,0,1);lamW=lam_LC_CMP(1,0,-1);}
		
		    j=IB(sat1[i],&rtk->opt);
            k=IB(sat2[i],&rtk->opt);

		    /* narrow-line ambiguity transformation matrix */
		    D[j + i * rtk->nx] = 1.0 / lamN;
		    D[k + i * rtk->nx] = -1.0 / lamN;
		    B1[i] = z[i];
		    N1[i] = ROUND(z[i]);
		}
		 /* covariance of narrow-lane ambiguities */
         matmul("TN",p,rtk->nx,rtk->nx,1.0,D,rtk->P,0.0,E);
         matmul("NN",p,p,rtk->nx,1.0,E,D,0.0,Q_);

		 //if((info=lambda(p,2,B1,Q_,N1,s))||(s[0]<=0.0)){
			// del++;
		 //}
         if ((info = mlambda(p, 2, B1, Q_, N1, s, &boot)) || (s[0] <= 0.0)) {
             del++;
         }
		 else {
			 rtk->sol.ratio=(float)(MIN(s[1]/s[0],999.9));
             if (rtk->sol.ratio < rtk->opt.thresar[0] * WL_ratio) {
                 trace(2, "ratio test error: %f %f ratio= %f ,threshold= %f \n", s[1],s[0],rtk->sol.ratio, rtk->opt.thresar[0] * WL_ratio);
                 del++;
             }
             else {
                 
                 free(D);free(E);free(B1);free(Q_);
                 return del;
             }  
		 }
		 free(D);free(E);free(B1);free(Q_);
    }
	return 0;
}
/* ppp partial ambiguity resolutions ------------------------------------------*/
static int ppp_amb_IFLAM_ELEPAR(rtk_t *rtk, const nav_t *nav, int m, int n, int *sat1, int *sat2, int *Nw, double *z, double *N1, const double *Q) {
	const double *lam;
	double lamN, lamW;
	double *D, *E, *B1, *Q_, s[2];
	int i, j, k, l, info, del, p;
	double ele[MAXOBS];
	char sys[32];

	/* init */
	del=1;p=m-del;

	/* extract elevation angle */
	for (i=0;i<m; i++) {
		ele[i]= rtk->ssat[sat1[i]-1].azel[1];
	}

	/* sort by variance */
	for (i=0;i<m;i++)for (j=1; j<m-i;j++) {
		if (ele[j]<=ele[j-1])continue;
		SWAP_I(sat1[j], sat1[j-1]);
		SWAP_I(sat2[j], sat2[j-1]);
		SWAP_I(Nw[j], Nw[j-1]);
		SWAP_D(N1[j], N1[j-1]);
		SWAP_D(z[j], z[j-1]);
		SWAP_D(ele[j], ele[j-1]);
	}

	/* iteration */
	while (p > 4) {
		D = zeros(rtk->nx, m);
		E = mat(m, rtk->nx);
		B1 = zeros(p, 1);
		//N1 = zeros(p, 2); //must remove
		Q_ = mat(p, p);

		p = m-del;
		for (i = 0; i < p; i++) {
		    if(satsys(sat1[i], NULL)==SYS_GPS)
		    {lamN=lam_LC(1,1,0);lamW=lam_LC(1,-1,0);}
		    else if(satsys(sat1[i],NULL)==SYS_GAL) // support galileo
		    {lamN=lam_LC(1,0,1);lamW=lam_LC(1,0,-1);}
            else if(satsys(sat1[i],NULL)==SYS_CMP) // support BDS
		    {lamN=lam_LC_CMP(1,0,1);lamW=lam_LC_CMP(1,0,-1);}
		
		    j=IB(sat1[i],&rtk->opt);
            k=IB(sat2[i],&rtk->opt);

		    /* narrow-line ambiguity transformation matrix */
		    D[j + i * rtk->nx] = 1.0 / lamN;
		    D[k + i * rtk->nx] = -1.0 / lamN;
		    B1[i] = z[i];
		    N1[i] = ROUND(z[i]);
		}
		 /* covariance of narrow-lane ambiguities */
         matmul("TN",p,rtk->nx,rtk->nx,1.0,D,rtk->P,0.0,E);
         matmul("NN",p,p,rtk->nx,1.0,E,D,0.0,Q_);

		 if((info=lambda(p,2,B1,Q_,N1,s))||(s[0]<=0.0)){
			 del++;
		 }
		 else {
			 rtk->sol.ratio=(float)(MIN(s[1]/s[0],999.9));
			 if(rtk->sol.ratio<rtk->opt.thresar[0]) del++;
             else {
                 free(D);free(E);free(B1);free(Q_);
                 return del;
             }  
		 }
		 free(D);free(E);free(B1);free(Q_);
    }
	return 0;
}
/* fix narrow-lane ambiguity by ILS ------------------------------------------*/
static int fix_amb_ILS(rtk_t *rtk, const obsd_t *obs, int nn,const nav_t *nav, int *sat1, int *sat2, int *NW, int n)
{
    double C1,C2,*B1,*N1,*NC,*D,*E,*Q,s[2],*NN1,lam_NL,lam1,lam2,var,v1,std[3],ep[6],lam_WL,L1fcb,L2fcb,ion,nlfcb,wlfcb,L4;
    //double sat1_nlfcb, sat2_nlfcb;
    int i,j,k,m=0,info,stat=0,flgs[MAXSAT]={0},max_flg=0,del=0,satG=0,satE=0,satC=0;
	ambc_t *amb1, *amb2;
	aug_t AUG_S[MAXSAT] = {0.0};
	char s0[4];
    double boot=0.0;
    
    B1=zeros(n,1); N1=zeros(n,2); D=zeros(rtk->nx,n); E=mat(n,rtk->nx);NN1=zeros(n,2);
    Q=mat(n,n); NC=mat(n,1);
    
    for (i=0;i<n;i++) {

         lam1 = lam_carr[0];
		 if(satsys(sat1[i], NULL)==SYS_GPS)
		 {lam2=lam_carr[1];lam_NL=lam_LC(1,1,0);}
		 else if(satsys(sat1[i],NULL)==SYS_GAL) 
		 {lam2=lam_carr[2];lam_NL=lam_LC(1,0,1);}
         else if(satsys(sat1[i],NULL)==SYS_CMP)
		 {lam1=lam_carr_CMP[0];lam2=lam_carr_CMP[2];lam_NL=lam_LC_CMP(1,0,1);}
    
         C1= SQR(lam2)/(SQR(lam2)-SQR(lam1));
         C2=-SQR(lam1)/(SQR(lam2)-SQR(lam1));
        
        /* check linear independency */
        if (!is_depend(sat1[i],sat2[i],flgs,&max_flg)) continue;
        
        j=IB(sat1[i],&rtk->opt);
        k=IB(sat2[i],&rtk->opt);
       

        /* float narrow-lane ambiguity (cycle) */
        B1[m]=(rtk->x[j]-rtk->x[k]+C2*lam2*NW[i])/lam_NL;
        N1[m]=ROUND(B1[m]);
        
        /* validation of narrow-lane ambiguity */
        if (fabs(N1[m]-B1[m])>rtk->opt.thresar[2]) continue;
        
        /* narrow-lane ambiguity transformation matrix */
        D[j+m*rtk->nx]= 1.0/lam_NL;
        D[k+m*rtk->nx]=-1.0/lam_NL;
        
        sat1[m]=sat1[i];
        sat2[m]=sat2[i];
        NW[m++]=NW[i];
    }
    if ((m<3 && rtk->opt.navsys==1)||(m<4 && rtk->opt.navsys>1)) {
        free(B1); free(N1); free(D); free(E); free(Q); free(NC); free(NN1);
        return 0;//m��СΪ3,support GAL+BDS
    }

    /* covariance of narrow-lane ambiguities */
    matmul("TN",m,rtk->nx,rtk->nx,1.0,D,rtk->P,0.0,E);
    matmul("NN",m,m,rtk->nx,1.0,E,D,0.0,Q);
    
    /* integer least square */
    if ((info=mlambda(m,2,B1,Q,N1,s,&boot))) {
        trace(2,"lambda error: info=%d\n",info);
        free(B1); free(N1); free(D); free(E); free(Q); free(NC); free(NN1);
        return 0;
    }
    if (s[0]<=0.0||boot<0.9) {//0.999
        trace(2, "validation error: s0=%f, s1=%f,ratio=%f,boot= %f",s[0],s[1],s[1]/s[0],boot);
        free(B1); free(N1); free(D); free(E); free(Q); free(NC); free(NN1);
        return 0;
    }
   
    rtk->sol.ratio=(float)(MIN(s[1]/s[0],999.9));
    
    /* varidation by ratio-test */
    if (rtk->opt.thresar[0]>0.0&&rtk->sol.ratio<rtk->opt.thresar[0]*WL_ratio) {
        trace(2, "do PAR\n");
		/* ppp partial ambiguity resolutions add by Xiao Yin*/
		del=ppp_amb_IFLAM_PAR(rtk,nav,m,m,sat1,sat2,NW,B1,N1,Q);
		if(del==0) {
            free(B1); free(N1); free(D); free(E); free(Q); free(NC); free(NN1);
            return 0;
        }
    }
    if ((m-del<3&&rtk->opt.navsys==1)||(m-del<4&&rtk->opt.navsys>1) ||((double)del/m>=0.5)) {
        free(B1); free(N1); free(D); free(E); free(Q); free(NC); free(NN1);
        return 0;
    }
    time2epoch(obs[0].time, ep);

    //trace(2, "fix success! *%5.0f%3.0f%3.0f%3.0f%3.0f%6.2f%4d%4d%4d\n", ep[0], ep[1], ep[2], ep[3], ep[4], ep[5], m - del , m , n);
    /* narrow-lane to iono-free ambiguity */
    for (i=0;i<m-del;i++) {

		NN1[i]=N1[i];

		 lam1=lam_carr[0]; 
		 if(satsys(sat1[i], NULL)==SYS_GPS)
		 {lam2=lam_carr[1];lam_NL=lam_LC(1,1,0);}
		 else if(satsys(sat1[i],NULL)==SYS_GAL) // support GAL
		 {lam2=lam_carr[2];lam_NL=lam_LC(1,0,1);}
         else if(satsys(sat1[i],NULL)==SYS_CMP) // support BDS
		 {lam1= lam_carr_CMP[0];lam2=lam_carr_CMP[2];lam_NL=lam_LC_CMP(1,0,1);}

		 C1= SQR(lam2)/(SQR(lam2)-SQR(lam1));
         C2=-SQR(lam1)/(SQR(lam2)-SQR(lam1));

        NC[i]=C1*lam1*N1[i]+C2*lam2*(N1[i]-NW[i]);
    }
    /* fixed solution */
    stat=fix_sol(rtk,sat1,sat2,NC,m-del);
    trace(2, "fix success! *%d%5.0f%3.0f%3.0f%3.0f%3.0f%6.2f%4d%4d%4d\n", stat,ep[0], ep[1], ep[2], ep[3], ep[4], ep[5], m - del, m, n);
    
	free(B1); free(N1); free(D); free(E); free(Q); free(NC); free(NN1);
    
    return stat;
}
int SeRefSearch(double *Bias,int *sat,int m){
    int i,j,k;
    double sum[MAXOBS]={0.0},temp=0.0;

    for (i=0;i<m;i++){
        for(j=0;j<m;j++){
            if(j!=i){
                temp=Bias[j]-Bias[i];
                sum[i]+=fabs(ROUND(temp)-temp);
            }
        }
    }
    /* Sat ����*/
    for (j=0;j<m-1;j++)for(k=j+1;k<m;k++) {
			if (sum[j]<=sum[k]) continue;
			SWAP_I(sat[j], sat[k]);
			SWAP_D(sum[j], sum[k]);
		}
    }
/* resolve integer ambiguity for ppp -----------------------------------------*/
extern int pppamb(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav,
                  const double *azel,int *exc)
{
#if 0
		const int sys[] = {SYS_GPS | SYS_QZS,SYS_GAL,SYS_CMP,0};
#else
	    const int sys[] = {SYS_GPS,SYS_GAL,SYS_CMP,0};
#endif
    double elmask,el[MAXOBS],fcb[MAXOBS]={0.0},lam_WL=0.0,BW[MAXOBS]={0.0};
    double sat_nlfcb;
    int i, j, k, m = 0, ns = 0, stat = 0, * NW, * sat1, * sat2, sat[MAXOBS];
    
    if (n<=0||rtk->opt.ionoopt!=IONOOPT_IFLC||rtk->opt.nf<2) return 0;
    
    trace(3,"pppamb: time=%s n=%d\n",time_str(obs[0].time,0),n);
    
    elmask=rtk->opt.elmaskar>0.0?rtk->opt.elmaskar:rtk->opt.elmin;
    
    sat1=imat(n*n,1); sat2=imat(n*n,1); NW=imat(n*n,1);
    
    /* average LC */
    average_LC(rtk,obs,n,nav,azel);
    

	/* Select Reference Satellite*/
	for (i=0; sys[i]; i++) {
		for (j=m=0; j<n; j++) {
			if (exc[j]||!(satsys(obs[j].sat,NULL)&sys[i])||
				!rtk->ssat[obs[j].sat-1].vsat[0]||azel[1+j*2]<elmask) continue;

			if (obs[j].L[0]==0) continue;

            /* BDS3 delete as ISB */
            if (obs[j].sat>125) continue;

			if (sys[i]==SYS_GPS){
                lam_WL=lam_LC(1,-1,0);
                BW[m]=(rtk->ambc+obs[j].sat-1)->LC[0]/lam_WL;
                fcb[m]=fabs(ROUND(BW[m])-BW[m]);
            }
            else if(sys[i]==SYS_GAL){
                lam_WL=lam_LC(1,0,-1);
                BW[m]=(rtk->ambc+obs[j].sat-1)->LC[1]/lam_WL;
                fcb[m]=fabs(ROUND(BW[m])-BW[m]);
            }
            else {
                lam_WL=lam_LC_CMP(1,0,-1);
                BW[m]=(rtk->ambc+obs[j].sat-1)->LC[1]/lam_WL;
                fcb[m]=fabs(ROUND(BW[m])-BW[m]);
            }
            
			sat[m]=obs[j].sat;
            el[m++]=azel[1+j*2]; 
		}
        /* Elevation Angle Method */
		for (j=0; j<m-1; j++)for (k=j+1; k<m; k++) {
			if (el[j]>=el[k]) continue;
			SWAP_I(sat[j], sat[k]);
			SWAP_D(el[j], el[k]);
		}
        /* FCB Method */
		//for (j=0; j<m-1; j++)for (k=j+1; k<m; k++) {
		//	if (fcb[j]<=fcb[k]) continue;
		//	SWAP_I(sat[j], sat[k]);
		//	SWAP_D(fcb[j], fcb[k]);
		//}

        /* FCB Search Method */
        //SeRefSearch(BW,sat,m);

		for (j=1;j<m; j++) {
			sat1[ns]=sat[j];
			sat2[ns]=sat[0];
            /* fix wide-lane ambiguity */
			if (fix_amb_WL(rtk,nav,sat1[ns],sat2[ns],NW+ns)) ns++;
		}
	}

    if (ns > 0) {
        int zyf = 1;
    }
	/* Select Reference Satellite for BDS3 */
	for (i=2; sys[i]; i++) {
		for (j=m=0; j<n; j++) {
			if (exc[j]||!(satsys(obs[j].sat,NULL)&sys[i])||
				!rtk->ssat[obs[j].sat-1].vsat[0]||azel[1+j*2]<elmask) continue;

            // sat_nlfcb=selfcb(rtk->sol.time,obs[j].sat,nav);
			if (obs[j].L[0]==0) continue;
		    // if(sat_nlfcb==0 ||sat_nlfcb==10 ) continue; 
             
            /* BDS3 delete as ISB */
            if (obs[j].sat<=125) continue;

			if (sys[i]==SYS_GPS){
                lam_WL=lam_LC(1,-1,0);
                BW[m]=(rtk->ambc+obs[j].sat-1)->LC[0]/lam_WL;
                fcb[m]=fabs(ROUND(BW[m])-BW[m]);
            }
            else if(sys[i]==SYS_GAL){
                lam_WL=lam_LC(1,0,-1);
                BW[m]=(rtk->ambc+obs[j].sat-1)->LC[1]/lam_WL;
                fcb[m]=fabs(ROUND(BW[m])-BW[m]);
            }
            else {
                lam_WL=lam_LC_CMP(1,0,-1);
                BW[m]=(rtk->ambc+obs[j].sat-1)->LC[1]/lam_WL;
                fcb[m]=fabs(ROUND(BW[m])-BW[m]);
            }
            
			sat[m]=obs[j].sat;
            el[m++]=azel[1+j*2]; 
		}
        /* Elevation Angle Method */
		for (j=0; j<m-1; j++)for (k=j+1; k<m; k++) {
			if (el[j]>=el[k]) continue;
			SWAP_I(sat[j], sat[k]);
			SWAP_D(el[j], el[k]);
		}

		for (j=1;j<m; j++) {
			sat1[ns]=sat[j];
			sat2[ns]=sat[0];
            /* fix wide-lane ambiguity */
			if (fix_amb_WL(rtk,nav,sat1[ns],sat2[ns],NW+ns)) ns++;
		}
	}

    if (ns > 0) {
        trace(2, "WL fix ns= %d \n)", ns);
    }
    /* fix narrow-lane ambiguity */
    if (rtk->opt.modear==ARMODE_PPPAR) {
        stat=fix_amb_ROUND(rtk,sat1,sat2,NW,ns);
    }
    else if (rtk->opt.modear==ARMODE_PPPAR_ILS) {
        stat=fix_amb_ILS(rtk,obs,n,nav,sat1,sat2,NW,ns);
    }
    free(sat1); free(sat2); free(NW);
    
    return stat;
}
