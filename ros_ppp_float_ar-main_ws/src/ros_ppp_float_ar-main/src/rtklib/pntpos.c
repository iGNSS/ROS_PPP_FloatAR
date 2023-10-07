/*------------------------------------------------------------------------------
* pntpos.c : standard positioning
*
*          Copyright (C) 2007-2018 by T.TAKASU, All rights reserved.
*
* version : $Revision:$ $Date:$
* history : 2010/07/28 1.0  moved from rtkcmn.c
*                           changed api:
*                               pntpos()
*                           deleted api:
*                               pntvel()
*           2011/01/12 1.1  add option to include unhealthy satellite
*                           reject duplicated observation data
*                           changed api: ionocorr()
*           2011/11/08 1.2  enable snr mask for single-mode (rtklib_2.4.1_p3)
*           2012/12/25 1.3  add variable snr mask
*           2014/05/26 1.4  support galileo and beidou
*           2015/03/19 1.5  fix bug on ionosphere correction for GLO and BDS
*           2018/10/10 1.6  support api change of satexclude()
*-----------------------------------------------------------------------------*/
#include "rtklib.h"

/* constants -----------------------------------------------------------------*/

#define SQR(x)      ((x)*(x))

#define NX          (4+3)       /* # of estimated parameters */

#define MAXITR      10          /* max number of iteration for point pos */
#define ERR_ION     5.0         /* ionospheric delay std (m) */
#define ERR_TROP    3.0         /* tropspheric delay std (m) */
#define ERR_SAAS    0.3         /* saastamoinen model error std (m) */
#define ERR_BRDCI   0.5         /* broadcast iono model error factor */
#define ERR_CBIAS   0.3         /* code bias error std (m) */
#define REL_HUMI    0.7         /* relative humidity for saastamoinen model */

static int eph_sel[]={ /* GPS,GLO,GAL,QZS,BDS,SBS */
	0,0,1,0,0,0
};
/* select ephememeris --------------------------------------------------------*/
static eph_t *seleph(gtime_t time, int sat, int iode, const nav_t *nav)
{
	double t,tmax,tmin;
	int i,j=-1,sys,sel=0;

	trace(4,"seleph  : time=%s sat=%2d iode=%d\n",time_str(time,3),sat,iode);

	sys=satsys(sat,NULL);
	switch (sys) {
	case SYS_GPS: tmax=MAXDTOE+1.0    ; sel=eph_sel[0]; break;
	case SYS_GAL: tmax=MAXDTOE_GAL    ; sel=eph_sel[2]; break;
	case SYS_QZS: tmax=MAXDTOE_QZS+1.0; sel=eph_sel[3]; break;
	case SYS_BDS: tmax=MAXDTOE_BDS+1.0; sel=eph_sel[4]; break;
	default: tmax=MAXDTOE+1.0; break;
	}
	tmin=tmax+1.0;

	for (i=0;i<nav->n;i++) {
		if (nav->eph[i].sat!=sat) continue;
		if (iode>=0&&nav->eph[i].iode!=iode) continue;
		if (sys==SYS_GAL&&sel) {
			if (sel==1&&!(nav->eph[i].code&(1<<9))) continue; /* I/NAV */
			if (sel==2&&!(nav->eph[i].code&(1<<8))) continue; /* F/NAV */
		}
		if ((t=fabs(timediff(nav->eph[i].toe,time)))>tmax) continue;
		if (iode>=0) return nav->eph+i;
		if (t<=tmin) {j=i; tmin=t;} /* toe closest to time */
	}
	if (iode>=0||j<0) {
		trace(3,"no broadcast ephemeris: %s sat=%2d iode=%3d\n",time_str(time,0),
			sat,iode);
		return NULL;
	}
	return nav->eph+j;
}

/* pseudorange measurement error variance ------------------------------------*/
static double varerr(const prcopt_t *opt, double el, int sys)
{
    double fact,varr;
    fact=sys==SYS_GLO?EFACT_GLO:(sys==SYS_SBS?EFACT_SBS:EFACT_GPS);
    varr=SQR(opt->err[0])*(SQR(opt->err[1])+SQR(opt->err[2])/sin(el));
    if (opt->ionoopt==IONOOPT_IFLC) varr*=SQR(3.0); /* iono-free */
    return SQR(fact)*varr;
}
/* get tgd parameter (m) -----------------------------------------------------*/
static double gettgd(int sat, const nav_t *nav)
{
    int i;
    for (i=0;i<nav->n;i++) {
        if (nav->eph[i].sat!=sat) continue;
        return CLIGHT*nav->eph[i].tgd[0];
    }
    return 0.0;
}
static double gettgd1(int sat, const nav_t *nav)
{
	int i;
	for (i=0;i<nav->n;i++) {
		if (nav->eph[i].sat!=sat) continue;
		return CLIGHT*nav->eph[i].tgd[1];
	}
	return 0.0;
}

/* mutipath correct-------------------------------------------------------------
* BeiDou satellite-induced code pseudorange variations correct
* args   : rtk_t *rtk       IO  rtk control/result struct
obsd_t *obs      IO  observation data
int    n         I   number of observation data
nav_t  *nav      I   navigation messages
* note   :  
* -----------------------------------------------------------------------------*/
extern double BDmultipath(obsd_t *obs, double el, double *dmp)
{
    int i, j, sat, b;
    double elev, a;
	char s0[32], prn[2];

    const static double IGSOCOEF[3][10] =  		/* m */
    {
        { -0.55, -0.40, -0.34, -0.23, -0.15, -0.04, 0.09, 0.19, 0.27, 0.35},	//B1
        { -0.71, -0.36, -0.33, -0.19, -0.14, -0.03, 0.08, 0.17, 0.24, 0.33},	//B2
        { -0.27, -0.23, -0.21, -0.15, -0.11, -0.04, 0.05, 0.14, 0.19, 0.32},	//B3
    };
	const static double MEOCOEF[3][10] =  		/* m */
	{
		{ -0.47, -0.38, -0.32, -0.23, -0.11, 0.06, 0.34, 0.69, 0.97, 1.05},	//B1
		{ -0.40, -0.31, -0.26, -0.18, -0.06, 0.09, 0.28, 0.48, 0.64, 0.69},	//B2
		{ -0.22, -0.15, -0.13, -0.10, -0.04, 0.05, 0.14, 0.27, 0.36, 0.47},	//B3
	};

	satno2id(obs->sat, s0); strncpy(prn, s0 + 1, 2); sat = atoi(prn);

	if (sat <= 5) return 0.0;// exclude GEO

	elev = el * R2D;// elevation

	for (j=0;j<3;j++) dmp[j] = 0.0;

	a = elev * 0.1;
	b = (int)(a);

	int bIGSO = (sat >= 6 && sat < 11) || 13 == sat || 16 == sat;

	if (bIGSO)   // IGSO(C06, C07, C08, C09, C10, C13, C16)
	{
		if (b < 0)
		{
			for (j=0; j < 3; j++) dmp[j] = IGSOCOEF[j][0];
		}
		else if (b >= 9)
		{
			for (j=0; j < 3; j++) dmp[j] = IGSOCOEF[j][9];
		}
		else
		{
			for (j=0; j < 3; j++) dmp[j] = IGSOCOEF[j][b] * (1.0 - a + b) + IGSOCOEF[j][b + 1] * (a - b);
		}
	}
	else if (sat == 11|| sat==12|| sat==14)   // MEO(C11, C12, C14)
	{
		if (b < 0)
		{
			for (j = 0; j < 3; j++) dmp[j] = MEOCOEF[j][0];
		}
		else if (b >= 9)
		{
			for (j = 0; j < 3; j++) dmp[j] = MEOCOEF[j][9];
		}
		else
		{
			for (j = 0; j < 3; j++) dmp[j] = MEOCOEF[j][b] * (1.0 - a + b) + MEOCOEF[j][b + 1] * (a - b);
		}
	}
	return 1;
}

/* psendorange with code bias correction -------------------------------------*/
static double prange(const obsd_t *obs, const nav_t *nav, const double *azel,
                     int iter, const prcopt_t *opt, double *var)
{
    const double *lam=nav->lam[obs->sat-1];
    double PC,P1,P2,P1_P2,P1_C1,P2_C2,gamma,dmp[3]={0.0};
    int i=0,j=1,sys;
    
    *var=0.0;
    
    if (!(sys=satsys(obs->sat,NULL))) return 0.0;
    
    /* L1-L2 for GPS/GLO/QZS, L1-L5 for GAL/SBS */
    if (NFREQ>=3&&(sys&(SYS_GAL|SYS_SBS|SYS_BDS))) j=2;
    
    if (NFREQ<2||lam[i]==0.0||lam[j]==0.0) return 0.0;

    /* single-frequency/P1-P2 not used */
    if (!obs->P[0]||!obs->P[j]) return 0.0;
    if(fabs(obs->P[0]-obs->P[j])>25)return 0.0;
    
    /* test snr mask */
    if (iter>0) {
        if (testsnr(0,i,azel[1],obs->SNR[i]*0.25,&opt->snrmask)) {
            trace(4,"snr mask: %s sat=%2d el=%.1f snr=%.1f\n",
                  time_str(obs->time,0),obs->sat,azel[1]*R2D,obs->SNR[i]*0.25);
            return 0.0;
        }
        if (opt->ionoopt==IONOOPT_IFLC) {
            if (testsnr(0,j,azel[1],obs->SNR[j]*0.25,&opt->snrmask)) return 0.0;
        }
    }
    gamma=SQR(lam[j])/SQR(lam[i]); /* f1^2/f2^2 */
    P1=obs->P[i];
    P2=obs->P[j];
    P1_P2=nav->cbias[obs->sat-1][0];
    P1_C1=nav->cbias[obs->sat-1][1];
    P2_C2=nav->cbias[obs->sat-1][2];
    
    /* if no P1-P2 DCB, use TGD instead */
    if (P1_P2==0.0&&(sys&(SYS_GPS|SYS_GAL|SYS_QZS))) {
        P1_P2=(1.0-gamma)*gettgd(obs->sat,nav);
    }

    if (opt->ionoopt==IONOOPT_IFLC) { /* dual-frequency */
        
        if (P1==0.0||P2==0.0) return 0.0;
        if (obs->code[i]==CODE_L1C) P1+=P1_C1; /* C1->P1 */
        if (obs->code[j]==CODE_L2C) P2+=P2_C2; /* C2->P2 */
        
        if (sys&(SYS_BDS)) {/* Code Pseudorange Multipath Correct*/
            BDmultipath(obs,azel[1],dmp);
            P1+= dmp[0];
            P2+= dmp[j];
        };
        /* iono-free combination */
        PC=(gamma*P1-P2)/(gamma-1.0);
    }
    else { /* single-frequency */
        
        if (P1==0.0) return 0.0;
        if (obs->code[i]==CODE_L1C) P1+=P1_C1; /* C1->P1 */
        if (sys&(SYS_BDS)) {
            BDmultipath(obs,azel[1],dmp);
            P1+= dmp[0];
        };
        PC=P1-P1_P2/(1.0-gamma);
    }
    if (opt->sateph==EPHOPT_SBAS) PC-=P1_C1; /* sbas clock based C1 */

	/*BDS���ݴ���ʱ��B1I��B2I��ҪTGD������B3I*/
	if (P1_P2==0.0&&(sys&(SYS_BDS))) {
		if(opt->mode==0){
			if(opt->ionoopt==IONOOPT_IFLC){
                P1_P2=gamma*gettgd(obs->sat,nav)/(gamma-1);
            }
            else
            {
                P1_P2=gettgd(obs->sat,nav);
            }
			PC -= P1_P2;
		}
		else
		{
			P1 = obs->P[0];
			P2 = obs->P[2];
			if (P1==0.0||P2==0.0) return 0.0;
			if (!BDmultipath(obs, azel[1], dmp)) return 0.0;
			P1 += dmp[0];
			P2 += dmp[2];
			gamma = SQR(lam[2]) / SQR(lam[0]);
			/* iono-free combination */
			PC = (gamma*P1 - P2) / (gamma - 1.0);
		}
	}
    
    *var=SQR(ERR_CBIAS);
    
    return PC;
}
/* ionospheric correction ------------------------------------------------------
* compute ionospheric correction
* args   : gtime_t time     I   time
*          nav_t  *nav      I   navigation data
*          int    sat       I   satellite number
*          double *pos      I   receiver position {lat,lon,h} (rad|m)
*          double *azel     I   azimuth/elevation angle {az,el} (rad)
*          int    ionoopt   I   ionospheric correction option (IONOOPT_???)
*          double *ion      O   ionospheric delay (L1) (m)
*          double *var      O   ionospheric delay (L1) variance (m^2)
* return : status(1:ok,0:error)
*-----------------------------------------------------------------------------*/
extern int ionocorr(gtime_t time, const nav_t *nav, int sat, const double *pos,
                    const double *azel, int ionoopt, double *ion, double *var)
{
    trace(4,"ionocorr: time=%s opt=%d sat=%2d pos=%.3f %.3f azel=%.3f %.3f\n",
          time_str(time,3),ionoopt,sat,pos[0]*R2D,pos[1]*R2D,azel[0]*R2D,
          azel[1]*R2D);
    
    /* broadcast model */
    if (ionoopt==IONOOPT_BRDC) {
        *ion=ionmodel(time,nav->ion_gps,pos,azel);
        *var=SQR(*ion*ERR_BRDCI);
        return 1;
    }
    /* sbas ionosphere model */
    if (ionoopt==IONOOPT_SBAS) {
        return sbsioncorr(time,nav,pos,azel,ion,var);
    }
    /* ionex tec model */
    if (ionoopt==IONOOPT_TEC) {
        return iontec(time,nav,pos,azel,1,ion,var);
    }
    /* qzss broadcast model */
    if (ionoopt==IONOOPT_QZS&&norm(nav->ion_qzs,8)>0.0) {
        *ion=ionmodel(time,nav->ion_qzs,pos,azel);
        *var=SQR(*ion*ERR_BRDCI);
        return 1;
    }
    /* lex ionosphere model */
    if (ionoopt==IONOOPT_LEX) {
        return lexioncorr(time,nav,pos,azel,ion,var);
    }
    *ion=0.0;
    *var=ionoopt==IONOOPT_OFF?SQR(ERR_ION):0.0;
    return 1;
}
/* tropospheric correction -----------------------------------------------------
* compute tropospheric correction
* args   : gtime_t time     I   time
*          nav_t  *nav      I   navigation data
*          double *pos      I   receiver position {lat,lon,h} (rad|m)
*          double *azel     I   azimuth/elevation angle {az,el} (rad)
*          int    tropopt   I   tropospheric correction option (TROPOPT_???)
*          double *trp      O   tropospheric delay (m)
*          double *var      O   tropospheric delay variance (m^2)
* return : status(1:ok,0:error)
*-----------------------------------------------------------------------------*/
extern int tropcorr(gtime_t time, const nav_t *nav, const double *pos,
                    const double *azel, int tropopt, double *trp, double *var)
{
    trace(4,"tropcorr: time=%s opt=%d pos=%.3f %.3f azel=%.3f %.3f\n",
          time_str(time,3),tropopt,pos[0]*R2D,pos[1]*R2D,azel[0]*R2D,
          azel[1]*R2D);
    
    /* saastamoinen model */
    if (tropopt==TROPOPT_SAAS||tropopt==TROPOPT_EST||tropopt==TROPOPT_ESTG) {
        *trp=tropmodel(time,pos,azel,REL_HUMI);
        *var=SQR(ERR_SAAS/(sin(azel[1])+0.1));
        return 1;
    }
    /* sbas troposphere model */
    if (tropopt==TROPOPT_SBAS) {
        *trp=sbstropcorr(time,pos,azel,var);
        return 1;
    }
    /* no correction */
    *trp=0.0;
    *var=tropopt==TROPOPT_OFF?SQR(ERR_TROP):0.0;
    return 1;
}
/* pseudorange residuals -----------------------------------------------------*/
static int rescode(int iter, const obsd_t *obs, int n, const double *rs,
                   const double *dts, const double *vare, const int *svh,
                   const nav_t *nav, const double *x, const prcopt_t *opt,
                   double *v, double *H, double *var, double *azel, int *vsat,
                   double *resp, int *ns, int flag)
{
    double r,dion,dtrp,vmeas,vion,vtrp,rr[3],pos[3],dtr,e[3],P,lam_L1,temp,snr=0;
    int i,j,nv=0,sys,mask[4]={0};
    
    trace(3,"resprng : n=%d\n",n);
    
    for (i=0;i<3;i++) rr[i]=x[i];
    dtr=x[3];
    
    ecef2pos(rr,pos);
    
    for (i=*ns=0;i<n&&i<MAXOBS;i++) {
        vsat[i]=0; azel[i*2]=azel[1+i*2]=resp[i]=0.0;
        
        if (!(sys=satsys(obs[i].sat,NULL))) continue;
        
        /* reject duplicated observation data */
        if (i<n-1&&i<MAXOBS-1&&obs[i].sat==obs[i+1].sat) {
            //trace(2,"duplicated observation data %s sat=%2d\n",
            //      time_str(obs[i].time,3),obs[i].sat);
            i++;
            continue;
        }
        /* geometric distance/azimuth/elevation angle */
        if ((r=geodist(rs+i*6,rr,e))<=0.0||
            satazel(pos,e,azel+i*2)<opt->elmin) continue;
        
        /* psudorange with code bias correction */
        if ((P=prange(obs+i,nav,azel+i*2,iter,opt,&vmeas))==0.0) continue;
        
        /* excluded satellite? */
        if (satexclude(obs[i].sat,vare[i],svh[i],opt)) continue;
        
        /* ionospheric corrections */
        if (!ionocorr(obs[i].time,nav,obs[i].sat,pos,azel+i*2,
                      iter>0?opt->ionoopt:IONOOPT_BRDC,&dion,&vion)) continue;
       
        if(opt->ionoopt==IONOOPT_IFLC) {
            dion=0.0;
        }
       
        /* GPS-L1 -> L1/B1 */
        if ((lam_L1=nav->lam[obs[i].sat-1][0])>0.0) {
            dion*=SQR(lam_L1/lam_carr[0]);
        }
        /* tropospheric corrections */
        if (!tropcorr(obs[i].time,nav,pos,azel+i*2,
                      iter>0?opt->tropopt:TROPOPT_SAAS,&dtrp,&vtrp)) {
            continue;
        }

		/* BDS adop LC */
		//if (sys==SYS_BDS&&opt->mode!=0) dion=0.0;

        /* pseudorange residual */
        v[nv]=P-(r+dtr-CLIGHT*dts[i*2]+dion+dtrp);
		
		temp=OMGE*((rs+i*6)[0]*rr[1]-(rs+i*6)[1]*rr[0])/CLIGHT; 
		snr=obs[i].SNR[0]/4;

        /* design matrix */
        for (j=0;j<NX;j++) H[j+nv*NX]=j<3?-e[j]:(j==3?1.0:0.0);
        
        /* time system and receiver bias offset correction */
        if      (sys==SYS_GLO) {v[nv]-=x[4]; H[4+nv*NX]=1.0; mask[1]=1;}
        else if (sys==SYS_GAL) {v[nv]-=x[5]; H[5+nv*NX]=1.0; mask[2]=1;}
        else if (sys==SYS_BDS) {v[nv]-=x[6]; H[6+nv*NX]=1.0; mask[3]=1;}
        else mask[0]=1;
        
        vsat[i]=1; resp[i]=v[nv]; (*ns)++;
        
        /* error variance */
        var[nv++]=varerr(opt,azel[1+i*2],sys)+vare[i]+vmeas+vion+vtrp;
        
        trace(4,"sat=%2d azel=%5.1f %4.1f res=%7.3f sig=%5.3f\n",obs[i].sat,
              azel[i*2]*R2D,azel[1+i*2]*R2D,resp[i],sqrt(var[nv-1]));
    }
    /* constraint to avoid rank-deficient */
    for (i=0;i<4;i++) {
        if (mask[i]) continue;
        v[nv]=0.0;
        for (j=0;j<NX;j++) H[j+nv*NX]=j==i+3?1.0:0.0;
        var[nv++]=0.01;
    }
    return nv;
}
/* validate solution ---------------------------------------------------------*/
static int valsol(const double *azel, const int *vsat, int n,
                  const prcopt_t *opt, const double *v, int nv, int nx,
                  char *msg)
{
    double azels[MAXOBS*2],dop[4],vv;
    int i,ns;
    
    trace(3,"valsol  : n=%d nv=%d\n",n,nv);
    
    /* chi-square validation of residuals  vv>chisqr[nv-nx-1]*/
    vv=dot(v,v,nv);
    //if (nv>nx&&vv>chisqr[nv-nx-1]) {
    //    sprintf(msg,"chi-square error nv=%d vv=%.1f cs=%.1f",nv,vv,chisqr[nv-nx-1]);
    //    return 0;
    //}
    /* large gdop check */
    for (i=ns=0;i<n;i++) {
        if (!vsat[i]) continue;
        azels[  ns*2]=azel[  i*2];
        azels[1+ns*2]=azel[1+i*2];
        ns++;
    }
    dops(ns,azels,opt->elmin,dop);
    if (dop[0]<=0.0||dop[0]>opt->maxgdop) {
        sprintf(msg,"gdop error nv=%d gdop=%.1f",nv,dop[0]);
        return 0;
    }
    return 1;
}
/* estimate receiver position ------------------------------------------------*/
static int estpos(const obsd_t *obs, int n, const double *rs, const double *dts,
                  const double *vare, const int *svh, const nav_t *nav,
                  const prcopt_t *opt, sol_t *sol, double *azel, int *vsat,
                  double *resp, char *msg)
{
    double x[NX]={0},dx[NX]={0.0},Q[NX*NX]={0.0},*v,*H,*var,sig;
    int i,j,k,info,stat,nv,ns,iter_flag;
    
    trace(3,"estpos  : n=%d\n",n);
    
    v=mat(n+4,1); H=mat(NX,n+4); var=mat(n+4,1);
    
    for (i=0;i<3;i++) x[i]=sol->rr[i];
	if(x[1]==0) iter_flag=0;
	else  iter_flag=1;

    for (i=0;i<MAXITR;i++) {
        
        /* pseudorange residuals */
        nv=rescode(i,obs,n,rs,dts,vare,svh,nav,x,opt,v,H,var,azel,vsat,resp,
                   &ns,iter_flag);
        
        if (nv<NX) {
            sprintf(msg,"lack of valid sats ns=%d",nv);
            break;
        }
		/* weight by variance */
		for (j=0;j<nv;j++) {
			sig=sqrt(var[j]);
			v[j]/=sig;
			for (k=0;k<NX;k++) H[k+j*NX]/=sig;
		}
        /* least square estimation */
        if ((info=lsq(H,v,NX,nv,dx,Q))) {
            sprintf(msg,"lsq error info=%d",info);
            break;
        }
        for (j=0;j<NX;j++) x[j]+=dx[j];
        
        if (norm(dx,NX)<1E-4) {
            sol->type=0;
            sol->time=timeadd(obs[0].time,-x[3]/CLIGHT);
            sol->dtr[0]=x[3]/CLIGHT; /* receiver clock bias (s) */
            sol->dtr[1]=x[4]/CLIGHT; /* glo-gps time offset (s) */
            sol->dtr[2]=x[5]/CLIGHT; /* gal-gps time offset (s) */
            sol->dtr[3]=x[6]/CLIGHT; /* bds-gps time offset (s) */
            /*for (j=0;j<6;j++) sol->rr[j]=j<3?x[j]:0.0;*/
			for (j=0;j<3;j++) sol->rr[j]=x[j];
            for (j=0;j<3;j++) sol->qr[j]=(float)Q[j+j*NX];
            sol->qr[3]=(float)Q[1];    /* cov xy */
            sol->qr[4]=(float)Q[2+NX]; /* cov yz */
            sol->qr[5]=(float)Q[2];    /* cov zx */
            sol->ns=(unsigned char)ns;
            sol->age=sol->ratio=0.0;
            
            /* validate solution */
            if ((stat=valsol(azel,vsat,n,opt,v,nv,NX,msg))) {
                sol->stat=opt->sateph==EPHOPT_SBAS?SOLQ_SBAS:SOLQ_SINGLE;
            }
            free(v); free(H); free(var);
            
            return stat;
        }
    }
    if (i>=MAXITR) sprintf(msg,"iteration divergent i=%d",i);
    
    free(v); free(H); free(var);
    
    return 0;
}
/* raim fde (failure detection and exclution) -------------------------------*/
static int raim_fde(const obsd_t *obs, int n, const double *rs,
                    const double *dts, const double *vare, const int *svh,
                    const nav_t *nav, const prcopt_t *opt, sol_t *sol,
                    double *azel, int *vsat, double *resp, char *msg)
{
    obsd_t *obs_e;
    sol_t sol_e={{0}};
    char tstr[32],name[16],msg_e[128];
    double *rs_e,*dts_e,*vare_e,*azel_e,*resp_e,rms_e,rms=100.0;
    int i,j,k,nvsat,stat=0,*svh_e,*vsat_e,sat=0;
    
    trace(3,"raim_fde: %s n=%2d\n",time_str(obs[0].time,0),n);
    
    if (!(obs_e=(obsd_t *)malloc(sizeof(obsd_t)*n))) return 0;
    rs_e = mat(6,n); dts_e = mat(2,n); vare_e=mat(1,n); azel_e=zeros(2,n);
    svh_e=imat(1,n); vsat_e=imat(1,n); resp_e=mat(1,n); 
    
    for (i=0;i<n;i++) {
        
        /* satellite exclution */
        for (j=k=0;j<n;j++) {
            if (j==i) continue;
            obs_e[k]=obs[j];
            matcpy(rs_e +6*k,rs +6*j,6,1);
            matcpy(dts_e+2*k,dts+2*j,2,1);
            vare_e[k]=vare[j];
            svh_e[k++]=svh[j];
        }
        /* estimate receiver position without a satellite */
        if (!estpos(obs_e,n-1,rs_e,dts_e,vare_e,svh_e,nav,opt,&sol_e,azel_e,
                    vsat_e,resp_e,msg_e)) {
            trace(3,"raim_fde: exsat=%2d (%s)\n",obs[i].sat,msg);
            continue;
        }
        for (j=nvsat=0,rms_e=0.0;j<n-1;j++) {
            if (!vsat_e[j]) continue;
            rms_e+=SQR(resp_e[j]);
            nvsat++;
        }
        if (nvsat<5) {
            trace(3,"raim_fde: exsat=%2d lack of satellites nvsat=%2d\n",
                  obs[i].sat,nvsat);
            continue;
        }
        rms_e=sqrt(rms_e/nvsat);
        
        trace(3,"raim_fde: exsat=%2d rms=%8.3f\n",obs[i].sat,rms_e);
        
        if (rms_e>rms) continue;
        
        /* save result */
        for (j=k=0;j<n;j++) {
            if (j==i) continue;
            matcpy(azel+2*j,azel_e+2*k,2,1);
            vsat[j]=vsat_e[k];
            resp[j]=resp_e[k++];
        }
        stat=1;
        *sol=sol_e;
        sat=obs[i].sat;
        rms=rms_e;
        vsat[i]=0;
      //strcpy(msg,msg_e);
    }
    if (stat) {
        time2str(obs[0].time,tstr,2); satno2id(sat,name);
       // trace(2,"%s: %s excluded by raim\n",tstr+11,name);
    }
    free(obs_e);
    free(rs_e ); free(dts_e ); free(vare_e); free(azel_e);
    free(svh_e); free(vsat_e); free(resp_e);
    return stat;
}
/* doppler residuals ---------------------------------------------------------*/
static int resdop(int iter,const obsd_t *obs, int n, const double *rs, const double *dts,
                  const nav_t *nav, const double *rr, const double *x,
                  const double *azel, const int *vsat, double *v, double *H)
{
    double lam,rate,pos[3],E[9],a[3]={0.0},e[3],vs[3],cosel=0.0,Vs,C1=0.0,dtrp=0.0,vtrp=0.0,decc=0.0,r;
    int i,j,nv=0,m;

    trace(3,"resdop  : n=%d\n",n);
    

    ecef2pos(rr,pos); xyz2enu(pos,E);
    
    for (i=0;i<n&&i<MAXOBS;i++) {
        
        lam=nav->lam[obs[i].sat-1][0];
        
        if (obs[i].D[0]==0.0||lam==0.0||!vsat[i]||norm(rs+3+i*6,3)<=0.0/*||obs[i].sat==73*/) 
		{
            continue;
        }
        /* line-of-sight vector in ecef */
/*------------------------------------------------------------------------------------*/
		//cosel=cos(azel[1+i*2]);
        //a[0]=sin(azel[i*2])*cosel;
        //a[1]=cos(azel[i*2])*cosel;
        //a[2]=sin(azel[1+i*2]);
        //matmul("TN",3,1,3,1.0,E,a,0.0,e);

		for (m=0;m<3;m++) e[m]=rs[m+i*6]-rr[m];
		r=norm(e,3);
		for (m=0;m<3;m++) e[m]/=r;
/*------------------------------------------------------------------------------------*/
        /* satellite velocity relative to receiver in ecef */
        for (j=0;j<3;j++) vs[j]=rs[j+3+i*6]-x[j];
        
        /* range rate with earth rotation correction */
        rate=dot(vs,e,3)+OMGE/CLIGHT*(rs[4+i*6]*rr[0]+rs[1+i*6]*x[0]-
                                      rs[3+i*6]*rr[1]-rs[  i*6]*x[1]);
        
/*------------------------------------------------------------------------------------*/
		/* doppler residual�������κθ�����ʱ�� */
        /*v[nv]=-lam*obs[i].D[0]-(rate+x[3]-CLIGHT*dts[1+i*2]);*/

/*------------------------------------------------------------------------------------*/
		/*���Ǽ���չ���Ķ�����Ӱ��İ취*/
		Vs=0.5*(SQR(rs[3+i*6])+SQR(rs[4+i*6])+SQR(rs[5+i*6]))/CLIGHT;/*�����ٶ�V^2*/
        v[nv]=-lam*obs[i].D[0]-(rate+x[3]-CLIGHT*dts[1+i*2]+0.5*Vs/CLIGHT);
	   //if(obs[i].sat<=35)v[nv]=-lam*obs[i].D[0]-(rate+x[3]-CLIGHT*dts[1+i*2]+0.5*Vs/CLIGHT);
	   //else v[nv]=-lam*obs[i].D[0]-(rate+x[4]-CLIGHT*dts[1+i*2]+0.5*Vs/CLIGHT);
/*------------------------------------------------------------------------------------*/
		/*����˫Ƶ���������Ӱ��İ취*/
		//C1=1.64694;//f1^2/f2^2	
		//v[nv]=-(C1*lam*obs[i].D[0]-sqrt(C1*lam*lam)*obs[i].D[1])/(C1-1)-(rate+x[3]-CLIGHT*dts[1+i*2]+0.5*Vs/CLIGHT);
/*------------------------------------------------------------------------------------*/
		/*���Ƕ������ӳٱ仯Ӱ��İ취*/
		//if (tropcorr(obs[i].time,nav,pos,azel+i*2,TROPOPT_SAAS,&dtrp,&vtrp))
		//{
		//	dtrp=dtrp*cos(PI-azel[1+2*i]);
		//	dtrp=dtrp*cos(azel[1+2*i])*(1.0E-4)/SQR(sin(azel[1+2*i])+0.00143/(tan(azel[1+2*i])+0.0445));
		//}
		//v[nv]=-lam*obs[i].D[0]-(rate+x[3]-CLIGHT*dts[1+i*2]+Vs+abs(dtrp));
/*------------------------------------------------------------------------------------*/
		/*���Ƿ�Բ�������Ӱ��İ취������취���Ǻܺã�ֱ�����Ӳ�����ȽϺ�*/
		//if (!(eph=seleph(obs[0].time,obs[i].sat,-1,nav))) return 0;
		//if (obs[i].sat == eph->sat)
		//{decc=(2*3.986004415E+14)/CLIGHT*( 1/eph->A-1/(SQR(rs[0+i*6])+SQR(rs[1+i*6])+SQR(rs[2+i*6])));}
		//v[nv]=-lam*obs[i].D[0]-(rate+x[3]-CLIGHT*dts[1+i*2]+Vs-decc);
/*------------------------------------------------------------------------------------*/
       // if (iter==0&&x[0]!=0&&x[1]!=0) {
	      //      trace(2,"%s %3d %4.1f %7.3f\n",time_str(obs[i].time,3),obs[i].sat,
		     //azel[1+i*2]*R2D,v[i]);
       //  }
        /* design matrix */
        for (j=0;j<4;j++) H[j+nv*4]=j<3?-e[j]:1.0;
		//for (j=0;j<3;j++) H[j+nv*5]=-e[j];
		//if(obs[i].sat<=35) H[3+nv*5]=1,H[4+nv*5]=0;
		//else  H[3+nv*5]=0,H[4+nv*5]=1;
        nv++;
    }
    return nv;
}
/* doppler residuals by variation-----------------------------------------------------*/
static int resdopp(int iter,const obsd_t *obs, int n, const double *rs, const double *dts,
	const nav_t *nav, const double *rr, const double *x,
	const double *azel, const int *vsat, double *v, double *H,double *var)
{
	double lam,rate,pos[3],E[9],a[3]={0.0},e[3],vs[3],cosel=0.0,Vs,C1=0.0,dtrp=0.0,vtrp=0.0,decc=0.0,r,snr,ep[6];
	int i,j,nv=0,m;

	trace(3,"resdop  : n=%d\n",n);


	ecef2pos(rr,pos); xyz2enu(pos,E);
	time2epoch(obs[0].time,ep);
	for (i=0;i<n&&i<MAXOBS;i++) {

		lam=nav->lam[obs[i].sat-1][0];

		if (obs[i].D[0]==0.0||lam==0.0||norm(rs+3+i*6,3)<=0.0||azel[1+i*2]<10*D2R/*||obs[i].sat==73||obs[i].sat==77*/) 
		{
			continue;
		}
		//if (obs[i].sat!=13&&obs[i].sat!=5&&obs[i].sat!=66&&obs[i].sat!=67&&obs[i].sat!=72&&obs[i].sat!=2) 
		//{
		//	continue;
		//}
		//if(ep[4]>=33&&ep[4]<42&&(obs[i].sat==5||obs[i].sat==2)){
		//	continue;
		//}
		/* line-of-sight vector in ecef */
		/*------------------------------------------------------------------------------------*/
		for (m=0;m<3;m++) e[m]=rs[m+i*6]-rr[m];
		r=norm(e,3);
		for (m=0;m<3;m++) e[m]/=r;
		/*------------------------------------------------------------------------------------*/
		/* satellite velocity relative to receiver in ecef */
		for (j=0;j<3;j++) vs[j]=rs[j+3+i*6]-x[j];

		/* range rate with earth rotation correction */
		rate=dot(vs,e,3)+OMGE/CLIGHT*(rs[4+i*6]*rr[0]+rs[1+i*6]*x[0]-
			rs[3+i*6]*rr[1]-rs[  i*6]*x[1]);

		/*------------------------------------------------------------------------------------*/
		/* doppler residual�������κθ�����ʱ�� */
		/*v[nv]=-lam*obs[i].D[0]-(rate+x[3]-CLIGHT*dts[1+i*2]);*/

		/*------------------------------------------------------------------------------------*/
		/*���Ǽ���չ���Ķ�����Ӱ��İ취*/
		Vs=0.5*(SQR(rs[3+i*6])+SQR(rs[4+i*6])+SQR(rs[5+i*6]))/CLIGHT;/*�����ٶ�V^2*/
		v[nv]=-lam*obs[i].D[0]-(rate+x[3]-CLIGHT*dts[1+i*2]+0.5*Vs/CLIGHT);

		snr=obs[i].SNR[0]/4;
		//���ʱ�� ���� ����������� ������ά�ٶ�
		//if (iter==1&&snr>30) {
		//	trace(2,"%s %3d %12.3f %12.3f %12.3f %12.3f\n",time_str(obs[i].time,3),obs[i].sat,-lam*obs[i].D[0]+CLIGHT*dts[1+i*2]
		//		,(rs+i*6)[3],(rs+i*6)[4],(rs+i*6)[5]);//����۲�ֵ���Ա�Matlab����
		//}


		//if(obs[i].sat<=35)v[nv]=-lam*obs[i].D[0]-(rate+x[3]-CLIGHT*dts[1+i*2]+0.5*Vs/CLIGHT);
		//else v[nv]=-lam*obs[i].D[0]-(rate+x[3]+x[4]-CLIGHT*dts[1+i*2]+0.5*Vs/CLIGHT);

		/*------------------------------------------------------------------------------------*/
		/* error variance */
		if(azel[1+i*2]>30*D2R) var[nv]=SQR(1);
	    else var[nv]=SQR(1)/SQR(2*sin(azel[1+i*2]));//SQR(0.002)/(2*sin(azel[1+i*2]))+SQR(0.0005);
		if(obs[i].sat==73||obs[i].sat==77){
			var[nv]=var[nv]+SQR(1);
			}
	
		//--------�������NLOS/MP-------
		//if (iter==0&&(fabs(x[0])+fabs(x[1])+fabs(x[2]))>0) {
		//	trace(2,"%s %3d %4.1f %7.3f\n",time_str(obs[i].time,3),obs[i].sat,
		//		azel[1+i*2]*R2D,v[nv]);
		//}

		/* design matrix */
		for (j=0;j<4;j++) H[j+nv*4]=j<3?-e[j]:1.0;
		//for (j=0;j<3;j++) H[j+nv*5]=-e[j];
		//if(obs[i].sat<=35) H[3+nv*5]=1,H[4+nv*5]=0;
		//else  H[3+nv*5]=1,H[4+nv*5]=1;
		nv++;
	}
	return nv;
}
/* doppler residuals by post-----------------------------------------------------*/
static int resdopp_post(const obsd_t *obs, int n, const double *rs, const double *dts,
	const nav_t *nav, const double *rr, const double *x,
	const double *azel, const int *vsat, double *v, double *H,double *var)
{
	double lam,rate,pos[3],E[9],a[3]={0.0},e[3],vs[3],cosel=0.0,Vs,C1=0.0,dtrp=0.0,vtrp=0.0,decc=0.0,r;
	int i,j,nv=0,m;

	trace(3,"resdop  : n=%d\n",n);


	ecef2pos(rr,pos); xyz2enu(pos,E);

	for (i=0;i<n&&i<MAXOBS;i++) {

		lam=nav->lam[obs[i].sat-1][0];

		if (obs[i].D[0]==0.0||lam==0.0/*||!vsat[i]*/||norm(rs+3+i*6,3)<=0.0||azel[1+i*2]<10*D2R/*||obs[i].sat==60*/) 
		{
			continue;
		}
		/* line-of-sight vector in ecef */
		/*------------------------------------------------------------------------------------*/
		for (m=0;m<3;m++) e[m]=rs[m+i*6]-rr[m];
		r=norm(e,3);
		for (m=0;m<3;m++) e[m]/=r;
		/*------------------------------------------------------------------------------------*/
		/* satellite velocity relative to receiver in ecef */
		for (j=0;j<3;j++) vs[j]=rs[j+3+i*6]-x[j];

		/* range rate with earth rotation correction */
		rate=dot(vs,e,3)+OMGE/CLIGHT*(rs[4+i*6]*rr[0]+rs[1+i*6]*x[0]-
			rs[3+i*6]*rr[1]-rs[  i*6]*x[1]);

		/*------------------------------------------------------------------------------------*/
		/* doppler residual�������κθ�����ʱ�� */
		/*v[nv]=-lam*obs[i].D[0]-(rate+x[3]-CLIGHT*dts[1+i*2]);*/

		/*------------------------------------------------------------------------------------*/
		/*�����ٶ�V^2*/
		Vs=0.5*(SQR(rs[3+i*6])+SQR(rs[4+i*6])+SQR(rs[5+i*6]))/CLIGHT;
		/*����˫Ƶ���������Ӱ��*/
	 //   C1=1.79327;//E1��E5a��C1=1.64694;//f1^2/f2^2	
		//v[nv]=-(C1*lam*obs[i].D[0]-sqrt(C1*lam*lam)*obs[i].D[2])/(C1-1)-(rate+x[3]-CLIGHT*dts[1+i*2]+0.5*Vs/CLIGHT);
		
		v[nv]=-lam*obs[i].D[0]-(rate+x[3]-CLIGHT*dts[1+i*2]+0.5*Vs/CLIGHT);
		//if(obs[i].sat==83)trace(2,"%s Nsat=%2d  %6.4f  %4.2f %7.6f\n",time_str(obs[i].time,3),obs[i].sat,v[nv],azel[1+i*2]*R2D,
		//	lam*obs[i].D[0]-nav->lam[obs[i].sat-1][2]*obs[i].D[2]);

		/*------------------------------------------------------------------------------------*/
		/* error variance */
		if((obs[i].sat==73||obs[i].sat==77)&&vsat[i]!=0) var[nv]=SQR(0.002)/(2*sin(azel[1+i*2]))+SQR(0.004);
		else var[nv]=SQR(0.002)/(2*sin(azel[1+i*2]))+SQR(0.0005);

		/* design matrix */
		//for (j=0;j<4;j++) H[j+nv*4]=j<3?-e[j]:1.0;
		//for (j=0;j<3;j++) H[j+nv*5]=-e[j];
		//if(obs[i].sat<=35) H[3+nv*5]=1,H[4+nv*5]=0;
		//else  H[3+nv*5]=0,H[4+nv*5]=1;
		nv++;
	}
	return nv;
}
/* estimate receiver velocity ------------------------------------------------*/
static void estvel(obsd_t *obs, int n, const double *rs, const double *dts,
                   const nav_t *nav, const prcopt_t *opt, sol_t *sol,
                   const double *azel, const int *vsat,ssat_t *ssat)
{
    double x[4]={0},dx[4],Q[16],*v,*H,*var,sig,tt=0.0,tropR,ionoR;
    int i,j,nv,k;
    
    trace(3,"estvel  : n=%d\n",n);
	for(i=0;i<3;i++)x[i]=sol->rr[3+i];
/*------------------------------------------------------------------------------------*/
	/*ǰ��������Ԫ���쵼�������չ۲�ֵ*/
	for(i=0;i<n;i++)for(j=0;j<1;j++){//NFREQ
		tt=timediff(obs[i].time,ssat[obs[i].sat-1].pt[1][j]);
		if (ssat[obs[i].sat-1].slip[j]==1||ssat[obs[i].sat-1].ph[1][j]==0||fabs(tt)> 1|| obs[i].L[j] == 0 )
		{
			obs[i].D[j] = 0;
		}
		else{
			tropR=(ssat[obs[i].sat-1].trop-ssat[obs[i].sat-1].trop0)/(nav->lam[obs[i].sat-1][0]);
			ionoR=(ssat[obs[i].sat-1].iono-ssat[obs[i].sat-1].iono0)/(nav->lam[obs[i].sat-1][0]);
			if(fabs(tropR)>0.1||fabs(ionoR)>0.1)tropR=0,ionoR=0;
			obs[i].D[j]=(-obs[i].L[j]+ssat[obs[i].sat-1].ph[1][j]-tropR-ionoR)/tt;
		}
	}
/*------------------------------------------------------------------------------------*/
    v=mat(n,1); H=mat(4,n);var=mat(n,1);
    
    for (i=0;i<MAXITR;i++) {
        
        /* doppler residuals */
        if ((nv=resdopp(i,obs,n,rs,dts,nav,sol->rr,x,azel,vsat,v,H,var))<4) {
            break;
        }

		/* weight by variance */
		for (j=0;j<nv;j++) {
			sig=sqrt(var[j]);
			v[j]/=sig;
			for (k=0;k<4;k++) H[k+j*4]/=sig;
		}
        /* least square estimation */
        if (lsq(H,v,4,nv,dx,Q)) break;
        
        for (j=0;j<4;j++) x[j]+=dx[j];
        
        if (norm(dx,4)<1E-6) {
            for (j=0;j<3;j++) sol->rr[j+3]=x[j];
		/*	for(i=0;i<nv;i++)  trace(2,"  %.5f",v[i]);*/
			//trace(2,"  %.8f",x[3]);
			//trace(2," \n");
			//x[0]=0.0;x[1]=0.0;x[2]=0.0;x[3]=dx[3];
			//nv=resdopp_post(obs,n,rs,dts,nav,sol->rr,x,azel,vsat,v,H,var);
			//trace(2,"%s %12.4f\n",time_str(obs[i].time,3),x[3]);
            break;
        }
    }
    free(v); free(H); free(var);
//	double x[5]={0},dx[5],Q[25],QA[16],xa[4]={0.0},*v,*H,*var,*HA,*HB,*HC,tt=0.0,CoeX=0,CoeY=0,CoeZ=0,CoeT=0,sig=0,tropR,ionoR;
//	double dVxyz[3]={0},true_Pos[3]={0.0},pos[3]={0},dVenu[3]={0.0};extern double TruePos[3];
//    int i,j,nv,k,info;
//    
//    trace(3,"estvel  : n=%d\n",n);
///*------------------------------------------------------------------------------------*/
//	/*ǰ��������Ԫ���쵼�������չ۲�ֵ*/
//	for (i = 0; i < n; i++)for (j = 0; j < NFREQ; j++) {
//		tt = timediff(obs[i].time, ssat[obs[i].sat - 1].pt[1][j]);
//		if (ssat[obs[i].sat - 1].slip[j] == 1 || ssat[obs[i].sat - 1].ph[1][j] == 0 || fabs(tt) > 1 || obs[i].L[j] == 0)
//		{
//			obs[i].D[j] = 0;
//		}
//		else {
//			tropR = (ssat[obs[i].sat - 1].trop - ssat[obs[i].sat - 1].trop0) / (nav->lam[obs[i].sat - 1][0]);
//			ionoR = (ssat[obs[i].sat - 1].iono - ssat[obs[i].sat - 1].iono0) / (nav->lam[obs[i].sat - 1][0]);
//			if (fabs(tropR) > 0.1 || fabs(ionoR) > 0.1)tropR = 0, ionoR = 0;
//			obs[i].D[j] = (-obs[i].L[j] + ssat[obs[i].sat - 1].ph[1][j] - tropR - ionoR) / tt;
//		}
	//}
///*------------------------------------------------------------------------------------*/
//    v=mat(n,1); H=mat(5,n);var=mat(n,1);HA=mat(4,n);HB=mat(n,1);HC=mat(4,n);
//    
//    for (i=0;i<MAXITR;i++) {
//        
//        /* doppler residuals */
//        if ((nv=resdopp(i,obs,n,rs,dts,nav,sol->rr,x,azel,vsat,v,H,var))<5) {
//            break;
//        }
//	   /* weight by variance */
//		for (j=0;j<nv;j++) {
//			sig=sqrt(var[j]);
//			v[j]/=sig;
//			for (k=0;k<5;k++) H[k+j*5]/=sig;
//		}
//        /* least square estimation */
//        if (lsq(H,v,5,nv,dx,Q)) break;
//        
//        for (j=0;j<5;j++) x[j]+=dx[j];
//        
//        if (norm(dx,5)<1E-6) {
//            for (j=0;j<3;j++) sol->rr[j+3]=x[j];
//				trace(2,"%s %12.4f %12.4f\n",time_str(obs[i].time,3),x[3],x[4]);
//			break;
//        }
//    }
//    free(v); free(H);free(var);free(HA);free(HB);free(HC);
}
/* single-point positioning ----------------------------------------------------
* compute receiver position, velocity, clock bias by single-point positioning
* with pseudorange and doppler observables
* args   : obsd_t *obs      I   observation data
*          int    n         I   number of observation data
*          nav_t  *nav      I   navigation data
*          prcopt_t *opt    I   processing options
*          sol_t  *sol      IO  solution
*          double *azel     IO  azimuth/elevation angle (rad) (NULL: no output)
*          ssat_t *ssat     IO  satellite status              (NULL: no output)
*          char   *msg      O   error message for error exit
* return : status(1:ok,0:error)
* notes  : assuming sbas-gps, galileo-gps, qzss-gps, compass-gps time offset and
*          receiver bias are negligible (only involving glonass-gps time offset
*          and receiver bias)
*-----------------------------------------------------------------------------*/
extern int pntpos(const obsd_t *obs, int n, const nav_t *nav,
                  const prcopt_t *opt, sol_t *sol, double *azel, ssat_t *ssat,
                  char *msg)
{
	prcopt_t opt_ = *opt; 
    double *rs,*dts,*var,*azel_,*resp,ep[6]={0.0};
    int i,stat,vsat[MAXOBS]={0},svh[MAXOBS];
    
    sol->stat=SOLQ_NONE;
    
    if (n<=0) {strcpy(msg,"no observation data"); return 0;}
    
    trace(3,"pntpos  : tobs=%s n=%d\n",time_str(obs[0].time,3),n);
    
    sol->time=obs[0].time; msg[0]='\0';	

	//����ʱ����ʾ��ǰ��������ʱ��
	time2epoch(sol->time,ep); 
	if (ep[3]==5&&ep[4]==2&&fabs(ep[5]-2)<=0.1)
	{
		trace(3,"����ʱ����UTC�����ձ�ʾ�Ա�Ա�");
	}

	rs=mat(6,n); dts=mat(2,n); var=mat(1,n); azel_=zeros(2,n); resp=mat(1,n);
    
    if (opt_.mode!=PMODE_SINGLE||(opt_.mode==PMODE_SINGLE&&opt_.sateph==1)) { /* for precise positioning */
#if 0
        opt_.sateph =EPHOPT_BRDC;
#endif
        //opt_.ionoopt=IONOOPT_BRDC;
        opt_.ionoopt=IONOOPT_IFLC;
        opt_.tropopt=TROPOPT_SAAS;
    }
    /* satellite positons, velocities and clocks */
    satposs(sol->time,obs,n,nav,opt_.sateph,rs,dts,var,svh);
    
    /* estimate receiver position with pseudorange */
    stat=estpos(obs,n,rs,dts,var,svh,nav,&opt_,sol,azel_,vsat,resp,msg);

    /* raim fde */
    if (!stat&&n>=6&&opt->posopt[4]) {
        stat=raim_fde(obs,n,rs,dts,var,svh,nav,&opt_,sol,azel_,vsat,resp,msg);
    }

/*------------------------------------------------------------------------------------*/
    /* estimate receiver velocity with doppler��ԭʼ������ԭʼ�����չ۲�ֵ���ĳ��� */
    //if (stat) estvel(obs,n,rs,dts,nav,&opt_,sol,azel_,vsat);

   /*estimate receiver velocity with differential carrier phase,�ز���λ��ַ�����*/
	satpossDopp(sol->time,obs,n,nav,opt_.sateph,rs,dts,var,svh);
	if (stat) estvel(obs,n,rs,dts,nav,&opt_,sol,azel_,svh,ssat);//�޸�vsat->svh
/*------------------------------------------------------------------------------------*/

    if (azel) {
        for (i=0;i<n*2;i++) azel[i]=azel_[i];
    }
    if (ssat) {
        for (i=0;i<MAXSAT;i++) {
            ssat[i].vs=0;
            ssat[i].azel[0]=ssat[i].azel[1]=0.0;
            ssat[i].resp[0]=ssat[i].resc[0]=0.0;
            ssat[i].snr_rover[0]=0;
        }
        for (i=0;i<n;i++) {
            ssat[obs[i].sat-1].azel[0]=azel_[  i*2];
            ssat[obs[i].sat-1].azel[1]=azel_[1+i*2];
            ssat[obs[i].sat-1].snr_rover[0]=obs[i].SNR[0];
            if (!vsat[i]) continue;
            ssat[obs[i].sat-1].vs=1;
            ssat[obs[i].sat-1].resp[0]=resp[i];
        }
    }
    free(rs); free(dts); free(var); free(azel_); free(resp);
    return stat;
}
