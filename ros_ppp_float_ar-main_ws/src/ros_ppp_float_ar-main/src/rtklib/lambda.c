/*------------------------------------------------------------------------------
* lambda.c : integer ambiguity resolution
*
*          Copyright (C) 2007-2008 by T.TAKASU, All rights reserved.
*
* reference :
*     [1] P.J.G.Teunissen, The least-square ambiguity decorrelation adjustment:
*         a method for fast GPS ambiguity estimation, J.Geodesy, Vol.70, 65-82,
*         1995
*     [2] X.-W.Chang, X.Yang, T.Zhou, MLAMBDA: A modified LAMBDA method for
*         integer least-squares estimation, J.Geodesy, Vol.79, 552-565, 2005
*
* version : $Revision: 1.1 $ $Date: 2008/07/17 21:48:06 $
* history : 2007/01/13 1.0 new
*           2015/05/31 1.1 add api lambda_reduction(), lambda_search()
*-----------------------------------------------------------------------------*/
#include "rtklib.h"

/* constants/macros ----------------------------------------------------------*/

#define LOOPMAX     10000           /* maximum count of search loop */

#define SGN(x)      ((x)<=0.0?-1.0:1.0)
#define ROUND(x)    (floor((x)+0.5))
#define SWAP(x,y)   do {double tmp_; tmp_=x; x=y; y=tmp_;} while (0)

/* LD factorization (Q=L'*diag(D)*L) -----------------------------------------*/
static int LD(int n, const double *Q, double *L, double *D)
{
    int i,j,k,info=0;
    double a,*A=mat(n,n);
    
    memcpy(A,Q,sizeof(double)*n*n);
    for (i=n-1;i>=0;i--) {
        if ((D[i]=A[i+i*n])<=0.0) {info=-1; break;}
        a=sqrt(D[i]);
        for (j=0;j<=i;j++) L[i+j*n]=A[i+j*n]/a;
        for (j=0;j<=i-1;j++) for (k=0;k<=j;k++) A[j+k*n]-=L[i+k*n]*L[i+j*n];
        for (j=0;j<=i;j++) L[i+j*n]/=L[i+i*n];
    }
    free(A);
    if (info) fprintf(stderr,"%s : LD factorization error\n",__FILE__);
    return info;
}
/* integer gauss transformation ----------------------------------------------*/
static void gauss(int n, double *L, double *Z, int i, int j)
{
    int k,mu;
    
    if ((mu=(int)ROUND(L[i+j*n]))!=0) {
        for (k=i;k<n;k++) L[k+n*j]-=(double)mu*L[k+i*n];
        for (k=0;k<n;k++) Z[k+n*j]-=(double)mu*Z[k+i*n];
    }
}
/* permutations --------------------------------------------------------------*/
static void perm(int n, double *L, double *D, int j, double del, double *Z)
{
    int k;
    double eta,lam,a0,a1;
    
    eta=D[j]/del;
    lam=D[j+1]*L[j+1+j*n]/del;
    D[j]=eta*D[j+1]; D[j+1]=del;
    for (k=0;k<=j-1;k++) {
        a0=L[j+k*n]; a1=L[j+1+k*n];
        L[j+k*n]=-L[j+1+j*n]*a0+a1;
        L[j+1+k*n]=eta*a0+lam*a1;
    }
    L[j+1+j*n]=lam;
    for (k=j+2;k<n;k++) SWAP(L[k+j*n],L[k+(j+1)*n]);
    for (k=0;k<n;k++) SWAP(Z[k+j*n],Z[k+(j+1)*n]);
}
/* lambda reduction (z=Z'*a, Qz=Z'*Q*Z=L'*diag(D)*L) (ref.[1]) ---------------*/
static void reduction(int n, double *L, double *D, double *Z)
{
    int i,j,k;
    double del;
    
    j=n-2; k=n-2;
    while (j>=0) {
        if (j<=k) for (i=j+1;i<n;i++) gauss(n,L,Z,i,j);
        del=D[j]+L[j+1+j*n]*L[j+1+j*n]*D[j+1];
        if (del+1E-6<D[j+1]) { /* compared considering numerical error */
            perm(n,L,D,j,del,Z);
            k=j; j=n-2;
        }
        else j--;
    }
}
/* modified lambda (mlambda) search (ref. [2]) -------------------------------*/
static int search(int n, int m, const double *L, const double *D,
                  const double *zs, double *zn, double *s)
{
    int i,j,k,c,nn=0,imax=0;
    double newdist,maxdist=1E99,y;
    double *S=zeros(n,n),*dist=mat(n,1),*zb=mat(n,1),*z=mat(n,1),*step=mat(n,1);
    
    k=n-1; dist[k]=0.0;
    zb[k]=zs[k];
    z[k]=ROUND(zb[k]); y=zb[k]-z[k]; step[k]=SGN(y);
    for (c=0;c<LOOPMAX;c++) {
        newdist=dist[k]+y*y/D[k];
        if (newdist<maxdist) {
            if (k!=0) {
                dist[--k]=newdist;
                for (i=0;i<=k;i++)
                    S[k+i*n]=S[k+1+i*n]+(z[k+1]-zb[k+1])*L[k+1+i*n];
                zb[k]=zs[k]+S[k+k*n];
                z[k]=ROUND(zb[k]); y=zb[k]-z[k]; step[k]=SGN(y);
            }
            else {
                if (nn<m) {
                    if (nn==0||newdist>s[imax]) imax=nn;
                    for (i=0;i<n;i++) zn[i+nn*n]=z[i];
                    s[nn++]=newdist;
                }
                else {
                    if (newdist<s[imax]) {
                        for (i=0;i<n;i++) zn[i+imax*n]=z[i];
                        s[imax]=newdist;
                        for (i=imax=0;i<m;i++) if (s[imax]<s[i]) imax=i;
                    }
                    maxdist=s[imax];
                }
                z[0]+=step[0]; y=zb[0]-z[0]; step[0]=-step[0]-SGN(step[0]);
            }
        }
        else {
            if (k==n-1) break;
            else {
                k++;
                z[k]+=step[k]; y=zb[k]-z[k]; step[k]=-step[k]-SGN(step[k]);
            }
        }
    }
    for (i=0;i<m-1;i++) { /* sort by s */
        for (j=i+1;j<m;j++) {
            if (s[i]<s[j]) continue;
            SWAP(s[i],s[j]);
            for (k=0;k<n;k++) SWAP(zn[k+i*n],zn[k+j*n]);
        }
    }
    free(S); free(dist); free(zb); free(z); free(step);
    
    if (c>=LOOPMAX) {
        fprintf(stderr,"%s : search loop count overflow\n",__FILE__);
        return -1;
    }
    return 0;
}
/* lambda/mlambda integer least-square estimation ------------------------------
* integer least-square estimation. reduction is performed by lambda (ref.[1]),
* and search by mlambda (ref.[2]).
* args   : int    n      I  number of float parameters
*          int    m      I  number of fixed solutions
*          double *a     I  float parameters (n x 1)
*          double *Q     I  covariance matrix of float parameters (n x n)
*          double *F     O  fixed solutions (n x m)
*          double *s     O  sum of squared residulas of fixed solutions (1 x m)
* return : status (0:ok,other:error)
* notes  : matrix stored by column-major order (fortran convension)
*-----------------------------------------------------------------------------*/
extern int lambda(int n, int m, const double *a, const double *Q, double *F,
                  double *s)
{
    int info;
    double *L,*D,*Z,*z,*E;
    
    if (n<=0||m<=0) return -1;
    L=zeros(n,n); D=mat(n,1); Z=eye(n); z=mat(n,1); E=mat(n,m);
    
    /* LD factorization */
    if (!(info=LD(n,Q,L,D))) {
        
        /* lambda reduction */
        reduction(n,L,D,Z);
        matmul("TN",n,1,n,1.0,Z,a,0.0,z); /* z=Z'*a */

        /* mlambda search */
        if (!(info=search(n,m,L,D,z,E,s))) {
            
            info=solve("T",Z,E,n,m,F); /* F=Z'\E */
        }
    }
    free(L); free(D); free(Z); free(z); free(E);
    return info;
}
const double BootSt_2Phi_sigRecip_1[461] =
{
    1.00000000000, 1.00000000000, 1.00000000000, 1.00000000000, 1.00000000000, 0.99999999999, 0.99999999999, 0.99999999999, 0.99999999999, 0.99999999998,
    0.99999999997, 0.99999999996, 0.99999999995, 0.99999999994, 0.99999999992, 0.99999999989, 0.99999999985, 0.99999999981, 0.99999999975, 0.99999999968,
    0.99999999959, 0.99999999947, 0.99999999933, 0.99999999915, 0.99999999892, 0.99999999864, 0.99999999830, 0.99999999788, 0.99999999736, 0.99999999672,
    0.99999999596, 0.99999999502, 0.99999999390, 0.99999999255, 0.99999999092, 0.99999998898, 0.99999998667, 0.99999998393, 0.99999998068, 0.99999997684,
    0.99999997233, 0.99999996703, 0.99999996082, 0.99999995357, 0.99999994514, 0.99999993534, 0.99999992399, 0.99999991087, 0.99999989575, 0.99999987835,
    0.99999985840, 0.99999983555, 0.99999980946, 0.99999977971, 0.99999974588, 0.99999970748, 0.99999966399, 0.99999961482, 0.99999955936, 0.99999949690,
    0.99999942670, 0.99999934794, 0.99999925975, 0.99999916117, 0.99999905115, 0.99999892860, 0.99999879229, 0.99999864095, 0.99999847317, 0.99999828748,
    0.99999808226, 0.99999785581, 0.99999760630, 0.99999733179, 0.99999703020, 0.99999669931, 0.99999633680, 0.99999594015, 0.99999550675, 0.99999503379,
    0.99999451832, 0.99999395722, 0.99999334720, 0.99999268480, 0.99999196637, 0.99999118807, 0.99999034588, 0.99998943558, 0.99998845273, 0.99998739270,
    0.99998625064, 0.99998502150, 0.99998369998, 0.99998228058, 0.99998075754, 0.99997912490, 0.99997737643, 0.99997550567, 0.99997350593, 0.99997137025,
    0.99996909141, 0.99996666195, 0.99996407415, 0.99996132003, 0.99995839134, 0.99995527956, 0.99995197591, 0.99994847134, 0.99994475652, 0.99994082187,
    0.99993665752, 0.99993225331, 0.99992759885, 0.99992268343, 0.99991749610, 0.99991202562, 0.99990626047, 0.99990018888, 0.99989379879, 0.99988707788,
    0.99988001356, 0.99987259297, 0.99986480299, 0.99985663023, 0.99984806106, 0.99983908158, 0.99982967764, 0.99981983483, 0.99980953852, 0.99979877381,
    0.99978752559, 0.99977577847, 0.99976351689, 0.99975072503, 0.99973738684, 0.99972348608, 0.99970900630, 0.99969393082, 0.99967824278, 0.99966192513,
    0.99964496062, 0.99962733182, 0.99960902114, 0.99959001080, 0.99957028286, 0.99954981923, 0.99952860168, 0.99950661182, 0.99948383112, 0.99946024094,
    0.99943582250, 0.99941055690, 0.99938442516, 0.99935740815, 0.99932948670, 0.99930064151, 0.99927085322, 0.99924010239, 0.99920836952, 0.99917563504,
    0.99914187933, 0.99910708275, 0.99907122558, 0.99903428811, 0.99899625058, 0.99895709323, 0.99891679629, 0.99887533998, 0.99883270454, 0.99878887021,
    0.99874381725, 0.99869752596, 0.99864997666, 0.99860114972, 0.99855102555, 0.99849958461, 0.99844680744, 0.99839267462, 0.99833716683, 0.99828026482,
    0.99822194940, 0.99816220152, 0.99810100219, 0.99803833254, 0.99797417380, 0.99790850734, 0.99784131462, 0.99777257724, 0.99770227694, 0.99763039558,
    0.99755691517, 0.99748181789, 0.99740508603, 0.99732670206, 0.99724664863, 0.99716490852, 0.99708146471, 0.99699630033, 0.99690939872, 0.99682074337,
    0.99673031799, 0.99663810645, 0.99654409284, 0.99644826144, 0.99635059672, 0.99625108338, 0.99614970631, 0.99604645060, 0.99594130160, 0.99583424482,
    0.99572526604, 0.99561435123, 0.99550148659, 0.99538665856, 0.99526985380, 0.99515105920, 0.99503026190, 0.99490744925, 0.99478260885, 0.99465572855,
    0.99452679643, 0.99439580080, 0.99426273023, 0.99412757353, 0.99399031975, 0.99385095819, 0.99370947841, 0.99356587020, 0.99342012360, 0.99327222891,
    0.99312217668, 0.99296995769, 0.99281556299, 0.99265898388, 0.99250021191, 0.99233923886, 0.99217605681, 0.99201065803, 0.99184303508, 0.99167318077,
    0.99150108813, 0.99132675048, 0.99115016135, 0.99097131455, 0.99079020413, 0.99060682436, 0.99042116979, 0.99023323520, 0.99004301562, 0.98985050631,
    0.98965570280, 0.98945860082, 0.98925919638, 0.98905748571, 0.98885346527, 0.98864713177, 0.98843848214, 0.98822751357, 0.98801422345, 0.98779860942,
    0.98758066935, 0.98736040131, 0.98713780363, 0.98691287485, 0.98668561372, 0.98645601923, 0.98622409057, 0.98598982715, 0.98575322862, 0.98551429480,
    0.98527302576, 0.98502942175, 0.98478348325, 0.98453521091, 0.98428460563, 0.98403166847, 0.98377640072, 0.98351880383, 0.98325887948, 0.98299662952,
    0.98273205600, 0.98246516116, 0.98219594742, 0.98192441738, 0.98165057382, 0.98137441971, 0.98109595820, 0.98081519259, 0.98053212639, 0.98024676323,
    0.97995910696, 0.97966916156, 0.97937693119, 0.97908242016, 0.97878563295, 0.97848657419, 0.97818524867, 0.97788166132, 0.97757581723, 0.97726772163,
    0.97695737991, 0.97664479759, 0.97632998033, 0.97601293393, 0.97569366433, 0.97537217760, 0.97504847994, 0.97472257770, 0.97439447732, 0.97406418539,
    0.97373170862, 0.97339705384, 0.97306022799, 0.97272123813, 0.97238009145, 0.97203679524, 0.97169135688, 0.97134378388, 0.97099408386, 0.97064226453,
    0.97028833371, 0.96993229931, 0.96957416934, 0.96921395192, 0.96885165523, 0.96848728757, 0.96812085733, 0.96775237296, 0.96738184302, 0.96700927615,
    0.96663468106, 0.96625806654, 0.96587944147, 0.96549881479, 0.96511619553, 0.96473159277, 0.96434501569, 0.96395647350, 0.96356597551, 0.96317353108,
    0.96277914962, 0.96238284062, 0.96198461363, 0.96158447823, 0.96118244410, 0.96077852092, 0.96037271847, 0.95996504656, 0.95955551504, 0.95914413383,
    0.95873091286, 0.95831586215, 0.95789899172, 0.95748031167, 0.95705983210, 0.95663756317, 0.95621351508, 0.95578769805, 0.95536012234, 0.95493079825,
    0.95449973610, 0.95406694624, 0.95363243906, 0.95319622495, 0.95275831434, 0.95231871770, 0.95187744550, 0.95143450824, 0.95098991643, 0.95054368061,
    0.95009581134, 0.94964631918, 0.94919521471, 0.94874250853, 0.94828821126, 0.94783233351, 0.94737488592, 0.94691587911, 0.94645532375, 0.94599323049,
    0.94552960997, 0.94506447288, 0.94459782987, 0.94412969162, 0.94366006879, 0.94318897207, 0.94271641211, 0.94224239960, 0.94176694519, 0.94129005956,
    0.94081175336, 0.94033203725, 0.93985092188, 0.93936841788, 0.93888453589, 0.93839928654, 0.93791268045, 0.93742472821, 0.93693544042, 0.93644482767,
    0.93595290052, 0.93545966953, 0.93496514524, 0.93446933818, 0.93397225886, 0.93347391778, 0.93297432541, 0.93247349222, 0.93197142864, 0.93146814510,
    0.93096365201, 0.93045795974, 0.92995107865, 0.92944301910, 0.92893379140, 0.92842340583, 0.92791187269, 0.92739920221, 0.92688540462, 0.92637049012,
    0.92585446889, 0.92533735107, 0.92481914678, 0.92429986613, 0.92377951918, 0.92325811596, 0.92273566650, 0.92221218077, 0.92168766874, 0.92116214032,
    0.92063560541, 0.92010807387, 0.91957955554, 0.91905006022, 0.91851959769, 0.91798817767, 0.91745580989, 0.91692250401, 0.91638826967, 0.91585311649,
    0.91531705404, 0.91478009186, 0.91424223946, 0.91370350631, 0.91316390185, 0.91262343548, 0.91208211658, 0.91153995446, 0.91099695844, 0.91045313777,
    0.90990850167, 0.90936305934, 0.90881681992, 0.90826979253, 0.90772198625, 0.90717341012, 0.90662407314, 0.90607398427, 0.90552315245, 0.90497158656,
    0.90441929545
};
double pBootStrapping(double *sig)
{
    int i,j;
    double dt,dRes=0.0;

    if (*sig>=0.3 ) dRes=0.9;
    else if (*sig<=0.07) dRes=1.0;
    else
    {
        dt=(*sig-0.07)*2000.0; 
        i = round(dt);
        j = i+1;

        dRes=(j-dt)*BootSt_2Phi_sigRecip_1[i]+(dt - i)*BootSt_2Phi_sigRecip_1[j];
    }
    return dRes;
}
/* mlambda integer least-square estimation ------------------------------
* integer least-square estimation. reduction is performed by lambda (ref.[1]),
* and search by mlambda (ref.[2]).
* args   : int    n      I  number of float parameters
*          int    m      I  number of fixed solutions
*          double *a     I  float parameters (n x 1)
*          double *Q     I  covariance matrix of float parameters (n x n)
*          double *F     O  fixed solutions (n x m)
*          double *s     O  sum of squared residulas of fixed solutions (1 x m)
* return : status (0:ok,other:error)
* notes  : matrix stored by column-major order (fortran convension)
*-----------------------------------------------------------------------------*/
extern int mlambda(int n, int m, const double *a, const double *Q, double *F,
                  double *s,double *boot)
{
    int info,i;
    double *L,*D,*Z,*z,*E,dt=0.0,sr=1.0;
    
    if (n<=0||m<=0) return -1;
    L=zeros(n,n); D=mat(n,1); Z=eye(n); z=mat(n,1); E=mat(n,m);
    
    /* LD factorization */
    if (!(info=LD(n,Q,L,D))) {
        
        /* lambda reduction */
        reduction(n,L,D,Z);
        matmul("TN",n,1,n,1.0,Z,a,0.0,z); /* z=Z'*a */

        for(i=0;i<n;i++){
            dt = sqrt(D[i]);
			sr *= pBootStrapping(&dt);
        }
        *boot=sr;
        /* mlambda search */
        if (!(info=search(n,m,L,D,z,E,s))) {
            
            info=solve("T",Z,E,n,m,F); /* F=Z'\E */
        }
    }
    free(L); free(D); free(Z); free(z); free(E);
    return info;
}
/* lambda reduction ------------------------------------------------------------
* reduction by lambda (ref [1]) for integer least square
* args   : int    n      I  number of float parameters
*          double *Q     I  covariance matrix of float parameters (n x n)
*          double *Z     O  lambda reduction matrix (n x n)
* return : status (0:ok,other:error)
*-----------------------------------------------------------------------------*/
extern int lambda_reduction(int n, const double *Q, double *Z)
{
    double *L,*D;
    int i,j,info;
    
    if (n<=0) return -1;
    
    L=zeros(n,n); D=mat(n,1);
    
    for (i=0;i<n;i++) for (j=0;j<n;j++) {
        Z[i+j*n]=i==j?1.0:0.0;
    }
    /* LD factorization */
    if ((info=LD(n,Q,L,D))) {
        free(L); free(D);
        return info;
    }
    /* lambda reduction */
    reduction(n,L,D,Z);
     
    free(L); free(D);
    return 0;
}
/* mlambda search --------------------------------------------------------------
* search by  mlambda (ref [2]) for integer least square
* args   : int    n      I  number of float parameters
*          int    m      I  number of fixed solutions
*          double *a     I  float parameters (n x 1)
*          double *Q     I  covariance matrix of float parameters (n x n)
*          double *F     O  fixed solutions (n x m)
*          double *s     O  sum of squared residulas of fixed solutions (1 x m)
* return : status (0:ok,other:error)
*-----------------------------------------------------------------------------*/
extern int lambda_search(int n, int m, const double *a, const double *Q,
                         double *F, double *s)
{
    double *L,*D;
    int info;
    
    if (n<=0||m<=0) return -1;
    
    L=zeros(n,n); D=mat(n,1);
    
    /* LD factorization */
    if ((info=LD(n,Q,L,D))) {
        free(L); free(D);
        return info;
    }
    /* mlambda search */
    info=search(n,m,L,D,a,F,s);
    
    free(L); free(D);
    return info;
}
