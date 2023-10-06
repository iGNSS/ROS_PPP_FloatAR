/*------------------------------------------------------------------------------
* tides.c : tidal displacement corrections
*
*          Copyright (C) 2015 by T.TAKASU, All rights reserved.
*
* options : -DIERS_MODEL use IERS tide model
*
* references :
*     [1] D.D.McCarthy, IERS Technical Note 21, IERS Conventions 1996, July 1996
*     [2] D.D.McCarthy and G.Petit, IERS Technical Note 32, IERS Conventions
*         2003, November 2003
*     [3] D.A.Vallado, Fundamentals of Astrodynamics and Applications 2nd ed,
*         Space Technology Library, 2004
*     [4] J.Kouba, A Guide to using International GNSS Service (IGS) products,
*         May 2009
*     [5] G.Petit and B.Luzum (eds), IERS Technical Note No. 36, IERS
*         Conventions (2010), 2010
*
* version : $Revision:$ $Date:$
* history : 2015/05/10 1.0  separated from ppp.c
*           2015/06/11 1.1  fix bug on computing days in tide_oload() (#128)
*-----------------------------------------------------------------------------*/
#include "rtklib.h"

static const char rcsid[] = "$Id:$";

#define SQR(x)      ((x)*(x))

#define AS2R        (D2R/3600.0)    /* arc sec to radian */
#define GME         3.986004415E+14 /* earth gravitational constant */
#define GMS         1.327124E+20    /* sun gravitational constant */
#define GMM         4.902801E+12    /* moon gravitational constant */

/* function prototypes -------------------------------------------------------*/
#ifdef IERS_MODEL
extern int dehanttideinel_(double *xsta, int *year, int *mon, int *day,
	double *fhr, double *xsun, double *xmon,
	double *DXTIDE);

/* ------------------------------------------------------------------------
* compute the out-of-phase corrections induced by
*  mantle anelasticity in the diurnal band
* args   : double *RecCorXYZ      I   receiver position {x,y,z} (m)
*		   double *SunCoorXYZ     I   sun position {x,y,z} (m)
*		   double *MoonCoorXYZ    I   moon position {x,y,z} (m)
*          double FAC2SUN	      I   Degree 2 TGP factor for the Sun
*          double FAC2MON         I   Degree 2 TGP factor for the Moon
* return : XCORSTA				  IO  Out of phase station corrections for diurnal band
* note   : see ref [5] ftp://tai.bipm.org/iers/conv2010/chapter7/dehanttideinel/
*-----------------------------------------------------------------------------*/
void ST1IDIU(double *RecCorXYZ, double *SunCoorXYZ, double *MoonCoorXYZ, double FAC2SUN, double FAC2MON, double *XCORSTA)
{
	/* THIS SUBROUTINE GIVES THE OUT-OF-PHASE CORRECTIONS INDUCED BY
	 MANTLE INELASTICITY IN THE DIURNAL BAND
	
	       INPUT : XSTA,XSUN,XMON,FAC2SUN,FAC2MON
	      OUTPUT : XCORSTA
	*/
	double DHI = -0.0025;
	double DLI = -0.0007;
	double RSTA = sqrt(pow(RecCorXYZ[0], 2) + pow(RecCorXYZ[1], 2) + pow(RecCorXYZ[3], 2));
	double SINPHI = RecCorXYZ[2] / RSTA;
	double COSPHI = sqrt(pow(RecCorXYZ[0], 2) + pow(RecCorXYZ[1], 2)) / RSTA;
	double COS2PHI = pow(COSPHI, 2) - pow(SINPHI, 2);
	double SINLA = RecCorXYZ[1] / SINPHI / RSTA;
	double COSLA = RecCorXYZ[0] / COSPHI / RSTA;

	double
		RSUN = sqrt(pow(SunCoorXYZ[0], 2) + pow(SunCoorXYZ[1], 2) + pow(SunCoorXYZ[2], 2));
	double
		RMON = sqrt(pow(MoonCoorXYZ[0], 2) + pow(MoonCoorXYZ[1], 2) + pow(MoonCoorXYZ[2], 2));

	double XSUN[3] = { SunCoorXYZ[0], SunCoorXYZ[1], SunCoorXYZ[2] };
	double XMON[3] = { MoonCoorXYZ[0], MoonCoorXYZ[1], MoonCoorXYZ[2] };
	double
		DRSUN = -3.0*DHI*SINPHI*COSPHI*FAC2SUN*XSUN[2] * (XSUN[0] * SINLA - XSUN[1] * COSLA) / pow(RSUN, 2);
	double
		DRMON = -3.0*DHI*SINPHI*COSPHI*FAC2MON*XMON[2] * (XMON[0] * SINLA - XMON[1] * COSLA) / pow(RMON, 2);
	double
		DNSUN = -3.0*DLI*COS2PHI*FAC2SUN*XSUN[2] * (XSUN[0] * SINLA - XSUN[1] * COSLA) / pow(RSUN, 2);
	double
		DNMON = -3.0*DLI*COS2PHI*FAC2MON*XMON[2] * (XMON[0] * SINLA - XMON[1] * COSLA) / pow(RMON, 2);
	double
		DESUN = -3.0*DLI*SINPHI*FAC2SUN*XSUN[2] * (XSUN[0] * COSLA + XSUN[1] * SINLA) / pow(RSUN, 2);
	double
		DEMON = -3.0*DLI*SINPHI*FAC2MON*XMON[2] * (XMON[0] * COSLA + XMON[1] * SINLA) / pow(RMON, 2);

	double DR = DRSUN + DRMON;
	double DN = DNSUN + DNMON;
	double DE = DESUN + DEMON;
	XCORSTA[0] = DR*COSLA*COSPHI - DE*SINLA - DN*SINPHI*COSLA;
	XCORSTA[1] = DR*SINLA*COSPHI + DE*COSLA - DN*SINPHI*SINLA;
	XCORSTA[2] = DR*SINPHI + DN*COSPHI;
}

/* ------------------------------------------------------------------------
* compute the out-of-phase corrections induced by
*  mantle anelasticity in the semi-diurnal band
* args   : double *RecCorXYZ      I   receiver position {x,y,z} (m)
*		   double *SunCoorXYZ     I   sun position {x,y,z} (m)
*		   double *MoonCoorXYZ    I   moon position {x,y,z} (m)
*          double FAC2SUN	      I   Degree 2 TGP factor for the Sun
*          double FAC2MON         I   Degree 2 TGP factor for the Moon
* return : XCORSTA				  IO  Out of phase station corrections for semi-diurnal band
* note   : see ref [5] ftp://tai.bipm.org/iers/conv2010/chapter7/dehanttideinel/
*-----------------------------------------------------------------------------*/
void ST1ISEM(double *RecCorXYZ, double *SunCoorXYZ, double *MoonCoorXYZ, double FAC2SUN, double FAC2MON, double *XCORSTA)
{
	/*
	 THIS SUBROUTINE GIVES THE OUT-OF-PHASE CORRECTIONS INDUCED BY
	 MANTLE INELASTICITY IN THE DIURNAL BAND
	
	       INPUT : XSTA,XSUN,XMON,FAC2SUN,FAC2MON
	      OUTPUT : XCORSTA
	*/
	double DHI = -0.0022;
	double DLI = -0.0007;
	double RSTA = sqrt(pow(RecCorXYZ[0], 2) + pow(RecCorXYZ[1], 2) + pow(RecCorXYZ[2], 2));
	double SINPHI = RecCorXYZ[2] / RSTA;
	double COSPHI = sqrt(pow(RecCorXYZ[0], 2) + pow(RecCorXYZ[1], 2)) / RSTA;
	double COS2PHI = pow(COSPHI, 2) - pow(SINPHI, 2);
	double SINLA = RecCorXYZ[1] / COSPHI / RSTA;
	double COSLA = RecCorXYZ[0] / COSPHI / RSTA;
	double COSTWOLA = pow(COSLA, 2) - pow(SINLA, 2);
	double SINTWOLA = 2.0*COSLA*SINLA;
	double
		RSUN = sqrt(pow(SunCoorXYZ[0], 2) + pow(SunCoorXYZ[1], 2) + pow(SunCoorXYZ[2], 2));
	double
		RMON = sqrt(pow(MoonCoorXYZ[0], 2) + pow(MoonCoorXYZ[1], 2) + pow(MoonCoorXYZ[2], 2));


	double XSUN[3] = { SunCoorXYZ[0], SunCoorXYZ[1], SunCoorXYZ[2] };
	double XMON[3] = { MoonCoorXYZ[0], MoonCoorXYZ[1], MoonCoorXYZ[2] };

	double
		DRSUN = -3.0 / 4.0*DHI*pow(COSPHI, 2)*FAC2SUN*((pow(XSUN[0], 2) - pow(XSUN[1], 2))*SINTWOLA - 2.*XSUN[0] * XSUN[1] * COSTWOLA) / pow(RSUN, 2);
	double
		DRMON = -3.0 / 4.0*DHI*pow(COSPHI, 2)*FAC2MON*((pow(XMON[0], 2) - pow(XMON[1], 2))*
		SINTWOLA - 2.*XMON[0] * XMON[1] * COSTWOLA) / pow(RMON, 2);
	double
		DNSUN = 3.0 / 2.0*DLI*SINPHI*COSPHI*FAC2SUN*((pow(XSUN[0], 2) - pow(XSUN[1], 2))*SINTWOLA - 2.*XSUN[0] * XSUN[1] * COSTWOLA) / pow(RSUN, 2);
	double
		DNMON = 3.0 / 2.0*DLI*SINPHI*COSPHI*FAC2MON*((pow(XMON[0], 2) - pow(XMON[1], 2))*SINTWOLA - 2.*XMON[0] * XMON[1] * COSTWOLA) / pow(RMON, 2);
	double
		DESUN = -3.0 / 2.0*DLI*COSPHI*FAC2SUN*((pow(XSUN[0], 2) - pow(XSUN[1], 2))*COSTWOLA + 2.*XSUN[0] * XSUN[1] * SINTWOLA) / pow(RSUN, 2);
	double
		DEMON = -3.0 / 2.0*DLI*COSPHI*FAC2MON*((pow(XMON[0], 2) - pow(XMON[1], 2))*COSTWOLA + 2.*XMON[0] * XMON[1] * SINTWOLA) / pow(RMON, 2);
	double DR = DRSUN + DRMON;
	double DN = DNSUN + DNMON;
	double DE = DESUN + DEMON;
	XCORSTA[0] = DR*COSLA*COSPHI - DE*SINLA - DN*SINPHI*COSLA;
	XCORSTA[1] = DR*SINLA*COSPHI + DE*COSLA - DN*SINPHI*SINLA;
	XCORSTA[2] = DR*SINPHI + DN*COSPHI;
}

/* ------------------------------------------------------------------------
* compute the corrections induced by the latitude
*  dependence given by L^1 in Mathews et al. 1991 (See References).
* args   : double *RecCorXYZ      I   receiver position {x,y,z} (m)
*		   double *SunCoorXYZ     I   sun position {x,y,z} (m)
*		   double *MoonCoorXYZ    I   moon position {x,y,z} (m)
*          double FAC2SUN	      I   Degree 2 TGP factor for the Sun
*          double FAC2MON         I   Degree 2 TGP factor for the Moon
* return : XCORSTA				  IO  Out of phase station corrections for semi-diurnal band
* note   : see ref [5] ftp://tai.bipm.org/iers/conv2010/chapter7/dehanttideinel/
*-----------------------------------------------------------------------------*/
void ST1L1(double *RecCorXYZ, double *SunCoorXYZ, double *MoonCoorXYZ, double FAC2SUN, double FAC2MON, double *XCORSTA)
{
	/*
	THIS SUBROUTINE GIVES THE CORRECTIONS INDUCED BY THE LATITUDE DEPENDENCE
	GIVEN BY L^(1) IN MAHTEWS ET AL (1991)
	
	      INPUT : XSTA,XSUN,XMON,FAC3SUN,FAC3MON
	     OUTPUT : XCORSTA
	*/
	double L1D = 0.0012;
	double L1SD = 0.0024;
	double RSTA = sqrt(pow(RecCorXYZ[0], 2) + pow(RecCorXYZ[1], 2) + pow(RecCorXYZ[2], 2));
	double SINPHI = RecCorXYZ[2] / RSTA;
	double COSPHI = sqrt(pow(RecCorXYZ[0], 2) + pow(RecCorXYZ[1], 2)) / RSTA;
	double SINLA = RecCorXYZ[1] / COSPHI / RSTA;
	double COSLA = RecCorXYZ[0] / COSPHI / RSTA;
	double
		RSUN = sqrt(pow(SunCoorXYZ[0], 2) + pow(SunCoorXYZ[1], 2) + pow(SunCoorXYZ[2], 2));
	double
		RMON = sqrt(pow(MoonCoorXYZ[0], 2) + pow(MoonCoorXYZ[1], 2) + pow(MoonCoorXYZ[2], 2));

	double XSUN[3] = { SunCoorXYZ[0], SunCoorXYZ[1], SunCoorXYZ[2] };
	double XMON[3] = { MoonCoorXYZ[0], MoonCoorXYZ[1], MoonCoorXYZ[2] };
	/*
	 FOR THE DIURNAL BAND
	*/
	double  L1 = L1D;
	double
		DNSUN = -L1*pow(SINPHI, 2)*FAC2SUN*XSUN[2] * (XSUN[0] * COSLA + XSUN[1] * SINLA) / pow(RSUN, 2);
	double
		DNMON = -L1*pow(SINPHI, 2)*FAC2MON*XMON[2] * (XMON[0] * COSLA + XMON[1] * SINLA) / pow(RMON, 2);
	double
		DESUN = L1*SINPHI*(pow(COSPHI, 2) - pow(SINPHI, 2))*FAC2SUN*XSUN[2] * (XSUN[0] * SINLA - XSUN[1] * COSLA) / pow(RSUN, 2);
	double
		DEMON = L1*SINPHI*(pow(COSPHI, 2) - pow(SINPHI, 2))*FAC2MON*XMON[2] * (XMON[0] * SINLA - XMON[1] * COSLA) / pow(RMON, 2);
	double DE = 3.0*(DESUN + DEMON);
	double DN = 3.0*(DNSUN + DNMON);
	XCORSTA[0] = -DE*SINLA - DN*SINPHI*COSLA;
	XCORSTA[1] = DE*COSLA - DN*SINPHI*SINLA;
	XCORSTA[2] = DN*COSPHI;
	/*
	 FOR THE SEMI-DIURNAL BAND
	*/
	L1 = L1SD;
	double COSTWOLA = pow(COSLA, 2) - pow(SINLA, 2);
	double SINTWOLA = 2.0*COSLA*SINLA;
	DNSUN = -L1 / 2.0*SINPHI*COSPHI*FAC2SUN*((pow(XSUN[0], 2) - pow(XSUN[1], 2))*COSTWOLA + 2.0*XSUN[0] * XSUN[1] * SINTWOLA) / pow(RSUN, 2);
	DNMON = -L1 / 2.0*SINPHI*COSPHI*FAC2MON*((pow(XMON[0], 2) - pow(XMON[1], 2))*COSTWOLA + 2.0*XMON[0] * XMON[1] * SINTWOLA) / pow(RMON, 2);
	DESUN = -L1 / 2.0*pow(SINPHI, 2)*COSPHI*FAC2SUN*((pow(XSUN[0], 2) - pow(XSUN[1], 2))*SINTWOLA - 2.0*XSUN[0] * XSUN[1] * COSTWOLA) / pow(RSUN, 2);
	DEMON = -L1 / 2.0*pow(SINPHI, 2)*COSPHI*FAC2MON*((pow(XMON[0], 2) - pow(XMON[1], 2))*SINTWOLA - 2.0*XMON[0] * XMON[1] * COSTWOLA) / pow(RMON, 2);
	DE = 3.0*(DESUN + DEMON);
	DN = 3.0*(DNSUN + DNMON);
	XCORSTA[0] = XCORSTA[0] - DE*SINLA - DN*SINPHI*COSLA;
	XCORSTA[1] = XCORSTA[1] + DE*COSLA - DN*SINPHI*SINLA;
	XCORSTA[2] = XCORSTA[2] + DN*COSPHI;
}

/*----------------------------------------------------------------------*/
void STEP2LON(double *RecCorXYZ, double fhr, double t, double *xcorsta)
{
	double deg2rad = 0.0174532925199432957690;

	/* *** cf. table 7.5b of IERS conventions 2003 (TN.32, pg.82)
	 *** columns are s,h,p,N',ps, dR(ip),dT(ip),dR(op),dT(op)
	 *** IERS cols.= s,h,p,N',ps, dR(ip),dR(op),dT(ip),dT(op)
	 *** units of mm*/
	int i, j;
	double
		datdi[9][5] = { { 0, 0, 1, 2, 2 }, { 0, 2, 0, 0, 0 }, { 0, 0, -1, 0, 0 }, { 1, 0, 0, 0, 1 }, { 0, 0, 0, 0, 0 },
		{ 0.47, -0.20, -0.11, -0.13, -0.05 }, { 0.23, -0.12, -0.08, -0.11, -0.05 },
		{ 0.16, -0.11, -0.09, -0.15, -0.06 }, { 0.07, -0.05, -0.04, -0.07, -0.03 } };

	double s = 218.316645630 + 481267.881940*t
		- 0.00146638890*t*t + 0.000001851390*t*t*t;
	double  pr = 1.396971278*t + 0.000308889*t*t
		+ 0.000000021*t*t*t + 0.000000007*t*t*t*t;
	s = s + pr;
	double h = 280.466450 + 36000.76974890*t
		+ 0.000303222220*t*t + 0.000000020*t*t*t - 0.00000000654*t*t*t*t;
	double p = 83.353243120
		+ 4069.013635250*t - 0.010321722220*t*t - 0.00001249910*t*t*t
		+ 0.000000052630*t*t*t*t;
	double zns = 234.955444990 + 1934.136261970*t
		- 0.002075611110*t*t - 0.000002139440*t*t*t + 0.000000016500*t*t*t*t;
	double ps = 282.937340980 + 1.719457666670*t
		+ 0.000456888890*t*t - 0.000000017780*t*t*t - 0.000000003340*t*t*t*t;

	double rsta = sqrt(pow(RecCorXYZ[0], 2) + pow(RecCorXYZ[1], 2) + pow(RecCorXYZ[2], 2));
	double sinphi = RecCorXYZ[2] / rsta;
	double cosphi = sqrt(pow(RecCorXYZ[0], 2) + pow(RecCorXYZ[1], 2)) / rsta;
	double sinla = RecCorXYZ[1] / sinphi / rsta;
	double cosla = RecCorXYZ[0] / cosphi / rsta;
	/*  reduce angles to between 0 and 360*/
	s = fmod(s, 360.0);
	h = fmod(h, 360.0);
	p = fmod(p, 360.0);
	zns = fmod(zns, 360.0);
	ps = fmod(ps, 360.0);
	double dr_tot = 0.0;
	double dn_tot = 0.0;

	for (i = 1; i<3; i++)
	{
		xcorsta[i] = 0.0;
	}


	/*             1 2 3 4   5   6      7      8      9
	  columns are s,h,p,N',ps, dR(ip),dT(ip),dR(op),dT(op)*/
	for (j = 0; j<5; j++)
	{
		double thetaf =
			(datdi[0][j] * s + datdi[1][j] * h + datdi[2][j] * p + datdi[3][j] * zns + datdi[4][j] * ps)*deg2rad;
		double dr = datdi[5][j] * (3.0*pow(sinphi, 2) - 1.0) / 2.*cos(thetaf) +
			datdi[7][j] * (3.0*pow(sinphi, 2) - 1.0) / 2.*sin(thetaf);
		double dn =
			datdi[6][j] * (cosphi*sinphi*2.0)*cos(thetaf) + datdi[8][j] * (cosphi*sinphi*2.0)*sin(thetaf);
		double de = 0.0;
		dr_tot = dr_tot + dr;
		dn_tot = dn_tot + dn;

		xcorsta[0] = xcorsta[0] + dr*cosla*cosphi - de*sinla - dn*sinphi*cosla;
		xcorsta[1] = xcorsta[1] + dr*sinla*cosphi + de*cosla - dn*sinphi*sinla;
		xcorsta[2] = xcorsta[2] + dr*sinphi + dn*cosphi;
	}

	for (i = 0; i<3; i++)
	{
		xcorsta[i] = xcorsta[i] / 1000.0;
	}
}
/*------------------------------------------------------------------------*/
void STEP2DIU(double *RecCorXYZ, double fhr, double t, double *XCORSTA)
{
	int i;
	double deg2rad = 0.0174532925199432957690;

	double datdi[9][31] = {

		{ -3.000, -3.000, -2.000, -2.000, -2.000, -1.000, -1.000, -1.000, 0.000, 0.000, 0.000, 0.000, 0.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 2.000, 2.000, 3.000, 3.000 },

		{ 0.000, 2.000, 0.000, 0.000, 2.000, 0.000, 0.000, 2.000, -2.000, 0.000, 0.000, 0.000, 2.000, -3.000, -2.000, -2.000, -1.000, -1.000, 0.000, 0.000, 0.000, 0.000, 1.000, 1.000, 1.000, 2.000, 2.000, -2.000, 0.000, 0.000, 0.000 },

		{ 2.000, 0.000, 1.000, 1.000, -1.000, 0.000, 0.000, 0.000, 1.000, -1.000, 1.000, 1.000, -1.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, -2.000, 0.000, 1.000, -1.000, 0.000, 0.000 },

		{ 0.000, 0.000, -1.000, 0.000, 0.000, -1.000, 0.000, 0.000, 0.000, 0.000, 0.000, 1.000, 0.000, 0.000, 1.000, 0.000, 0.000, 0.000, -1.000, 0.000, 1.000, 2.000, 0.000, 0.000, 1.000, 0.000, 0.000, 0.000, 0.000, 0.000, 1.000 },

		{ 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 1.000, 0.000, 0.000, -1.000, 1.000, 0.000, 0.000, 0.000, 0.000, -1.000, 1.000, -1.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 },

		{ -0.010, -0.010, -0.020, -0.080, -0.020, -0.100, -0.510, 0.010, 0.010, 0.020, 0.060, 0.010, 0.010, -0.060, 0.010, -1.230, 0.020, 0.040, -0.220, 12.000, 1.730, -0.040, -0.500, 0.010, -0.010, -0.010, -0.110, -0.010, -0.020, 0.000, 0.000 },

		{ -0.010, -0.010, -0.010, 0.000, -0.010, 0.000, 0.000, 0.000, 0.000, 0.010, 0.000, 0.000, 0.000, 0.000, 0.000, -0.070, 0.000, 0.000, 0.010, -0.780, -0.120, 0.000, -0.010, 0.000, 0.000, 0.000, 0.010, 0.000, 0.020, 0.010, 0.010 },

		{ 0.000, 0.000, 0.000, 0.010, 0.000, 0.000, -0.020, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.060, 0.000, 0.000, 0.010, -0.670, -0.100, 0.000, 0.030, 0.000, 0.000, 0.000, 0.010, 0.000, 0.000, 0.000, 0.000 },

		{ 0.000, 0.000, 0.000, 0.010, 0.000, 0.000, 0.030, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.010, 0.000, 0.000, 0.000, -0.030, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.010, 0.010, 0.000 } };


	double s = 218.316645630 + 481267.881940*t
		- 0.00146638890*t*t + 0.00000185139*t*t*t;
	double tau = fhr*15.0 + 280.46061840 + 36000.77005360*t
		+ 0.000387930*t*t - 0.00000002580*t*t*t - s;
	double pr = 1.396971278*t + 0.000308889*t*t
		+ 0.000000021*t*t*t + 0.000000007*t*t*t*t;
	s = s + pr;
	double  h = 280.466450 + 36000.76974890*t
		+ 0.000303222220*t*t + 0.000000020*t*t*t - 0.00000000654*t*t*t*t;
	double   p = 83.353243120 + 4069.013635250*t
		- 0.010321722220*t*t - 0.00001249910*t*t*t + 0.000000052630*t*t*t*t;
	double zns = 234.955444990 + 1934.136261970*t
		- 0.002075611110*t*t - 0.000002139440*t*t*t + 0.000000016500*t*t*t*t;
	double  ps = 282.937340980 + 1.719457666670*t
		+ 0.00045688889*t*t - 0.000000017780*t*t*t - 0.000000003340*t*t*t*t;

	/**** reduce angles to between 0 and 360*/

	s = fmod(s, 360.0);
	tau = fmod(tau, 360.0);
	h = fmod(h, 360.0);
	p = fmod(p, 360.0);
	zns = fmod(zns, 360.0);
	ps = fmod(ps, 360.0);

	double rsta = sqrt(pow(RecCorXYZ[0], 2) + pow(RecCorXYZ[1], 2) + pow(RecCorXYZ[2], 2));
	double sinphi = RecCorXYZ[2] / rsta;
	double cosphi = sqrt(pow(RecCorXYZ[0], 2) + pow(RecCorXYZ[1], 2)) / rsta;

	double sinla = RecCorXYZ[1] / sinphi / rsta;
	double cosla = RecCorXYZ[0] / cosphi / rsta;
	double zla = atan2(RecCorXYZ[1], RecCorXYZ[0]);
	XCORSTA[0] = 0; XCORSTA[1] = 0; XCORSTA[2] = 0;
	for ( i = 0; i<31; i++)
	{
		double thetaf = (tau + datdi[0][i] * s + datdi[1][i] * h
			+ datdi[2][i] * p + datdi[3][i] * zns + datdi[4][i] * ps) *deg2rad;
		double dr =
			datdi[5][i] * 2.0*sinphi*cosphi*sin(thetaf + zla) + datdi[6][i] * 2.0*sinphi*cosphi*cos(thetaf + zla);
		double dn =
			datdi[7][i] * (pow(cosphi, 2) - pow(sinphi, 2))*sin(thetaf + zla) + datdi[8][i] * (pow(cosphi, 2) - pow(sinphi, 2))*cos(thetaf + zla);

		double de =
			datdi[7][i] * sinphi*cos(thetaf + zla) - datdi[8][i] * sinphi*sin(thetaf + zla);

		XCORSTA[0] = XCORSTA[0] + dr*cosla*cosphi - de*sinla - dn*sinphi*cosla;
		XCORSTA[1] = XCORSTA[1] + dr*sinla*cosphi + de*cosla - dn*sinphi*sinla;
		XCORSTA[2] = XCORSTA[2] + dr*sinphi + dn*cosphi;
	}

	for (i = 0; i<3; i++)
	{
		XCORSTA[i] = XCORSTA[i] / 1000.0;
	}
}

void iau_CAL2JD(int IY, int IM, int ID, double *DJM0, double *DJM)
{
	int J, MY, IYPMY;

	/*  Earliest year allowed (4800BC)*/
	int IYMIN = -4799;

	/*  Month lengths in days*/
	int MTAB[12] = { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };

	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

	/*  Preset status.*/
	J = 0;

	/*  Validate year.*/
	if (IY<IYMIN)
	{
		J = -1;
	}
	else
	{
		/*     Validate month.*/
		if (IM >= 1 && IM <= 12)
		{

			/*       Allow for leap year.*/
			if (IY % 4 == 0)
			{
				MTAB[1] = 29;
			}
			else
			{
				MTAB[1] = 28;
			}
			if ((IY % 100) == 0 && (IY % 400) != 0)
			{
				MTAB[1] = 28;
			}

			/*        Validate day.*/
			if (ID<1 || ID> MTAB[IM])
			{
				J = -3;
			}
			/*       Result.*/
			MY = (IM - 14) / 12;
			IYPMY = IY + MY;
			*DJM0 = 2400000.5;
			*DJM = ((1461.0 * (IYPMY + 4800)) / 4 + (367 * (IM - 2 - 12 * MY)) / 12
				- (3 * ((IYPMY + 4900) / 100)) / 4 + ID - 2432076);

		}
		else
		{
			/*        Bad month*/
			J = -2;
		}
	}
}

void iau_DAT(int IY, int IM, int ID, double FD, double* DELTAT)
{
	int J;

	/* Release year for this version of iau_DAT*/
	int IYV = 2009;

	/* Number of Delta(AT) changes (increase by 1 for each new leap second)*/
	const int NDAT = 39;

	/*  Number of Delta(AT) expressions before leap seconds were introduced*/
	const int NERA1 = 14;

	/*  Dates (year, month) on which new Delta(AT) came into force*/
	int
		IDATE[2][39] = { { 1960, 1961, 1961, 1962, 1963, 1964, 1964, 1964, 1965, 1965, 1965, 1965,
		1966, 1968, 1972, 1972, 1973, 1974, 1975, 1976, 1977, 1978, 1979, 1980, 1981, 1982,

		1983, 1985, 1988, 1990, 1991, 1992, 1993, 1994, 1996, 1997, 1999, 2006, 2009 }, { 1, 1, 8, 1, 11, 1, 4, 9, 1, 3, 7, 9, 1, 2, 1, 7, 1, 1, 1, 1, 1, 1, 1, 1, 7, 7, 7, 7, 1, 1, 1, 7, 7, 7, 1, 7, 1, 1, 1 } };

	/*  New Delta(AT) which came into force on the given dates*/
	double
		DATS[39] = { 1.4178180, 1.4228180, 1.3728180, 1.8458580, 1.9458580, 3.2401300, 3.3401300, 3.4401300,

		3.5401300, 3.6401300, 3.7401300, 3.8401300, 4.3131700, 4.2131700, 10, 11, 12, 13, 14, 15, 16, 17,
		18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34 };

	/*  Reference dates (MJD) and drift rates (s/day), pre leap seconds*/
	double DRIFT[2][14] = { { 37300, 37300, 37300, 37665, 37665, 38761,

		38761, 38761, 38761, 38761, 38761, 38761, 39126, 39126 }, { 0.001296, 0.001296, 0.001296, 0.0011232, 0.0011232,

		0.001296, 0.001296, 0.001296, 0.001296, 0.001296, 0.001296, 0.001296, 0.002592, 0.002592 } };

	/* Miscellaneous local variables*/
	int MORE;
	int JS, M, IS;   /*I, */
	double DAT, DJM0, DJM;

	/* Dates and Delta(AT)s*/

	/* Reference dates and drift rates*/


	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

	/* Initialize the result to zero and the status to OK.*/
	DAT = 0.0;
	JS = 0;

	/*  If invalid fraction of a day, set error status and give up.*/
	if (FD<0.0 || FD >= 1.0)
	{
		JS = -4;
		*DELTAT = DAT;
		J = JS;
		return;
	}

	/*  Convert the date into an MJD.*/
	iau_CAL2JD(IY, IM, ID, &DJM0, &DJM);

	/* If invalid year, month, or day, give up.*/
	if (JS < 0)
	{
		*DELTAT = DAT;
		J = JS;
		return;
	}

	/*  If pre-UTC year, set warning status and give up.*/
	if (IY < IDATE[0][0])
	{
		JS = 1;
		*DELTAT = DAT;
		J = JS;
		return;
	}

	/*  If suspiciously late year, set warning status but proceed.*/
	if (IY > IYV + 5)
	{
		JS = 1;
	}
	/*  Combine year and month.*/
	M = 12 * IY + IM;

	/*  Prepare to search the tables.*/
	MORE = 1;

	/*  Find the most recent table entry.*/
	int i = 0;
	for ( i = NDAT; i>1; i--)
	{
		if (MORE)
		{
			IS = i;
			MORE = M < (12 * IDATE[0][i] + IDATE[1][i]);
		}
	}

	/* Get the Delta(AT).*/
	DAT = DATS[IS - 1];

	/*  If pre-1972, adjust for drift.*/
	if (IS < NERA1)
	{
		DAT = DAT + (DJM + FD - DRIFT[0][IS]) * DRIFT[1][IS];
	}
	*DELTAT = DAT;
}

/* ------------------------------------------------------------------------
* compute the tidal corrections of station displacements
*  caused by lunar and solar gravitational attraction (see References).
*  The computations are calculated by the following steps:
*
*  Step 1): General degree 2 and degree 3 corrections + CALL ST1IDIU
*  + CALL ST1ISEM + CALL ST1L1.
*
*  Step 2): CALL STEP2DIU + CALL STEP2LON
* args   : double *xsta      I		Geocentric position of the IGS station (ITRF)
*		   double *XSUN      I      Geocentric position of the Sun (ITRF)
*		   double *XMON      I      Geocentric position of the Moon (ITRF)
*		   int    *YR        I      Year(UTC)
*          int    *MONTH     I      Month(UTC)
*          int    *DAY       I      Day of Month(UTC)
*          double *FHR       I      Hour in the day(UTC)
*
* return :double *DXTIDE     IO     Displacement vector(ITRF)
*All components are expressed in meters.
*---------------------------------------------------------------------------- - */
extern int dehanttideinel_(double *xsta, int *year, int *mon, int *day,
	double *fhr, double *xsun, double *xmon,
	double *DXTIDE)
{
	int i;
	double JJM0, JJM1;
	/* Nominal 2nd and 3rd degree load love and shida numbers*/
	double H20 = 0.6078;
	double L20 = 0.0847;
	double  H3 = 0.292;
	double  L3 = 0.015;
	double
		distSun = sqrt(pow(xsun[0], 2) + pow(xsun[1], 2) + pow(xsun[2], 2));
	double
		distMoon = sqrt(pow(xmon[0], 2) + pow(xmon[1], 2) + pow(xmon[2], 2));
	double
		distStn = sqrt(pow(xsta[0], 2) + pow(xsta[1], 2) + pow(xsta[2], 2));

	double scsun = dot(xsta, xsun, 3) / distStn / distSun;
	double scmon = dot(xsta, xmon, 3) / distStn / distMoon;

	/*/ COMPUTATION OF NEW H2 AND L2*/
	/*----------------------------*/
	double  cosphi = sqrt(pow(xsta[0], 2) + pow(xsta[1], 2)) / distStn;
	double  H2 = H20 - 0.0006*(1 - 3. / 2 * pow(cosphi, 2));
	double  L2 = L20 + 0.0002*(1 - 3. / 2 * pow(cosphi, 2));

	/* P2-TERM	 -------*/
	double  P2SUN = 3 * (H2 / 2 - L2)*pow(scsun, 2) - H2 / 2;
	double  P2MON = 3 * (H2 / 2 - L2)*pow(scmon, 2) - H2 / 2;

	/* P3-TERM	 -------*/
	double  P3SUN = 5.0 / 2.0*(H3 - 3 * L3)*pow(scsun, 3) + 3.0 / 2.0*(L3 - H3)*scsun;
	double  P3MON = 5.0 / 2.0*(H3 - 3 * L3)*pow(scmon, 3) + 3.0 / 2.0*(L3 - H3)*scmon;


	/* TERM IN DIRECTION OF SUN/MOON VECTOR
	 ------------------------------------*/
	double  X2SUN = 3 * L2*scsun;
	double  X2MON = 3 * L2*scmon;
	double  X3SUN = 3.0*L3 / 2.0*(5 * pow(scsun, 2) - 1);
	double  X3MON = 3.0*L3 / 2.0*(5 * pow(scmon, 2) - 1);


	/* FACTORS FOR SUN/MOON
	 --------------------*/
	double  MASS_RATIO_SUN = 332945.943062;
	double  MASS_RATIO_MOON = 0.012300034;
	double  RE = 6378136.55;
	double  FAC2SUN = MASS_RATIO_SUN*RE*pow(RE / distSun, 3);
	double  FAC2MON = MASS_RATIO_MOON*RE*pow(RE / distMoon, 3);
	double  FAC3SUN = FAC2SUN*(RE / distSun);
	double  FAC3MON = FAC2MON*(RE / distMoon);


	/* TOTAL DISPLACEMENT
	 ------------------*/
	/*DXTIDE= 0.0;*/
	for (i = 0; i < 3; i++){
		DXTIDE[i] = FAC2SUN*(X2SUN*xsun[i] / distSun + P2SUN*xsta[i] / distStn) +
			FAC2MON*(X2MON*xmon[i] / distMoon + P2MON*xsta[i] / distStn) +
			FAC3SUN*(X3SUN*xsun[i] / distSun + P3SUN*xsta[i] / distStn) +
			FAC3MON*(X3MON*xmon[i] / distMoon + P3MON*xsta[i] / distStn);
	}


	/* CORRECTIONS FOR THE OUT-OF-PHASE PART OF LOVE NUMBERS (PART H_2^(0)I
	            AND L_2^(0)I )
	 FIRST, FOR THE DIURNAL BAND*/
	double XCORSTA[3];
	ST1IDIU(xsta, xsun, xmon, FAC2SUN, FAC2MON, XCORSTA);
	for (i = 0; i<3; i++)
	{
		DXTIDE[i] = DXTIDE[i] + XCORSTA[i];
	}

	/* SECOND, FOR THE SEMI-DIURNAL BAND*/
	
	ST1ISEM(xsta, xsun, xmon, FAC2SUN, FAC2MON, XCORSTA);

	for (i = 0; i<3; i++)
	{
		DXTIDE[i] = DXTIDE[i] + XCORSTA[i];
	}

	/* CORRECTIONS FOR THE LATITUDE DEPENDENCE OF LOVE NUMBERS (PART L^(1) )*/
	
	ST1L1(xsta, xsun, xmon, FAC2SUN, FAC2MON, XCORSTA);
	for (i = 0; i<3; i++)
	{
		DXTIDE[i] = DXTIDE[i] + XCORSTA[i];
	}
	/*CONSIDER CORRECTIONS FOR STEP 2*/
	/*
	 CORRECTIONS FOR THE DIURNAL BAND:
	
	 FIRST, WE NEED TO KNOW THE DATE CONVERTED IN JULIAN CENTURIES*/
	iau_CAL2JD(*year, *mon, *day, &JJM0, &JJM1);
	double FHRD = *fhr;
	double T = ((JJM0 - 2451545.0) + JJM1 + FHRD) / 36525.0;
	double DTT;
	iau_DAT(*year, *mon, *day, FHRD, &DTT);
	DTT = DTT + 32.184;
	T = T + DTT / (3600.0*24.0*36525.0);
	/*
	  SECOND, WE CAN CALL THE SUBROUTINE STEP2DIU, FOR THE DIURNAL BAND CORRECTIONS,
	   (in-phase and out-of-phase frequency dependence):
	*/
	STEP2DIU(xsta, FHRD, T, XCORSTA);
	for (i = 0; i<3; i++)
	{
		DXTIDE[i] = DXTIDE[i] + XCORSTA[i];
	}
	/*  CORRECTIONS FOR THE LONG-PERIOD BAND,
	   (in-phase and out-of-phase frequency dependence):
	*/
	STEP2LON(xsta, FHRD, T, XCORSTA);
	for (i = 0; i<3; i++)
	{
		DXTIDE[i] = DXTIDE[i] + XCORSTA[i];
	}

}
#endif

/* solar/lunar tides (ref [2] 7) ---------------------------------------------*/
/*#ifndef IERS_MODEL*/
static void tide_pl(const double *eu, const double *rp, double GMp,
	const double *pos, double *dr)
{
	const double H3 = 0.292, L3 = 0.015;
	double r, ep[3], latp, lonp, p, K2, K3, a, H2, L2, dp, du, cosp, sinl, cosl;
	int i;

	trace(4, "tide_pl : pos=%.3f %.3f\n", pos[0] * R2D, pos[1] * R2D);

	if ((r = norm(rp, 3)) <= 0.0) return;

	for (i = 0; i<3; i++) ep[i] = rp[i] / r;

	K2 = GMp / GME*SQR(RE_WGS84)*SQR(RE_WGS84) / (r*r*r);
	K3 = K2*RE_WGS84 / r;
	latp = asin(ep[2]); lonp = atan2(ep[1], ep[0]);
	cosp = cos(latp); sinl = sin(pos[0]); cosl = cos(pos[0]);

	/* step1 in phase (degree 2) */
	p = (3.0*sinl*sinl - 1.0) / 2.0;
	H2 = 0.6078 - 0.0006*p;
	L2 = 0.0847 + 0.0002*p;
	a = dot(ep, eu, 3);
	dp = K2*3.0*L2*a;
	du = K2*(H2*(1.5*a*a - 0.5) - 3.0*L2*a*a);

	/* step1 in phase (degree 3) */
	dp += K3*L3*(7.5*a*a - 1.5);
	du += K3*(H3*(2.5*a*a*a - 1.5*a) - L3*(7.5*a*a - 1.5)*a);

	/* step1 out-of-phase (only radial) */
	du += 3.0 / 4.0*0.0025*K2*sin(2.0*latp)*sin(2.0*pos[0])*sin(pos[1] - lonp);
	du += 3.0 / 4.0*0.0022*K2*cosp*cosp*cosl*cosl*sin(2.0*(pos[1] - lonp));

	dr[0] = dp*ep[0] + du*eu[0];
	dr[1] = dp*ep[1] + du*eu[1];
	dr[2] = dp*ep[2] + du*eu[2];

	trace(5, "tide_pl : dr=%.3f %.3f %.3f\n", dr[0], dr[1], dr[2]);
}
/* displacement by solid earth tide (ref [2] 7) ------------------------------*/
static void tide_solid(const double *rsun, const double *rmoon,
	const double *pos, const double *E, double gmst, int opt,
	double *dr)
{
	double dr1[3], dr2[3], eu[3], du, dn, sinl, sin2l;

	trace(3, "tide_solid: pos=%.3f %.3f opt=%d\n", pos[0] * R2D, pos[1] * R2D, opt);

	/* step1: time domain */
	eu[0] = E[2]; eu[1] = E[5]; eu[2] = E[8];
	tide_pl(eu, rsun, GMS, pos, dr1);
	tide_pl(eu, rmoon, GMM, pos, dr2);

	/* step2: frequency domain, only K1 radial */
	sin2l = sin(2.0*pos[0]);
	du = -0.012*sin2l*sin(gmst + pos[1]);

	dr[0] = dr1[0] + dr2[0] + du*E[2];
	dr[1] = dr1[1] + dr2[1] + du*E[5];
	dr[2] = dr1[2] + dr2[2] + du*E[8];

	/* eliminate permanent deformation */
	if (opt & 8) {
		sinl = sin(pos[0]);
		du = 0.1196*(1.5*sinl*sinl - 0.5);
		dn = 0.0247*sin2l;
		dr[0] += du*E[2] + dn*E[1];
		dr[1] += du*E[5] + dn*E[4];
		dr[2] += du*E[8] + dn*E[7];
	}
	trace(5, "tide_solid: dr=%.3f %.3f %.3f\n", dr[0], dr[1], dr[2]);
}
/*#endif*/ /* !IERS_MODEL */

/* displacement by ocean tide loading (ref [2] 7) ----------------------------*/
static void tide_oload(gtime_t tut, const double *odisp, double *denu)
{
	const double args[][5] = {
		{ 1.40519E-4, 2.0, -2.0, 0.0, 0.00 },  /* M2 */
		{ 1.45444E-4, 0.0, 0.0, 0.0, 0.00 },  /* S2 */
		{ 1.37880E-4, 2.0, -3.0, 1.0, 0.00 },  /* N2 */
		{ 1.45842E-4, 2.0, 0.0, 0.0, 0.00 },  /* K2 */
		{ 0.72921E-4, 1.0, 0.0, 0.0, 0.25 },  /* K1 */
		{ 0.67598E-4, 1.0, -2.0, 0.0, -0.25 },  /* O1 */
		{ 0.72523E-4, -1.0, 0.0, 0.0, -0.25 },  /* P1 */
		{ 0.64959E-4, 1.0, -3.0, 1.0, -0.25 },  /* Q1 */
		{ 0.53234E-5, 0.0, 2.0, 0.0, 0.00 },  /* Mf */
		{ 0.26392E-5, 0.0, 1.0, -1.0, 0.00 },  /* Mm */
		{ 0.03982E-5, 2.0, 0.0, 0.0, 0.00 }   /* Ssa */
	};
	const double ep1975[] = { 1975, 1, 1, 0, 0, 0 };
	double ep[6], fday, days, t, t2, t3, a[5], ang, dp[3] = { 0 };
	int i, j;

	trace(3, "tide_oload:\n");

	/* angular argument: see subroutine arg.f for reference [1] */
	time2epoch(tut, ep);
	fday = ep[3] * 3600.0 + ep[4] * 60.0 + ep[5];
	ep[3] = ep[4] = ep[5] = 0.0;
	days = timediff(epoch2time(ep), epoch2time(ep1975)) / 86400.0 + 1.0;
	t = (27392.500528 + 1.000000035*days) / 36525.0;
	t2 = t*t; t3 = t2*t;

	a[0] = fday;
	a[1] = (279.69668 + 36000.768930485*t + 3.03E-4*t2)*D2R; /* H0 */
	a[2] = (270.434358 + 481267.88314137*t - 0.001133*t2 + 1.9E-6*t3)*D2R; /* S0 */
	a[3] = (334.329653 + 4069.0340329577*t - 0.010325*t2 - 1.2E-5*t3)*D2R; /* P0 */
	a[4] = 2.0*PI;

	/* displacements by 11 constituents */
	for (i = 0; i<11; i++) {
		ang = 0.0;
		for (j = 0; j<5; j++) ang += a[j] * args[i][j];
		for (j = 0; j<3; j++) dp[j] += odisp[j + i * 6] * cos(ang - odisp[j + 3 + i * 6] * D2R);
	}
	denu[0] = -dp[1];
	denu[1] = -dp[2];
	denu[2] = dp[0];

	trace(5, "tide_oload: denu=%.3f %.3f %.3f\n", denu[0], denu[1], denu[2]);
}
/* iers mean pole (ref [7] eq.7.25) ------------------------------------------*/
static void iers_mean_pole(gtime_t tut, double *xp_bar, double *yp_bar)
{
	const double ep2000[] = { 2000, 1, 1, 0, 0, 0 };
	double y, y2, y3;

	y = timediff(tut, epoch2time(ep2000)) / 86400.0 / 365.25;

	if (y<3653.0 / 365.25) { /* until 2010.0 */
		y2 = y*y; y3 = y2*y;
		*xp_bar = 55.974 + 1.8243*y + 0.18413*y2 + 0.007024*y3; /* (mas) */
		*yp_bar = 346.346 + 1.7896*y - 0.10729*y2 - 0.000908*y3;
	}
	else { /* after 2010.0 */
		*xp_bar = 23.513 + 7.6141*y; /* (mas) */
		*yp_bar = 358.891 - 0.6287*y;
	}
}
/* displacement by pole tide (ref [7] eq.7.26) --------------------------------*/
static void tide_pole(gtime_t tut, const double *pos, const double *erpv,
	double *denu)
{
	double xp_bar, yp_bar, m1, m2, cosl, sinl;

	trace(3, "tide_pole: pos=%.3f %.3f\n", pos[0] * R2D, pos[1] * R2D);

	/* iers mean pole (mas) */
	iers_mean_pole(tut, &xp_bar, &yp_bar);

	/* ref [7] eq.7.24 */
	m1 = erpv[0] / AS2R - xp_bar*1E-3; /* (as) */
	m2 = -erpv[1] / AS2R + yp_bar*1E-3;

	/* sin(2*theta) = sin(2*phi), cos(2*theta)=-cos(2*phi) */
	cosl = cos(pos[1]);
	sinl = sin(pos[1]);
	denu[0] = 9E-3*sin(pos[0])    *(m1*sinl - m2*cosl); /* de= Slambda (m) */
	denu[1] = -9E-3*cos(2.0*pos[0])*(m1*cosl + m2*sinl); /* dn=-Stheta  (m) */
	denu[2] = -33E-3*sin(2.0*pos[0])*(m1*cosl + m2*sinl); /* du= Sr      (m) */

	trace(5, "tide_pole : denu=%.3f %.3f %.3f\n", denu[0], denu[1], denu[2]);
}
/* tidal displacement ----------------------------------------------------------
* displacements by earth tides
* args   : gtime_t tutc     I   time in utc
*          double *rr       I   site position (ecef) (m)
*          int    opt       I   options (or of the followings)
*                                 1: solid earth tide
*                                 2: ocean tide loading
*                                 4: pole tide
*                                 8: elimate permanent deformation
*          double *erp      I   earth rotation parameters (NULL: not used)
*          double *odisp    I   ocean loading parameters  (NULL: not used)
*                                 odisp[0+i*6]: consituent i amplitude radial(m)
*                                 odisp[1+i*6]: consituent i amplitude west  (m)
*                                 odisp[2+i*6]: consituent i amplitude south (m)
*                                 odisp[3+i*6]: consituent i phase radial  (deg)
*                                 odisp[4+i*6]: consituent i phase west    (deg)
*                                 odisp[5+i*6]: consituent i phase south   (deg)
*                                (i=0:M2,1:S2,2:N2,3:K2,4:K1,5:O1,6:P1,7:Q1,
*                                   8:Mf,9:Mm,10:Ssa)
*          double *dr       O   displacement by earth tides (ecef) (m)
* return : none
* notes  : see ref [1], [2] chap 7
*          see ref [4] 5.2.1, 5.2.2, 5.2.3
*          ver.2.4.0 does not use ocean loading and pole tide corrections
*-----------------------------------------------------------------------------*/
extern void tidedisp(gtime_t tutc, const double *rr, int opt, const erp_t *erp,
	const double *odisp, double *dr)
{
	gtime_t tut;
	double pos[2], E[9], drt[3] = { 0 }, denu[3], rs[3], rm[3], gmst, erpv[5] = { 0 };
	int i;
#ifdef IERS_MODEL
	double ep[6], fhr;
	int year, mon, day;
#endif

	trace(3, "tidedisp: tutc=%s\n", time_str(tutc, 0));

	if (erp) geterp(erp, tutc, erpv);

	tut = timeadd(tutc, erpv[2]);

	dr[0] = dr[1] = dr[2] = 0.0;

	if (norm(rr, 3) <= 0.0) return;

	pos[0] = asin(rr[2] / norm(rr, 3));
	pos[1] = atan2(rr[1], rr[0]);
	xyz2enu(pos, E);

	if (opt & 1) { /* solid earth tides */

		/* sun and moon position in ecef */
		sunmoonpos(tutc, erpv, rs, rm, &gmst);

#ifdef IERS_MODEL
		time2epoch(tutc, ep);
		year = (int)ep[0];
		mon = (int)ep[1];
		day = (int)ep[2];
		fhr = ep[3] + ep[4] / 60.0 + ep[5] / 3600.0;

		/* --------------call DEHANTTIDEINEL------------------ */
		/**  Test case:
		*     given input: XSTA(1) = 4075578.385D0 meters
		*                  XSTA(2) =  931852.890D0 meters
		*                  XSTA(3) = 4801570.154D0 meters
		*                  XSUN(1) = 137859926952.015D0 meters
		*                  XSUN(2) = 54228127881.4350D0 meters
		*                  XSUN(3) = 23509422341.6960D0 meters
		*                  XMON(1) = -179996231.920342D0 meters
		*                  XMON(2) = -312468450.131567D0 meters
		*                  XMON(3) = -169288918.592160D0 meters
		*                  YR      = 2009
		*                  MONTH   = 4
		*                  DAY     = 13
		*                  FHR     = 0.00D0 seconds
		*
		*     expected output:  DXTIDE(1) = 0.7700420357108125891D-01 meters
		*                       DXTIDE(2) = 0.6304056321824967613D-01 meters
		*                       DXTIDE(3) = 0.5516568152597246810D-01 meters
		**/
		/*rr[0] = 4075578.385; rr[1] = 931852.890; rr[2] = 4801570.154;
		rs[0] = 137859926952.015; rs[1] = 54228127881.4350; rs[2] = 23509422341.6960;
		rm[0] = -179996231.920342; rm[1] = -312468450.131567; rm[2] = -169288918.592160;

		year = 2009; mon = 4; day = 13; fhr = 0.0;*/

		dehanttideinel_((double *)rr, &year, &mon, &day, &fhr, rs, rm, drt);

#else
		tide_solid(rs, rm, pos, E, gmst, opt, drt);
#endif
		for (i = 0; i<3; i++) dr[i] += drt[i];
	}
	if ((opt & 2) && odisp) { /* ocean tide loading */
		tide_oload(tut, odisp, denu);
		matmul("TN", 3, 1, 3, 1.0, E, denu, 0.0, drt);
		for (i = 0; i<3; i++) dr[i] += drt[i];
	}
	if ((opt & 4) && erp) { /* pole tide */
		tide_pole(tut, pos, erpv, denu);
		matmul("TN", 3, 1, 3, 1.0, E, denu, 0.0, drt);
		for (i = 0; i<3; i++) dr[i] += drt[i];
	}
	trace(5, "tidedisp: dr=%.3f %.3f %.3f\n", dr[0], dr[1], dr[2]);
}
