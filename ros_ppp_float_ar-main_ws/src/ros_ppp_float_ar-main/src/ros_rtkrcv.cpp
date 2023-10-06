/*------------------------------------------------------------------------------
 * rtkrcv.c : rtk-gps/gnss receiver console ap
 *
 *          Copyright (C) 2009-2015 by T.TAKASU, All rights reserved.
 *
 * notes   :
 *     current version does not support win32 without pthread library
 *
 * version : $Revision:$ $Date:$
 * history : 2009/12/13 1.0  new
 *           2010/07/18 1.1  add option -m
 *           2010/08/12 1.2  fix bug on ftp/http
 *           2011/01/22 1.3  add option misc-proxyaddr,misc-fswapmargin
 *           2011/08/19 1.4  fix bug on size of arg solopt arg for rtksvrstart()
 *           2012/11/03 1.5  fix bug on setting output format
 *           2013/06/30 1.6  add "nvs" option for inpstr*-format
 *           2014/02/10 1.7  fix bug on printing obs data
 *                           add print of status, glonass nav data
 *                           ignore SIGHUP
 *           2014/04/27 1.8  add "binex" option for inpstr*-format
 *           2014/08/10 1.9  fix cpu overload with abnormal telnet shutdown
 *           2014/08/26 1.10 support input format "rt17"
 *                           change file paths of solution status and debug trace
 *           2015/01/10 1.11 add line editting and command history
 *                           separate codes for virtual console to vt.c
 *           2015/05/22 1.12 fix bug on sp3 id in inpstr*-format options
 *           2015/07/31 1.13 accept 4:stat for outstr1-format or outstr2-format
 *                           add reading satellite dcb
 *           2015/12/14 1.14 add option -sta for station name (#339)
 *           2015/12/25 1.15 fix bug on -sta option (#339)
 *           2015/01/26 1.16 support septentrio
 *           2016/07/01 1.17 support CMR/CMR+
 *           2016/08/20 1.18 add output of patch level with version
 *           2016/09/05 1.19 support ntrip caster for output stream
 *           2016/09/19 1.20 support multiple remote console connections
 *                           add option -w
 *           2017/09/01 1.21 add command ssr
 *-----------------------------------------------------------------------------*/
#include <stdlib.h>
#include <signal.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <netinet/tcp.h>
#include <arpa/inet.h>
#include <netdb.h>
#include <errno.h>
#include "rtklib.h"
#include "vt.h"
#include "ros/ros.h"
#include <gnss_comm/gnss_constant.hpp>
#include "gnss_comm/gnss_ros.hpp"
#include "gnss_comm/gnss_utility.hpp"
#include "rtkrcv/corr_msg.h"

#define PRGNAME "rtkrcv"                    /* program name */
#define CMDPROMPT "rtkrcv> "                /* command prompt */
#define MAXCON 32                           /* max number of consoles */
#define MAXARG 10                           /* max number of args in a command */
#define MAXCMD 256                          /* max length of a command */
#define MAXSTR 1024                         /* max length of a stream */
#define OPTSDIR "."                         /* default config directory */
#define OPTSFILE "rtkrcv.conf"              /* default config file */
#define NAVIFILE "rtkrcv.nav"               /* navigation save file */
#define STATFILE "rtkrcv_%Y%m%d%h%M.stat"   /* solution status file */
#define TRACEFILE "rtkrcv_%Y%m%d%h%M.trace" /* debug trace file */
#define INTKEEPALIVE 1000                   /* keep alive interval (ms) */

#define ESC_CLEAR "\033[H\033[2J" /* ansi/vt100 escape: erase screen */
#define ESC_RESET "\033[0m"       /* ansi/vt100: reset attribute */
#define ESC_BOLD "\033[1m"        /* ansi/vt100: bold */

#define SQRT(x) ((x) <= 0.0 || (x) != (x) ? 0.0 : sqrt(x))

/* type defintions -----------------------------------------------------------*/

typedef struct
{                    /* console type */
   int state;        /* state (0:stop,1:run) */
   vt_t *vt;         /* virtual terminal */
   pthread_t thread; /* console thread */
} con_t;

/* function prototypes -------------------------------------------------------*/
extern FILE *popen(const char *, const char *);
extern int pclose(FILE *);

/* global variables ----------------------------------------------------------*/
static rtksvr_t svr;  /* rtk server struct */
static stream_t moni; /* monitor stream */

static int intflg = 0; /* interrupt flag (2:shtdown) */

static char passwd[MAXSTR] = "admin"; /* login password */
static int timetype = 0;              /* time format (0:gpst,1:utc,2:jst,3:tow) */
static int soltype = 0;               /* sol format (0:dms,1:deg,2:xyz,3:enu,4:pyl) */
static int solflag = 2;               /* sol flag (1:std+2:age/ratio/ns) */
static int strtype[] = {              /* stream types */
                        STR_SERIAL, STR_NONE, STR_NONE, STR_NONE, STR_NONE, STR_NONE, STR_NONE, STR_NONE};
static char strpath[8][MAXSTR] = {"", "", "", "", "", "", "", ""}; /* stream paths */
static int strfmt[] = {                                            /* stream formats */
                       STRFMT_UBX, STRFMT_RTCM3, STRFMT_SP3, SOLF_LLH, SOLF_NMEA};
static int svrcycle = 10;              /* server cycle (ms) */
static int timeout = 10000;            /* timeout time (ms) */
static int reconnect = 10000;          /* reconnect interval (ms) */
static int nmeacycle = 5000;           /* nmea request cycle (ms) */
static int buffsize = 32768;           /* input buffer size (bytes) */
static int navmsgsel = 0;              /* navigation mesaage select */
static char proxyaddr[256] = "";       /* http/ntrip proxy */
static int nmeareq = 0;                /* nmea request type (0:off,1:lat/lon,2:single) */
static double nmeapos[] = {0, 0, 0};   /* nmea position (lat/lon/height) (deg,m) */
static char rcvcmds[3][MAXSTR] = {""}; /* receiver commands files */
static char startcmd[MAXSTR] = "";     /* start command */
static char stopcmd[MAXSTR] = "";      /* stop command */
static int modflgr[256] = {0};         /* modified flags of receiver options */
static int modflgs[256] = {0};         /* modified flags of system options */
static int moniport = 0;               /* monitor port */
static int keepalive = 0;              /* keep alive flag */
static int start = 0;                  /* auto start */
static int fswapmargin = 30;           /* file swap margin (s) */
static char sta_name[256] = "";        /* station name */

static prcopt_t prcopt;            /* processing options */
static solopt_t solopt[2] = {{0}}; /* solution options */
static filopt_t filopt = {""};     /* file options */

/* receiver options table ----------------------------------------------------*/
#define TIMOPT "0:gpst,1:utc,2:jst,3:tow"
#define CONOPT "0:dms,1:deg,2:xyz,3:enu,4:pyl"
#define FLGOPT "0:off,1:std+2:age/ratio/ns"
#define ISTOPT "0:off,1:serial,2:file,3:tcpsvr,4:tcpcli,7:ntripcli,8:ftp,9:http"
#define OSTOPT "0:off,1:serial,2:file,3:tcpsvr,4:tcpcli,6:ntripsvr,11:ntripc_c"
#define FMTOPT "0:rtcm2,1:rtcm3,2:oem4,3:oem3,4:ubx,5:swiftnav,6:hemis,7:skytraq,8:gw10,9:javad,10:nvs,11:binex,12:rt17,13:sbf,14:cmr,15:tersus,18:sp3"
#define NMEOPT "0:off,1:latlon,2:single"
#define SOLOPT "0:llh,1:xyz,2:enu,3:nmea,4:stat"
#define MSGOPT "0:all,1:rover,2:base,3:corr"

static opt_t rcvopts[] = {
    {"console-passwd", 2, (void *)passwd, ""},
    {"console-timetype", 3, (void *)&timetype, TIMOPT},
    {"console-soltype", 3, (void *)&soltype, CONOPT},
    {"console-solflag", 0, (void *)&solflag, FLGOPT},

    {"inpstr1-type", 3, (void *)&strtype[0], ISTOPT},
    {"inpstr2-type", 3, (void *)&strtype[1], ISTOPT},
    {"inpstr3-type", 3, (void *)&strtype[2], ISTOPT},
    {"inpstr1-path", 2, (void *)strpath[0], ""},
    {"inpstr2-path", 2, (void *)strpath[1], ""},
    {"inpstr3-path", 2, (void *)strpath[2], ""},
    {"inpstr1-format", 3, (void *)&strfmt[0], FMTOPT},
    {"inpstr2-format", 3, (void *)&strfmt[1], FMTOPT},
    {"inpstr3-format", 3, (void *)&strfmt[2], FMTOPT},
    {"inpstr2-nmeareq", 3, (void *)&nmeareq, NMEOPT},
    {"inpstr2-nmealat", 1, (void *)&nmeapos[0], "deg"},
    {"inpstr2-nmealon", 1, (void *)&nmeapos[1], "deg"},
    {"inpstr2-nmeahgt", 1, (void *)&nmeapos[2], "m"},
    {"outstr1-type", 3, (void *)&strtype[3], OSTOPT},
    {"outstr2-type", 3, (void *)&strtype[4], OSTOPT},
    {"outstr1-path", 2, (void *)strpath[3], ""},
    {"outstr2-path", 2, (void *)strpath[4], ""},
    {"outstr1-format", 3, (void *)&strfmt[3], SOLOPT},
    {"outstr2-format", 3, (void *)&strfmt[4], SOLOPT},
    {"logstr1-type", 3, (void *)&strtype[5], OSTOPT},
    {"logstr2-type", 3, (void *)&strtype[6], OSTOPT},
    {"logstr3-type", 3, (void *)&strtype[7], OSTOPT},
    {"logstr1-path", 2, (void *)strpath[5], ""},
    {"logstr2-path", 2, (void *)strpath[6], ""},
    {"logstr3-path", 2, (void *)strpath[7], ""},

    {"misc-svrcycle", 0, (void *)&svrcycle, "ms"},
    {"misc-timeout", 0, (void *)&timeout, "ms"},
    {"misc-reconnect", 0, (void *)&reconnect, "ms"},
    {"misc-nmeacycle", 0, (void *)&nmeacycle, "ms"},
    {"misc-buffsize", 0, (void *)&buffsize, "bytes"},
    {"misc-navmsgsel", 3, (void *)&navmsgsel, MSGOPT},
    {"misc-proxyaddr", 2, (void *)proxyaddr, ""},
    {"misc-fswapmargin", 0, (void *)&fswapmargin, "s"},

    {"misc-startcmd", 2, (void *)startcmd, ""},
    {"misc-stopcmd", 2, (void *)stopcmd, ""},

    {"file-cmdfile1", 2, (void *)rcvcmds[0], ""},
    {"file-cmdfile2", 2, (void *)rcvcmds[1], ""},
    {"file-cmdfile3", 2, (void *)rcvcmds[2], ""},

    {"", 0, NULL, ""}};

constexpr char opt_file[] = "rtkrcv.conf";

#define TIMOPT "0:gpst,1:utc,2:jst,3:tow"
#define CONOPT "0:dms,1:deg,2:xyz,3:enu,4:pyl"
#define FLGOPT "0:off,1:std+2:age/ratio/ns"
#define ISTOPT "0:off,1:serial,2:file,3:tcpsvr,4:tcpcli,6:ntripcli,7:ftp,8:http"
#define OSTOPT "0:off,1:serial,2:file,3:tcpsvr,4:tcpcli,6:ntripsvr,11:ntripc_c"
#define FMTOPT "0:rtcm2,1:rtcm3,2:oem4,3:oem3,4:ubx,5:ss2,6:hemis,7:skytraq,8:gw10,9:javad,10:nvs,11:binex,12:rt17,13:sbf,14:cmr,15:tersus,18:sp3"
#define NMEOPT "0:off,1:latlon,2:single"
#define SOLOPT "0:llh,1:xyz,2:enu,3:nmea,4:stat"
#define MSGOPT "0:all,1:rover,2:base,3:corr"

static ros::Publisher pos_pub;
static ros::Publisher meas_pub;
static ros::Publisher ephem_pub;
static ros::Publisher glo_ephem_pub;
static ros::Publisher iono_pub;
static ros::Publisher ssr_pub;

extern "C" void publish_nmea(rtksvr_t *svr)
{

   if (svr->raw->cstrParsedGPGGA[0])
   {
      printf("CALL ON: %s", svr->raw->cstrParsedGPGGA);
      // xxx
      svr->raw->cstrParsedGPGGA[0] = 0x00;
   }

   if (svr->raw->cstrParsedGPRMC[0])
   {
      printf("CALL ON: %s", svr->raw->cstrParsedGPRMC);
      // xxx
      svr->raw->cstrParsedGPRMC[0] = 0x00;
   }
}

extern "C" void publish_position(const sol_t *sol)
{
   gnss_comm::GnssPVTSolnMsg pvt_msg;

   int gps_week;
   pvt_msg.time.tow = time2gpst(sol->time, &gps_week);
   pvt_msg.time.week = gps_week;

   switch (sol->stat)
   {
   case SOLQ_SINGLE:
      pvt_msg.fix_type = 3;
      break;
   case SOLQ_PPP:
      pvt_msg.fix_type = 6;
      break;
   case SOLQ_NONE:
   default:
      pvt_msg.fix_type = 0;
      break;
   }

   pvt_msg.valid_fix = pvt_msg.fix_type != 0 ? 1 : 0;
   pvt_msg.diff_soln = 0;
   pvt_msg.carr_soln = 0;
   pvt_msg.num_sv = sol->ns;

   double pos[3];
   ecef2pos(sol->rr, pos);

   pvt_msg.latitude = pos[0] * R2D;
   pvt_msg.longitude = pos[1] * R2D;
   pvt_msg.altitude = pos[2];
   pvt_msg.height_msl = pos[2] - geoidh(pos);

   pvt_msg.h_acc = 0;
   pvt_msg.v_acc = 0;
   pvt_msg.p_dop = 0;

   double vel[3];
   ecef2enu(pos, sol->rr + 3, vel);

   pvt_msg.vel_n = vel[1];
   pvt_msg.vel_e = vel[0];
   pvt_msg.vel_d = -vel[2];
   pvt_msg.vel_acc = 0;

   pos_pub.publish(pvt_msg);
}

static uint32_t satno_cvt(int satno)
{
   int prn;
   int sys = satsys(satno, &prn);

   return gnss_comm::sat_no(sys, prn);
}

extern "C" void publish_measurement(obsd_t *obs, int n)
{
   // std::cout<<"test_range_measure"<<std::endl;
   gnss_comm::GnssMeasMsg gnss_meas_msg;

   // HACK

   for (int i = 0; i < n; i++)
   {
      gnss_comm::GnssObsMsg obs_msg;

      int week = 0;
      double tow = time2gpst(obs[i].time, &week);
      obs_msg.time.week = week;
      obs_msg.time.tow = tow;
      obs_msg.sat = satno_cvt(obs[i].sat);

      int prn;
      int sys = satsys(obs[i].sat, &prn);

      for (int n = 0; n < NFREQ; n++)
      {
         obs_msg.freqs.push_back(code2freq(sys, obs[i].code[n], 0));
         obs_msg.CN0.push_back(obs[i].SNR[n] * 0.001);
         obs_msg.LLI.push_back(obs[i].LLI[n] & LLI_SLIP);
         obs_msg.code.push_back(obs[i].code[n]);
         obs_msg.psr.push_back(obs[i].P[n]);
         obs_msg.psr_std.push_back(0);
         obs_msg.cp.push_back(obs[i].L[n]);
         obs_msg.cp_std.push_back(0);
         obs_msg.dopp.push_back(obs[i].D[n]);
         obs_msg.dopp_std.push_back(0);

         uint8_t status = (obs[i].LLI[n] & LLI_HALFC ? 4 : 0) | (obs[i].LLI[n] & LLI_HALFS ? 8 : 0) | 3;
         obs_msg.status.push_back(status);
      }

      gnss_meas_msg.meas.push_back(obs_msg);
   }

   meas_pub.publish(gnss_meas_msg);
}

extern "C" void publish_ephemeris(eph_t *eph)
{
   std::cout << "ephem\n"
             << std::endl;
   gnss_comm::GnssEphemMsg ephem_msg;

   int week, sys;
   double tow;

   ephem_msg.sat = satno_cvt(eph->sat);
   sys = satsys(eph->sat, NULL);

   if (sys == SYS_GLO)
      return;

   tow = time2gpst(eph->ttr, &week);
   ephem_msg.ttr.week = week;
   ephem_msg.ttr.tow = tow;

   tow = time2gpst(eph->toe, &week);
   ephem_msg.toe.week = week;
   ephem_msg.toe.tow = tow;

   tow = time2gpst(eph->toc, &week);
   ephem_msg.toc.week = week;
   ephem_msg.toc.tow = tow;

   ephem_msg.toe_tow = eph->toes;
   ephem_msg.week = eph->week;
   ephem_msg.iode = eph->iode;
   ephem_msg.iodc = eph->iodc;
   ephem_msg.health = eph->svh;
   ephem_msg.code = eph->code;
   ephem_msg.ura = uravalue(sys, eph->sva);
   ephem_msg.A = eph->A;
   ephem_msg.e = eph->e;
   ephem_msg.i0 = eph->i0;
   ephem_msg.omg = eph->omg;
   ephem_msg.OMG0 = eph->OMG0;
   ephem_msg.M0 = eph->M0;
   ephem_msg.delta_n = eph->deln;
   ephem_msg.OMG_dot = eph->OMGd;
   ephem_msg.i_dot = eph->idot;
   ephem_msg.cuc = eph->cuc;
   ephem_msg.cus = eph->cus;
   ephem_msg.crc = eph->crc;
   ephem_msg.crs = eph->crs;
   ephem_msg.cic = eph->cic;
   ephem_msg.cis = eph->cis;
   ephem_msg.af0 = eph->f0;
   ephem_msg.af1 = eph->f1;
   ephem_msg.af2 = eph->f2;
   ephem_msg.tgd0 = eph->tgd[0];
   ephem_msg.tgd1 = eph->tgd[1];
   ephem_msg.A_dot = eph->Adot;
   ephem_msg.n_dot = eph->ndot;

   ephem_pub.publish(ephem_msg);
}

extern "C" void publish_glo_ephemeris(geph_t *geph)
{
   std::cout << "glo_ephem\n"
             << std::endl;
   gnss_comm::GnssGloEphemMsg glo_ephem_msg;

   int week;
   double tow;

   glo_ephem_msg.sat = satno_cvt(geph->sat);
   tow = time2gpst(geph->tof, &week);
   glo_ephem_msg.ttr.week = week;
   glo_ephem_msg.ttr.tow = tow;
   tow = time2gpst(geph->toe, &week);
   glo_ephem_msg.toe.week = week;
   glo_ephem_msg.toe.tow = tow;
   glo_ephem_msg.freqo = geph->frq;
   glo_ephem_msg.iode = geph->iode;
   glo_ephem_msg.health = geph->svh;
   glo_ephem_msg.age = geph->age;
   /* mapping GLONASS sva to ura*/
   if (geph->sva == 2)        glo_ephem_msg.ura = 2.5;   
   else if (geph->sva < 5)    glo_ephem_msg.ura = geph->sva+1.0;
   else if (geph->sva == 5)   glo_ephem_msg.ura = 7.0;
   else if (geph->sva < 10)   glo_ephem_msg.ura = 10.0+(geph->sva-6)*2.0;
   else if (geph->sva < 15)   glo_ephem_msg.ura = 1 << (geph->sva-5);
   else                glo_ephem_msg.ura = -1.0;

   glo_ephem_msg.pos_x = geph->pos[0];
   glo_ephem_msg.pos_y = geph->pos[1];
   glo_ephem_msg.pos_z = geph->pos[2];
   glo_ephem_msg.vel_x = geph->vel[0];
   glo_ephem_msg.vel_y = geph->vel[1];
   glo_ephem_msg.vel_z = geph->vel[2];
   glo_ephem_msg.acc_x = geph->acc[0];
   glo_ephem_msg.acc_y = geph->acc[1];
   glo_ephem_msg.acc_z = geph->acc[2];
   glo_ephem_msg.tau_n = geph->taun;
   glo_ephem_msg.gamma = geph->gamn;
   glo_ephem_msg.delta_tau_n = geph->dtaun;

   glo_ephem_pub.publish(glo_ephem_msg);
}

extern "C" void publish_iono(rtksvr_t *svr)
{
   gnss_comm::StampedFloat64Array iono_msg;

   int sys;

   if (svr->nav.eph->ttr.time != 0)
      iono_msg.header.stamp = ros::Time(static_cast<double>(svr->nav.eph->ttr.time) + svr->nav.eph->ttr.sec);

   sys = satsys(svr->nav.eph->sat, NULL);
   switch (sys)
   {
   case SYS_GPS:
      memcpy(&iono_msg.data, svr->nav.ion_gps, sizeof(svr->nav.ion_gps));
      break;
   case SYS_GLO:
      break;
   case SYS_GAL:
      memcpy(&iono_msg.data, svr->nav.ion_gal, sizeof(svr->nav.ion_gal));
      break;
   case SYS_BDS:
      memcpy(&iono_msg.data, svr->nav.ion_bds, sizeof(svr->nav.ion_bds));
      break;
   }

   iono_pub.publish(iono_msg);
}

extern "C" void publish_corrections(ssr_t *ssr_data)
{
   rtkrcv::corr_msg corrs_ros_msg;

   for (int i = 0; i < MAXSAT; i++)
   {
      ssr_t *ssr = ssr_data + i;
      if (!ssr->t0[0].time)
         continue;

      rtkrcv::ssr_msg ssr_ros_msg;

      char satid[32];
      satno2id(i + 1, satid);

      ssr_ros_msg.sat = std::string(satid);

      for (int n = 0; n < 6; n++)
      {
         int week;
         double tow;
         tow = time2gpst(ssr->t0[n], &week);

         ssr_ros_msg.t0[n].week = week;
         ssr_ros_msg.t0[n].tow = tow;
      }

      for (int n = 0; n < 6; n++)
      {
         ssr_ros_msg.udi[n] = ssr->udi[n];
      }

      for (int n = 0; n < 6; n++)
      {
         ssr_ros_msg.iod_ssr[n] = ssr->iod[n];
      }

      ssr_ros_msg.iode = ssr->iode;
      ssr_ros_msg.ura = ssr->ura;

      for (int n = 0; n < 3; n++)
      {
         ssr_ros_msg.orbit_corr[n] = ssr->deph[n];
      }

      for (int n = 0; n < 3; n++)
      {
         ssr_ros_msg.orbit_corr[n + 3] = ssr->ddeph[n];
      }

      for (int n = 0; n < 3; n++)
      {
         ssr_ros_msg.clock_corr[n] = ssr->dclk[n];
      }

      ssr_ros_msg.hrclk = ssr->hrclk;

      for (int n = 0; n < MAXCODE; n++)
      {
         ssr_ros_msg.code_bias[n] = ssr->cbias[n];
      }

      for (int n = 0; n < MAXCODE; n++)
      {
         ssr_ros_msg.phase_bias[n] = ssr->pbias[n];
      }

      corrs_ros_msg.corr_data.push_back(ssr_ros_msg);
   }

   ssr_pub.publish(corrs_ros_msg);
}

static int start_rtklib()
{
   double npos[3] = {0, 0, 0};
   char *ropts[] = {"-EPHALL", "", ""};
   char *paths[] = {
       strpath[0], strpath[1], strpath[2], strpath[3], strpath[4], strpath[5],
       strpath[6], strpath[7]};
   char s1[3][MAXRCVCMD] = {"", "", ""}, *cmds[] = {NULL, NULL, NULL};
   char s2[3][MAXRCVCMD] = {"", "", ""}, *cmds_periodic[] = {NULL, NULL, NULL};
   char errmsg[2048] = "";
   int stropt[8] = {0};

   stropt[0] = timeout;
   stropt[1] = reconnect;
   stropt[2] = 1000;
   stropt[3] = buffsize;
   stropt[4] = fswapmargin;
   strsetopt(stropt);

   solopt[0].posf = strfmt[3];
   solopt[1].posf = strfmt[4];

   if (!rtksvrstart(&svr, svrcycle, buffsize, strtype, paths, strfmt, navmsgsel,
                    cmds, cmds_periodic, ropts, nmeacycle, nmeareq, npos, &prcopt,
                    solopt, nullptr, errmsg))
   {
      ROS_ERROR("rtklib start error (%s)\n", errmsg);
      return 0;
   }

   return 1;
}

int main(int argc, char **argv)
{
   ros::init(argc, argv, "rtkrcv");

   ros::NodeHandle node;

   pos_pub = node.advertise<gnss_comm::GnssPVTSolnMsg>("pvt", 100);
   meas_pub = node.advertise<gnss_comm::GnssMeasMsg>("range_meas", 100);
   ephem_pub = node.advertise<gnss_comm::GnssEphemMsg>("ephem", 100);
   glo_ephem_pub = node.advertise<gnss_comm::GnssGloEphemMsg>("glo_ephem", 100);
   iono_pub = node.advertise<gnss_comm::StampedFloat64Array>("iono_params", 100);
   ssr_pub = node.advertise<rtkrcv::corr_msg>("corrections", 100);

   std::string config_file, output_dir;
   node.getParam("rtkrcv/config_file", config_file);
   node.getParam("rtkrcv/output_dir", output_dir);

   ROS_INFO("config file: %s", config_file.c_str());

   rtksvrinit(&svr);

   resetsysopts();
   if (!loadopts(config_file.c_str(), rcvopts) || !loadopts(config_file.c_str(), sysopts))
   {
      ROS_ERROR("no options file: %s.", opt_file);
      exit(-1);
   }
   getsysopts(&prcopt, solopt, &filopt);

   std::string output_file = output_dir + "/" + strpath[4];
   strcpy(strpath[4], output_file.c_str());
   ROS_INFO("output path: %s", strpath[4]);

   if (!start_rtklib())
      return -1;

   ros::spin();

   return 0;
}
