/* Forward-only light matching. No measured masses, pole seeds, constituent
 * mass proxies, alpha(0), or hadronic vacuum-polarization data are consumed.
 * SMDR source remains unchanged. OPH_Decouple_thZW is a generated, minimally
 * corrected adapter to Martin, Phys.Rev.D99 (2019) 033007, arXiv:1812.04100.
 */
#include "smdr_internal.h"
void OPH_Decouple_thZW(int);
SMDR_REAL OPH_Fthresh(SMDR_REAL, SMDR_REAL, SMDR_REAL);
SMDR_REAL SMDR_Fthresh(SMDR_REAL, SMDR_REAL, SMDR_REAL);
SMDR_REAL SMDR_deltam1f(SMDR_REAL, SMDR_REAL);

static long double number(const char *s) { return strtold(s, NULL); }
static void poison(void) {
#include "source_ew_poison.h"
  SMDR_Delta_alpha_had_5_MZ_in=NAN;
}
static void output53(void) {
  printf("{\"Q\":%.18Lg,\"alphaS\":%.18Lg,\"alpha\":%.18Lg,"
    "\"masses\":{\"u\":%.18Lg,\"d\":%.18Lg,\"s\":%.18Lg,"
    "\"c\":%.18Lg,\"b\":%.18Lg,\"e\":%.18Lg,"
    "\"mu\":%.18Lg,\"tau\":%.18Lg}}\n",
    Q_53,alphaS_53,alpha_53,mu_53,md_53,ms_53,mc_53,mb_53,
    melectron_53,mmuon_53,mtau_53);
}
static void one_loop_decouple(void) {
  long double a=g2*gp2/(4*PI*g2pgp2),as=g32/(4*PI);
  long double lt=logl(T/Q2),lw=logl(W/Q2);
  long double de=ONELOOPFACTOR*SMDR_deltam1f(-1,-0.5);
  long double du=ONELOOPFACTOR*SMDR_deltam1f(2.L/3,0.5);
  long double dd=ONELOOPFACTOR*SMDR_deltam1f(-1.L/3,-0.5);
  long double db=dd+ONELOOPFACTOR*.75L*yt2*(W2*logl(T/W)/((T-W)*(T-W))
    -W/(T-W)+5.L/6-lt);
  Q_53=Q;
  alphaS_53=as*(1+ONELOOPFACTOR*g32*2.L/3*lt);
  alpha_53=a*(1+a/(4*PI)*(2.L/3-7*lw+16.L/9*lt));
  mu_53=yu*v/SQRT2*(1+du);mc_53=yc*v/SQRT2*(1+du);
  md_53=yd*v/SQRT2*(1+dd);ms_53=ys*v/SQRT2*(1+dd);
  mb_53=yb*v/SQRT2*(1+db);
  melectron_53=ye*v/SQRT2*(1+de);mmuon_53=ymu*v/SQRT2*(1+de);
  mtau_53=ytau*v/SQRT2*(1+de);
}
int main(int argc,char **argv) {
  if(argc<2) return 2;
  poison();
  if(!strcmp(argv[1],"run")) {
    if(argc!=10) return 3;
    long double as,a,cu,cd,ce;
    int code=SMDR_RGrun_QCDQED(number(argv[2]),number(argv[3]),atoi(argv[7]),
      atoi(argv[4]),atoi(argv[5]),atoi(argv[6]),number(argv[8]),number(argv[9]),
      &as,&a,&cu,&cd,&ce);
    printf("{\"code\":%d,\"alphaS\":%.18Lg,\"alpha\":%.18Lg,"
      "\"cu\":%.18Lg,\"cd\":%.18Lg,\"ce\":%.18Lg}\n",code,as,a,cu,cd,ce);
    return 0;
  }
  if(!strcmp(argv[1],"threshold")) {
    if(argc!=9) return 4;
    long double as,a,zu,zd,ze;
    SMDR_QCDQEDmatch(argv[2],number(argv[3]),number(argv[4]),number(argv[5]),
      number(argv[6]),atoi(argv[7]),atoi(argv[8]),&as,&a,&zu,&zd,&ze);
    printf("{\"alphaS\":%.18Lg,\"alpha\":%.18Lg,\"zu\":%.18Lg,"
      "\"zd\":%.18Lg,\"ze\":%.18Lg}\n",as,a,zu,zd,ze);
    return 0;
  }
  if(!strcmp(argv[1],"Fcontrol")) {
    printf("{\"fixed_base\":%.18Lg,\"fixed_rescaled\":%.18Lg,"
      "\"fixed_changed_Q\":%.18Lg,\"original_base\":%.18Lg,"
      "\"original_rescaled\":%.18Lg}\n",
      OPH_Fthresh(2,3,5),OPH_Fthresh(14,21,35),OPH_Fthresh(2,3,11),
      SMDR_Fthresh(2,3,5),SMDR_Fthresh(14,21,35));
    return 0;
  }
  if(argc!=17) return 5;
  SMDR_Q_in=number(argv[2]);SMDR_v_in=number(argv[3]);
  SMDR_gp_in=number(argv[4]);SMDR_g_in=number(argv[5]);SMDR_g3_in=number(argv[6]);
  SMDR_yu_in=number(argv[7]);SMDR_yc_in=number(argv[8]);SMDR_yt_in=number(argv[9]);
  SMDR_yd_in=number(argv[10]);SMDR_ys_in=number(argv[11]);SMDR_yb_in=number(argv[12]);
  SMDR_ye_in=number(argv[13]);SMDR_ymu_in=number(argv[14]);SMDR_ytau_in=number(argv[15]);
  SMDR_lambda_in=number(argv[16]);
  SMDR_m2_in=-SMDR_lambda_in*SMDR_v_in*SMDR_v_in;SMDR_Lambda_in=0;
  SMDR_Load_Inputs();SMDR_Update();
  if(!strcmp(argv[1],"fixed")) OPH_Decouple_thZW(4);
  else if(!strcmp(argv[1],"original")) SMDR_Decouple_thZW(4);
  else if(!strcmp(argv[1],"one")) one_loop_decouple();
  else return 6;
  output53();return 0;
}
