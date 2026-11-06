/* Exact additional one-loop light-fermion W/Z self-energy pieces.
 * The generic functions are the same fermion terms used in pinned SMDR
 * src/MW.c:244-270 and src/MZ.c:232-267, now summed over the first two
 * generations. Their already-present massless contributions are subtracted.
 * CKM=I and massless neutrinos are the declared running/action branch.
 */
#include "smdr_internal.h"
SMDR_COMPLEX SMDR_MW_1loop(void);
SMDR_COMPLEX SMDR_MZ_1loop(void);
int SMDR_MW_DoTSILW(float);
int SMDR_MZ_DoTSILZ(float);

static SMDR_COMPLEX wferm(SMDR_REAL x, SMDR_REAL y) {
  return (TSIL_A(x,Q2)*(2*x-2*y-4*W)
          +TSIL_A(y,Q2)*(-2*x+2*y-4*W)-4*W*(x+y)+4*W*W/3
          +TSIL_B(x,y,W,Q2)*(2*x*x-4*x*y+2*y*y+2*W*(x+y)-4*W*W))/(3*v2);
}
static SMDR_COMPLEX zferm(SMDR_REAL x,SMDR_REAL left,SMDR_REAL right) {
  SMDR_COMPLEX loop=TSIL_B(x,x,Z,Q2);
  return (left*left+right*right)*(2*(x-Z)*loop/3-4*TSIL_A(x,Q2)/3+2*(Z-6*x)/9)
         -4*left*right*x*loop;
}
static SMDR_COMPLEX light_w(SMDR_REAL factor) {
  SMDR_REAL mass[6]={yu*yu*v2/2,yc*yc*v2/2,yd*yd*v2/2,ys*ys*v2/2,ye*ye*v2/2,ymu*ymu*v2/2};
  for(int j=0;j<6;j++) mass[j]*=factor*factor;
  SMDR_COMPLEX zero=wferm(0,0);
  return 3*(wferm(mass[0],mass[2])+wferm(mass[1],mass[3])-2*zero)
          +wferm(0,mass[4])+wferm(0,mass[5])-2*zero;
}
static SMDR_COMPLEX light_z(SMDR_REAL factor) {
  SMDR_REAL up[2]={yu,yc},down[2]={yd,ys},lepton[2]={ye,ymu};
  SMDR_COMPLEX answer=0;
  for(int j=0;j<2;j++) {
    answer+=3*(zferm(factor*factor*up[j]*up[j]*v2/2,auL,auR)-zferm(0,auL,auR));
    answer+=3*(zferm(factor*factor*down[j]*down[j]*v2/2,adL,adR)-zferm(0,adL,adR));
    answer+=zferm(factor*factor*lepton[j]*lepton[j]*v2/2,aeL,aeR)-zferm(0,aeL,aeR);
  }
  return answer;
}
static void complex_json(SMDR_COMPLEX z) {printf("[%.18Lg,%.18Lg]",creall(z),cimagl(z));}
int main(int argc,char **argv) {
  if(argc!=17) return 2;
  SMDR_Q_in=strtold(argv[2],0);SMDR_v_in=strtold(argv[3],0);
  SMDR_gp_in=strtold(argv[4],0);SMDR_g_in=strtold(argv[5],0);SMDR_g3_in=strtold(argv[6],0);
  SMDR_yu_in=strtold(argv[7],0);SMDR_yc_in=strtold(argv[8],0);SMDR_yt_in=strtold(argv[9],0);
  SMDR_yd_in=strtold(argv[10],0);SMDR_ys_in=strtold(argv[11],0);SMDR_yb_in=strtold(argv[12],0);
  SMDR_ye_in=strtold(argv[13],0);SMDR_ymu_in=strtold(argv[14],0);SMDR_ytau_in=strtold(argv[15],0);
  SMDR_lambda_in=strtold(argv[16],0);SMDR_m2_in=-SMDR_lambda_in*SMDR_v_in*SMDR_v_in;
  SMDR_Lambda_in=0;SMDR_Delta_alpha_had_5_MZ_in=NAN;
#include "source_ew_poison.h"
  SMDR_Load_Inputs();SMDR_Update();
  SMDR_COMPLEX dw=light_w(1),dz=light_z(1);
  SMDR_COMPLEX mwzero=wferm(0,0)-4*W*W*(1-3*TSIL_B(0,0,W,Q2))/(9*v2);
  SMDR_COMPLEX mzzero=zferm(0,auL,auR)-(auL2+auR2)*Z*(2.L/9-2*TSIL_B(0,0,Z,Q2)/3);
  SMDR_COMPLEX zw=light_w(0),zz=light_z(0);
  SMDR_COMPLEX smallw=light_w(.0001L),smallz=light_z(.0001L);

  /* Independent source implementation control: vary only the three massive
   * SMDR flavors. Bosonic and already-massless terms cancel in the difference.
   * The difference must equal the generic functions above, with no fit. */
  SMDR_MW_DoTSILW(1);SMDR_COMPLEX sw0=SMDR_MW_1loop();
  SMDR_MZ_DoTSILZ(1);SMDR_COMPLEX sz0=SMDR_MZ_1loop();
  SMDR_COMPLEX fw0=3*wferm(T,b)+wferm(0,tau);
  SMDR_COMPLEX fz0=3*zferm(T,auL,auR)+3*zferm(b,adL,adR)+zferm(tau,aeL,aeR);
  yt*=.93L;yb*=1.1L;ytau*=.8L;SMDR_Update();
  SMDR_MW_DoTSILW(1);SMDR_COMPLEX sw1=SMDR_MW_1loop();
  SMDR_MZ_DoTSILZ(1);SMDR_COMPLEX sz1=SMDR_MZ_1loop();
  SMDR_COMPLEX fw1=3*wferm(T,b)+wferm(0,tau);
  SMDR_COMPLEX fz1=3*zferm(T,auL,auR)+3*zferm(b,adL,adR)+zferm(tau,aeL,aeR);
  printf("{\"delta_complex_pole_squared_GeV2\":{\"W\":");complex_json(ONELOOPFACTOR*dw);
  printf(",\"Z\":");complex_json(ONELOOPFACTOR*dz);
  printf("},\"raw_self_energy_coefficients_GeV2\":{\"W\":");complex_json(dw);
  printf(",\"Z\":");complex_json(dz);
  printf("},\"controls\":{\"analytic_massless_W_abs\":%.18Lg,\"analytic_massless_Z_abs\":%.18Lg,"
         "\"all_light_zero_W_abs\":%.18Lg,\"all_light_zero_Z_abs\":%.18Lg,"
         "\"light_Yukawas_times_1e_minus4_W_abs\":%.18Lg,\"light_Yukawas_times_1e_minus4_Z_abs\":%.18Lg,"
         "\"independent_SMDR_massive_difference_W_abs\":%.18Lg,\"independent_SMDR_massive_difference_Z_abs\":%.18Lg}}\n",
         cabsl(mwzero),cabsl(mzzero),cabsl(zw),cabsl(zz),cabsl(smallw),cabsl(smallz),
         cabsl(sw1-sw0-fw1+fw0),cabsl(sz1-sz0-fz1+fz0));
  return 0;
}
