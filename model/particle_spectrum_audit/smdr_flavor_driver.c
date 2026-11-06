/* Forward-only audit driver for the pinned SMDR library.
 * Every input is supplied by the caller; no experimental values are loaded.
 * Modes beta1/beta2 expose polynomial-control evaluations. Mode match runs
 * the same W/Z/H/top/GF routines as the frozen comparison, at input Q0 only.
 * All nine Yukawas can be nonzero. SMDR's own documented approximations
 * remain; this wrapper does not make its matching fully massive/all-flavor.
 */
#include "smdr_internal.h"

static void reload(void) {SMDR_Load_Inputs(); SMDR_Update();}
static void obs(const char *name, long double m, long double w,
                long double bw, long double bww) {
  printf("\"%s\":{\"mass_GeV\":%.18Lg,\"width_GeV\":%.18Lg,"
         "\"BW_mass_GeV\":%.18Lg,\"BW_width_GeV\":%.18Lg}",name,m,w,bw,bww);
}
int main(int argc, char **argv) {
  if(argc!=17) {
    fprintf(stderr,"Expected mode Q0 v gY g2 g3 yu yc yt yd ys yb ye ymu ytau lambda\n");
    return 2;
  }
  SMDR_Q_in=strtold(argv[2],0); SMDR_v_in=strtold(argv[3],0);
  SMDR_gp_in=strtold(argv[4],0); SMDR_g_in=strtold(argv[5],0);
  SMDR_g3_in=strtold(argv[6],0);
  SMDR_yu_in=strtold(argv[7],0); SMDR_yc_in=strtold(argv[8],0);
  SMDR_yt_in=strtold(argv[9],0); SMDR_yd_in=strtold(argv[10],0);
  SMDR_ys_in=strtold(argv[11],0); SMDR_yb_in=strtold(argv[12],0);
  SMDR_ye_in=strtold(argv[13],0); SMDR_ymu_in=strtold(argv[14],0);
  SMDR_ytau_in=strtold(argv[15],0); SMDR_lambda_in=strtold(argv[16],0);
  SMDR_m2_in=-SMDR_lambda_in*SMDR_v_in*SMDR_v_in; SMDR_Lambda_in=0;
  SMDR_Delta_alpha_had_5_MZ_in=NAN;
  SMDR_MTPOLE_TOLERANCE=1e-12L;
#include "source_ew_poison.h"
  reload();
  if(!strcmp(argv[1],"beta1") || !strcmp(argv[1],"beta2")) {
    SMDR_Betas(argv[1][4]-'0');
    printf("[%.18Lg,%.18Lg,%.18Lg,%.18Lg,%.18Lg,%.18Lg,%.18Lg,"
           "%.18Lg,%.18Lg,%.18Lg,%.18Lg,%.18Lg,%.18Lg]\n",
           sqrtl(5.L/3.L)*beta_gp,beta_g,beta_g3,beta_yu,beta_yc,beta_yt,
           beta_yd,beta_ys,beta_yb,beta_ye,beta_ymu,beta_ytau,beta_k);
    return 0;
  }
  if(strcmp(argv[1],"match")) return 3;
  long double m20=SMDR_Eval_m2(-1,0),m21=SMDR_Eval_m2(-1,1),m22=SMDR_Eval_m2(-1,2);
  SMDR_m2_in=m22;
  printf("{\"schema\":\"oph.round2.SMDR_nonzero_flavor_forward.v1\","
         "\"m2_minimum_controls\":{\"0\":%.18Lg,\"1\":%.18Lg,\"2\":%.18Lg},\"rows\":[",m20,m21,m22);
  for(int order=0;order<=2;order++) {
    reload(); if(order) printf(",");
    printf("{\"scale_multiplier\":1,\"matching_order\":%d,\"working\":{"
           "\"Q\":%.18Lg,\"gY\":%.18Lg,\"g2\":%.18Lg,\"g3\":%.18Lg,"
           "\"yu\":%.18Lg,\"yc\":%.18Lg,\"yt\":%.18Lg,"
           "\"yd\":%.18Lg,\"ys\":%.18Lg,\"yb\":%.18Lg,"
           "\"ye\":%.18Lg,\"ymu\":%.18Lg,\"ytau\":%.18Lg,"
           "\"lambda\":%.18Lg,\"v\":%.18Lg,\"m2\":%.18Lg},",
           order,Q,gp,g,g3,yu,yc,yt,yd,ys,yb,ye,ymu,ytau,k,v,m2);
    long double mt,gt,mh,gh,mw,gw,mwb,gwb,mz,gz,mzb,gzb;
    SMDR_Eval_Mt_pole(-1,1,order,order,&mt,&gt);obs("top_method1",mt,gt,0,0);printf(",");
    if(order==2) {reload();SMDR_Eval_Mt_pole(-1,1,4,2,&mt,&gt);obs("top_QCD4_EW2",mt,gt,0,0);printf(",");}
    reload(); SMDR_Eval_Mh_pole(-1,order,&mh,&gh);obs("higgs",mh,gh,0,0);printf(",");
    reload(); SMDR_Eval_MW_pole(-1,order,&mw,&gw,&mwb,&gwb);obs("W",mw,gw,mwb,gwb);printf(",");
    reload(); SMDR_Eval_MZ_pole(-1,order,&mz,&gz,&mzb,&gzb);obs("Z",mz,gz,mzb,gzb);printf(",");
    reload(); long double gfraw=SMDR_Eval_GFermi(-1,order);
    long double gfmu=0.00000051862L/(sqrtl(2.L)*v*v);
    printf("\"GF_raw_GeVm2\":%.18Lg,\"GF_removed_finite_muon_term_GeVm2\":%.18Lg,"
           "\"GF_local_chiral_external_legs_GeVm2\":%.18Lg}",gfraw,gfmu,gfraw-gfmu);
  }
  printf("]}\n"); return 0;
}
