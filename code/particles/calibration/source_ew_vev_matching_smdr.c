/* Conditional full-SM matching bridge. GPL-2.0-or-later SMDR is linked,
 * not vendored. No reference model or experimental value is loaded.
 * All Yukawas except yt vanish. The empirical finite-muon phase-space
 * term is algebraically removed to report the local charged-current
 * Wilson coefficient in this chiral truncation. Every scale is transported
 * from the SAME Q0 boundary with two-loop RG, including the running VEV.
 */
#include "smdr.h"
static void reset_at(long double q) { SMDR_RGeval_SM(q,2); SMDR_Update(); }
static void obs(const char*name,long double m,long double w,long double bw,long double bww){printf("\"%s\":{\"mass_GeV\":%.18Lg,\"width_GeV\":%.18Lg,\"BW_mass_GeV\":%.18Lg,\"BW_width_GeV\":%.18Lg}",name,m,w,bw,bww);}
int main(int argc, char **argv){
 if(argc != 8) {fprintf(stderr,"Expected Q0 v gY g2 g3 yt lambda\n"); return 2;}
 SMDR_Q_in=strtold(argv[1],0); SMDR_v_in=strtold(argv[2],0);
 SMDR_gp_in=strtold(argv[3],0); SMDR_g_in=strtold(argv[4],0);
 SMDR_g3_in=strtold(argv[5],0); SMDR_yt_in=strtold(argv[6],0);
 SMDR_lambda_in=strtold(argv[7],0);
 SMDR_yb_in=SMDR_yc_in=SMDR_ys_in=SMDR_yu_in=SMDR_yd_in=0;
 SMDR_ytau_in=SMDR_ymu_in=SMDR_ye_in=0;
 SMDR_Lambda_in=0; SMDR_m2_in=-SMDR_lambda_in*SMDR_v_in*SMDR_v_in;
 SMDR_Delta_alpha_had_5_MZ_in=NAN;
 SMDR_MTPOLE_TOLERANCE=1e-12L;
#include "source_ew_poison.h"
 SMDR_Load_Inputs(); SMDR_Update();
 long double m20=SMDR_Eval_m2(-1,0),m21=SMDR_Eval_m2(-1,1),m22=SMDR_Eval_m2(-1,2);
 SMDR_m2_in=m22;
 printf("{\"schema\":\"oph.smdr_forward_matching.v1\",\"m2_minimum_controls\":{\"0\":%.18Lg,\"1\":%.18Lg,\"2\":%.18Lg},\"rows\":[",m20,m21,m22);
 int first=1;
 for(int mult=1;mult<=4;mult*=2){
  long double q=SMDR_Q_in*mult;
  for(int order=0;order<=2;order++){
   reset_at(q);
   if(!first)printf(","); first=0;
   printf("{\"scale_multiplier\":%d,\"matching_order\":%d,\"working\":{\"Q\":%.18Lg,\"gY\":%.18Lg,\"g2\":%.18Lg,\"g3\":%.18Lg,\"yt\":%.18Lg,\"lambda\":%.18Lg,\"v\":%.18Lg,\"m2\":%.18Lg,\"mt_MSbar\":%.18Lg},",mult,order,SMDR_Q,SMDR_gp,SMDR_g,SMDR_g3,SMDR_yt,SMDR_lambda,SMDR_v,SMDR_m2,SMDR_yt*SMDR_v/sqrtl(2));
   long double mt,gt,mh,gh,mw,gw,mwb,gwb,mz,gz,mzb,gzb;
   SMDR_Eval_Mt_pole(-1,0,order,order,&mt,&gt);obs("top_method0",mt,gt,0,0);printf(",");
   reset_at(q); SMDR_Eval_Mt_pole(-1,1,order,order,&mt,&gt);obs("top_method1",mt,gt,0,0);printf(",");
   if(order==2){reset_at(q);SMDR_Eval_Mt_pole(-1,1,4,2,&mt,&gt);obs("top_QCD4_EW2",mt,gt,0,0);printf(",");}
   reset_at(q);SMDR_Eval_Mh_pole(-1,order,&mh,&gh);obs("higgs",mh,gh,0,0);printf(",");
   reset_at(q);SMDR_Eval_MW_pole(-1,order,&mw,&gw,&mwb,&gwb);obs("W",mw,gw,mwb,gwb);printf(",");
   reset_at(q);SMDR_Eval_MZ_pole(-1,order,&mz,&gz,&mzb,&gzb);obs("Z",mz,gz,mzb,gzb);printf(",");
   reset_at(q);long double gfraw=SMDR_Eval_GFermi(-1,order),gfmu=0.00000051862L/(sqrtl(2)*SMDR_v*SMDR_v),gf=gfraw-gfmu;
   printf("\"GF_raw_GeVm2\":%.18Lg,\"GF_removed_finite_muon_term_GeVm2\":%.18Lg,\"GF_local_zero_muon_GeVm2\":%.18Lg}",gfraw,gfmu,gf);
  }
 }
 printf("]}\n"); return 0;
}
