#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../../CalcHEP_src/include/extern.h"
#include "../../CalcHEP_src/include/VandP.h"
#include "autoprot.h"
extern int  FError;
/*  Special model functions  */

int nModelParticles=20;
static ModelPrtclsStr ModelPrtcls_[20]=
{
  {"A","A", 22, "0","0",2,1,0}
, {"Z","Z", 23, "MZ","wZ",2,1,0}
, {"G","G", 21, "0","0",2,8,0}
, {"W+","W-", 24, "MW","wW",2,1,3}
, {"ne","Ne", 12, "0","0",1,1,0}
, {"e","E", 11, "0","0",1,1,-3}
, {"nm","Nm", 14, "0","0",1,1,0}
, {"m","M", 13, "Mm","0",1,1,-3}
, {"nl","Nl", 16, "0","0",1,1,0}
, {"l","L", 15, "Ml","0",1,1,-3}
, {"u","U", 2, "Mu","0",1,3,2}
, {"d","D", 1, "Md","0",1,3,-1}
, {"c","C", 4, "Mc","0",1,3,2}
, {"s","S", 3, "Ms","0",1,3,-1}
, {"t","T", 6, "Mt","wtop",1,3,2}
, {"b","B", 5, "Mb","0",1,3,-1}
, {"h","h", 25, "Mh","wh",0,1,0}
, {"~H+","~H-", 37, "MHC","wHC",0,1,3}
, {"~x1","~X1", 35, "Mdm1","wX1",0,1,0}
, {"~x2","~X2", 45, "Mdm2","wX2",0,1,0}
};
ModelPrtclsStr *ModelPrtcls=ModelPrtcls_; 
int nModelVars=28;
int nModelFunc=27;
static char*varNames_[55]={
 "EE","alfSMZ","SW","MZ","Q","Mtp","MbMb","McMc","wZ","wW"
,"Mm","Ml","Mu","Md","Ms","wtop","Mh","la3","la2","laS"
,"laS1","laS2","laS21","Mdm1","Mdm2","sinDm","MHC","muppS","CW","MW"
,"Lqcd","Mb","Mt","Mc","PI","V","cosDm","muSq","mu_sh","mu2q"
,"la4","ahF_c","ahF_b","ahF_t","ahF_l","a_hV_W","a_hS_Hc","aQCD_h","Rqcd_h","Mbp"
,"Mcp","Quq","Qdq","LGGh","LAAh"};
char**varNames=varNames_;
static REAL varValues_[55]={
   3.122300E-01,  1.184000E-01,  4.810000E-01,  9.118700E+01,  1.000000E+02,  1.735000E+02,  4.230000E+00,  1.270000E+00,  2.502000E+00,  2.094000E+00
,  1.057000E-01,  1.777000E+00,  1.000000E-02,  1.000000E-02,  2.000000E-01,  1.442000E+00,  1.250000E+02,  2.000000E-01,  1.000000E-02,  1.000000E-01
,  1.000000E-01,  1.000000E-01,  1.000000E-01,  5.000000E+01,  2.000000E+02,  1.000000E-01,  2.000000E+02,  1.000000E+02};
REAL*varValues=varValues_;
int calcMainFunc(void)
{
   int i;
   static REAL * VV=NULL;
   static int iQ=-1;
   static int cErr=1;
   REAL *V=varValues;
   FError=0;
   if(VV && cErr==0)
   { for(i=0;i<nModelVars;i++) if(i!=iQ && VV[i]!=V[i]) break;
     if(i==nModelVars)      {if(iQ>=0 && VV[iQ]!=V[iQ]) goto FirstQ; else return 0;} 
   }
  cErr=1;
   V[28]=sqrt(1-pow(V[2],2));
   if(!isfinite(V[28]) || FError) return 28;
   V[29]=V[3]*V[28];
   if(!isfinite(V[29]) || FError) return 29;
   V[30]=initQCD5(V[1],V[7],V[6],V[5]);
   if(!isfinite(V[30]) || FError) return 30;
 FirstQ:
 cErr=1;
   V[31]=MbEff(V[4]);
   if(!isfinite(V[31]) || FError) return 31;
   V[32]=MtEff(V[4]);
   if(!isfinite(V[32]) || FError) return 32;
   V[33]=McEff(V[4]);
   if(!isfinite(V[33]) || FError) return 33;
   V[34]=acos(-1);
   if(!isfinite(V[34]) || FError) return 34;
   V[35]=2*V[29]/(V[0])*V[2];
   if(!isfinite(V[35]) || FError) return 35;
   V[36]=sqrt(1-pow(V[25],2));
   if(!isfinite(V[36]) || FError) return 36;
   V[37]=pow(V[24],2)*pow(V[25],2)+pow(V[23],2)*pow(V[36],2)-pow(V[35],2)*V[20]/(2);
   if(!isfinite(V[37]) || FError) return 37;
   V[38]=-4*(pow(V[24],2)-pow(V[23],2))*V[25]*V[36]/(V[35])/(M_SQRT2);
   if(!isfinite(V[38]) || FError) return 38;
   V[39]=pow(V[26],2)-V[17]*pow(V[35],2)/(2);
   if(!isfinite(V[39]) || FError) return 39;
   V[40]=2*(pow(V[23],2)*pow(V[25],2)+pow(V[24],2)*pow(V[36],2)-pow(V[26],2))/(pow(V[35],2));
   if(!isfinite(V[40]) || FError) return 40;
   V[41]=-V[0]/(V[29])*V[33]/(V[2])/(2)/(V[33]);
   if(!isfinite(V[41]) || FError) return 41;
   V[42]=-V[0]/(V[29])*V[31]/(V[2])/(2)/(V[31]);
   if(!isfinite(V[42]) || FError) return 42;
   V[43]=-V[0]/(V[29])*V[32]/(V[2])/(2)/(V[32]);
   if(!isfinite(V[43]) || FError) return 43;
   V[44]=-V[0]/(V[29])*V[11]/(V[2])/(2)/(V[11]);
   if(!isfinite(V[44]) || FError) return 44;
   V[45]=V[0]*V[29]/(V[2])/(pow(V[29],2));
   if(!isfinite(V[45]) || FError) return 45;
   V[46]=-2/(V[0])*V[29]*V[2]*V[17]/(pow(V[26],2));
   if(!isfinite(V[46]) || FError) return 46;
   V[47]=alphaQCD(V[16])/(V[34]);
   if(!isfinite(V[47]) || FError) return 47;
   V[48]=sqrt(1+V[47]*(149/(double)((12))+V[47]*(68.6482-V[47]*212.447)));
   if(!isfinite(V[48]) || FError) return 48;
   V[49]=V[6]*(1+4/(double)((3))*alphaQCD(V[6])/(V[34]));
   if(!isfinite(V[49]) || FError) return 49;
   V[50]=V[7]*(1+4/(double)((3))*alphaQCD(V[7])/(V[34]));
   if(!isfinite(V[50]) || FError) return 50;
   V[51]=4/(double)((9));
   if(!isfinite(V[51]) || FError) return 51;
   V[52]=1/(double)((9));
   if(!isfinite(V[52]) || FError) return 52;
   V[53]=-cabs(hGGeven(V[16],V[47],3,1,3,V[49],V[42],1,3,V[50],V[41],1,3,V[5],V[43]));
   if(!isfinite(V[53]) || FError) return 53;
   V[54]=-cabs(V[52]*hAAeven(V[16],V[47],1,1,3,V[49],V[42])+V[51]*hAAeven(V[16],V[47],2,1,3,V[5],V[43],1,3,V[50],V[41])+hAAeven(V[16],V[47],3,1,1,V[11],V[44],2,1,V[29],V[45],0,1,V[26],V[46]));
   if(!isfinite(V[54]) || FError) return 54;
   if(VV==NULL) 
   {  VV=malloc(sizeof(REAL)*nModelVars);
      for(i=0;i<nModelVars;i++) if(strcmp(varNames[i],"Q")==0) iQ=i;
   }
   for(i=0;i<nModelVars;i++) VV[i]=V[i];
   cErr=0;
   return 0;
}
