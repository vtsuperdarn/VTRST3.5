/* fit_raw.h
   ========
   Author: AJR
*/

#ifndef _FIT_RAW_H
#define _FIT_RAW_H

void FitACFFree(struct FitBlock *fptr); 
struct FitBlock *FitACFMake(struct RadarSite *hd,int year);
int fit_raw(struct RadarParm *prm,struct RawData *ptr,struct FitBlock *input,
            struct FitData *fit, int *fittype);
void setup_fblk(struct RadarParm *prm, struct RawData *raw,struct FitBlock *input)


#endif
