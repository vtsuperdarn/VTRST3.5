#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <unistd.h>
#include <zlib.h>
#include "rmath.h"
#include "rtypes.h"
#include "dmap.h"
#include "rprm.h"
#include "radar.h"
#include "rawdata.h"
#include "fitdata.h"
#include "fitblk.h"
#include "fit_acf.h"
#include "fitacfversion.h"

#define GOOSEBAY 1

/*function to free memory allocated for a fitblk object*/
void FitACFFree(struct FitBlock *fptr){
    if(fptr->prm.pulse !=NULL) free(fptr->prm.pulse);
    if(fptr->prm.lag[0] !=NULL) free(fptr->prm.lag[0]);
    if(fptr->prm.lag[1] !=NULL) free(fptr->prm.lag[1]);
    if(fptr->acfd !=NULL) free(fptr->acfd);
    if(fptr->xcfd !=NULL) free(fptr->xcfd);
}

 /*function to create a fitblk object*/
struct FitBlock *FitACFMake(struct RadarSite *hd, int year){
    int i;
    struct FitBlock *fptr;
    fptr=malloc(sizeof(struct FitBlock));
    if (fptr==NULL) return NULL;
    if (year < 1993) fptr->prm.old=1;
    for (i=0;i<3;i++) fptr->prm.interfer[i]=hd->interfer[i];
    fptr->prm.bmsep=hd->bmsep;
    fptr->prm.phidiff=hd->phidiff;
    fptr->prm.tdiff=hd->tdiff;
    fptr->prm.vdir=hd->vdir;
    fptr->prm.maxbeam=hd->maxbeam;
    fptr->prm.pulse=NULL;
    fptr->prm.lag[0]=NULL;
    fptr->prm.lag[1]=NULL;
    fptr->prm.pwr0=NULL;
    fptr->acfd=NULL;
    fptr->xcfd=NULL;
    return fptr;
}

/*a function to set up a fitblk object for use in fitting*/
void setup_fblk(struct RadarParm *prm, struct RawData *raw,struct FitBlock *fblk){
    int i,j,n;
    void *tmp=NULL;

    if (prm->time.yr < 1993) fblk->prm.old=1;

    /*initialize some params*/
    fit->revision.major=FITACF_MAJOR_REVISION;
    fit->revision.minor=FITACF_MINOR_REVISION;
    fblk->prm.xcf=prm->xcf;
    fblk->prm.tfreq=prm->tfreq;
    fblk->prm.noise=prm->noise.search;
    fblk->prm.nrang=prm->nrang;
    fblk->prm.smsep=prm->smsep;
    fblk->prm.nave=prm->nave;
    fblk->prm.mplgs=prm->mplgs;
    fblk->prm.mpinc=prm->mpinc;
    fblk->prm.txpl=prm->txpl;
    fblk->prm.lagfr=prm->lagfr;
    fblk->prm.mppul=prm->mppul;
    fblk->prm.bmnum=prm->bmnum;
    fblk->prm.cp=prm->cp;
    fblk->prm.channel=prm->channel;
    fblk->prm.offset=prm->offset; /* stereo offset */

    /* need to incorporate Sessai's code for setting the offset
        for legacy data here.*/

    if(fblk->prm.pulse==NULL) tmp=malloc(sizeof(int)*fblk->prm.mppul);
    else tmp=realloc(fblk->prm.pulse,sizeof(int)*fblk->prm.mppul);
    if(tmp==NULL) return -1;
    fblk->prm.pulse=tmp;
    for(i=0;i<fblk->prm.mppul;i++) fblk->prm.pulse[i]=prm->pulse[i];

    for(n=0;n<2;n++){
        if (fblk->prm.lag[n]==NULL) tmp=malloc(sizeof(int)*(fblk->prm.mplgs+1));
        else tmp=realloc(fblk->prm.lag[n],sizeof(int)*(fblk->prm.mplgs+1));
        if (tmp==NULL) return -1;
        fblk->prm.lag[n]=tmp;
        for(i=0;i<=fblk->prm.mplgs;i++) fblk->prm.lag[n][i]=prm->lag[n][i];
    }



    if (fblk->prm.pwr0==NULL) tmp=malloc(sizeof(int)*fblk->prm.nrang);
    else tmp=realloc(fblk->prm.pwr0,sizeof(int)*fblk->prm.nrang); 
    if (tmp==NULL) return -1;
    fblk->prm.pwr0=tmp;

    if (fblk->acfd==NULL) tmp=malloc(sizeof(struct complex)*fblk->prm.nrang*
                                    fblk->prm.mplgs);
    else tmp=realloc(fblk->acfd,sizeof(struct complex)*fblk->prm.nrang*
                                   fblk->prm.mplgs); 
    if (tmp==NULL) return -1;
    fblk->acfd=tmp;

    if (fblk->xcfd==NULL) tmp=malloc(sizeof(struct complex)*fblk->prm.nrang*
                                        fblk->prm.mplgs);
    else tmp=realloc(fblk->xcfd,sizeof(struct complex)*fblk->prm.nrang*
                    fblk->prm.mplgs); 
    if (tmp==NULL) return -1;
    fblk->xcfd=tmp;

    memset(fblk->acfd,0,sizeof(struct complex)*fblk->prm.nrang*
            fblk->prm.mplgs);   
    memset(fblk->xcfd,0,sizeof(struct complex)*fblk->prm.nrang*
            fblk->prm.mplgs);   



    for (i=0;i<fblk->prm.nrang;i++) {
        fblk->prm.pwr0[i]=raw->pwr0[i];
        if (raw->acfd[0] !=NULL) {
            for (j=0;j<fblk->prm.mplgs;j++) {
                fblk->acfd[i*fblk->prm.mplgs+j].x=raw->acfd[0][i*fblk->prm.mplgs+j];
                fblk->acfd[i*fblk->prm.mplgs+j].y=raw->acfd[1][i*fblk->prm.mplgs+j];
            }
        }
        if (raw->xcfd[0] !=NULL) {
            for (j=0;j<fblk->prm.mplgs;j++) {
                fblk->xcfd[i*fblk->prm.mplgs+j].x=raw->xcfd[0][i*fblk->prm.mplgs+j];
                fblk->xcfd[i*fblk->prm.mplgs+j].y=raw->xcfd[1][i*fblk->prm.mplgs+j];
            }
        }
    }
}

/*top level function for doing the fitting, sets things up for fitacf, fitex2, or lmfit*/
int fit_raw(struct RadarParm *prm, struct RawData *raw,
            struct FitBlock *fblk, struct FitData *fit){

    int goose;
    int lag_lim = 5;
    struct FitACFBadSample badsmp;
    int *badlag=NULL;
    int i=0,k;
    double *pwrd=NULL,*pwrt=NULL;
    double mnpwr, skylog, freq_to_vel, range;
    double xomega=0.0;
    double noise_pwr=0.0; 
    int ni;

    /*set up the fblk object, initialze range gates*/
    setup_fblk(prm, raw, fblk);
    FitSetRng(fit,fblk->prm.nrang);
    if(fblk->prm.xcf){
        FitSetXrng(fit,fblk->prm.nrang);
        FitSetElv(fit,fblk->prm.nrang);
    }
  
    /*goose bay has intereferometer in front of main array*/
    goose=(prm->stid==GOOSEBAY);

    /*initialize noise levels*/
    fit.noise->skynoise=0.0;
    fit.noise->lag0=0.0;
    fit.noise->vel=0.0;

    if (fblk->prm.nave <= 1) return 0;

    freq_to_vel = C/(4*PI)/(fblk->prm.tfreq * 1000.0);

    badlag=malloc(sizeof(int)*fblk->prm.nrang*fblk->prm.mplgs);
    if (badlag==NULL) return -1;

    pwrd=malloc(sizeof(double)*fblk->prm.nrang);
    if (pwrd==NULL){
        free(badlag);
        return -1;
    }
    pwrt=malloc(sizeof(double)*fblk->prm.nrang);
    if (pwrt==NULL){
        free(badlag);
        free(pwrd);
        return -1;
    }

    if (fblk->prm.channel==0) FitACFBadlags(&fblk->prm,&badsmp);  
    else FitACFBadlagsStereo(&fblk->prm,&badsmp);  


    /* Determine the lag_0 noise level (0 dB reference) and the noise level at 
    which fit_acf is to quit (average power in the 
    fluctuations of the acfs which are pure noise) */

    for(i=0; i < fblk->prm.nrang; i++){
        pwrd[i] = (double)fblk->prm.pwr0[i]; 
        /* transfer powers into local array */
        pwrt[i] = pwrd[i];
    }
    qsort(pwrt, fblk->prm.nrang, sizeof(double), dbl_cmp);
    /* determine the average lag0 power of the 10 lowest power acfs */
    
    ni = 0;
    i = 0;
    mnpwr = 0.0;
    
    /*  look for the lowest 10 values of lag0 power and average to 
      get the noise level.  Ignore values that are exactly 0.  If
      you can't find 10 useable values within the first 1/3 of the
      sorted power list, then just use whatever you got in that 
      first 1/3.  If you didn't get any useable values, then use
      the NOISE parameter */
        
    while ((ni < 10) && (i < fblk->prm.nrang/3)){
        if(pwrt[i]) ++ni;
        mnpwr += pwrt[i++];  
    }

    ni = (ni > 0) ? ni :  1;
    mnpwr = mnpwr/ni;
    if (mnpwr < 1.0) mnpwr = fblk->prm.noise;
    fit.noise->skynoise = mnpwr;

    /* Now determine the level which will be used as the cut-off power 
        for fit_acf.  This is the average power at all non-zero lags of all
        acfs which have lag0 power < 1.6*mnpwr + 1 stnd. deviation from that
        average power level */

    noise_pwr = noise_stat(mnpwr,&fblk->prm,&badsmp,fblk->acfd); 

    /* convert the lag0 powers to dB */

    if(&fit->noise->skynoise > 0.0) skylog = 10.0 * log10(&fit->noise->skynoise);
    else skylog = 0.0;

    for(i=0; i<fblk->prm.nrang; i++){   

        pwrd[i] = pwrd[i] - &fit->noise->skynoise;
        if(pwrd[i] <= 0.0) fit->rng[i].p_0 = -50.0;
        else fit->rng[i].p_0 = 10.0*log10(pwrd[i]) - skylog;
    }

    /*    reset the output arrays */

    for(i=0; i<fblk->prm.nrang; i++){
        fit->rng[i].p_l = -50.0;
        fit->rng[i].p_s = -50.0;
        fit->rng[i].p_l_err= 0.0;
        fit->rng[i].p_s_err= 0.0;
        fit->rng[i].w_l = 0.0;
        fit->rng[i].w_s = 0.0;
        fit->rng[i].w_l_err = 0.0;
        fit->rng[i].w_s_err = 0.0;
        fit->rng[i].v = 0.0;
        fit->rng[i].v_err = 0.0;
        fit->rng[i].phi0 = 0.0;
        fit->rng[i].phi0_err=0.0;
        fit->rng[i].sdev_l = 0.0;
        fit->rng[i].sdev_s = 0.0;
        fit->rng[i].sdev_phi = 0.0;
        fit->rng[i].gsct = 0.0;
        fit->rng[i].qflg = 0;
        fit->rng[i].nump=0;
        if(fit->xrng !=NULL){
            fit->xrng[i].p_l = -50.0;
            fit->xrng[i].p_s = -50.0;
            fit->xrng[i].p_l_err= 0.0;
            fit->xrng[i].p_s_err= 0.0;
            fit->xrng[i].w_l = 0.0;
            fit->xrng[i].w_s = 0.0;
            fit->xrng[i].w_l_err = 0.0;
            fit->xrng[i].w_s_err = 0.0;
            fit->xrng[i].v = 0.0;
            fit->xrng[i].v_err = 0.0;
            fit->xrng[i].phi0 = 0.0;
            fit->xrng[i].phi0_err=0.0;
            fit->xrng[i].sdev_l = 0.0;
            fit->xrng[i].sdev_s = 0.0;
            fit->xrng[i].sdev_phi = 0.0;
            fit->xrng[i].gsct = 0.0;
            fit->xrng[i].qflg = 0;
            fit->xrng[i].nump=0;

            fit->elv[i].normal= 0.0;
            fit->elv[i].low = 0.0;
            fit->elv[i].high = 0.0;
        }
    }

    /* ----------------------------------------------------------------------*/
    /*    Now do the fits for each acf */

    for (k=0, i=0; k<fblk->prm.nrang;k++){

        fit->rng[k].qflg = fit_acf(&fblk->acfd[k*fblk->prm.mplgs], k+1,
                                    &badlag[k*fblk->prm.mplgs],&badsmp,
                                    lag_lim,&fblk->prm,noise_pwr,0,0.0,&fit->rng[k]);
        xomega=fit->rng[k].v;
        if(fit->rng[k].qflg == 1){
            /* several changes have been made here to 
                fix an apparent problem in handling HUGE_VAL.
             
                If there are too few points in an ACF to allow
                the error on a parameter to be calculated then
                the subroutine fit_acf sets the value to HUGE_VAL.

                However, in this routine the error values are converted
                to natural units (e.g. velocity instead of frequency).
                It appears that multiplying HUGE_VAL by something causes
                a floating point exception that then sets the result of
                the calculation to 0.  Thus the error values that were being
                stored in the file would be zero instead of HUGE_VAL.

                The code now checks to see if the value is set to
                HUGE_VAL before doing the conversion.  If it is then
                instead of a converted version the error value is
                reset to HUGE_VAL.
            */

            /* convert power from natural log to dB */

            fit->rng[k].p_l = fit->rng[k].p_l*LN_TO_LOG - skylog;
            fit->rng[k].p_s = fit->rng[k].p_s*LN_TO_LOG - skylog;

            it->rng[k].p_l_err = (fit->rng[k].p_l_err == HUGE_VAL) ?
                                    HUGE_VAL:fit->rng[k].p_l_err*LN_TO_LOG;

            fit->rng[k].p_s_err = (fit->rng[k].p_s_err == HUGE_VAL) ?
                                    HUGE_VAL:fit->rng[k].p_s_err*LN_TO_LOG;

            /* convert Doppler frequency to velocity */
            fit->rng[k].v = fblk->prm.vdir*freq_to_vel*fit->rng[k].v;

            /* flag absurdly high velocities with qflg of 8 */
            if(fit->rng[k].v > (freq_to_vel* (PI* 1000.0* 1000.0)/ fblk->prm.mpinc))
                fit->rng[k].qflg= 8;     
      
            fit->rng[k].v_err = (fit->rng[k].v_err == HUGE_VAL) ?
                                HUGE_VAL:freq_to_vel*fit->rng[k].v_err;

            /* convert decay parameters to spectral widths */

            fit->rng[k].w_l = freq_to_vel*2*fit->rng[k].w_l;
            fit->rng[k].w_l_err = (fit->rng[k].w_l_err == HUGE_VAL) ?
                                    HUGE_VAL:freq_to_vel*2*fit->rng[k].w_l_err;

            /* sigma is returned as sigma**2 so check the sign for validity
            if sigma**2 is negative take sqrt of the abs and transfer the sign */

            fit->rng[k].w_s = (fit->rng[k].w_s >= 0) ? 
                                sqrt(fit->rng[k].w_s) : -sqrt(-fit->rng[k].w_s);


            if ((fit->rng[k].w_s !=0.0) && (fit->rng[k].w_s_err != HUGE_VAL))  
                fit->rng[k].w_s_err = 0.5*fit->rng[k].w_s_err/fabs(fit->rng[k].w_s);
            else fit->rng[k].w_s_err=HUGE_VAL;

            fit->rng[k].w_s = 3.33*freq_to_vel*fit->rng[k].w_s;
            fit->rng[k].w_s_err = (fit->rng[k].w_s_err == HUGE_VAL) ?
                        HUGE_VAL :
                        3.33*freq_to_vel*fit->rng[k].w_s_err;
       

            /*  Now check the values of power, velocity and width
            to see if this should be flagged as ground-scatter */
            
            if (fit->rng[k].gsct == 0) fit->rng[k].gsct=ground_scatter(&fit->rng[k]); 
        }
    
        if((fblk->prm.xcf==0) || (fit->rng[k].qflg !=1)){
            if((fit->rng[k].qflg == 1)) i++;
            continue;
        }

        if(fit->rng[k].qflg == 1) i++;
    }

    free(badlag);
    free(pwrd);
    free(pwrt);

    return 0;
}
