/* make_fit.c
     ==========
     Author: R.J.Barnes
*/

/*
 (c) 2010 JHU/APL & Others - Please Consult LICENSE.superdarn-rst.3.2-beta-4-g32f7302.txt for more information.
 
 
 
*/

#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <string.h>
#include <time.h>
#include <zlib.h>
 
#include "rtypes.h"
#include "option.h"

#include "dmap.h"
#include "rprm.h"
#include "rawdata.h"
#include "fitblk.h"
#include "fitdata.h"
#include "radar.h"

#include "fit_raw.h"

#include "fitacf.h"
#include "rawread.h"
#include "fitwrite.h"

#include "oldrawread.h"
#include "oldfitwrite.h"

#include "errstr.h"
#include "hlpstr.h"


/*initialize objects*/
struct RadarParm *prm;
struct RawData *raw;
struct FitData *fit;
struct FitBlock *fblk;

struct RadarNetwork *network;  
struct Radar *radar;
struct RadarSite *site;

struct OptionData opt;

int main(int argc,char *argv[]){

    /* File format transistion
     * ------------------------
     * 
     * When we switch to the new file format remove any reference
     * to "new". Change the command line option "new" to "old" and
     * remove "old=!new".
     */


    unsigned char old=0;
    unsigned char new=0;

    char *envstr;
    int status;
    int arg;

    unsigned char help=0;
    unsigned char option=0;

    unsigned char vb=0;

    FILE *fp=NULL;
    struct OldRawFp *rawfp=NULL;
    FILE *fitfp=NULL;
    FILE *inxfp=NULL;  
    int irec=1;
    int drec=2;
    int dnum=0;
    int fitex2flg=0,lmfitflg=0;
    int fittype=0;
    time_t ctime;
    int c,n;
    char command[128];
    char tmstr[40];
 
    prm=RadarParmMake();
    raw=RawMake();
    fit=FitMake();

    OptionAdd(&opt,"-help",'x',&help);
    OptionAdd(&opt,"-option",'x',&option);

    OptionAdd(&opt,"vb",'x',&vb);

    OptionAdd(&opt,"new",'x',&new);

    OptionAdd(&opt,"fitex2",'x',&fitex2flg);
    OptionAdd(&opt,"lmfit",'x',&lmfitflg);

    arg=OptionProcess(1,argc,argv,&opt,NULL);

    old=!new;
    /*if fitex2 or lmfit is set, turn off fitacf*/
    fitacfflg = !(fitex2flg || lmfitflg);

    /*check for fit type, 0=fitacf (default), 1=fitex2, 2=lmfit*/
    if(fitex2flg) fittype = 1;
    else if(lmfitflg) fittype = 2;

    if (help==1){
        OptionPrintInfo(stdout,hlpstr);
        exit(0);
    }

    if (option==1){
        OptionDump(stdout,&opt);
        exit(0);
    }


    if ((old) && (argc-arg<2)){
        OptionPrintInfo(stdout,hlpstr);
        exit(-1);
    }

    envstr=getenv("SD_RADAR");
    if (envstr==NULL) {
        fprintf(stderr,"Environment variable 'SD_RADAR' must be defined.\n");
        exit(-1);
    }

    fp=fopen(envstr,"r");

    if (fp==NULL) {
        fprintf(stderr,"Could not locate radar information file.\n");
        exit(-1);
    }

    network=RadarLoad(fp);
    fclose(fp); 
    if (network==NULL) {
        fprintf(stderr,"Failed to read radar information.\n");
        exit(-1);
    }

    envstr=getenv("SD_HDWPATH");
    if (envstr==NULL) {
        fprintf(stderr,"Environment variable 'SD_HDWPATH' must be defined.\n");
        exit(-1);
    }

    /*load radar hardware tables*/
    RadarLoadHardware(envstr,network);
    

    /*read the rawacf data*/
    if(old){
        rawfp=OldRawOpen(argv[arg],NULL);
        if (rawfp==NULL) {
            fprintf(stderr,"File not found.\n");
            exit(-1);
        }
        status=OldRawRead(rawfp,prm,raw);  
    } 
    else{ 
        if (arg==argc) fp=stdin;
        else fp=fopen(argv[arg],"r");
        if (fp==NULL){
            fprintf(stderr,"File not found.\n");
            exit(-1);
        }
        status=RawFread(fp,prm,raw);
    }

    /*load radar data*/
    radar=RadarGetRadar(network,prm->stid);
    if(radar==NULL){
        fprintf(stderr,"Failed to get radar information.\n");
        exit(-1);
    }

    /*get radar site information*/
    site=RadarYMDHMSGetSite(radar,prm->time.yr,prm->time.mo,
                            prm->time.dy,prm->time.hr,prm->time.mt,
                            prm->time.sc);
    if(site==NULL){
        fprintf(stderr,"Failed to get site information.\n");
        exit(-1);
    }


    command[0]=0;
    n=0;
    for (c=0;c<argc;c++){
        n+=strlen(argv[c])+1;
        if (n>127) break;
        if (c !=0) strcat(command," ");
        strcat(command,argv[c]);
    }


    /*check for verbose operation*/
    if(vb) 
        fprintf(stderr,"%d-%d-%d %d:%d:%d beam=%d\n",prm->time.yr,prm->time.mo,
                prm->time.dy,prm->time.hr,prm->time.mt,prm->time.sc,prm->bmnum);

    /*set up fblk object*/
    fblk=FitACFMake(site,prm->time.yr); 

    /*do the fitting*/
    fit_raw(prm,raw,fblk,fit,fittype);
    /*FitACF(prm,raw,fblk,fit);*/
    
    /*open the file(s)*/
    if(old){
        char vstr[256];
        fitfp=fopen(argv[arg+1],"w");
        if(fitfp==NULL){
            fprintf(stderr,"Could not create fit file.\n");
            exit(-1);
        }
        if(argc-arg>2){
            inxfp=fopen(argv[arg+2],"w");
            if(inxfp==NULL){
                fprintf(stderr,"Could not create inx file.\n");
                exit(-1);
            }
        }
        sprintf(vstr,"%d.%d",fit->revision.major,fit->revision.minor);
        OldFitHeaderFwrite(fitfp,"make_fit","fitacf",vstr);
        if (inxfp !=NULL) OldFitInxHeaderFwrite(inxfp,prm);
    }


    /*loop through all records*/
    do{
        /*set some params*/
        ctime = time((time_t) 0);
        RadarParmSetOriginCommand(prm,command);
        strcpy(tmstr,asctime(gmtime(&ctime)));
        tmstr[24]=0;
        RadarParmSetOriginTime(prm,tmstr);  

        /*write the fit file*/
        if(old){
            dnum=OldFitFwrite(fitfp,prm,fit,NULL);
            if (inxfp !=NULL) OldFitInxFwrite(inxfp,drec,dnum,prm);
            drec+=dnum;
            irec++;
        } 
        else status=FitFwrite(stdout,prm,fit);
        
        /*read the next record*/
        if (old) status=OldRawRead(rawfp,prm,raw);
        else status=RawFread(fp,prm,raw);

        /*verbose output*/
        if (vb) 
            fprintf(stderr,"%d-%d-%d %d:%d:%d beam=%d\n",prm->time.yr,prm->time.mo,
                    prm->time.dy,prm->time.hr,prm->time.mt,prm->time.sc,prm->bmnum);


        /*do the fitting*/
        if (status==0) fit_raw(prm,raw,fblk,fit,&fittype);
        /*if (status==0) FitACF(prm,raw,fblk,fit);*/

    
    } while(status==0);
    
    FitACFFree(fblk);
    if (old) OldRawClose(rawfp);
    return 0;
}













