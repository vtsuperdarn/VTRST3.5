# Makefile for hmb_plot
# =====================
# Author: R.J.Barnes
# by R.J.Barnes
#
#
include $(MAKECFG).$(SYSTEM)


INCLUDE=-I$(IPATH)/base -I$(IPATH)/general -I$(IPATH)/analysis \
        -I$(IPATH)/superdarn 
SRC=hlpstr.h errstr.h hmb_plot.c
OBJS = hmb_plot.o
DSTPATH = $(BINPATH)
OUTPUT = hmb_plot

SLIB=-lpng -lm -lz
LIBS=-lcnvmap.1 -loldcnvmap.1 -lgrd.1 -loldgrd.1 -lrfile.1 -ldmap.1 -lrtime.1 -lopt.1 -lrcnv.1 -lgrplot.1 -liplot.1 -lsplot.1 -lfontdb.1 -lrps.1 -lrfbuf.1 -lrxml.1 
XLIB=-lxwin.1

include $(MAKEBIN).$(SYSTEM)
