DEBUG          = no

# set  fortran compiler
FC		= gfortran

# set cc compiler
CC		= gcc


CCOPT   = -O3 -m64 -D_LARGEFILE64_SOURCE -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64
FCOPT   = -finit-local-zero -cpp -O3 -D_LARGEFILE64_SOURCE -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64

ifeq ($(DEBUG),yes)
FFLAGS  = -g -Wall -fbacktrace -fbounds-check
CFLAGS	= -g -Wall
LDFLAGS	= -g
endif

FFLAGS += $(FCOPT)
CFLAGS += $(CCOPT)
LDFLAGS+=

LIBS 	= -ldl -lc -lpthread -lm -L./fftpack5/ -lfftpack5


OBJS= module_init_params.o cinv.o eispack.o module_data.o module_antenna.o quickselect.o sort1.o sort2.o select_raies.o sourc_cour.o  \
			vp_noise.o module_azim.o smooth.o pass1.o pass2.o multi_courant.o un_courant_music.o rad2ascii.o radial_LatLon_WGS84_r16.o  \
			read_alcov2.o sort2radial.o

HDRS= 

.SUFFIXES : .o .F90 .f90 .c .a .F

.F90.o  :
	$(FC) -c $(FFLAGS) $<

.F.o  :
	$(FC) -c $(FFLAGS) $<

.f90.o  :
	$(FC) -c $(FFLAGS) $<

.c.o    :
	$(CC) -c $(CFLAGS) $<


default all: libfftpack5.a sort2radial

libfftpack5.a: 
	-cd fftpack5 && make

sort2radial: $(OBJS) $(HDRS)
	$(FC) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

clean:
	-rm -f sort2radial
	-rm -f *.o
	-rm -f *.mod
	-rm -f *genmod.f90
	-cd fftpack5 && make clean
