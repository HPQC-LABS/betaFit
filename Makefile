# option to select compiler (intel 'ifort', or SUN 'f77', 'f90 or f95)
#   as in :   make FC=ifort
ifndef $FC
   FC = ifort                    # default compiler on this machine
endif
#
#option to choose level of optimization:  make DEBUG=1
ifndef DEBUG
    FFLAGS =  -O3 -u            # DEBUG not specified, so optimize
  else                           
    FFLAGS =  -g -u             # DEBUG=1  (basic DEBUG option)
    ifeq ($(DEBUG),2)           
        FFLAGS =  -C -g -u       # DEBUG=2  (higher-level DEBUG option)
    endif
#                 can add more options for  DEBUG=3 , etc., as desired
endif
#
# as usual, list the objects
#
OBJECTS = betaFITret.o nllssrr.o llsqf.o AF3x3ret.o

betafit: $(OBJECTS)
	$(FC) $(FFLAGS) $(OBJECTS) -o betafRet.x

#--------------------------------------------------------------
# Form from Sean Mcleod,  13 February 2008
#--------------------------------------------------------------





