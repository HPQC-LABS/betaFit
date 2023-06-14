# option to select compiler (intel 'ifort', or SUN 'f77', 'f90 or f95)
#   as in :   make FC=ifort
ifndef $FC
   FC = nagfor
#  FC = f90                    # default compiler on this machine
#  FC = gfortran               # a compiler for testing purposes
endif
#
#option to choose level of optimization:  make debug=1
ifndef debug
    FFLAGS =  -O2             # debug not specified, so optimize
  else                           
#   FFLAGS =  -q  -C          # debug=1  (basic gtortran debug option)
#    FFLAGS =  -g     -u       # debug=1  (basic debug option)
    ifeq ($(debug),2)
        FFLAGS =  -C -g        # debug=2  (higher-level debug option)
    endif
#                 can add more options for  debug=3 , etc., as desired
  endif
#
# as usual, list the objects
#
OBJECTS = 1-s2.0-S0022407316300966-mmc2.o

fit: $(OBJECTS)
	$(FC) $(FFLAGS) $(OBJECTS) -o betaf.x

#--------------------------------------------------------------
# Shell Form from Sean Mcleod,  13 February 2008
#--------------------------------------------------------------
# To run the code, execute:
# ./betaf.x < input.5 > out.6

# This will create out.6 and fort.7. 

# If you are using the 1-s2.0-S0022407316300966-mmc2.txt from LeRoy's paper as your .f file then MAKE SURE TO COMMENT the two lines starting with = signs or the code will likely not compile!!!

# Ther's also an error in mm3.txt from LeRoy's paper: the PE-MLR (polynomial-exponent MLR, vs spline-exponent MLR) section needs to have non-zero p-values. See the MLR folder that was created in this commit.
