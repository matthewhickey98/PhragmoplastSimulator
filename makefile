MAGICK_DIR = opt/homebrew/Cellar/graphicsmagick/1.3.36/

CXXFLAG = -Wall -std=c++17 -mmacosx-version-min=10.15 

INCLUDES += -I/usr/local/include 
INCLUDES += -I/$(MAGICK_DIR)/include/GraphicsMagick/
INCLUDES += -I/Applications/MATLAB_R2021a.app/extern/include



L += -L/$(LDIR) 

L += -L/$(MAGICK_DIR)/lib 

L += -lGraphicsMagick++
L += -lGraphicsMagick
L += -lGraphicsMagickWand

CXX = g++
FILES = Bitmap readpars RandomGen FrameName Microtubule Phragmoplast 
OBJ = Microtubule.o Phragmoplast.o Bitmap.o readpars.o RandomGen.o FrameName.o
APPS = php_simulator_3s php_mkfig_3s

.cpp:
	@$(CXX) -c $< -o $*.o $(CFXXLAG) $(INCLUDES)

.cpp.o: 
	@$(CXX) -c $< -o $*.o $(CXXFLAG) $(INCLUDES)

all: $(APPS)

php_simulator_3s: php_simulator_3s.o $(OBJ)
	@$(CXX) -o $@ $(CXXFLAG) $(INCLUDES) $@.o $(OBJ) $(L)

php_mkfig_3s: php_mkfig_3s.o $(OBJ)
	@$(CXX) -o $@ $(CXXFLAG) $(INCLUDES) $@.o $(OBJ) $(L)

clean:
	/bin/rm -f *.o core.* $(APPS)

	

