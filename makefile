CXXFLAG = -Wall -g  -arch x86_64 #-D_THREAD_SAFE
LDIR = usr/local/lib/
#IDIR = $(LDIR)/include
INCLUDES += -I/usr/local/include \
-I/opt/homebrew/Cellar/graphicsmagick/1.3.36/include/GraphicsMagick

#-I/usr/X11R6/include/X11/magick/ -I/usr/X11R6/lib 
#-I/opt/apps/imagemagick/6.9.11-23/include/ImageMagick-6 -I/opt/apps/imagemagick/6.9.11-23/lib  -I/opt/apps/imagemagick/6.9.11-23/
L += -L/$(LDIR) -L/opt/homebrew/Cellar/graphicsmagick/1.3.36/lib -lGraphicsMagick++ -lGraphicsMagick -lGraphicsMagickWand  #-lfreetype -lbz2 -lz  -lm -lpthread -llcms2 #-lboost_python27 -lboost_python27-mt#-lltdl
CXX = g++

OBJ = Phragmoplast.o Microtuble.o Bitmap.o readpars.o RandomGen.o FrameName.o  #t_except.o
#str_buff.o ics.o
APPS = php_simulator_3s php_mkfig_3s

.cpp:
	$(CXX) -c $< -o $*.o $(CFXXLAG) $(INCLUDES) $(L)
.cpp.o: 
	$(CXX) -c $< -o $*.o $(CXXFLAG) $(INCLUDES) $(L)

all: $(APPS)


php_simulator_3s: php_simulator_3s.o $(OBJ)
	$(CXX) -o $@ $(CXXFLAG) $(INCLUDES) $@.o $(OBJ) $(L)

php_mkfig_3s: php_mkfig_3s.o $(OBJ)
	$(CXX) -o $@ $(CXXFLAG) $(INCLUDES) $@.o $(OBJ) $(L)

clean:
	/bin/rm -f *.o core.* $(APPS)

