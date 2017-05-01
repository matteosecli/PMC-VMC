TEMPLATE = app
CONFIG  += console c++11
CONFIG  -= app_bundle
CONFIG  -= qt

QMAKE_CXXFLAGS += -std=c++11

# Not sure if this is needed. Should give extra speed to Armadillo even for simple operations.
QMAKE_CXXFLAGS += -fopenmp
QMAKE_LFLAGS   += -fopenmp

# Manually include a more recent version of Armadillo with 'regspace'.
#QMAKE_CXXFLAGS += -I ../VariationalMC/armadillo/usr/local/include

# It's not a nice thing, it's better to properly set the environment.
#INCLUDEPATH += /opt/intel/compilers_and_libraries/linux/lib/intel64

# Not really needed if you include Armadillo as above. I think you can also skip LAPACK and BLAS.
LIBS += -larmadillo -llapack -lblas

SOURCES += main.cpp \
    lib.cpp

HEADERS += \
    lib.h

