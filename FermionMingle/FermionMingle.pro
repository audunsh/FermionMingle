TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    basis/basis.cpp \
    basis/basisbank.cpp \
    basis/contracted.cpp \
    basis/primitive.cpp \
    integrator/boysfunction.cpp \
    integrator/integrator.cpp \
    interface/fmingle.cpp \
    solvers/ccsolve.cpp \
    solvers/hfsolve.cpp \
    solvers/rhfsolve.cpp \
    solvers/uhfsolve.cpp

include(deployment.pri)
qtcAddDeployment()

HEADERS += \
    basis/basis.h \
    basis/basisbank.h \
    basis/contracted.h \
    basis/primitive.h \
    integrator/boysfunction.h \
    integrator/integrator.h \
    interface/fmingle.h \
    solvers/ccsolve.h \
    solvers/hfsolve.h \
    solvers/rhfsolve.h \
    solvers/uhfsolve.h


LIBS += -larmadillo -lblas -llapack
