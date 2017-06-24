## Process this file with automake to produce Makefile.in

AM_CFLAGS = -ansi -Wall -pedantic
lib_LTLIBRARIES = libnwalign.la
libnwalign_la_SOURCES = nwalign.c nwalign.h blosum.h

bin_PROGRAMS = nwalign
nwalign_SOURCES = main.c
nwalign_LDADD = libnwalign.la

