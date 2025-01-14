{\rtf1\ansi\ansicpg1252\cocoartf2706
\cocoatextscaling0\cocoaplatform0{\fonttbl\f0\froman\fcharset0 TimesNewRomanPS-BoldMT;\f1\froman\fcharset0 TimesNewRomanPSMT;}
{\colortbl;\red255\green255\blue255;\red0\green0\blue0;}
{\*\expandedcolortbl;;\cssrgb\c0\c0\c0;}
\margl1440\margr1440\vieww11520\viewh8400\viewkind0
\deftab720
\pard\pardeftab720\partightenfactor0

\f0\b\fs48 \cf0 \expnd0\expndtw0\kerning0
\outl0\strokewidth0 \strokec2 Quantifying the influence of biophysical factors in shaping brain communication through remnant functional networks
\fs37\fsmilli18667 \
\
\

\f1\b0\fs28 This repository contains at the code and data created to the reproduce the main figures in the paper \'93\outl0\strokewidth0 Quantifying the influence of biophysical factors in shaping brain communication through remnant functional networks\'94\outl0\strokewidth0 \strokec2  published in Network Neuroscience. \
\
\
Analysis_Figures_X: are scripts for running and generating the figures associated with the main text. \
\
\pard\pardeftab720\partightenfactor0
\cf0 \outl0\strokewidth0 Data.m: contained the FC from the HCP dataset, estimated biophysical networks (sc, ec, gc, and rc) and data from the LA5c dataset. \
\pard\pardeftab720\partightenfactor0
\cf0 \outl0\strokewidth0 \strokec2 \
create_RFN.m: is the main function for creating the Remnant Functional Networks (RFN). \
\
my_delta.m: estimates the percent difference between features estimated from fully connected FC network and RFN. \
\
percent_overlap.m: calculates the percent of connections in FC that have an underlying biophysical connection. \
\
(un)wrapMat.m: (un)vectorized the brain networks. \
\
\
\
\pard\pardeftab720\partightenfactor0

\f0\b\fs36 \cf0 \outl0\strokewidth0 \

\fs37\fsmilli18667 \outl0\strokewidth0 \strokec2 \
}