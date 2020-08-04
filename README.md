# JAABAvAnnotation
This repository contains code that analyzes the variability in human annotation of fly behavior from video recordings, and those that characterizes the correlation between human annotation confidence and confidence output from JAABA, a supervised machine learning program for automated animal behavior analysis. 

The code in this repository was used to conduct analysis and draw figures in the paper

Xubo Leng, Margot Wohl, Kenichi Ishii, Pavan Nayak, Kenta Asahina (2020). Quantitative comparison of Drosophila behavior annotations by human observers and a machine learning algorithm. Submitted.  

Details on functions of the code are described in Materials and Methods section of the paper.

## System requirements
The code is developed in MATLAB R2018a and tested on that version. It should work with higher versions of MATLAB. The operating system used in development is Windows 10. 

## Software dependencies
To perform analysis on existing JAABA outputs from one JAABA classifier, the code does not depend on other software; in cases where re-running JAABA analysis is necessary (for example, when testing multiple JAABA classifiers), please install [JAABA](https://github.com/kristinbranson/JAABA). 

## Installation
To install, simply clone this repository to your computer. 

## Instruction for use
The primary program entrypoint is run_analyze_human_jaaba_annot_corr_all.m. Given the JAABA classifiers and the corresponding detection thresholds, the program first runs JAABA to produce JAABA outputs, then uses the OrgData function to produce a FLYMAT MAT-file summarizing JAABA predictions. Based on FLYMAT files, the program finally calls the analyze_human_jaaba_annot_corr_v3 function to produce bout_matches MAT-files and to analyze and visualize correlation between annotation and prediction confidence values (shown as main figures of the paper). 

If JAABA, desired classifiers and the JAABA output folder structure are set up correctly, a user should be able to run run_analyze_human_jaaba_annot_corr_all directly. Otherwise, with the FLYMAT files, a user should be able to run analyze_human_jaaba_annot_corr_v3 separately to reproduce the same results.

## Contributions
Code in this repository is mostly written by Xubo Leng. It also contains several scripts modified from original work by Margot Wohl: specifically, they are MANvsJAABAduration_xubo.m, OrgData020816_XuboCopy.m and smoothing.m. 
