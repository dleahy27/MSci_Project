#! /usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import lhapdf
import patchpdfplotter as ppdf

# for LaTeX font on plots
# first need to install latex packages
# [in terminal] sudo apt install -q cm-super dvipng texlive-latex-extra texlive-latex-recommended
plt.rcParams['font.size'] = 15
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['text.usetex'] = True

# get PDF member object from set data file
# need to lhapdf install file first
pdf = lhapdf.mkPDF("NNPDF40_nlo_pch_as_01180")

# make 1D plot of PDFs of all flavours
ppdf.plotPatchPDF1DAllFlavours(pdf, pdf.xfxQ2, num=100, axis=0, index_q=10, fname='PPEplots/patchpdfplots/patchpdf_1D_lowQ', axes=['$x$', '$Q^2$', '$xf(x, Q^2)$'])
ppdf.plotPatchPDF1DAllFlavours(pdf, pdf.xfxQ2, num=100, axis=0, index_q=40, fname='PPEplots/patchpdfplots/patchpdf_1D_highQ', axes=['$x$', '$Q^2$', '$xf(x, Q^2)$'])

# plot PDFs
ppdf.plotPatchPDF(21, pdf.xfxQ2, xstart=-3, xend=0, qstart=0, qend=8, num=50, fname='PPEplots/patchpdfplots/patchpdf_g', axes=['$x$', '$Q^2$', '$xf(x, Q^2)$'])
ppdf.plotPatchPDF(2, pdf.xfxQ2, xstart=-3, xend=0, qstart=0, qend=8, num=50, fname='PPEplots/patchpdfplots/patchpdf_u', axes=['$x$', '$Q^2$', '$xf(x, Q^2)$'])

# plot PDF derivatives
ppdf.plotPatchPDFDerivative(21, pdf.xfxQ2, order=1, num=100, fname='PPEplots/patchpdfplots/patchpdf_g', axes=['$x$', '$Q^2$', '$xf(x, Q^2)$'], plot1d=True, index_x=10, index_y=40)
ppdf.plotPatchPDFDerivative(21, pdf.xfxQ2, order=2, num=100, fname='PPEplots/patchpdfplots/patchpdf_g', axes=['$x$', '$Q^2$', '$xf(x, Q^2)$'], plot1d=True, index_x=10, index_y=40)
ppdf.plotPatchPDFDerivative(2, pdf.xfxQ2, order=1, num=100, fname='PPEplots/patchpdfplots/patchpdf_u', axes=['$x$', '$Q^2$', '$xf(x, Q^2)$'], plot1d=True, index_x=30, index_y=25)
ppdf.plotPatchPDFDerivative(2, pdf.xfxQ2, order=2, num=100, fname='PPEplots/patchpdfplots/patchpdf_u', axes=['$x$', '$Q^2$', '$xf(x, Q^2)$'], plot1d=True, index_x=30, index_y=25)