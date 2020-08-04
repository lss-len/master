#!/usr/bin/python
# -*- coding: utf-8 -*-
# Author:lss
# @Date: 20-6-17

import matplotlib.pyplot as plt
import numpy as np
import scipy.interpolate as itp
from matplotlib.backends.backend_pdf import PdfPages
import glob,os,sys
from astropy.io import fits

def read_rebintxt(name):
    line = open(name).readline()
    jd = line.split()[1]
    syserr = line.split()[2]
    data = np.loadtxt(name)
    wave = data[:,0]
    flux = data[:,1]
    err = data[:,2]
    return jd,wave,flux,err,syserr


def readlist(listname):
    """ read fits file in the list """

    lst = [i.split()[0] for i in file(listname)]

    jd_tot = []
    wave_tot = []
    flux_tot = []
    err_tot = []
    syserr_tot = []

    for i in xrange(len(lst)):

        jd, wave, flux, err, syserr = read_rebintxt(lst[i])
        jd_tot.append(float(jd))
        wave_tot.append(wave)
        flux_tot.append(flux)
        err_tot.append(err)
        syserr_tot.append(float(syserr))

    jd_tot = np.array(jd_tot)
    wave_tot = np.array(wave_tot)
    flux_tot = np.array(flux_tot)
    err_tot = np.array(err_tot)
    syserr_tot = np.array(syserr_tot)

    return jd_tot, wave_tot, flux_tot, err_tot, syserr_tot, lst


def determine_systematic_uncertainty(xc, yc, eyc):
    """ determine systematic uncertainty of light curve """

    nfilter = 5
    critec = 3.0
    niter = 50

    ic = np.arange(len(xc))

    # print 'mask points:'
    tmp_xc = xc * 1.0
    tmp_yc = yc * 1.0
    tmp_eyc = eyc * 1.0
    tmp_ic = ic * 1
    # print 'num:', len(tmp_xc)
    tmp_num = len(tmp_xc)
    for i in xrange(niter):
        tmp_yc0 = median_filter(tmp_yc, nfilter)
        tmp_dyc = tmp_yc - tmp_yc0
        tmp_mean = tmp_dyc.mean()
        tmp_std = tmp_dyc.std()
        index = np.where(np.abs(tmp_dyc - tmp_mean) <= critec * tmp_std)
        tmp_xc = tmp_xc[index[0]]
        tmp_yc = tmp_yc[index[0]]
        tmp_eyc = tmp_eyc[index[0]]
        tmp_ic = tmp_ic[index[0]]
        # print 'num:', len(tmp_xc)
        if tmp_num != len(tmp_xc):
            tmp_num = len(tmp_xc)
        else:
            break
    tmp_yc0 = median_filter(tmp_yc, nfilter)
    tmp_dyc = tmp_yc - tmp_yc0
    meanc = tmp_dyc.mean()
    stdc = tmp_dyc.std()
    # print 'sys err:', stdc, 'fraction sys err:', stdc / tmp_yc.mean()

    # === plot result ===
    xlim1, xlim2 = plotrange(xc)
    # fig = plt.figure(figsize = (8, 8))

    # ax1 = fig.add_subplot(211)
    # ax1.errorbar(outindex(xc, tmp_ic), outindex(yc, tmp_ic)
    #        , yerr = outindex(eyc, tmp_ic), fmt = '.', color = 'r')
    # ax1.errorbar(tmp_xc, tmp_yc, yerr = tmp_eyc, fmt = '.', color = 'b')
    # ax1.plot(tmp_xc, tmp_yc0, 'g')
    # ax1.plot(tmp_xc, tmp_yc0 - stdc, 'g--')
    # ax1.plot(tmp_xc, tmp_yc0 + stdc, 'g--')
    # ax1.set_xlim(xlim1, xlim2)
    # [i.set_visible(False) for i in ax1.get_xticklabels()]
    # ax1.set_ylabel('$F$', size = 14)

    # ax2 = fig.add_subplot(212)
    # ax2.errorbar(tmp_xc, tmp_dyc, yerr = tmp_eyc, fmt = '.', color = 'b')
    # f = interpolate.interp1d(tmp_xc, tmp_yc0)
    # ax2.errorbar(outindex(xc, tmp_ic), outindex(yc, tmp_ic) - f(outindex(xc, tmp_ic))
    #        , yerr = outindex(eyc, tmp_ic), fmt = '.', color = 'r')
    # ax2.plot(tmp_xc, meanc + np.zeros(len(tmp_xc)), 'g')
    # ax2.plot(tmp_xc, meanc - stdc + np.zeros(len(tmp_xc)), 'g--')
    # ax2.plot(tmp_xc, meanc + stdc + np.zeros(len(tmp_xc)), 'g--')
    # ax2.set_xlim(xlim1, xlim2)
    # [i.set_visible(False) for i in ax2.get_xticklabels()]
    # ax2.set_ylabel('$F - F_{med}$', size = 14)

    # plt.show()

    # new_err = []
    # for i in xrange(len(eyc)):
    #    if eyc[i] > stdc:
    #        new_err.append(eyc[i])
    #    else:
    #        new_err.append(stdc)
    # new_err = np.array(new_err)

    if stdc ** 2 > np.sum(eyc ** 2) / float(len(eyc)):
        scatter = (stdc ** 2 - np.sum(eyc ** 2) / float(len(eyc))) ** 0.5
    else:
        scatter = 0.0
    print 'intrinsic scatter:', scatter
    print eyc
    new_err = []
    for i in xrange(len(eyc)):
        new_err.append((scatter ** 2 + eyc[i] ** 2) ** 0.5)
    new_err = np.array(new_err)

    return new_err

def plotrange(x):
    x1 = np.min(x)
    x2 = np.max(x)
    dx = x2 - x1
    x1 = x1 - 0.1 * dx
    x2 = x2 + 0.1 * dx
    return x1, x2

def plotrangeerr(x, ex):
    x1 = np.min(x - ex)
    x2 = np.max(x + ex)
    dx = x2 - x1
    x1 = x1 - 0.1 * dx
    x2 = x2 + 0.1 * dx
    return x1, x2


def plotrangeerr2(x, ex1, ex2):
    x1 = np.min(x - ex1)
    x2 = np.max(x + ex2)
    dx = x2 - x1
    x1 = x1 - 0.05 * dx
    x2 = x2 + 0.05 * dx
    return x1, x2


def median_filter(x, n):
    """
    n = 3, 5, 7, 9...
    """
    nedge = (n - 1) / 2
    x1 = np.ones(nedge) * x[0]
    x2 = np.ones(nedge) * x[-1]
    newx = np.append(x1, x)
    newx = np.append(newx, x2)
    median = []
    for i in xrange(len(x)):
        # print newx[i: i + n]
        y = np.median(newx[i: i + n])
        # print y
        median.append(y)
    median = np.array(median)
    return median

def inte_spec(jd_tot, wave_tot, flux_tot, err_tot, syserr_tot,contil1, contil2, contir1, contir2,hb1,hb2,lst,conti51001,conti51002):
    conti_right = []
    conti_right_err = []
    conti_left = []
    conti_left_err = []
    line = []
    line_err = []

    conti5100 = []
    conti5100_err = []

    fig = plt.figure()
    plt.ion()
    pdf = PdfPages('plt_windows_spec.pdf')
    for i in range(len(wave_tot)):
        wave = wave_tot[i]
        flux = flux_tot[i]
        err = err_tot[i]

        # === left continuum window ===
        index = np.where((wave >= contil1) & (wave < contil2))
        flux_con_left = np.median(flux[index[0]])
        wave_con_left = np.mean(wave[index[0]])
        conti_left.append(flux_con_left)
        conti_left_err.append(np.std(flux[index[0]]) / (float(len(index[0]))) ** 0.5)

        # === right continuum window ===
        index = np.where((wave >= contir1) & (wave < contir2))
        flux_con_right = np.median(flux[index[0]])
        wave_con_right = np.mean(wave[index[0]])
        conti_right.append(flux_con_right)
        conti_right_err.append(np.std(flux[index[0]]) / (float(len(index[0]))) ** 0.5)

        # === f5100 continuum  ===
        index = np.where((wave >= conti51001) & (wave < conti51002))
        flux5100 = np.median(flux[index[0]])
        wave5100 = np.mean(wave[index[0]])
        conti5100.append(flux5100)
        conti5100_err.append(np.std(flux[index[0]]) / (float(len(index[0]))) ** 0.5)


        # === create interpolation ===
        fcon = flux_con_left + (flux_con_right - flux_con_left) / \
               (wave_con_right - wave_con_left) * (wave - wave_con_left)

        flux = flux - fcon

        # === hbeta flux =======

        index = np.where((wave>=hb1) &(wave<=hb2))
        line.append(np.sum(flux[index[0]])*(wave[1]-wave[0]))
        line_err.append(np.sum((err[index[0]])**2)**0.5*(wave[1]-wave[0]))

        plt.errorbar(wave_tot[i], flux_tot[i], yerr=err_tot[i], ls='none', marker='o', ms=1, linewidth=0.4)
        arg = np.where((wave_tot[i] > 5075) & (wave_tot[i] < 5125))
        sn1 = np.mean(flux_tot[i][arg]) / np.mean(err_tot[i][arg])

        plt.axvline(4861, linestyle='--', color='r', linewidth=0.8)
        plt.axvline(contil1, linestyle='--', color='C2', linewidth=0.8)
        plt.axvline(contil2, linestyle='--', color='C2', linewidth=0.8)
        plt.axvline(contir1, linestyle='--', color='C2', linewidth=0.8)
        plt.axvline(contir2, linestyle='--', color='C2', linewidth=0.8)

        plt.axvline(hb1, linestyle='--', color='C3', linewidth=0.8)
        plt.axvline(hb2, linestyle='--', color='C3', linewidth=0.8)
        plt.plot([wave_con_left,wave_con_right],[flux_con_left,flux_con_right],linestyle='--',color='C4')
        plt.scatter(wave_con_left,flux_con_left,color='C5')
        plt.scatter(wave_con_right,flux_con_right,color = 'C5')
        plt.xlim(4500, 5200)
        plt.title(lst[i] + ' ' + str(sn1))

        pdf.savefig()
        plt.clf()
    pdf.close()

    (conti5100,conti5100_err,line,line_err) = np.array((conti5100,conti5100_err,line,line_err))

    conti_err_plussyserr = []
    for i in xrange(len(conti5100_err)):
        conti_err_plussyserr.append((conti5100_err[i] ** 2 + (syserr_tot[i] * conti5100[i]) ** 2) ** 0.5)
    conti_err_plussyserr = np.array(conti_err_plussyserr)

    conti_new_err = determine_systematic_uncertainty(jd_tot, conti5100, conti_err_plussyserr)
    hb_new_err = determine_systematic_uncertainty(jd_tot, line,line_err)

    f1 = open('lc.txt','w')
    f2 = open('lc-syserr.txt','w')
    for j in range(len(jd_tot)):
        text1 = '%s %s %s %s %s\n'%(jd_tot[j],line[j],line_err[j],conti5100[j],conti_err_plussyserr[j])
        text2 = '%s %s %s %s %s\n'%(jd_tot[j],line[j],hb_new_err[j],conti5100[j],conti_new_err[j])
        f1.write(text1)
        f2.write(text2)
    f1.close()
    f2.close()



def main():
    curdir = os.getcwd()
    paradir = os.path.split(curdir)[0]
    parameters = open(paradir + os.sep + 'sumline.para').readlines()
    print paradir
    z = float(parameters[1].split()[0])
    contil1 = float(parameters[2].split()[0])
    contil2 = float(parameters[2].split()[1])
    contir1 = float(parameters[3].split()[0])
    contir2 = float(parameters[3].split()[1])
    hb1 = float(parameters[4].split()[0])
    hb2 = float(parameters[4].split()[1])

    conti51001 = 5075  # left edge of right continuum window
    conti51002 = 5125  # right edge of right continuum window

    listname = sys.argv[1]
    jd_tot, wave_tot, flux_tot, err_tot, syserr_tot, lst = readlist(listname)
    wave_tot =wave_tot/(1+z)
    ### plot rebin spectrum ##
    fig = plt.figure()
    plt.ion()
    pdf = PdfPages('rebin_spec.pdf')
    for i in range(len(lst)):
        # fig = plt.figure()
        plt.errorbar(wave_tot[i],flux_tot[i],yerr=err_tot[i],ls='none', marker='o',ms=1,linewidth=0.4)
        arg = np.where((wave_tot[i] > 5075) & (wave_tot[i]< 5125))
        sn1 = np.mean(flux_tot[i][arg]) / np.mean(err_tot[i][arg])

        plt.axvline(4861,linestyle='--', color='r', linewidth=0.8)
        plt.axvline(contil1, linestyle='--', color='C2', linewidth=0.8)
        plt.axvline(contil2, linestyle='--', color='C2', linewidth=0.8)
        plt.axvline(contir1, linestyle='--', color='C2', linewidth=0.8)
        plt.axvline(contir2, linestyle='--', color='C2', linewidth=0.8)

        plt.axvline(hb1, linestyle='--', color='C3', linewidth=0.8)
        plt.axvline(hb2, linestyle='--', color='C3', linewidth=0.8)

        plt.axvline(5878/(1+z), linestyle='--', color='C4', linewidth=0.8)
        plt.axvline(5908/(1+z), linestyle='--', color='C4', linewidth=0.8)

        plt.xlim(4500,5200)


        plt.title(lst[i]+' '+str(sn1))
        # plt.title(lst[i]+' '+str(sn1))
        pdf.savefig()
        plt.clf()
    pdf.close()

    ###



    #
    inte_spec(jd_tot, wave_tot, flux_tot, err_tot, syserr_tot,contil1, contil2, contir1, contir2,hb1,hb2,lst,conti51001,conti51002)




if __name__ == '__main__':
    main()

