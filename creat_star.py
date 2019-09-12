#!/usr/bin/python

import sys, math,os
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits as pyfits



def readcagn(fitsname):
    fits = pyfits.open(fitsname)
    jd = fits[0].header['JD']
    wave0 = fits[0].header['CRVAL1']
    dwave = fits[0].header['CD1_1']
    nwave = fits[0].header['NAXIS1']
    wave = wave0 + dwave * np.arange(nwave)
    data = fits[0].data
    flux = data[0]
    err = data[3]
    fits.close()
    return jd, wave, flux, err

def readfits(fn):
    """fits file to txt"""
    fit = pyfits.open(fn)
    beg = fit[0].header['CRVAL1']
    step = fit[0].header['CD1_1']
    size = fit[0].header['NAXIS1']
    jd = fit[0].header['JD']
    wave = np.arange(size) * step + beg
    flux = fit[0].data[0]
    if len(fit[0].data) ==4:
        err = fit[0].data[3]
    else:
        err = fit[0].data[1]
    return jd,wave,flux,err

def continuum_subtraction(wave, flux,err, contil1, contil2, contir1, contir2,conti51001,conti51002):

    """ continuum subtraction """
    # === left continuum window ===
    index = np.where((wave >= contil1) & (wave < contil2))
    flux_con_left = np.median(flux[index[0]])
    wave_con_left = np.mean(wave[index[0]])

    # === right continuum window ===
    index = np.where((wave >= contir1) & (wave < contir2))
    flux_con_right = np.median(flux[index[0]])
    wave_con_right = np.mean(wave[index[0]])
    # conti.append(flux_con_right)
    # conti_err.append(np.std(flux[index[0]]) / (float(len(index[0]))) ** 0.5)
    # fig = plt.figure()
    # ax = fig.add_subplot(211)
    # index =np.where((wave>4700)&(wave<5500))
    #
    # ax.errorbar(wave[index[0]], flux[index[0]], yerr=err[index[0]], fmt='-', capsize=0.0, color='b')
    # ax.plot(np.array([wave_con_left, wave_con_right]), np.array([flux_con_left, flux_con_right]),
    #         'ro')
    # ax.plot(np.array([wave_con_left, wave_con_right]), np.array([flux_con_left, flux_con_right]),
    #         'r-')
    # ax.axvline(contil1, linestyle='--', color='g')
    # ax.axvline(contil2, linestyle='--', color='g')
    # ax.axvline(contir1, linestyle='--', color='g')
    # ax.axvline(contir2, linestyle='--', color='g')
    # ax.axvline(conti51002, linestyle='--', color='g')
    # ax.axvline(conti51001, linestyle='--', color='g')
    #
    # ax.axvline(hb2, linestyle='--', color='r')
    # ax.axvline(hb1, linestyle='--', color='r')
    # fcon = flux_con_left + (flux_con_right - flux_con_left) / \
    #        (wave_con_right - wave_con_left) * (wave - wave_con_left)
    # hbeta = flux - fcon
    # ax = fig.add_subplot(212)
    # ax.errorbar(wave[index[0]], hbeta[index[0]], yerr=err[index[0]], fmt='-', capsize=0.0, color='b')
    # ax.axvline(contil1, linestyle='--', color='g')
    # ax.axvline(contil2, linestyle='--', color='g')
    # ax.axvline(contir1, linestyle='--', color='g')
    # ax.axvline(contir2, linestyle='--', color='g')
    # ax.axvline(hb1, linestyle='--', color='r')
    # ax.axvline(hb2, linestyle='--', color='r')
    # ax.axhline(0.0, color='r')
    # plt.show()

    # === create interpolation ===
    fcon = flux_con_left + (flux_con_right - flux_con_left) / \
           (wave_con_right - wave_con_left) * (wave - wave_con_left)
    hbeta = flux - fcon

    index = np.where((wave>=conti51001)&(wave<conti51002))
    f5100=np.median(flux[index[0]])
    f5100_err = (np.std(flux[index[0]])/(float(len(index[0])))**0.5)
    return hbeta,f5100,f5100_err


def main():
    lstname = 'corlist'  # list of input
    parameters = open('sumline.para').readlines()
    z = float(parameters[1].split()[0])
    con_l_1 = float(parameters[2].split()[0])*(1+z)
    con_l_2 = float(parameters[2].split()[1])*(1+z)
    con_r_1 = float(parameters[3].split()[0])*(1+z)
    con_r_2 = float(parameters[3].split()[1])*(1+z)
    hb1 = float(parameters[4].split()[0])*(1+z)
    hb2 = float(parameters[4].split()[1])*(1+z)
    conti51001 = 5075*(1+z)  # left edge of right continuum window
    conti51002 = 5125*(1+z)

    lst = [i.split('\n')[0] for i in file(lstname)]
    clst = ['c' + i for i in lst]
    # cclstname ='ccclist'
    # cclst = [i.split()[0] for i in file(cclstname)]
    print hb1,hb2

    for k in xrange(len(clst)):
        calib = 'fits_out/' + clst[k]
        jd_calib, wave_calib, flux_calib, err_calib = readcagn(calib)
        # calib = 'fits_com/' + cclst[k]
        # jd_calib, wave_calib, flux_calib, err_calib = readfits(calib)
        hbeta, f5100, f5100_err=continuum_subtraction(wave_calib,flux_calib,err_calib,con_l_1,con_l_2,con_r_1,con_r_2,conti51001,conti51002)

        output = open(calib.split('/')[-1] + '.agn', 'w')
        text = '# JD %.8f  f5100 %.6e  f5100_err %.6e\n' % (jd_calib,f5100, f5100_err)
        output.write(text)
        index = np.where((wave_calib>= hb1)&(wave_calib<=hb2))

        for j in xrange(len(index[0])):
            text = '%.6f  %.6e  %.6e\n' % (wave_calib[index[0]][j], hbeta[index[0]][j], err_calib[index[0]][j])
            # text = '%.6f  %.6e  %.6e\n' % (wave_calib[index[0]][j], flux_calib[index[0]][j], err_calib[index[0]][j])
            output.write(text)
        output.close()
if __name__ == "__main__":
    main()
