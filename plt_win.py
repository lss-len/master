#!/usr/bin/python

import os,glob,shutil
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits as pyfits
from matplotlib.backends.backend_pdf import PdfPages

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

def continuum_subtraction(wave, flux):

    """ continuum subtraction """
    # === left continuum window ===
    index = np.where((wave >= contil1) & (wave < contil2))
    flux_con_left = np.median(flux[index[0]])
    wave_con_left = np.mean(wave[index[0]])

    # === right continuum window ===
    index = np.where((wave >= contir1) & (wave < contir2))
    flux_con_right = np.median(flux[index[0]])
    wave_con_right = np.mean(wave[index[0]])
    # === create interpolation ===
    fcon = flux_con_left + (flux_con_right - flux_con_left) / \
           (wave_con_right - wave_con_left) * (wave - wave_con_left)
    hbeta = flux - fcon

    index = np.where((wave>=conti51001)&(wave<conti51002))
    f5100=np.median(flux[index[0]])
    f5100_err = (np.std(flux[index[0]])/(float(len(index[0])))**0.5)
    return hbeta,f5100,f5100_err

def plt_win(dir,lst):
    pdf = PdfPages(dir+'.pdf')
    for k in range(len(lst)):
        calib = dir+os.sep+lst[k]
        jd,wave,flux,err=readfits(calib)
        """ continuum subtraction """
        # === left continuum window ===
        index = np.where((wave >= contil1) & (wave < contil2))
        flux_con_left = np.median(flux[index[0]])
        wave_con_left = np.mean(wave[index[0]])
        # === right continuum window ===
        index = np.where((wave >= contir1) & (wave < contir2))
        flux_con_right = np.median(flux[index[0]])
        wave_con_right = np.mean(wave[index[0]])

        fig = plt.figure()
        ax = fig.add_subplot(211)
        index =np.where((wave>4700)&(wave<5500))
        ax.errorbar(wave[index[0]], flux[index[0]], yerr=err[index[0]],linewidth=0.9, color='C1')
        ax.plot(np.array([wave_con_left, wave_con_right]), np.array([flux_con_left, flux_con_right]),
                'ro')
        ax.plot(np.array([wave_con_left, wave_con_right]), np.array([flux_con_left, flux_con_right]),
                'r-')
        ax.axvline(contil1, linestyle='--', color='g')
        ax.axvline(contil2, linestyle='--', color='g')
        ax.axvline(contir1, linestyle='--', color='g')
        ax.axvline(contir2, linestyle='--', color='g')
        ax.axvline(conti51002, linestyle='--', color='g')
        ax.axvline(conti51001, linestyle='--', color='g')
        ax.axvline(hb2, linestyle='--', color='r')
        ax.axvline(hb1, linestyle='--', color='r')
        ax.set_title(lst[k])
        fcon = flux_con_left + (flux_con_right - flux_con_left) / \
               (wave_con_right - wave_con_left) * (wave - wave_con_left)
        hbeta = flux - fcon
        ax = fig.add_subplot(212)
        ax.errorbar(wave[index[0]], hbeta[index[0]], yerr=err[index[0]],linewidth=0.9 , color='C2')
        ax.axvline(contil1, linestyle='--', color='g')
        ax.axvline(contil2, linestyle='--', color='g')
        ax.axvline(contir1, linestyle='--', color='g')
        ax.axvline(contir2, linestyle='--', color='g')
        ax.axvline(hb1, linestyle='--', color='r')
        ax.axvline(hb2, linestyle='--', color='r')
        ax.axhline(0.0, color='r')
        pdf.savefig()
        plt.clf()
    pdf.close()

parameters = open('sumline.para').readlines()
z = float(parameters[1].split()[0])
contil1 = float(parameters[2].split()[0])*(1+z)
contil2 = float(parameters[2].split()[1])*(1+z)
contir1 = float(parameters[3].split()[0])*(1+z)
contir2 = float(parameters[3].split()[1])*(1+z)
hb1 = float(parameters[4].split()[0])*(1+z)
hb2 = float(parameters[4].split()[1])*(1+z)
conti51001 = 5075*(1+z)  # left edge of right continuum window
conti51002 = 5125*(1+z)

def creat_txt(dir,lst):
    for k in range(len(lst)):
        calib = dir+os.sep+lst[k]

        jd,wave,flux,err=readfits(calib)
        hbeta,f5100,f5100_err = continuum_subtraction(wave,flux)
        output = open(calib.split('/')[-1]+'.agn','w')
        text = '# JD %.8f f5100 %.6e f5100_err %.6e\n'%(jd,f5100,f5100_err)

        output.write(text)
        index = np.where((wave>= hb1)&(wave<=hb2))
        # print wave[720],flux[720],hbeta[720]
        for j in range(len(index[0])):
            text = '%.6f  %.6e  %.6e\n' % (wave[index[0]][j],hbeta[index[0]][j],err[index[0]][j])
            output.write(text)
        output.close()
    print wave[700],hbeta[700]
def move(dir,lst):
    if not os.path.isdir(dir):
        print dir
        os.mkdir(dir)
    for i in lst:
        shutil.move(i,dir+os.sep+i)



def main():
    lstname = 'corlist'  # list of input
    lst = [i.split()[0] for i in file(lstname)]
    clst = ['c' + i for i in lst]
    cclstname ='ccclist'
    cclst = [i.split()[0] for i in file(cclstname)]

    plt_win('fits_out',clst)
    plt_win('fits_com',cclst)
    creat_txt('fits_out',clst)
    creat_txt('fits_com',cclst)
    caw = glob.glob('caw*.agn')
    move('spec_caw',caw)
    ccaw = glob.glob('ccaw*.agn')
    move('spec_ccaw',ccaw)


if __name__ == "__main__":
    main()
