#!/usr/bin/env /proj/sot/ska/bin/python

#########################################################################################
#                                                                                       #
#       sim_twist_trend_plot.py: create trend plots for sim twist etc                   #
#                                                                                       #
#               author: t. isobe (tisobe@cfa.harvard.edu)                               #
#                                                                                       #
#               last update: Jun 05, 2017                                               #
#                                                                                       #
#########################################################################################

import os
import sys
import re
import string
import random
import math
import time
import numpy
import astropy.io.fits  as pyfits
import Ska.engarchive.fetch as fetch
import Chandra.Time
#
#--- interactive plotting module
#
import mpld3
from mpld3 import plugins, utils
#
#--- pylab plotting routine related modules
#
import matplotlib as mpl

if __name__ == '__main__':
    mpl.use('Agg')
from pylab import *
import matplotlib.pyplot       as plt
import matplotlib.font_manager as font_manager
import matplotlib.lines        as lines

#
#--- reading directory list
#
path = '/data/mta/Script/ALIGNMENT/Sim_twist/house_keeping/dir_list_py'

f= open(path, 'r')
data = [line.strip() for line in f.readlines()]
f.close()

for ent in data:
    atemp = re.split(':', ent)
    var  = atemp[1].strip()
    line = atemp[0].strip()
    exec "%s = %s" %(var, line)
#
#--- append  pathes to private folders to a python directory
#
sys.path.append(bin_dir)
sys.path.append(mta_dir)
#
#--- import several functions
#
import convertTimeFormat          as tcnv       #---- contains MTA time conversion routines
import mta_common_functions       as mcf        #---- contains other functions commonly used in MTA scripts
import robust_linear              as rfit       #---- robust fit rountine

#
#--- temp writing file name
#
rtail  = int(time.time())
zspace = '/tmp/zspace' + str(rtail)

mon_list = [0, 31, 59, 90, 120, 151, 181, 212, 234, 373, 304, 334]
#
#--- break points of the line fitting
#
b_period = [1999.0, 2003.42, 2006.94]

#------------------------------------------------------------------------------------------
#-- sim_twist_trend_plot: create trend plots for sim twist etc                           --
#------------------------------------------------------------------------------------------

def sim_twist_trend_plot(inyear):
    """
    create trend plots for sim twist etc
    input:  inyear  --- if it is given, create the plots for the year
    output: <web_dir>/Plots>/sim_plot_<year>.png
            <web_dir>/Plots>/twist_plot_<year>.png
            <web_dir>/Plots>/dtheta_plot_<year>.png
            if year is not given, it creates each year plot for all year from 1999 to present and
            the full range plots
    """

    if inyear != '':
        syear = inyear
        pyear = inyear
        eyear = inyear + 1
        lchk  = 0
    else:
        syear = 1999
        out   = time.gmtime()
        cyear = int(float(out.tm_year))
        cmon  = int(float(out.tm_mon))
        eyear = cyear + 1
        if cmon < 3:
            pyear = cyear - 1
        else:
            pyear = cyear
        lchk  = 1

    full_data_info = []
    full_data_extr = []
    full_acis_i    = []
    full_acis_s    = []
    full_hrc_i     = []
    full_hrc_s     = []
    for year in range(syear, eyear):

        data_info      = read_data_info(year)
        data_extr      = read_data_extracted(year)
        [acis_i, acis_s, hrc_i, hrc_s] = dtheta_on_inst(data_extr)
#
#--- update plots for this year and possibly the last year (for the first two months of the year)
#
        if year >= pyear:
            plot_data_info(data_info, year=year)
            plot_data_extr(data_extr, year=year)
            plot_dtheta(acis_i, acis_s, hrc_i, hrc_s, year=year)

        if lchk == 1:
            if year == syear:
                full_data_info = data_info
                full_data_extr = data_extr
                full_acis_i    = acis_i
                full_acis_s    = acis_s
                full_hrc_i     = hrc_i
                full_hrc_s     = hrc_s
            else:
                for k in range(0, len(data_info)):
                    full_data_info[k] = full_data_info[k] + data_info[k]
    
                for k in range(0, len(data_extr)):
                    full_data_extr[k] = full_data_extr[k] + data_extr[k]
    
                    full_acis_i[k]    = full_acis_i[k]    + acis_i[k]                
                    full_acis_s[k]    = full_acis_s[k]    + acis_s[k]                
                    full_hrc_i[k]     = full_hrc_i[k]     + hrc_i[k]
                    full_hrc_s[k]     = full_hrc_s[k]     + hrc_s[k]


    if lchk == 1:
        plot_data_info(full_data_info)
        plot_data_extr(full_data_extr)
        plot_dtheta(full_acis_i, full_acis_s, full_hrc_i, full_hrc_s)

#------------------------------------------------------------------------------------------
#-- read_data_info: read data info data for a given year                                 --
#------------------------------------------------------------------------------------------

def read_data_info(year):
    """
    read data info data for a given year
    input:  year    --- the year of the data set
    output: [time, sim_x, sim_y, sim_z, pitch, yaw]
    """

    infile = data_dir + 'data_info_' + str(year)
    data   = read_data(infile)

    time   = []
    sim_x  = []
    sim_y  = []
    sim_z  = []
    pitch  = []
    yaw    = []
    for ent in data:
        atemp = re.split('\s+', ent)
        start = convert_to_fyear(atemp[3])
        stop  = convert_to_fyear(atemp[4])
        mid   = 0.5 *(start + stop)
        try:
            val1 = float(atemp[5])
            val2 = float(atemp[6])
            val3 = float(atemp[7])
            val4 = float(atemp[8]) * 3600.0
            val5 = float(atemp[9]) * 3600.0

            time.append(mid)
            sim_x.append(val1)
            sim_y.append(val2)
            sim_z.append(val3)
            pitch.append(val4)
            yaw.append(val5)
        except:
            continue

    return [time, sim_x, sim_y, sim_z, pitch, yaw]

#------------------------------------------------------------------------------------------
#-- read_data_extracted: read extracted data for a given year                            --
#------------------------------------------------------------------------------------------

def read_data_extracted(year):
    """
    read extracted data for a given year
    input:  year    --- the year of the data set
    ouptput: [time, dy, dz, dtheta, inst]
    """
#
#--- sometime there are out of spot data in the file; so remove them
#
    xbot   = str(year) + ':001:00:00:00'
    xbot   =  out = Chandra.Time.DateTime(xbot).secs

    infile = data_dir + 'data_extracted_' + str(year)
    data   = read_data(infile)
    
    time   = []
    dy     = []
    dz     = []
    dtheta = []
    inst   = []
    for ent in data:
        atemp = re.split('\s+', ent)
        try:
            if float(atemp[0]) < xbot:
                continue
            tval  = convert_to_fyear(atemp[0])
            yval  = float(atemp[3])
            zval  = float(atemp[4])
            hval  = float(atemp[5]) * 3600.0

            time.append(tval)
            dy.append(yval)
            dz.append(zval)
            dtheta.append(hval)
            inst.append(atemp[2])
        except:
            continue

    return [time, dy, dz, dtheta, inst]

#------------------------------------------------------------------------------------------
#-- dtheta_on_inst: separate the dtheta data by instruments                             ---
#------------------------------------------------------------------------------------------

def dtheta_on_inst(data):
    """
    separate the dtheta data by instruments
    input:  data    --- data: [time, data1, data2,...]
    output: [acis_i, acis_s, hrc_i, hrc_s] each has the same stracture as data but only
            for that instrument
    """
    
    acis_i = [[], [], [], [], []]
    acis_s = [[], [], [], [], []]
    hrc_i  = [[], [], [], [], []]
    hrc_s  = [[], [], [], [], []]

    for k in range(0, len(data[0])):
        if data[4][k] == 'ACIS-I':
            for m in range(0, 5):
                acis_i[m].append(data[m][k])

        elif data[4][k] == 'ACIS-S':
            for m in range(0, 5):
                acis_s[m].append(data[m][k])

        elif data[4][k] == 'HRC-I':
            for m in range(0, 5):
                hrc_i[m].append(data[m][k])

        elif data[4][k] == 'HRC-S':
            for m in range(0, 5):
                hrc_s[m].append(data[m][k])

    return [acis_i, acis_s, hrc_i, hrc_s]
        

#------------------------------------------------------------------------------------------
#-- plot_data_info: plot trend of the information related to sim twist                   --
#------------------------------------------------------------------------------------------

def plot_data_info(data_info, year =''):
    """
    plot trend of the information related to sim twist
    input:  data_info   --- data: [time, sim_x, sim_y, sim_z, pitchamp, yawamp]
            year        --- year of the data set. if not given, plot the entire range
                            from the beginning
    output: <web_dir>/Plots/sim_plot_<year>.png   _<year> can be ""
    """
    
    ylims = [[-2.0, 0.5], [-0.05, 0.05] , [-300, 300], [0, 100], [0, 100]]
    ylab  = ['sim_x (mm)', 'sim_y (mm)', 'sim_z (mm)', 'pitchamp (sec)', 'yawamp (sec)']

    if year == "":
        outname = 'sim_plot.png'
    else:
        outname = 'sim_plot_' + str(year) + '.png'

    plot_panel(data_info, ylims, ylab, outname)


#------------------------------------------------------------------------------------------
#-- plot_data_extr: plot sim twist trend                                                 --
#------------------------------------------------------------------------------------------

def plot_data_extr(data_extr, year=''):
    """
    plot sim twist trend
    input:  data_extr   --- data: [time, dy, dz, dtheta]
            year        --- year of the data set. if not given, plot the entire range 
                            from the beginning
    output: <web_dir>/Plots/twist_plot_<year>.png   _<year> can be ""
    """
    
    ylims = [[-0.5, 2.0], [-0.1, 2.4], [-50, 50]]
    ylab  = ['dy (mm)', 'dz (mm)', 'dtheta (sec)']

    if year == '':
        outname = 'twist_plot.png'
    else:
        outname = 'twist_plot_' + str(year) + '.png'

    if year == "":
        plot_panel(data_extr, ylims, ylab, outname, lfit=2)
    else:
        plot_panel(data_extr, ylims, ylab, outname, lfit=1)

#------------------------------------------------------------------------------------------
#-- plot_dtheta: plot dthera trends                                                      --
#------------------------------------------------------------------------------------------

def plot_dtheta(acis_i, acis_s, hrc_i, hrc_s, year=''):
    """
    plot dthera trends
    input:  acis_i  --- acis i data set[[time],[data set1],...]
            acis_s  --- acis s data set
            hrc_i   --- hrc i data set
            hrc_s   --- hrc s data set
            year    --- year of the data set. if not given, plot the entire range 
                        from the beginning 
    output: <web_dir>/Plots/dtheta_plot_<year>.png _<year> can be ""
    """


    ylims = [[-30, 20], [-40, 10], [10, 60], [-40, 10]]
    ylab  = ['ACIS-I','ACIS-S','HRC-I','HRC-S']

    if year == '':
        outname = 'dtheta_plot.png'
    else:
        outname = 'dtheta_plot_' + str(year) + '.png'

    plot_panel2(acis_i, acis_s, hrc_i, hrc_s, 3, ylims, ylab, outname, lfit=1)



#------------------------------------------------------------------------------------------
#-- plot_panel: plot sim twist and sim information trends                               ---
#------------------------------------------------------------------------------------------

def plot_panel(data, ylims, ylab, outname, lfit=0):
    """
    plot sim twist and sim information trends
    input:  data    --- data: [time, dy, dz, dthera] for sim twist plots
                              [time, sim_x, sim_y, sim_z, pithcamp, yawamp] for sim infor plot
            ylims   --- a list of y limits [(ymin, ymax),....]
            ylab    ----a list of ylabel 
            outname --- output file name (without png)
            lfit    --- the indicator to show whether to fit a line. 1: yes
    output: <web_dir>/Plots/<outname>.png
    """
#
#--- data are separated into a few sections. here we set the end of the each
#--- section in fractional year. the last closing data is set to year 4000
#--- see b_period for the starting date (in fractional year)
#
    s_period = b_period
    r_len    = len(b_period)
    e_period = []
    for k in range(1, r_len):
        e_period.append(b_period[k])

    e_period.append(4000.0)


    plt.close('all')
    dlen = len(ylab)
    mlen = dlen - 1
#
#--- set panel parameters
#
    plt.subplots_adjust(hspace=0.08)
    props = font_manager.FontProperties(size=9)
#
#--- set xrange
#
    [atime, xmin, xmax, xlabel, xtext] = set_x_range(data[0])
    xrange= [xmin, xmax]
#
#--- plot each panel
#
    for k in range(0, dlen):
        j = k +  1
        line = str(dlen) + '1' + str(j)
        exec "ax%s = plt.subplot(%s)" % (str(k), line)
        exec "ax%s.set_xlim(xmin=xmin, xmax=xmax, auto=False)" % (str(k))
#
#--- if ymin and ymax are given, use them
#
        if len(ylims) == dlen:
            ymin = ylims[k][0]
            ymax = ylims[k][1]
            exec "ax%s.set_ylim(ymin=ymin, ymax=ymax, auto=False)" % (str(k))

        ydata = data[k+1]
        exec "ax%s.plot(atime, ydata, color='blue', marker='.', markersize=1, lw=0)" % (str(k))
#
#--- for the case that a fitting line is requested
        if lfit == 1:
            [a, b, e] = rfit.robust_fit(atime, ydata)
            y1        = a + b * xmin
            y2        = a + b * xmax
            yrange    = [y1, y2]
            exec "ax%s.plot(xrange, yrange, color='red', lw=2)" % (str(k))

            ydiff     = ymax - ymin
            ytext     = ymax -0.1 * ydiff
            line      = "Slope: " + '%3.3f' % (round(b, 3))
            plt.text(xtext, ytext, line)

        elif lfit > 1:
            for m in range(0, r_len):
                [xrange, yrange, ytext, line] \
                    = fit_line_period(atime, ydata, b_period[m], e_period[m],  m,  ymin, ymax, bot=0)
                exec "ax%s.plot(xrange, yrange, color='red', lw=2)" % (str(k))
                plt.text(xtext, ytext, line)
#
#--- y axis label
#
        exec "ax%s.set_ylabel('%s')" % (str(k), ylab[k])
#
#--- x axis label
#
    exec 'ax%s.set_xlabel("%s")' % (str(mlen), xlabel)
#
#--- add x ticks label only on the last panel
#
    for k in range(0, dlen):
        ax = 'ax' + str(k)

        if k != mlen:
            exec "line = %s.get_xticklabels()" % (ax)
            for label in  line:
                label.set_visible(False)
        else:
            pass

#
#--- save the plot in a file
#
    fig    = matplotlib.pyplot.gcf()
    height = 2.0 * dlen + 0.5
    fig.set_size_inches(10.0, height)

    outfile = web_dir + 'Plots/' + outname
    plt.savefig(outfile, format='png', dpi=100.0)



#------------------------------------------------------------------------------------------
#-- plot_panel2: plotting dtheta information                                            ---
#------------------------------------------------------------------------------------------

def plot_panel2(data1, data2, data3, data4, ydpos, ylims, ylab, outname, lfit=0):
    """
    plotting dtheta information
    input:  data1   --- data for acis-i data is a list of list and the first is a list of time
            data2   --- data for acis-s
            data3   --- data for hrc-i 
            data4   --- data for hrc-s 
            ydpos   --- the index of the data you want to use
            ylim    --- a list of lists of y limits [(ymin ymax),...]
            ylab    --- a list of labels for each data set to be plotted
            outname --- a name of output file
            lfit    --- indicator of whether to fit line lfit=1: yes
    output: <web_dir>/Plots/<outname>.png
    """

    plt.close('all')
    dlen = len(ylab)
    mlen = dlen - 1
#
#--- set panel parameters
#
    plt.subplots_adjust(hspace=0.08)
    props = font_manager.FontProperties(size=9)
#
#--- set xrange
#
    atime = []
    data  = []
    [xlist, xmin, xmax, xlabel, xtext] = set_x_range(data1[0])
    atime.append(xlist)
    data.append(data1[ydpos])

    [xlist, xmin, xmax, xlabel, xtext] = set_x_range(data2[0])
    atime.append(xlist)
    data.append(data2[ydpos])

    [xlist, xmin, xmax, xlabel, xtext] = set_x_range(data3[0])
    atime.append(xlist)
    data.append(data3[ydpos])

    [xlist, xmin, xmax, xlabel, xtext] = set_x_range(data4[0])
    atime.append(xlist)
    data.append(data4[ydpos])

    xrange= [xmin, xmax]
#
#--- plot each panel
#
    for k in range(0, dlen):
        j = k +  1
        line = str(dlen) + '1' + str(j)
        exec "ax%s = plt.subplot(%s)" % (str(k), line)
        exec "ax%s.set_xlim(xmin=xmin, xmax=xmax, auto=False)" % (str(k))
#
#--- if ymin and ymax are given, use them
#
        if len(ylims) == dlen:
            ymin = ylims[k][0]
            ymax = ylims[k][1]
            exec "ax%s.set_ylim(ymin=ymin, ymax=ymax, auto=False)" % (str(k))

        xdata = atime[k]
        ydata = data[k]
        exec "ax%s.plot(xdata, ydata, color='blue', marker='.', markersize=1, lw=0)" % (str(k))
#
#--- for the case that a fitting line is requested
        if lfit == 1:
            [xt, yt]  = remove_out_layer(xdata, ydata, ymin, ymax)
            [a, b, e] = rfit.robust_fit(xt, yt)
            y1        = a + b * xmin
            y2        = a + b * xmax
            yrange    = [y1, y2]
            exec "ax%s.plot(xrange, yrange, color='red', lw=2)" % (str(k))

            ydiff     = ymax - ymin
            ytext     = ymax -0.1 * ydiff
            line      = "Slope: " + '%3.3f' % (round(b, 3))
            plt.text(xtext, ytext, line)
#
#--- y axis label
#
        exec "ax%s.set_ylabel('%s')" % (str(k), ylab[k])
#
#--- x axis label
#
    exec 'ax%s.set_xlabel("%s")' % (str(mlen), xlabel)
#
#--- add x ticks label only on the last panel
#
    for k in range(0, dlen):
        ax = 'ax' + str(k)

        if k != mlen:
            exec "line = %s.get_xticklabels()" % (ax)
            for label in  line:
                label.set_visible(False)
        else:
            pass

#
#--- save the plot in a file
#
    fig    = matplotlib.pyplot.gcf()
    height = 2.0 * dlen + 0.5
    fig.set_size_inches(10.0, height)

    outfile = web_dir + 'Plots/' + outname
    plt.savefig(outfile, format='png', dpi=100.0)

 
#------------------------------------------------------------------------------------------
#-- set_x_range: setting x range. if less than 2 years, use yday, otherwise, fractional year 
#------------------------------------------------------------------------------------------

def set_x_range(xdata):
    """
    setting x range. if less than 2 years, use yday, otherwise, fractional year
    input:  xdata   --- x data
    output: ndata   --- modified data (for the case if it is in yday)
            xmin    --- x min
            xmax    --- x max
            xlabel  --- label for x axis 
            xtext   --- x location for the text display
    """
    xmin = min(xdata)
    year = int(xmin)
    xmax = max(xdata)

    xdiff = xmax - xmin
    if xdiff < 2.0:
        if tcnv.isLeapYear(year) == 1:
            base = 366
        else:
            base = 365

        ndata = []
        for ent in xdata:
            val = base * (ent - year)
            ndata.append(val)
 
            xlabel = 'Time (yday @ year=' + str(year) + ')'

            xmin   = 0
            xmax   = base

    else:
        ndata = xdata
        xmin  = int(xmin)
        xmax  = int(xmax) + 1
        xlabel = 'Time (year)'

    xdiff = xmax - xmin
    xtext = xmin + 0.1 * xdiff

    return [ndata, xmin, xmax, xlabel, xtext]
 
#------------------------------------------------------------------------------------------
#-- remove_out_layer: drop outlyer data points                                           --
#------------------------------------------------------------------------------------------

def remove_out_layer(x, y, ymin, ymax):
    """
    drop outlyer data points
    input:  x       --- x data
            y       --- y data
            ymin    --- y min of the plotting range
            ymax    --- y max of the plotting range
    output: nx      --- x data without outlayers
            ny      --- y data without outlayers
    """

    ya = []
    for m in range(0, len(y)):
        if y[m] >= ymin and y[m] <= ymax:
                ya.append(y[m])

    avg   = numpy.mean(ya)
    std   = numpy.std(ya)
    bot   = avg - 3.0 * std
    top   = avg + 3.0 * std

    nx    = []
    ny    = []
    for k in range(0, len(x)):
        if y[k] >= bot and y[k] <= top:
            nx.append(x[k])
            ny.append(y[k])

    return [nx, ny]

#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------

def read_data(infile, clean=0, emp=0):

    try:
        f    = open(infile, 'r')
        data = [line.strip() for line in f.readlines()]
        f.close()
        if len(data) > 0:
#
#--- if emp == 1, remove the empty elements
#
            if emp == 1:
                out = []
                for ent in data:
                    if ent != '':
                        out.append(ent)
                data = out
    
        if clean == 1:
            mcf.rm_file(infile)
    except:
        data = []
    
    return data

#------------------------------------------------------------------------------------------
#-- convert_to_fyear: convert time format to fractional year                            ---
#------------------------------------------------------------------------------------------

def convert_to_fyear(cdate):
    """
    convert time format to fractional year
    input:  cdate   --- time in either <yyyy>-<mm>-<dd>T<hh>:<mm>:<ss> or seconds from 1998.1.1
    output: fyear   --- time in fractional year
    """

    mc = re.search('T', cdate)
#
#--- for the case the time format is: <yyyy>-<mm>-<dd>T<hh>:<mm>:<ss>
#
    if mc is not None:
        atemp = re.split('T', cdate)
        btemp = re.split('-', atemp[0])
        ctemp = re.split(':', atemp[1])
        year  = float(btemp[0])
        mon   = float(btemp[1])
        day   = float(btemp[2])
        hh    = float(ctemp[0])
        mm    = float(ctemp[1])
        ss    = float(ctemp[2])

        ydate = mon_list[int(mon)-1] + day 
        if tcnv.isLeapYear(year) == 1:
            if mon > 2:
                ydate += 1
#
#---- for the case the time format is seconds from 1998.1.1
#
    elif mcf.chkNumeric(cdate):
        out   = Chandra.Time.DateTime(float(cdate)).date
        atemp = re.split(':', out)

        year  = float(atemp[0])
        ydate = float(atemp[1])
        hh    = float(atemp[2])
        mm    = float(atemp[3])
        ss    = float(atemp[4])

    else:
        atemp = re.split(':', cdate)
        year  = float(atemp[0])
        ydate = float(atemp[1])
        hh    = float(atemp[2])
        mm    = float(atemp[3])
        ss    = float(atemp[4])

    ydate = ydate + hh / 24.0 + mm / 1440.0 + ss / 86400.0

    if tcnv.isLeapYear(year) == 1:
        base = 366.0
    else:
        base = 365.0

    fyear = year + ydate / base

    return fyear

#------------------------------------------------------------------------------------------
#-- quickChandra_time: axTime3 replacement                                             ----
#------------------------------------------------------------------------------------------

def quickChandra_time(ent):
    """
    axTime3 replacement
    input:  ent --- either seconds from 1998.1.1 or date in <yyyy>:<ddd>:<hh>:<mm>:<ss>
    output: out --- either seconds from 1998.1.1 or date in <yyyy>:<ddd>:<hh>:<mm>:<ss>
    """

    if mcf.chkNumeric(ent):
        out = Chandra.Time.DateTime(float(ent)).date
    else:
        out = Chandra.Time.DateTime(str(ent)).secs

    return out

#------------------------------------------------------------------------------------------
#-- fit_line_period: fit a line on the given data section and return the results         --
#------------------------------------------------------------------------------------------

def fit_line_period(x, y, x_start, x_stop, m,  ymin, ymax, bot=0):
    """
    fit a line on the given data section and return the result
    input:  x       --- x data
            y       --- y data
            x_start --- a list of the starting points
            x_stop  --- a list of the stopping points
            m       --- section position of the lists above
            ymin    --- y min
            yax     --- y max
            bot     --- an indicator of where to print the text. bot =1 to bottom side
    output: xrange  --- a list of start and stop position of the line in x 
            yrange  --- a list of start and stop position of the line in y
            ytext   --- y position of the text
            line    --- a line to be printed
    """
    dlen = len(x)
    xn   = []
    yn   = []
#
#--- select data for the limited range
#
    for k in range(0, dlen):
        if x[k] >= x_start and x[k] < x_stop:
            if y[k] >= ymin and y[k] <= ymax:
                xn.append(x[k])
                yn.append(y[k])

    try:
#
#--- compute the fitted line with robust method
#
        [a, b, e] = rfit.robust_fit(list(xn), list(yn))
        y1        = a + b * x_start
        y2        = a + b * x_stop
        xrange    = [x_start, x_stop]
        yrange    = [y1, y2]
#
#--- set the text postion and create the text to be printed
#
        ydiff     = ymax - ymin
        if bot == 0:
            ytext     = ymax -0.1 * ydiff * (m + 1)
        else:
            ytext     = ymax -0.1 * ydiff * (m + 1) -0.5 * ydiff


        if x_start < 2000:
            line      = "Slope (%4.1f < year):  %3.3f" % (x_stop, round(b, 3))
        elif x_stop > 3000:
            line      = "Slope (year > %4.1f):  %3.3f" % (x_start, round(b, 3))
        else:
            line      = "Slope (%4.1f <  year < %4.1f): %3.3f" % (x_start, x_stop, round(b, 3))

#
#--- for the case the fitting failed
#
    except:
        ydiff     = ymax - ymin
        if bot == 0:
            ytext     = ymax -0.1 * ydiff * (m + 1)
        else:
            ytext     = ymax -0.1 * ydiff * (m + 1) - 0.5 * ydiff

        xrange = [0, 0]
        yrange = [0, 0]
        line   = 'NA'

    return [xrange, yrange, ytext, line]


#------------------------------------------------------------------------------------------

if __name__ == "__main__":

    if len(sys.argv) > 1:
        year   = int(float(sys.argv[1]))
    else:
        year   = ''

    sim_twist_trend_plot(year)






