#!/usr/bin/env /proj/sot/ska/bin/python

#########################################################################################
#                                                                                       #
#               fid_light_data_extract.py: extract fid light information                #
#                                                                                       #
#               author: t. isobe (tisobe@cfa.harvard.edu)                               #
#                                                                                       #
#               last update: Apr 19, 2017                                               #
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

#
#--- from ska
#
from Ska.Shell import getenv, bash
ascdsenv = getenv('source /home/ascds/.ascrc -r release; source /home/mta/bin/reset_param', shell='tcsh')
ascdsenv['MTA_REPORT_DIR'] = '/data/mta/Script/ACIS/CTI/Exc/Temp_comp_area/'

import Chandra.Time
import Ska.engarchive.fetch as fetch
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
from DBI import *

#
#--- temp writing file name
#
rtail  = int(10000 * random.random())       #---- put a romdom # tail so that it won't mix up with other scripts space
zspace = '/tmp/zspace' + str(rtail)

mon_list = [0, 31, 59, 90, 120, 151, 181, 212, 234, 373, 304, 334]

#------------------------------------------------------------------------------------------
#-- fid_light_data_extract: extract aca position i and j, sim postion, and creates a table 
#------------------------------------------------------------------------------------------

def fid_light_data_extract(tstart='', tstop='', year=''):
    """
    extract aca position i and j, sim postion, and creates a table
    input:  tstart  --- starting time
            tstop   --- stopping time
            year    --- year
            if these are not given, the the interval of the previous month are used
    output: udated data file, e.g., I-1, S-3, H-I-1, H-S-3 etc in <data_dir>
    """
#
#--- if the starting and stopping time are not given, set them
#
    if tstart == '':
#
#--- current time
#
        today = time.strftime("%Y:%j:%H:%M:00", time.gmtime())
        year  = int(float(time.gmtime().tm_year))
        mon   = int(float(time.gmtime().tm_mon))
#
#--- set the data extraction time span to the entire previous month
#
        pyear = year
        pmon  = mon - 1
        if pmon < 1:
            pmon   = 12
            pyear -= 1
            lydate = find_ydate(pyear, pmon, 1, string=1)
            tstart = str(pyear) + ':' + lydate + ':00:00:00'

            lydate = find_ydate(year, mon, 1, string=1)
            tstop  = str(year)  + ':' + lydate + ':00:00:00'

#
#--- get fid light information
#
    [file_id, fid_detect] = extract_fidlight_data(tstart, tstop)
#
#--- get oterh informatin related to the fid information and print out the results
#
    get_acen_data(tstart, tstop, file_id, fid_detect)

#------------------------------------------------------------------------------------------
#-- get_acen_data: for given fid light information, extract acen fits files and extract needed information
#------------------------------------------------------------------------------------------

def get_acen_data(tstart, tstop, file_id, fid_detect):
    """
    for given fid light information, extract acen fits files and extract needed information
    input:  tstart  --- starting time
            tstop   --- stopping time
            file_id --- id of the acen fits file
            fid_detect  --- a dictionary of information [<slot id>, <id string> <id number>]
    output: udated data file, e.g., I-1, S-3, H-I-1, H-S-3 etc in <data_dir>
    """
#
#--- first just get a list of potential acent fits files
#
    acent_list = call_arc5gl('browse', 'pcad', 1, tstart=tstart, tstop=tstop, filetype='acacent', sub='aca')
#
#--- compare the list with a file id, and if it is found, procceed farther
#
    for fname in file_id:
        for comp in acent_list:
            mc = re.search(fname, comp)
            if mc is not None:
                filename = comp
                break
#
#--- extract an acen fits file
#
        [fits] = call_arc5gl('retrieve', 'pcad', 1, tstart='', tstop='', filetype='acacent', filename=filename, sub='aca')

        ff     = pyfits.open(fits)
        data   = ff[1].data
        ff.close()
        mcf.rm_file(fits)
#
#--- extract needed information for each slot 
#
        for m in range(0, len(fid_detect[fname][0])):
            slot_id =  fid_detect[fname][0][m]
            mask    = data['slot'] == slot_id
            out     = data[mask]
            time    = out['time']
            cent_i  = out['cent_i']
            cent_j  = out['cent_j']
            ang_y   = out['ang_y']
            ang_z   = out['ang_z']
            alg     = out['alg']
     
            ofile  = data_dir + fid_detect[fname][1][m]
            if os.path.isfile(ofile):
                fo    = open(ofile, 'a')
            else:
                fo    = open(ofile, 'w')
#
#---- take 5 min average for the data
#
            begin  = time[0]
            end    = begin + 300.0
            m      = 0
            k_list = []
            for k in range(m, len(time)):
                if time[k] < begin:
                    continue
                elif time[k] >= begin and time[k] < end:
                    k_list.append(k)
                else:
                    try:
                        atime   = numpy.mean(time[k_list[0]:k_list[-1]])
                        acent_i = numpy.mean(cent_i[k_list[0]:k_list[-1]])
                        acent_j = numpy.mean(cent_j[k_list[0]:k_list[-1]])
                        aang_y  = numpy.mean(ang_y[k_list[0]:k_list[-1]])
                        aang_z  = numpy.mean(ang_z[k_list[0]:k_list[-1]])
                        aalg    = alg[k_list[-1]]
                    except:
                        continue
#
#--- find fapos and tscpos info near to the given time interval
#
                    flist   = fetch.MSID('3fapos',  time[k_list[0]], time[k_list[-1]])
                    tslist  = fetch.MSID('3tscpos', time[k_list[0]], time[k_list[-1]])
                    fapos   = numpy.mean(flist.vals)
                    tscpos  = numpy.mean(tslist.vals)
#
#--- print out the results
#
                    line    = str(atime) + '\t'
                    line    = line  + str(fid_detect[fname][0][m]) + '\t'
                    line    = line  + str(fid_detect[fname][2][m]) + '\t'
                    line    = line  + str(aalg)    + '\t'
                    line    = line  + str(format(acent_i, '.3f'))  + '\t'
                    line    = line  + str(format(acent_j, '.3f'))  + '\t'
                    line    = line  + str(format(aang_y,  '.6f'))  + '\t'
                    line    = line  + str(format(aang_z,  '.6f'))  + '\t'
                    line    = line  + str(fapos)   + '\t'
                    line    = line  + str(tscpos)  + '\n'
                    fo.write(line)
    
                    k_list = []
                    begin  = end
                    end   += 300.0

            fo.close()

#------------------------------------------------------------------------------------------
#-- extract_fidlight_data: get fid light information from fidpr fits file                --
#------------------------------------------------------------------------------------------

def extract_fidlight_data(tstart, tstop):
    """
    get fid light information from fidpr fits file
    input:  tstart      --- starting time
            tstop       --- stopping time
    output: file_id     --- acent fits file id
            fid_detect  --- a dictionary of information [<slot id>, <id string> <id number>]
                            id string is I-<#>/S-<#> for ACIS and H-I-<#>/H-S-<#> for HRC

        pcadf286408332N003_fidpr1.fits
        ROW    slot       id_string    id_num

        1          0 ACIS-I-1              1
        2          1 ACIS-I-5              5
        3          2 ACIS-I-6              6
    """
#
#--- retrieve fidpr fits file
#
    fid_list = call_arc5gl('retrieve', 'pcad', 1, tstart=tstart, tstop=tstop, filetype='fidprops', sub='aca')

    file_id    = []
    fid_detect = {}
    for fits in fid_list:
        atemp   = re.split('pcadf', fits)
        btemp   = re.split('N',  atemp[1])
        file_n  = btemp[0]
#
#--- get fit light infor. see above for the structure of the table data
#
        #fits    = fits + '.gz'
        ff      = pyfits.open(fits)
        data    = ff[1].data
        ff.close()
        slot    = data['slot']
        id      = get_id_name(data['id_string'])
        id_n    = data['id_num']
#
#--- data is saved as a dictionary form
#
        fid_detect[file_n] = [slot, id, id_n]
        file_id.append(file_n)

        mcf.rm_file(fits)

    return [file_id, fid_detect]

#------------------------------------------------------------------------------------------
#-- get_id_name: constract chip id                                                       --
#------------------------------------------------------------------------------------------

def get_id_name(i_list):
    """
    constract chip id
    input:  i_list  --- a list of id string from fidpr fits file
    output: out     --- a list of chip ids
    """

    out = []
    for ent in i_list:
        mc = re.search('ACIS', ent)
        if mc is not None:
            atemp = re.split('ACIS-', ent)
            out.append(atemp[1])
        else:
            atemp = re.split('HRC-',  ent)
            name = 'H-' + str(atemp[1])
            out.append(name)

    return out
    
#------------------------------------------------------------------------------------------
#-- find_ydate: find ydate for given year/mon/day                                        --
#------------------------------------------------------------------------------------------

def find_ydate(year, mon, day, string=0):
    """
    find ydate for given year/mon/day
    input:  year    --- year
            mon     --- month
            day     --- day of month
            string  --- indicator whether to convert into string. 0: no, 1: yes
    output: ydate   --- ydate
    """

    ydate = mon_list[mon-1] + day 
    if (tmcv.isLeapYear(year) == 1) and mon > 2:
        ydate += 1

    if string == 1:
        lydate = str(ydate)
        if ydate < 10:
            lydate = '00' + lydate
        elif ydate < 100:
            lydate = '0' + lydate

        return lydate
    else:
        return ydate


#------------------------------------------------------------------------------------------
#-- call_arc5gl: using arc5gl to extract a fits file list or a file itself              ---
#------------------------------------------------------------------------------------------

def call_arc5gl(op, detector, level, tstart='', tstop='', filetype='', sub='', filename=''):
    """
    using arc5gl to extract a fits file list or a file itself
    input:  op      ---- operation: retreive/browse
            detector    --- detector
            level       --- level
            tstart      --- starting time if file name is provided, it is ignored
            tstop       --- stopping time if file name is provided, it is ignored
            filetype    --- filetype if file name is provided, it is ignored
            sub         --- sub detector name; default ""
            filename    --- file name; default ""
    output: flist       --- a list of fits file, either results of browse or extracted file name
    """

    line = 'operation=' + op + '\n'
    line = line + 'dataset=flight\n'
    line = line + 'detector='     + detector    + '\n'

    if sub != '':
        line = line + 'subdetector=' + sub      + '\n'

    line = line + 'level='        + str(level)  + '\n'
    line = line + 'filetype='     + filetype    + '\n'

    if filename == '':
        line = line + 'tstart='   + str(tstart) + '\n'
        line = line + 'tstop='    + str(tstop)  + '\n'
    else:
        line = line + 'filename=' + filename    + '\n'

    line = line + 'go\n'
    fo = open(zspace, 'w')
    fo.write(line)
    fo.close()

    cmd = ' /proj/axaf/simul/bin/arc5gl -user isobe -script ' + zspace + ' >ztemp_out'
    run_ascds(cmd)
    mcf.rm_file(zspace)

    data = read_data('ztemp_out', clean=1)
    flist = []
    for ent in data:
        mc = re.search('fits', ent)
        if mc is not None:
            atemp = re.split('\s+', ent)
            if len(atemp) > 1:
                flist.append(str(atemp[0]))
            else:
                flist.append(ent)

    return flist


#------------------------------------------------------------------------------------------
#-- run_ascds: run the command in ascds environment                                      --
#------------------------------------------------------------------------------------------

def run_ascds(cmd):
    """
    run the command in ascds environment
    input:  cmd --- command line
    output: command results
    """
    acmd = '/usr/bin/env PERL5LIB=""  ' + cmd

    try:
        bash(acmd, env=ascdsenv)
    except:
        try:
            bash(acmd, env=ascdsenv)
        except:
            pass

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

if __name__ == "__main__":

    if len(sys.argv) > 3:
        tstart = sys.argv[1]
        tstop  = sys.argv[2]
        year   = int(float(sys.argv[3]))
    else:
        tstart = ''
        tstop  = ''
        year   = ''

    fid_light_data_extract(tstart, tstop, year)




