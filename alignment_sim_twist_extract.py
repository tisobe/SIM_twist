#!/usr/bin/env /proj/sot/ska/bin/python

#########################################################################################
#                                                                                       #
#           alignment_sim_twist_extract.py: extract sim related data                    #
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
#-- alignment_sim_twist_extract: extract sim twist data and update the databases         --
#------------------------------------------------------------------------------------------

def alignment_sim_twist_extract(tstart='', tstop='', year=''):
    """
    extract sim twist data and update the databases
    input:  tstart  --- starting time in the format of <yyyy>:<ddd>:<hh>:<mm>:<ss>
                        default: "" --- the script will find the starting time
            tstop   --- stopping time in the format of <yyyy>:<ddd>:<hh>:<mm>:<ss>
                        default: "" --- the script will find the stopping time
            year    --- year, default: "".
    output: <data_dir>/data_info_<year>
            <data_dir>/data_extracted_<year>
    """
#
#--- if the starting and stopping time are given, use them
#
    if tstart != '':
        ending = extract_data(tstart, tstop, year, ldir ='./')
#
#--- otherwise, find stating and stopping time
#
    else:
#
#--- current time
#
        today = time.strftime("%Y:%j:%H:%M:00", time.gmtime())
        year  = int(float(time.gmtime().tm_year))
#
#--- find the last entry time
#
        ifile  = data_dir + 'data_info_' + str(year)
        if os.path.isfile(ifile):
            data  = read_data(ifile, emp=1)
        else:
#
#--- if the current year data files are not there, just creates an empty list
#
            data = []
#
#--- if this is a new year, read the last year's data
#
        if len(data) == 0:
            ifile2  = data_dir + 'data_info_' + str(year-1)
            tstart = get_last_entry_time(ifile2)
    
            tstop  = str(year) + ':001:00:00:00'

            ending = extract_data(tstart, tstop, year-1)
#
#--- check the last entry again for the next round
#
            tstart = convert_time_format(ending)
#
#--- work on this year's data
#
        else: 
            tstart = get_last_entry_time(ifile)
    
        ending = extract_data(tstart, today, year)
    
#------------------------------------------------------------------------------------------
#-- get_last_entry_time: find the last entry date from the database                     ---
#------------------------------------------------------------------------------------------

def get_last_entry_time(infile, pos=4):
    """
    find the last entry date from the database
    input:  infile  --- input file, ususally data_infor_<year> and assume that
                        pos-th entry is the (ending) time
            pos     --- position of the time element in the data (usually 4th)
    """

    data   = read_data(infile)

    if len(data) == 0:
        print "something wrong in starting time; exist"
        exit(1)
    else:
        atemp  = re.split('\s+', data[-1])
        ltime  = atemp[pos]
#
#--- convert time format to be used
#
        ftime  = convert_time_format(ltime)

        return ftime


#------------------------------------------------------------------------------------------
#-- convert_time_format: convert data fromat                                             --
#------------------------------------------------------------------------------------------

def convert_time_format(ltime):
    """
    convert data fromat from <yyyy>-<mm>-<dd>T<hh>:<mm>:<ss> to <yyyy>:<ddd>:<hh>:<mm>:<ss>
    input:  ltime   time in the format of <yyyy>-<mm>-<dd>T<hh>:<mm>:<ss>
    output: date    time in the format of <yyyy>:<ddd>:<hh>:<mm>:<ss>
    """
    
    atemp = re.split('T', ltime)
    btemp = re.split('-', atemp[0])
    year  = int(float(btemp[0]))
    mon   = int(float(btemp[1]))
    day   = int(float(btemp[2]))

    yday  = mon_list[mon-1] + day
    if (tcnv.isLeapYear(year) == 1) and (mon > 2):
        yday += 1

    lyday = str(yday)
    if yday < 10:
        lyday = '00' + lyday
    elif yday < 100:
        lyday = '0'  + lyday

    date = str(year) + ':' + lyday + ':' + atemp[1]

    return date


#------------------------------------------------------------------------------------------
#-- extract_data: using arc5gl, extract fits file and extract data                       --
#------------------------------------------------------------------------------------------

def extract_data(tstart, tstop, year, ldir=''):
    """
    using arc5gl, extract fits file and extract data
    input:  tstart  --- starting time in <yyyy>:<ddd>:<hh>:<mm>:<ss>
            tstop   --- stopping time in <yyyy>:<ddd>:<hh>:<mm>:<ss>
            year    --- the year of the data
            ldir    --- output directory. default is "" and if that is the case,
                        it will use <data_dir>
    output: <ldir>/data_info_<year>
            <ldir>/data_extracted_<year>
            ending  --- the last date in the corrected data in <yyyy>:<ddd>:<hh>:<mm>:<ss> 
    """
#
#--- if the output directory is not specified, use <data_dir> 
#
    if ldir == '':
        ldir = data_dir
#
#--- extract the possible candidates for the data
#
    inlist   = call_arc5gl('browse', 'pcad', '1', tstart, tstop, sub='aca')
    aca_list = []
#
#---- we need only "asol" files
#
    for ent in inlist:
        mc = re.search('asol', ent)
        if mc is not None:
            atemp = re.split('\s+', ent)
            aca_list.append(atemp[0])

    info_list = []
    data_list = []
#
#--- now extract the "asol" files one by one
#
    for infits in aca_list:
        fout   = call_arc5gl('retrieve', 'pcad', '1', tstart='', tstop='', sub='aca', file=infits)
        try:
            fits   = fout[0]
        except:
            continue
#
#--- get header information
#
        [out, obsid, inst] = get_header_info(fits)
        info_list.append(out)
#
#--- get the data we need
#
        out = get_data(fits, obsid, inst)
        data_list = data_list + out

        mcf.rm_file(fits)
#
#--- write the extracted information
#
    outfile = ldir + 'data_info_' + str(year)
    if os.path.isfile(outfile):
        fo      = open(outfile, 'a')
    else:
        fo      = open(outfile, 'a')
    for line in info_list:
        fo.write(line)
    fo.close()
#
#--- the last date of the data extraction (closing time)
#
    atemp  = re.split('\s+', info_list[-1])
    ending = atemp[4]
#
#---- write the extracted data
#
    outfile = ldir + 'data_extracted_' + str(year)
    if os.path.isfile(outfile):
        fo      = open(outfile, 'a')
    else:
        fo      = open(outfile, 'w')
    for line in data_list:
        fo.write(line)
    fo.close()

    return ending

#------------------------------------------------------------------------------------------
#-- get_header_info: read the header of the fits file and return the information needed  --
#------------------------------------------------------------------------------------------

def get_header_info(fits):
    """
    read the header of the fits file and return the information needed
    input:  fits    --- fits file name
    output: line    --- a line contains the information needed:
                        <fits name> <obsid> <instrument> <date obs> <date end>
                        <sim x> <sim y> <sim z> <pitchamp> <yawamp>
            obsid   --- obsid   
            inst    --- instrument
                            ---- these are duplicate but easier for the calling function to use
    """
#
#--- read header of the fits file. assume that extension is 1.
#
    ff = pyfits.open(fits)
    head = ff[1].header
    date_obs = head['DATE-OBS']
    date_end = head['DATE-END']
    sim_x    = head['SIM_X']
    sim_y    = head['SIM_Y']
    sim_z    = head['SIM_Z']
    pitchamp = head['PITCHAMP']
    yawamp   = head['YAWAMP']
    obsid    = head['OBS_ID']
    ff.close()
#
#--- for the instrument, call sql to find it
#
    try:
        sqlcmd   = 'select instrument from target where obsid=' + str(obsid)
        [inst]   = readSQL(sqlcmd, 'axafocat')
    except:
        inst = 'NA'

    fits2    = fits.replace('.gz', '')
    line     = fits2 + '\t'
    line     = line + obsid         + '\t'
    line     = line + inst          + '\t'
    line     = line + date_obs      + '\t'
    line     = line + date_end      + '\t'
    line     = line + str(sim_x)    + '\t'
    line     = line + str(sim_y)    + '\t'
    line     = line + str(sim_z)    + '\t'
    line     = line + str(pitchamp) + '\t'
    line     = line + str(yawamp)   + '\n'

    return [line, obsid, inst]


#------------------------------------------------------------------------------------------
#-- get_data: extract data from the fits file                                           ---
#------------------------------------------------------------------------------------------

def get_data(fits, obsid, inst):
    """
    extract data from the fits file
    input:  fits    --- fits file name
            obsid   --- obsid
            inst    --- instrument
    output: out_list--- a list of lines containing the data 
                        <time> <obsid> <instrument> <dy> <dz> <dtheta>
    """

    ff = pyfits.open(fits)
    data = ff[1].data
    tin  = data.field('time')
    dyin = data.field('dy')
    dzin = data.field('dz')
    dtin = data.field('dtheta')
    ff.close()

    time   = []
    d_dict = {}

    dlen = len(tin)
    step = int(dlen/1170) + 2

    for k in range(0, step):
        begin = 1170 * k 
        end   = begin + 1170
        
        if begin >= dlen:
            break
        elif end > dlen:
#
#--- the last set of the data must have at least 500 data points; otherwise drop it
#
            if (dlen - end) > 500:
                try:
                    v1 = int(numpy.median(tin[begin:]))
                    v2 = numpy.mean(dyin[begin:])
                    v3 = numpy.mean(dzin[begin:])
                    v4 = numpy.mean(dtin[begin:])
                    time.append(v1)
                    d_dict[str(v1)] = [v2, v3, v4]
                except:
                    pass
        else:
            try:
                v1 = int(numpy.median(tin[begin:end]))
                v2 = numpy.mean(dyin[begin:end])
                v3 = numpy.mean(dzin[begin:end])
                v4 = numpy.mean(dtin[begin:end])
                time.append(v1)
                d_dict[str(v1)] = [v2, v3, v4]
            except:
                pass
#
#--- remove duplicate
#
    time = list(set(time))

    out_list = []
    for k in range(0, len(time)):
        stime = str(time[k])
        line = stime        + '\t'
        line = line + obsid + '\t'
        line = line + inst  + '\t'
        line = line + str(d_dict[stime][0]) + '\t'
        line = line + str(d_dict[stime][1]) + '\t'
        line = line + str(d_dict[stime][2]) + '\n'

        out_list.append(line)

    return out_list

#------------------------------------------------------------------------------------------
#-- call_arc5gl: using arc5gl to extract a fits file list or a file itself              ---
#------------------------------------------------------------------------------------------

def call_arc5gl(op, detector, level, tstart='', tstop='', sub='', file=''):
    """
    using arc5gl to extract a fits file list or a file itself
    input:  op      ---- operation: retreive/browse
            detector    --- detector
            level       --- level
            tstart      --- starting time if file name is provided, it is ignored
            tstop       --- stopping time if file name is provided, it is ignored
            sub         --- sub detector name; default ""
            file        --- file name; default ""
    output: flist       --- a list of fits file, either results of browse or extracted file name
    """

    line = 'operation=' + op + '\n'
    line = line + 'dataset=flight\n'
    line = line + 'detector='     + detector    + '\n'

    if sub != '':
        line = line + 'subdetector=' + sub      + '\n'

    line = line + 'level='        + str(level)  + '\n'

    if file == '':
        line = line + 'tstart='   + str(tstart) + '\n'
        line = line + 'tstop='    + str(tstop)  + '\n'
    else:
        line = line + 'filename=' + file        + '\n'

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
            flist.append(str(ent))

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
#-- readSQL: for a given sql command and database, return result in dictionary form.     --
#------------------------------------------------------------------------------------------

def readSQL(cmd, database):
    """
    for a given sql command and database, return result in dictionary form.
    input   cmd         --- sql command
            database    --- the database to be used
    output: outdat      --- results
    """
    db_user = 'mtaops_internal_web'
    server  = 'ocatsqlsrv'

    line      = pass_dir + '.targpass_internal'
    f         = open(line, 'r')
    db_passwd = f.readline().strip()
    f.close()

    db = DBI(dbi='sybase', server=server, user=db_user, passwd=db_passwd, database=database)

    outdata = db.fetchone(cmd)
    atemp   = re.split('select', cmd)
    btemp   = re.split('from', atemp[1])
    ctemp   = re.split('\s+|,', btemp[0])

    results = []
    for ent in ctemp:
        if ent == '':
            continue
        ent = ent.strip()
        try:
            results.append(outdata[ent])
        except:
            pass

    return results



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

    alignment_sim_twist_extract(tstart, tstop, year)




