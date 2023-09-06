#!/usr/bin/env python
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# This code is an analysis module for WAD-QC 2.0: a server for automated 
# analysis of medical images for quality control.
#
# The WAD-QC Software can be found on 
# https://bitbucket.org/MedPhysNL/wadqc/wiki/Home
# 
#
# Changelog:
#   20230906: remove deprecated warning; Pillow 10.0.0
#   20200508: dropping support for python2; dropping support for WAD-QC 1; toimage no longer exists in scipy.misc
#   20200421: Fix some frequent OCR problems
#   20190426: Fix for matplotlib>3
#   20180913: new format of config for OCR: ocr_regions = {name: {prefix:, suffix:, type:, xywh}}; 
#             tesseract wants black text on white
#   20170925: distinquish params for actions
#   20170907: Added uni_range_model parameter
#   20170830: Added rev_bbox and auto_suffix parameters
#   20170510: Added rgbchannel param, defaults to 'B'; 
#             added optional parameters cluster_model, uni_start, ocr_threshold, ocr_zoom; 
#             removed uni_low; removed reverb.jpg output (added to uniformity)
#   20170201: bugfix for parameter mixup uni_depth and uni_low
#   20161221: changes in default positions for uniformity, scaling of uniformity, extra config params (PvH)
#   20161220: removed class variables; removed testing stuff
#   20160825: added extra config parameters (PvH)
#   20160802: sync with pywad1.0
#   20160622: removed adding limits (now part of analyzer)
#   20160620: remove quantity and units
#
# mkdir -p TestSet/StudyCurve
# mkdir -p TestSet/Config
# cp ~/Downloads/1/us_philips_*.xml TestSet/Config/
# ln -s /home/nol/WAD/pyWADdemodata/US/US_AirReverberations/dicom_curve/ TestSet/StudyCurve/
# ./QCUS_wadwrapper.py -d TestSet/StudyEpiqCurve/ -c Config/us_philips_epiq_instance.json -r results_epiq.json
#

__version__ = '20230906'
__author__ = 'aschilham'

import os
# this will fail unless wad_qc is already installed
from wad_qc.module import pyWADinput
from wad_qc.modulelibs import wadwrapper_lib

import numpy as np

try:
    from scipy.misc import toimage
except (ImportError, AttributeError) as e:
    try:
        from wad_qc.modulelibs.wadwrapper_lib import toimage as toimage
    except (ImportError, AttributeError) as e:
        msg = "Function 'toimage' cannot be found. Either downgrade scipy or upgrade WAD-QC."
        raise AttributeError("{}: {}".format(msg, e))

# sanity check: we need at least scipy 0.10.1 to avoid problems mixing PIL and Pillow
import scipy
scipy_version = [int(v) for v in scipy.__version__ .split('.')]
if scipy_version[0] == 0:
    if scipy_version[1]<10 or (scipy_version[1] == 10 and scipy_version[1]<1):
        raise RuntimeError("scipy version too old. Upgrade scipy to at least 0.10.1")

if not 'MPLCONFIGDIR' in os.environ:
    try:
        # new method
        from importlib.metadata import version as pkg_version
    except:
        # deprecated method
        import pkg_resources
        def pkg_version(what):
            return pkg_resources.get_distribution(what).version
    try:
        #only for matplotlib < 3 should we use the tmp work around, but it should be applied before importing matplotlib
        matplotlib_version = [int(v) for v in pkg_version("matplotlib").split('.')]
        if matplotlib_version[0]<3:
            os.environ['MPLCONFIGDIR'] = "/tmp/.matplotlib" # if this folder already exists it must be accessible by the owner of WAD_Processor 
    except:
        os.environ['MPLCONFIGDIR'] = "/tmp/.matplotlib" # if this folder already exists it must be accessible by the owner of WAD_Processor 

import matplotlib
matplotlib.use('Agg') # Force matplotlib to not use any Xwindows backend.

import QCUS_lib
import ocr_lib

def logTag():
    return "[QCUS_wadwrapper] "

# MODULE EXPECTS PYQTGRAPH DATA: X AND Y ARE TRANSPOSED!

##### Series wrappers
def setup_series(inputfile, params, headers_only, for_action):
    """
    Shared routine to set runtime parameters and build structure
    Workflow:
        1. Set runtime parameters
        2. Check data format
        3. Build and populate qcstructure
    """
    # 1. Set runtime parameters
    # 2. Check data format
    rgbchannel = params.get('rgbchannel', 'B')
    dcmInfile, pixeldataIn, dicomMode = wadwrapper_lib.prepareInput(inputfile, headers_only=headers_only, logTag=logTag(), rgbchannel=rgbchannel)

    # 3. Build and populate qcstructure
    qclib = QCUS_lib.US_QC(guimode=False)
    cs = QCUS_lib.USStruct(dcmInfile,pixeldataIn,dicomMode)
    
    if for_action == 'ocr':
        optional_pars = [ # parname, attribname, type, default
            ('verbose', 'verbose', 'boolean', False),        # by default do not produce detailed logging
            ('auto_suffix', 'auto_suffix', 'boolean', True), # by default try to add a probename suffix to each result
            # ocr_threshold and ocr_zoom handled in OCR action
        ]

        required_pars = [
        ]

    elif for_action == 'qc':
        optional_pars = [ # parname, attribname, type, default
            ('verbose', 'verbose', 'boolean', False),        # by default do not produce detailed logging
            ('auto_suffix', 'auto_suffix', 'boolean', True), # by default try to add a probename suffix to each result
            ('signal_thresh', 'signal_thresh', 'int', None),  # no default value
            ('cluster_mode', 'cluster_mode', 'string', None), # no default value
            ('uni_start', 'uni_start', 'string', None), # no default value
            ('rev_bbox', 'rev_forcebbox', 'stringint4', None), # no default value
            ('uni_range_model', 'uni_range_model', 'string', None), # no default value
            ]
        required_pars = [ # parname, attribname, type
            ('uni_filter', 'uni_filter', 'int'),
            ('uni_delta', 'uni_delta', 'float'),
            ('uni_depth', 'uni_depth', 'float'),
            ('sen_filter', 'sen_filter', 'int'),
            ('sen_delta', 'sen_delta', 'float'),
            ('ver_offset', 'ver_offset', 'int'),
            ('hor_offset', 'hor_offset', 'int'),
            ('fitcircle_frac', 'fitcircle_frac', 'float'),
            ('cluster_fminsize', 'cluster_fminsize', 'float'),
        ]
        
    elif for_action == 'hdrs':
        optional_pars = [ # parname, attribname, type, default
            ('verbose', 'verbose', 'boolean', False),        # by default do not produce detailed logging
            ('auto_suffix', 'auto_suffix', 'boolean', True), # by default try to add a probename suffix to each result
            ]
        required_pars = [
        ]

    # optional parameters
    for parname, attname, partype, default in optional_pars:
        try:
            val = params[parname]
            if partype == 'boolean':
                if isinstance(val, str):
                    val = val.lower() in  ['1', 'true', 'y', 'yes']
            elif partype == 'int':
                val = int(val)
            elif partype == 'stringint4':
                val =  [int(v) for v in val.split(';')]

            setattr(cs, attname, val)
        except:
            if not default is None: # only if a default value is given
                setattr(cs, attname, default)

        # optional parameters
    
    #required
    for parname, attname, partype in required_pars:
        val = params[parname]
        if partype == 'boolean':
            if isinstance(val, str):
                val = val.lower() in  ['1', 'true', 'y', 'yes']
        elif partype == 'int':
            val = int(val)
        elif partype == 'stringint4':
            val =  [int(v) for v in val.split(';')]

        setattr(cs, attname, val)

    
        
    return qclib,cs

def get_idname_from_ocr(data, params):
    """
    Separate function to generate idname from ocr, because it is needed by several functions
    """
    try:
        dummy = params['OCR_probeID:xywh']
        # build a fake action
        idparams = {}
        base = 'OCR_probeID'
        for tag in ['xywh', 'type', 'prefix', 'suffix']:
            try:
                name = '%s:%s'%(base, tag)
                val = params[name]
                idparams[name] = val
            except:
                pass
        values, error, msg = OCR(data, None, {'params':idparams}, idname='')
        if error:
            raise ValueError("Cannot find values for %s: %s"%(base, msg))
        idname = '_'+ values[base].replace('/','-')
    except Exception as e:
        idname = None # OCR cannot be found
    
    return idname

def qc_series(data, results, action, idname):
    """
    US Reverberations in Air analysis:
        Check the uniformity of the reverberation patterns

    Params needs to define:
       nothing yet
       
    Workflow:
        1. Set runtime parameters
        2. Check data format
        3. Build and populate qcstructure
        4. Run tests
        5. Build xml output
        6. Build artefact picture thumbnail
    """

    # 1.-3. 
    try:
        params = action['params']
    except KeyError:
        params = {}

    inputfile = data.series_filelist[0]  # give me a [filename]
    qclib,cs = setup_series(inputfile, params, headers_only=False, for_action='qc')
    if cs.auto_suffix:
        if idname is None:
            idname = '_'+qclib.imageID(cs,probeonly=True)
    else:
        idname = ''
    cs.resultslabel = idname[1:]

    # 4. Run tests
    error = qclib.Analyse(cs)

    # 5. Build xml output
    labvals = qclib.reportEntries(cs)

    #    image_fnames = [] # filenames of generated images

    for key,val,lev in labvals:
        varname = key+str(idname)
        results.addFloat(varname, val)

    # 6. also store images as results
    for fn in cs.image_fnames:
        varname = os.path.splitext(fn)[0]
        results.addObject(varname, fn)

    return cs

def acqdatetime_series(data, results, action):
    """
    Read acqdatetime from dicomheaders and write to IQC database

    Workflow:
        1. Read only headers
    """
    try:
        import pydicom as dicom
    except ImportError:
        import dicom
    try:
        params = action['params']
    except KeyError:
        params = {}

    ## 1. read only headers
    dcmInfile = dicom.read_file(data.series_filelist[0][0], stop_before_pixels=True)

    dt = wadwrapper_lib.acqdatetime_series(dcmInfile)

    results.addDateTime('AcquisitionDateTime', dt) 

def header_series(data, results, action, idname=None):
    """
    Read selected dicomfields and write to IQC database

    Workflow:
        1. Set runtime parameters
        2. Check data format
        3. Build and populate qcstructure
        4. Run tests
        5. Build xml output
    """
    # 1.-3.
    try:
        params = action['params']
    except KeyError:
        params = {}
    inputfile = data.series_filelist[0]  # give me a [filename]
    qclib,cs = setup_series(inputfile, params, headers_only=True, for_action='hdrs')
    ## find probename
    if cs.auto_suffix:
        if idname is None:
            idname = '_'+qclib.imageID(cs,probeonly=True)
    else:
        idname = ''
    cs.resultslabel = idname[1:]

    # 4. Run tests
    dicominfo = qclib.DICOMInfo(cs)

        
    ## 2. Add results to 'result' object
    varname = 'pluginversion'+idname
    results.addString(varname, str(qclib.qcversion))
    for di in dicominfo:
        varname = di[0]+idname
        results.addString(varname, str(di[1])[:min(len(str(di[1])),100)])
    

def OCR(data, results, action, idname):
    """
    Use pyOCR which for OCR
    returns rect rois for plotting in overview
    If results is None, just return the OCR content, do not add to results
    """
    try:
        params = action['params']
    except KeyError:
        params = {}

    # optional parameters
    ocr_options = {}
    for lab in ['ocr_threshold', 'ocr_zoom', 'ocr_border']:
        if lab in params:
            ocr_options[lab] = int(params[lab])

    inputfile = data.series_filelist[0]  # give me a [filename]
    qclib,cs = setup_series(inputfile, params, headers_only=False, for_action='ocr')
    
    if cs.auto_suffix:
        if idname is None:
            idname = '_'+qclib.imageID(cs, probeonly=True)
    else:
        idname = ''
        
    rectrois = []
    error = False
    msg = ''
    values = {}
    # solve ocr params
    ocr_regions = params.get('ocr_regions',{}) # new format

    regions = {}
    for ocrname,ocrparams in ocr_regions.items():
        regions[ocrname] = {'prefix':'', 'suffix':''}
        for key,val in ocrparams.items():
            if key == 'xywh':
                regions[ocrname]['xywh'] = [int(p) for p in val.split(';')]
            elif key == 'prefix':
                regions[ocrname]['prefix'] = val
            elif key == 'suffix':
                regions[ocrname]['suffix'] = val
            elif key == 'type':
                regions[ocrname]['type'] = val

    for name, region in regions.items():
        rectrois.append([ (region['xywh'][0],region['xywh'][1]), 
                          (region['xywh'][0]+region['xywh'][2],region['xywh'][1]+region['xywh'][3])])

        txt, part = ocr_lib.OCR(cs.pixeldataIn, region['xywh'], **ocr_options)
        uname = name+str(idname)
        if region['type'] == 'object':
            im = toimage(part) 
            fn = '%s.jpg'%uname
            im.save(fn)
            results.addObject(uname, fn)
            
        else:
            try:
                value = ocr_lib.txt2type(txt, region['type'], region['prefix'], region['suffix'])
                if not results is None:
                    if region['type'] == 'float':
                        results.addFloat(uname, value)
                    elif region['type'] == 'string':
                        results.addString(uname, value)
                    elif region['type'] == 'bool':
                        results.addBool(uname, value)
                else:
                    values[uname] = value
            except:
                print("error", uname, value)
                error = True
                msg += uname + ' '
                im = toimage(part) 
                fn = '%s.jpg'%uname
                im.save(fn)
                

    if results is None:
        return values, error, msg

    return rectrois, error, msg

def writeimages(cs, ocr_rois, idname):
    # also run ocr_series; needed as part of qc because of the boxes it generates
    xtra= {'rectrois': ocr_rois }
    qclib = QCUS_lib.US_QC(guimode=False)

    if cs.auto_suffix:
        if idname is None:
            idname = '_'+qclib.imageID(cs,probeonly=True)
    else:
        idname = ''
    cs.resultslabel = idname[1:]

    fname = 'overview%s.jpg'%str(idname)
    qclib.saveAnnotatedImage(cs, fname, what='overview',xtra=xtra)
    results.addObject(os.path.splitext(fname)[0],fname)
    
if __name__ == "__main__":
    data, results, config = pyWADinput()

    instances = data.getAllInstances()
    if len(instances) != 1:
        print('%s Error! Number of instances not equal to 1 (%d). Exit.'%(logTag(),len(instances)))

    error = False
    msg = ''
    cs = None
    # read runtime parameters for module
    idname = None

    if 'ocr_series' in config['actions'].keys():
        idname = get_idname_from_ocr(data, config['actions']['ocr_series']['params'])

    for name,action in config['actions'].items():
        if name == 'acqdatetime':
            acqdatetime_series(data, results, action)

        elif name == 'header_series':
            header_series(data, results, action, idname)
        
        elif name == 'qc_series':
            cs = qc_series(data, results, action, idname)

        elif name == 'ocr_series':
            ocr_rois, error, msg = OCR(data, results, action, idname)

    #label = instance.DeviceSerialNumber+'__'+''.join(instance.TransducerData).strip()
    writeimages(cs, ocr_rois, idname)
    
    results.write()

    if error:
        raise ValueError('%s Cannot read OCR box for %s'%(logTag(),msg))
    