import os, shutil
from glob import glob
from copy import copy
import numpy as np
import nibabel as nib
import nipype.interfaces.fsl.maths as math
from scipy.io import loadmat, savemat


# path to directory where scans are located
# NOTE: In its current iteration, if directory has any files within
# that are not name using niak conventions, script will fail
scan_dir = '/Users/jakevogel/bellec_lab/asl_project/processed'

# search string for asl scans in directory.
scn1_str = '_sess2_asl1.nii.gz'
scn2_str = '_sess2_asl2.nii.gz'

# Path to directory you wish output directory to be. It will be
# created if it doesn't already exist.
# NOTE: ONLY SET OUTPUT DIR AS SCAN (INPUT) DIR IF YOU WANT YOUR SCAN
# AND MAT FILES TO BE OVERWRITTEN!!!!
output_dir = os.path.join(scan_dir,'cbf')

mask = '/Users/jakevogel/bellec_lab/asl_project/bin_no.nii.gz'

# change to true if you want script to create cbf image
do_cbf = False

# change to true if you want script to smooth data
do_smth = True

# enter smoothing kernel here. Must be a float
fwhm = 6.0

def make_control_constant_image(scan_lst):
    cmnd = 'fslmaths '
    scan_len = len(scan_lst)
    ctrl_mean = os.path.join(output_dir,'mean_ctrl_img')
    for scn in scan_lst:
        if scn == scan_lst[-1]:
            cmnd = cmnd+'%s -div %s %s'%(scn,scan_len,ctrl_mean)
        else:
            cmnd = cmnd+'%s -add '%(scn)

    os.system(cmnd)

    ctrl_mean = ctrl_mean+'.nii.gz'

    return ctrl_mean

def define_denominator(ctrl_mean):
    alpha = 0.85 # assumed inversion efficiency of pCASL
    T1b = 1.664 #T1 of blood, in s
    w = 0.900 # post labeling delay in s
    tau = 1.5 # labeling duration in s
    e = 2.71828 # euler's number

    term3 = e**((-w)/T1b) - e**((-(tau + w)) / T1b)
    factor = (2*alpha) * T1b * term3

    denom = os.path.join(output_dir,'denom')

    mul = math.BinaryMaths()
    mul.inputs.in_file = ctrl_mean
    mul.inputs.operand_value = factor
    mul.inputs.operation = 'mul'
    mul.inputs.out_file = denom
    mul.inputs.output_type = 'NIFTI_GZ'
    mul.inputs.ignore_exception = True
    mul.run()

    denom = denom+'.nii.gz'

    return denom

def surround_avg(scan1,scan2):

    jnk,flnm = os.path.split(scan1)
    subn,jnk = flnm.split('_sess')
    savgfl = os.path.join(output_dir,'%s_surroundavg'%(subn))

    divisor = nib.load(scan1).get_data().shape[3]

    immth = math.BinaryMaths()
    immth.inputs.in_file = scan1
    immth.inputs.operand_file = scan2
    immth.inputs.operation = 'sub'
    immth.inputs.out_file = savgfl
    immth.inputs.output_type = 'NIFTI_GZ'
    immth.inputs.ignore_exception = True
    immth.run()

    savgfl = savgfl+'.nii.gz'

    return savgfl

def convert_2_cbf(savgfl,denom):

    jnk,flnm = os.path.split(savgfl)
    subn,jnk = flnm.split('_sur')
    cbf = os.path.join(output_dir,'%s_cbf'%(subn))

    l = 0.9 #blood/tissue water partition coefficient, in g/ml
    num = os.path.join(output_dir,'d_num')
    mul = math.BinaryMaths()
    mul.inputs.in_file = savgfl
    mul.inputs.operand_value = l
    mul.inputs.operation = 'mul'
    mul.inputs.out_file = num
    mul.inputs.output_type = 'NIFTI_GZ'
    mul.inputs.ignore_exception = True
    mul.run()

    num = num+'.nii.gz'

    immth = math.BinaryMaths()
    immth.inputs.in_file = num
    immth.inputs.operand_file = denom
    immth.inputs.operation = 'div'
    immth.inputs.out_file = cbf
    immth.inputs.output_type = 'NIFTI_GZ'
    immth.inputs.ignore_exception = True
    immth.run()

    cbf = cbf+'.nii.gz'

    os.system('rm %s'%(os.path.join(output_dir,'d_*.nii.gz')))

    return cbf

def mask_img(img,mask):

    msk = math.ApplyMask()
    msk.inputs.in_file = img
    msk.inputs.mask_file = mask
    msk.inputs.out_file = img
    msk.inputs.output_type = 'NIFTI_GZ'
    msk.inputs.ignore_exception = True
    msk.run()

def smooth_img(img,fwhm):

    smth = math.IsotropicSmooth()
    smth.inputs.in_file = img
    smth.inputs.fwhm = fwhm
    smth.inputs.out_file = img
    smth.inputs.output_type = 'NIFTI_GZ'
    smth.inputs.ignore_exception = True
    smth.run()

if __name__ == '__main__':
    """This script will perform surround averaging on pairs of asl images (tag and control),
	assuming the first scan of the pair is the tag scan. In addition, if specified, a cbf
    image will be created which converts voxel values to physiological units: ml/100g/min

	CBF equation is a single-blood-compartment model, borrowed from Chen et al., 2011 Journal
	of Magnetic Resonance Imaging
	"""

    if not os.path.isdir(output_dir):
        os.system('mkdir %s'%(output_dir))

    allscnz = sorted(glob(os.path.join(scan_dir,'*%s'%(scn2_str))))
    os.system('mkdir %s'%(output_dir))

    if do_cbf:
        print 'creating equilibrium magnetization image...'
        m0 = make_control_constant_image(allscnz)

        print 'preparing model...'
        denom = define_denominator(m0)

    for scan2 in allscnz:
        pref,jnk = scan2.split(scn2_str[:5])
        scan1 = pref+scn1_str
        jnk,subn = os.path.split(pref)

        print 'performing surround averaging for subject %s'%(subn)
        savgfl = surround_avg(scan1,scan2)

        if do_cbf:
            print 'calibrating cbf can to physiological units for subject %s'%(subn)
            cbf = convert_2_cbf(savgfl,denom)

            print 'masking image'
            mask_img(cbf,mask)

            print 'created new file %s'%(cbf)

            os.system('rm %s'%(os.path.join(output_dir,'denom*')))

            if do_smth:
                print 'smoothing image %s'%(cbf)
                smooth_img(cbf,fwhm)

        else:
            print 'masking image'
            mask_img(savgfl,mask)
            print 'created new file %s'%(savgfl)

            if do_smth:
                print 'smoothing image %s'%(savgfl)
                smooth_img(savgfl,fwhm)
