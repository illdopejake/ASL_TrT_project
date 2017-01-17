import os, shutil
from glob import glob
from copy import copy
import numpy as np
import nibabel as nib
import nipype.interfaces.fsl.utils as fsl
from scipy.io import loadmat, savemat

#####PLEASE NOTE########
#--Script only works with nii images
#
#--Script relies somewhat on niak naming convention. May require some 
#slight modification for other naming conventions, or an update to a 
#more flexible system
########################


# path to directory where scans are located
# NOTE: In its current iteration, if directory has any files within
# that are not named using niak conventions, script will fail 
scan_dir = '/Users/jakevogel/bellec_lab/asl_project/scans_new_niak/nifti'

# search string for asl scans in directory.
scn1_str = '_sess2_asl1.nii.gz'
scn2_str = '_sess2_asl2.nii.gz'

# path to _extra.mat files. Can be the same as scan_dir
matdir = scan_dir

# Path to directory you wish output directory to be. It will be
# created if it doesn't already exist.
# NOTE: ONLY SET OUTPUT DIR AS SCAN (INPUT) DIR IF YOU WANT YOUR SCAN
# AND MAT FILES TO BE OVERWRITTEN!!!!
output_dir = os.path.join(scan_dir,'test_delete')

# When files are recompiled after removing frames, a TR must be specified
tr = 4.22

# If scans are preprocessed using any NIAK version before 0.17.0, set 
# this variable to 'old'. Otherwise, if using any newer verison of NIAK
# (for example, NIAK 1), set this variable to 'new'
version = 'new'

def find_true_index_of_delframes(mat,todel,version='old'):
    """ given a dict indicating which frames are missing (todel) and a
    .mat file matching the correct scan, will output a new dict that
    accounts for the fact that removed frames will change the position
    of frames removed later in sequence. Essentially repositions index
    of frames-to-remove to account for earlier-removed frames.

    Example:
    Scan 2 is missing frames 2-5. If del = {16-19} todel will change to
    {12:4}, because four frames were already removed, 12 is the correct
    index to match the volume to be removed.

    """

    fail = False

    # identified all missing frames
    if version == 'new':
        ms = mat['mask_scrubbing']
    else:
        ms = mat['mask_suppressed']
    s_frames = []
    for i,supp in enumerate(ms.tolist()):
        if supp[0] == 1:
            s_frames.append(i)

    # make sure this worked
    if len(mat['time_frames'][0]) != (len(ms) - len(s_frames)):
        print 'the indentification of missing frames failed somehow.'
        failed = True

    # figure out correct position of frames-to-remove
    new_del = {}
    s_dict = redictify_frames_list(s_frames)


    for k,v in todel.iteritems():
        subber = 0
        for x,y in s_dict.iteritems():
            if k > x:
                subber = subber + y
        new_del.update({k-subber: v})

    for chk in new_del.keys():
        if chk < 0:
            fail = True

    return new_del, fail

def compare_missing_frames(mat1,mat2,version='old'):
    """ given two _extra.mat files, will ascertain where there are
    discrepancies in frames between the two scans. Outputs two
    dictionaries containing pairs of start and stop indexes indicating
    where frames must be removed to harmonize scans


    Example: Scan 1 is missing frames 2-5 and Scan 2 is missing frames
    2-5 and 16-19.
        Output:
        del1 = {} -- because no frames missing in Scan 1 were missing in Scan 2
        del2 = {16:4} -- indicates four consecutive frames, starting from
                         position 16, are missing from Scan 2, but not Scan 1.

    NOTE: If using any NIAK version before 0.17.0, leave version set to 'old'.
    If using a later version of NIAK, like NIAK 1 for example, set version to
    'new'

    """

    fail = False

    if version == 'new':
        ms1 = mat1['mask_scrubbing']
        ms2 = mat2['mask_scrubbing']
    else:
        ms1 = mat1['mask_suppressed']
        ms2 = mat2['mask_suppressed']

    # browse mask_suppressed of both matfiles to determine which frames
    # 							to remove
    del1 = []
    del2 = []

    if len(ms1) != len(ms2):
        fail = True
    else:
        for i,supp in enumerate(ms1):
            if supp != ms2[i]:
                if supp[0] == 1:
                    del2.append(i)
                else:
                    del1.append(i)
        if len(del1) > 0: # turn frame lists into dicts for later
            del1 = redictify_frames_list(del1)
        if len(del2) > 0:
            del2 = redictify_frames_list(del2)

    return del1, del2, fail


def redictify_frames_list(frameslst):
    """turn frames list into a dictionary where key = index of frame
    to remove and value = # of consecutive frames to remove starting
    from key

    Example: Scan 1 is missing frames 2-5 and Scan 2 is missing frames
    2-5 and 16-19.
        Output:
        del1 = {} -- because no frames missing in Scan 1 were missing in Scan 2
        del2 = {16:4} -- indicates four consecutive frames, starting from
                         position 16, are missing from Scan 2, but not Scan 1.
    """

    frame_dict = {}
    used = []

    for ind in frameslst:
        cntr = copy(ind)
        if ind not in used:
            count = 1
            while cntr+1 in frameslst:
                used.append(cntr + 1)
                count = count + 1
                cntr = cntr + 1
            frame_dict.update({ind: count})

    return frame_dict

def remove_dupes(lst):
    for itm in lst:
        if lst.count(itm) > 1:
            for i in range(lst.count(itm)-1):
                lst.remove(itm)

    return lst

def remove_frames(scan,todel,matfile,tr,order=''):
    """ given a 4d scan, a list of frames to remove, and an _extra.mat file,
    will remove indicated frames from the time dimension of scan and update
    _extra.mat file to reflect changes.

    tr should reflect the actual TR of the acquisition data.

    order indicates which scan is being edited for recording purposes only
    and can be left blank

    """

    # turn todel dict back into a list
    to_remove = []
    for k in sorted(todel.keys()):
         for i in range(k,k+todel[k]):
             to_remove.append(i)

    to_remove = remove_dupes(to_remove) # remove duplicates in list of frames)    

    print 'removing the following frames from scan %s: %s'%(order,to_remove)
    print 'splitting scan into volumes and removing frames'

    # split 4D image into volumes
    split = fsl.Split()
    split.inputs.in_file = scan
    split.inputs.dimension = 't'
    split.inputs.output_type = 'NIFTI_GZ'
    split.run() # outputs frames in same directory as volxxxx.nii.gz 

    # remove desired frames
    for rem_ind in to_remove:
        to_del = os.path.join(scan_dir,'vol%04d.nii.gz'%(rem_ind))
        os.remove(to_del)


    # merge remaining frames back into single image
    jnk,flnme = os.path.split(scan)


    vols = sorted(glob(os.path.join(scan_dir,'vol*')))

    merge = fsl.Merge()
    merge.inputs.in_files = vols
    merge.inputs.dimension = 't'
    merge.inputs.tr = tr
    merge.inputs.merged_file = os.path.join(output_dir, flnme)
    merge.inputs.output_type = 'NIFTI_GZ'
    merge.inputs.ignore_exception

    print 'creating new file %s'%(os.path.join(output_dir,flnme))
    merge.run()

    print 'updating _extra.mat file'

    # update _extra.mat file
    update_extra_mat(matfile,to_remove)

    print 'cleaning up'

    # clean up
    vols = sorted(glob(os.path.join(scan_dir,'vol*')))
    for vol in vols:
        os.remove(vol)


def update_extra_mat(matfile,to_remove):
    """ updates the time_frames, confounds and mask_suppressed arrays to
    reflect the removed volumes. However, does not change other items in
    _extra.mat file

    """

    mat = loadmat(matfile)
    # update time_frames
    ntf = np.delete(mat['time_frames'][0],to_remove)
    mat.update({'time_frames': ntf})

    # update confounds
    ncon = np.delete(mat['confounds'],to_remove,axis = 0)
    mat.update({'confounds': ncon})

    # update mask_suppressed
    ms = mat['mask_suppressed']
    for supp in to_remove:
        ms[supp][0] = 1
    mat.update({'mask_suppressed': ms})

    # save updated mat file
    jnk, flnme = os.path.split(matfile)
    savemat(os.path.join(output_dir,flnme),mat)

def check(scan1,scan2):
    """ performs an extremely crude check on the whether the script
    actually harmonized the two scans, by comparing the shape of the
    two scans to see if they are equal. This is obviously not very
    conclusive.

    """

    fail = False

    # load scans
    jnk,nme1 = os.path.split(scan1)
    scn1 = nib.load(os.path.join(output_dir,nme1))

    jnk,nme2 = os.path.split(scan2)
    scn2 = nib.load(os.path.join(output_dir,nme2))

    if scn2.shape != scn1.shape:
        fail = True

    return fail

def harmonize_mask_scrubbing(matfile1,matfile2):

    fail = False
    mat1 = loadmat(matfile1)
    mat2 = loadmat(matfile2)
    ms1 = mat1['mask_scrubbing']
    ms2 = mat2['mask_scrubbing']
    new_ms = []
    if len(ms1) != len(ms2):
        fail = True
    for i in range(len(ms1)):
        if ms1[i][0] == 1. or ms2[i][0] == 1.:
            new_ms.append([1.])
        else:
            new_ms.append([1.])
    mat1.update({'mask_scrubbing': new_ms})
    mat2.update({'mask_scrubbing': new_ms})
    jnk, flnme = os.path.split(matfile1)
    savemat(os.path.join(output_dir,flnme),mat1)
    jnk, flnme = os.path.split(matfile2)
    savemat(os.path.join(output_dir,flnme),mat2)

    return fail

if __name__ == '__main__':
    """ This script is designed to harmonized the volumes of two paried
    (dual echo) asl scans, such that the volumes removed from each scan
    are identical. The pursose is to allow for operations such as
    surround averaging/subtraction.

    While it was designed for asl data, it should work fine with any type
    4D niak data

    The script will
        i) assess which absolute frames must be removed from each of the
        two scans by comparing 'mask_suppressed' from the _extra.mat file

        ii) figure out the true index of volumes to be removed with
        respect to the current (scrubbed) list of volumes

        iii) remove volumes from both scans so that they share the exact
        same list of volumes

        vi) update components of the _extra.mat file to reflect these
        changes

    #####PLEASE NOTE########
    --Script only works with nii images

    --Script relies somewhat on niak naming convention. May require
    some slight modification for other naming conventions,
    or an update to a more flexible system
    ########################

    """

    cwd = os.getcwd()
    os.chdir(scan_dir)

    fail_log = {} # to keep track of subjects that fail. should make it write to a log

    processed = [] # to avoid processing a subject twice    

    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

    # grab subjects and data
    allscnz = sorted(glob(os.path.join(scan_dir,'*.ni*')))
    for strng in allscnz:

        pth,nme = os.path.split(strng)
        jnk,sid,_ = nme.split('_s')

        if sid in processed:
            continue
        else:
            processed.append(sid)

        print 'processing subject %s'%(sid)

        scan1 = os.path.join(pth,'%s_s%s%s'%(jnk,sid,scn1_str))
        scan2 = os.path.join(pth,'%s_s%s%s'%(jnk,sid,scn2_str))

        jnk,txt = os.path.split(scan1)
        if scan1[-2:] == 'gz':
            premat,ex1,ex2 = txt.split('.')
        else:
            premat,ex1 = txt.split('.')
        matfile1 = os.path.join(matdir, '%s_extra.mat'%(premat))
        mat1 = loadmat(matfile1)

        jnk,txt = os.path.split(scan2)
        if scan2[-2:] == 'gz':
            premat,ex1,ex2 = txt.split('.')
        else:
            premat,ex1 = txt.split('.')
        matfile2 = os.path.join(matdir, '%s_extra.mat'%(premat))
        mat2 = loadmat(matfile2)


        if version == 'old':
            # determine which absolute frames must be removed from each scan
            print 'determining which frames to remove'
            del1,del2,fail = compare_missing_frames(mat1,mat2,version)
            if fail:
                print 'Something went wrong while determining absolute frames to remove'
                fail_log.update({sid: 'absolute_frames'})
                continue

            if len(del1) < 1 & len(del2) < 1:
                print 'Scans already harmonized, moving on'
                continue

            # determine true index of volumes to remove  
            if len(del1) > 0:
                del1, fail = find_true_index_of_delframes(mat1,del1,version)
                if fail:
                    print 'Something went wrong while determining true frames to remove'
                    fail_log.update({sid: 'true frames, scan1'})

                # remove frames from image
                remove_frames(scan1,del1,matfile1,tr,order=1)
            else:
                print 'no deletions necessary for scan 1'
                shutil.copy(scan1,output_dir)


            if len(del2) > 0:
                del2, fail = find_true_index_of_delframes(mat2,del2)
                if fail:
                    print 'Something went wrong while determining true frames to remove'
                    fail_log.update({sid: 'true_frames, scan2'})

                # remove frames from image
                remove_frames(scan2,del2,matfile2,tr,order=2)
            else:
                print 'no deletions necessary for scan 2'
                shutil.copy(scan2,output_dir)

            # loosely double-check if scans are harmonized by inquiring their shape
            print 'checking my work'
            fail = check(scan1,scan2)
            if fail:
                print 'Something went wrong with the check-test.'
                fail_log.update({sid: 'check test failed'})
            else:
                print 'everything looks good for subject %s, wrapping up'%(sid)

        elif version == 'new':
            fail = harmonize_mask_scrubbing(matfile1,matfile2)
            if fail:
                print 'Something went wrong with the harmonizing the extra_mat file'
            else:
                shutil.copy(scan1,output_dir)
                shutil.copy(scan2,output_dir)
                print 'finished harmonizing subject %s'%(sid)

    print 'here is a list of failed subjects and reasons'
    print fail_log
    os.chdir(cwd)
