# -*- coding: utf-8 -*-
"""
Author: Anjali Balagopal, University of Texas Southwestern Medical center
Created on 09-17-2018
"""

import numpy as np
import cv2
import glob
import os
import time
import json
import math
import pydicom as dicom

def get_dir(data_dir, S):
    file_dirs = list()
    if S == 'CT':
        for subject in glob.glob(os.path.join(data_dir, 'CT/')):
            file_dirs.append(subject)
    if S == 'RT':
        for subject in glob.glob(os.path.join(data_dir,'RT/')):
            file_dirs.append(subject)
    return file_dirs

def extract_study(dcm_path, S):
    ######################################
    # Description: Read all dicom files and extract the study and the patient id
    #
    # parameter:
    # 	dcm_path: the path to the dicom folder from where the files have to be extracted
    #
    # return:
    #	study: dict[]
    #   patient id
    ######################################
    study = {}
    dcm_files = os.listdir(dcm_path)
    # if S == 'RT':
        # print(dcm_files)
    for dcm_file in dcm_files:
        data = dicom.read_file(os.path.join(dcm_path, dcm_file))
        study_id = data.StudyID
        patient_id = data.PatientID

        if not study_id in study:
            if S == 'CT':
                study[study_id] = {'CT': []}
            if S == 'RT':
                study[study_id] = {'RT': []}

        if S == 'CT':
            study[study_id]['CT'].append(os.path.join(dcm_path, dcm_file))
        if S == 'RT':
            if (dcm_file.startswith('RT') or dcm_file.startswith('RS')):
                study[study_id]['RT'].append(os.path.join(dcm_path, dcm_file))
        else:
            continue

    count = 0
    del_list = []
    for idx in study:
        empty = False
        for item in study[idx]:
            if study[idx][item] == []:
                empty = True
                break
        if not empty:
            count += 1
        else:
            del_list.append(idx)

    for idx in del_list:
        del study[idx]

    return study, patient_id

def extract_ct(ct_dcm_path_list, dst_path, patient_ID_CT):
    ######################################
    # parameter:
    # 	ct_dcm_path_list: the list of CT .dcm files
    #	dst_path: dst path to save the .npy files
    # 	resize_: int, will resize the numpy array to [resize_, resize_]
    #
    # return:
    #	offset: [offset_start_slice, offset_last_slice]
    #	spacing: [float, float], spacing in x and y axis
    #	order: list, use order to match each slice in a pair of CT and mask
    ######################################

    ct_dcms = ct_dcm_path_list
    slice_num = len(ct_dcms)
    offset_start = [0, 0, -1.0]
    offset_end = [0, 0, -1.0]
    order_slice = {}
    spacing = [1.0, 1.0]
    # print('There are {} slices in the patient!'.format(slice_num))

    start_time = time.time()
    ct_slice_1 = dicom.read_file(os.path.join(ct_dcms[0]))
    ct_all_slices = np.zeros((slice_num, ct_slice_1.Rows, ct_slice_1.Columns), dtype=float)

    ori_size = [ct_slice_1.Rows, ct_slice_1.Columns]

    tmp_all_slices = {}
    start_slice = 1000
    end_slice = -1
    for i in range(slice_num):
        ct_slice_path = os.path.join(ct_dcms[i])
        ct_slice_data = dicom.read_file(ct_slice_path)
        current_slice = np.zeros((ct_slice_data.Rows, ct_slice_data.Columns), dtype=float)

        if i == 1 and (not (ct_slice_data.Columns == 512 and ct_slice_data.Rows == 512)):
            print('The rows&columns of\'{}\' are not 512 !'.format(ct_slice_path))
            ori_size = [ct_slice_data.Rows, ct_slice_data.Columns]

        current_idx = ct_slice_data.SliceLocation
        order_slice[current_idx] = ct_slice_data.ImagePositionPatient[2]
        if current_idx < start_slice:
            start_slice = current_idx
            offset_start = [float(item) for item in ct_slice_data.ImagePositionPatient]
        if current_idx > end_slice:
            end_slice = current_idx
            offset_end = [float(item) for item in ct_slice_data.ImagePositionPatient]

        if i == 0:
            spacing = [float(ct_slice_data.PixelSpacing[0]), float(ct_slice_data.PixelSpacing[1])]

        pixel_data = ct_slice_data.PixelData
        intercept = np.float32(ct_slice_data.RescaleIntercept)
        slop = np.float32(ct_slice_data.RescaleSlope)
        SliceThickness = np.float32(ct_slice_data.SliceThickness)
        byte_num = len(pixel_data) // (ct_slice_data.Rows * ct_slice_data.Columns)
        dt = np.uint8

        if byte_num == 2:
            dt = np.int16
        elif byte_num == 4:
            dt = np.int32
        elif byte_num == 8:
            dt = np.int64

        current_slice = np.copy(
            np.reshape(np.frombuffer(pixel_data, dtype=dt, count=ct_slice_data.Rows * ct_slice_data.Columns),
                       current_slice.shape).astype(np.float32))

        current_slice = current_slice * slop + intercept
        tmp_all_slices[current_idx] = current_slice
    # sort
    sorted_tmp_all_slices = sorted(tmp_all_slices)
    idx = 0
    for ii in sorted_tmp_all_slices:
        ct_all_slices[idx, :, :] = tmp_all_slices[ii]
        idx += 1
    # print('{}/{} slices have been processed!'.format(slice_num, slice_num))
    ct_all_slices = np.transpose(ct_all_slices, (1, 2, 0))
    # save ct slice
    np.save(os.path.join(dst_path, 'CT_{}_0.npy'.format(patient_ID_CT)), ct_all_slices)
    order = [float(order_slice[key]) for key in sorted(order_slice)]
    offset = {'start': offset_start, 'end': offset_end}

    # during = math.floor(time.time() - start_time)
    # print("[Time for saving CT ]: {} hrs {} mins {} secs".format(during // 3600, (during % 3600 // 60), (during % 3600 % 60)))

    return offset, spacing, order, ori_size, SliceThickness

def create_structures(rt_struct_rois,roi,contours_sequence, slice_num, filled=True):
    current_roi_id = rt_struct_rois.index(roi)
    contour_sequence = contours_sequence[current_roi_id].ContourSequence
    current_mask_volume = np.zeros((slice_num, ori_size[0], ori_size[1]), dtype=float)

    for i in range(len(contour_sequence)):
        points_num = int(contour_sequence[i].NumberOfContourPoints)
        contour_data = contour_sequence[i].ContourData
        tmp_contour = [[[int((contour_data[j * 3 + 0] - offset.get('start')[0]) / spacing[0]),
                         int((contour_data[j * 3 + 1] - offset.get('start')[1]) / spacing[1])] for j in
                        range(points_num)]]
        tmp_contour = np.asarray(tmp_contour)

        tmp_distance = [abs(float(contour_data[2]) - v) for v in order]
        tmp_idx = tmp_distance.index(min(tmp_distance))
        fill_flag = -1
        if not filled:
            fill_flag = 3
        result = cv2.drawContours(np.zeros((ori_size[0], ori_size[1])), tmp_contour, -1, 1.0, fill_flag)
        current_mask_volume[tmp_idx, :, :] = current_mask_volume[tmp_idx, :, :] + result
    current_mask_volume[current_mask_volume >= 0.5] = 1.0
    current_mask_volume[current_mask_volume < 0.5] = 0.0
    current_mask_volume = np.transpose(current_mask_volume, (1, 2, 0))
    return current_mask_volume

def extract_contours(rt_struct_path, dst_path, patient_id_RT, ctv_, bladder_, rectum_, order):
    ######################################
    # parameter:
    # 	rt_struct_path: RT .dcm files
    #	dst_path: dst path to save the .npy files
    # 	offset: [offset_start_slice, offset_last_slice]
    #	spacing: [float, float], spacing in x and y axis
    #	order: list, use order to match each slice in a pair of CT and mask
    #	ori_size: [int, int] generally [512*512], in some cases, the size of ROI will be [272*272]
    #	roi_lists ( ctv_, bladder_, rectum_) : list, the list of original ROIs that need to save as numpy arrays, e.g. ['Bladder', 'CTV']
    #
    # return:
    #	None
    ######################################
    # print('********************************************* Extracting ROI masks for Patient {} ****************************************************'.format(patient_id_RT))
    start_time = time.time()

    # extract contours and save as *.png
    rt_struct_data = dicom.read_file(rt_struct_path)
    rt_struct_rois = [rt_struct_data.StructureSetROISequence[ii].ROIName for ii in
                      range(len(rt_struct_data.StructureSetROISequence))]
    # print(rt_struct_rois)

    slice_num = len(order)

    contours_sequence = []
    if rt_struct_data.dir('contour')[0]:
        contours_sequence = rt_struct_data.ROIContourSequence
    else:
        print('ERROR: There is no contour in the file!')
        return

    roi_list = rt_struct_rois

    #### Finding Structures using keywords
    for idx in range(len(roi_list)):
        if not bladder_:
            if (roi_list[idx].lower() == 'bladder') or (roi_list[idx].lower() == 'bladderfull'):
                bladder_.append(roi_list[idx])
        if not bladder_:
            if (roi_list[idx].lower() == 'rectum'):
                rectum_.append(roi_list[idx])
        if not ctv_:
            if (roi_list[idx].lower() == 'CTV'):
                rectum_.append(roi_list[idx])

    current_mask_volume = create_structures(rt_struct_rois,ctv_list,contours_sequence, slice_num, filled=True)
    np.save(os.path.join(dst_path, 'ROI_{}_0_7.npy'.format(patient_id_RT)), current_mask_volume.astype(np.bool))   ## However you want to save it
    print('{}....... saved!'.format(ctv_list))

    if(len(bladder_list)== 0):
       print("Cannot find any Bladder structure")
    elif(len(bladder_list) > 1):
       print("Found muliple Bladder structures, find the right structure from :", bladder_list, " Line 231")
    else:
        current_mask_volume = create_structures(rt_struct_rois, bladder_list[0], contours_sequence, slice_num, filled=True)
        np.save(os.path.join(dst_path, 'ROI_{}_0_5.npy'.format(patient_id_RT)), current_mask_volume.astype(np.bool)) ## However you want to save it
        print('{}....... saved!'.format(bladder_list[0]))

    if(len(rectum_list)== 0):
       print("Cannot find any Rectum structure")
    elif(len(rectum_list) > 1):
       print("Found muliple Rectum structures, find the right structure from :", rectum_list, " Line 233")
    else:
       current_mask_volume = create_structures(rt_struct_rois, rectum_list[0], contours_sequence, slice_num,
                                               filled=True)
       np.save(os.path.join(dst_path, 'ROI_{}_0_4.npy'.format(patient_id_RT)), current_mask_volume.astype(np.bool)) ## However you want to save it
       print('{}....... saved!'.format(rectum_list[0]))


##################### dicom 2 numpy ##########################

if __name__ == '__main__':
    with open('config_dicomtonumpy.json', 'r') as f:
        cfg = json.load(f)
    numpy_dst_dir = cfg['numpy_dst_dir']
    data_dir = cfg['data_dir']


    ########## If you already know the contour names for each structure to be extracted use this section of the code to create a list.
    ########## I have created an excel sheet with all names and is reading from there.
    ########## Another option which is included in the code is to search using keywords which is coded in extract_contours() function

    from xlrd import open_workbook

    wb = open_workbook('pat_list.xlsx')
    pat_list = []
    alias_list = []
    ctv_list = []
    bladder_list = []
    rectum_list = []
    for s in wb.sheets():
        if s.name == 'contournames':
            num_cols = s.ncols  # Number of columns
            for row_idx in range(0, s.nrows):  # Iterate through rows
                pat_list.append(str.zfill(str(int(s.cell(row_idx, 0).value)),10))
                alias_list.append(str(int(s.cell(row_idx, 1).value)))
                ctv_list.append(str(s.cell(row_idx, 2).value))
                bladder_list.append(str(s.cell(row_idx, 3).value))
                rectum_list.append(str(s.cell(row_idx, 4).value))


    import shutil
    for p in range(len(pat_list)):
        print("***************************************** DICOM TO NUMPY CONVERSION FOR PATIENT {} ***********************************************".format(pat_list[p]))
        print(alias_list[p])
        CT_file_dirs = get_dir(dicom_dir, 'CT')
        RT_file_dirs = get_dir(dicom_dir, 'RT')

        for CT_file_path in CT_file_dirs:
            RT_file_path = RT_file_dirs[CT_file_dirs.index(CT_file_path)]

            study_CT, patient_id_CT = extract_study(CT_file_path, 'CT')
            study_RT, patient_id_RT = extract_study(RT_file_path, 'RT')

            if not study_CT:
                print('There are no dicoms in this folder!')

            # load the first study, or load all studies by using a loop
            study_id = list(study_CT.keys())[0]

            # CT dicom to numpy
            offset, spacing, order, ori_size, SliceThickness = extract_ct(study_CT[study_id]['CT'],numpy_dst_dir,alias_list[p])
            #################################################################
            # mask dicom to numpy
            rt_path = study_RT[study_id]['RT'][0]
            rt_data = dicom.read_file(rt_path)
            extract_contours(rt_path, numpy_dst_dir,alias_list[p], ctv_list[p], bladder_list[p], rectum_list[p], order)

