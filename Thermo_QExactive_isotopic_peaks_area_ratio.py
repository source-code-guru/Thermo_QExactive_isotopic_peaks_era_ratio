"""
A python script which extracts masses intensities and the mono/first isotopic peaks integrals ratios for QExactive Orbitrap LC/MS-MS Thermo Fisher .raw files
Author : Marc HAEGELIN
Date : 02-10-2025
"""

from __future__ import print_function

import clr
import sys
import time
import os
import numpy as np
import pandas as pd
from scipy.optimize import leastsq
import scipy.integrate as integrate
import itertools
import subprocess
import re
import pkg_resources

clr.AddReference('System')
clr.AddReference('System.Collections')

clr.AddReference('ThermoFisher.CommonCore.Data')
clr.AddReference('ThermoFisher.CommonCore.RawFileReader')
clr.AddReference('ThermoFisher.CommonCore.BackgroundSubtraction')
clr.AddReference('ThermoFisher.CommonCore.MassPrecisionEstimator')

from ThermoFisher.CommonCore.Data.Business import ChromatogramSignal, ChromatogramTraceSettings, DataUnits, Device, Scan, TraceType
from ThermoFisher.CommonCore.Data.FilterEnums import MSOrderType
from ThermoFisher.CommonCore.Data.Interfaces import IChromatogramSettings, IScanEventBase, IScanFilter
from ThermoFisher.CommonCore.RawFileReader import RawFileReaderAdapter


#Delta in mass for ion masses to search for (in ppm)
delta_mass_ppm = 10.0
#Detection intensity threshold
detection_threshold = 0      
#Delta retention time allowed in minutes
delta_rt = 10.0
#CSV file path containing the masses and retention times to search for
csv_file_mass = "C:\\Users\\Marc\\Documents\\Test\\list_mass_test.csv"
#Folder containing the raw files
working_dir="C:\\Users\\Marc\\Documents\\Test\\"

# Make sure .NET SDK >= 8.x is installed
def check_dotnet_version(min_major: int = 8) -> bool:
    try:
        # Run command: dotnet --list-sdks
        output_dotnet = subprocess.run(
            ["dotnet", "--list-sdks"],
            capture_output=True,
            text=True,
            check=True
        )

        # Typical lines look like: "8.0.100 [/usr/share/dotnet/sdk]"
        for line in output_dotnet.stdout.splitlines():
            m = re.match(r"\s*(\d+)\.(\d+)\.(\d+)", line)
            if not m:
                continue

            major = int(m.group(1))
            if major >= min_major:
                return True

        return False

    except FileNotFoundError:
        print("An exception was raised : 'dotnet' command not found (is .NET installed and in PATH?)")
        return False
    except subprocess.CalledProcessError as e:
        print(f"An exception was raised : dotnet returned non-zero exit code: {e}")
        return False
    except Exception as e:
        print(f"An exception was raised : {e}")
        return False

             
def check_pythonnet_version(min_major=2, min_minor=5):
    try:
        v = pkg_version("pythonnet")  # ex: "3.0.3"
        major, minor = (int(x) for x in v.split(".")[:2])
        return (major > min_major) or (major == min_major and minor >= min_minor)
    except Exception as e:
        print(f"An exception was raised : {e}")
        return False

def load_file_MS(rawFile):

    scanFilter = IScanFilter(rawFile.GetFilterForScanNumber(1)) # First scan is MS
    scanEvent = IScanEventBase(rawFile.GetScanEventForScanNumber(1)) #scanEvent for first MS

    #Enumerate scans that have the first scan's filter
    IEnumerator = rawFile.GetFilteredScanEnumerator(scanFilter)

    scanNumbers = []
    ret_times = []

    #Count number of scans matching the filter
    nscans=0
    for i in IEnumerator:
        scanEvent_i = IScanEventBase(rawFile.GetScanEventForScanNumber(i))
        if scanEvent_i.SourceFragmentation==scanEvent.SourceFragmentation:
            scanNumbers.append(i)
            ret_times.append(rawFile.RetentionTimeFromScanNumber(i))
            nscans+=1

    intensities = []
    masses = []
    scans = []
    charges = []
    for scanNumber in scanNumbers: # for each scan
        try:
            #Read the scan
            scan = Scan.FromFile(rawFile, scanNumber)
            scans.append(scanNumber)
            if scan.HasCentroidStream: # if centroided data is there
                size = len(scan.CentroidScan.Masses)
                #Extract intensities, masses and charges
                i = []
                m = []
                c = []
                for j in range(size):
                    i.append(scan.CentroidScan.Intensities[j])
                    m.append(scan.CentroidScan.Masses[j])
                    c.append(scan.CentroidScan.Charges[j])
                #Insert intensities, masses and charges to final output data
                intensities.append(i)
                masses.append(m)
                charges.append(c)
            else: # No centroided data
                intensities.append([])
                masses.append([])
                charges.append([])
        except Exception as ex:
            print('Error reading spectrum {} - {}'.format(scanNumber, str(ex)))

    return ret_times, scans, masses, intensities, charges

global center
global intensity

def peak_fitting(params, x, y):
    #h, w, z = params[0], params[1], params[2]
    w = params[0]
    z=center
    h=intensity
    residual = y-(np.array([h]*len(x))*np.exp(-pow(np.array(x)-np.array([z]*len(x)),2)/(2*(pow(w,2) if w!=0 else 0.0001))))
    return residual

def scan_files(working_dir):
    #List existing files
    files = os.listdir(working_dir)
    #Execution time data for each file
    total_file_execution_time=np.array([],dtype=np.float64)
    
    #Check if the directory path working_dir is correct
    if(not ( os.path.exists(working_dir) and os.path.isdir(working_dir) ) ):
        print("The directory path " + working_dir + " does not exist")
        exit(1)
        
    #Check if the file csv_file_mass exists
    if(not ( os.path.exists(os.path.join(csv_file_mass)) and os.path.isfile(os.path.join(csv_file_mass)) ) ):
        print("The file path " + os.path.join(csv_file_mass) + " does not exist")
        exit(1)
        
    #List of the files which contain the .raw extension
    raw_files_extensions = [file for file in files if file.endswith('.raw')]
    #Check if there is at least one raw file in the directory
    if(len(raw_files_extensions)<1):
        print("There no .raw file in " + os.path.join(working_dir))
        exit(1)
    
    file_name_list = []
    ions_masses_list = []
    charges_list = []
    first_iso_charges_list = []
    scans_list = []
    ret_times_list = []
    intensity_mass_list = []
    first_iso_intensity_mass_list = []
    first_over_mono_iso_ratio_list = []
    basepeak_list = []
    
    for i in range(len(files)): # For each file in working_dir
        if(files[i].endswith('.raw')): # If it is a raw file
            #Display file name
            print("file name=", files[i])
            begin_preprocessing=time.time()
            # Check to see if the specified RAW file exists
            if not os.path.isfile(os.path.join(working_dir, files[i])):
                print('The file doesn\'t exist in the specified location - {}'.format(files[i]))
                sys.exit()

            # Create the IRawDataPlus object for accessing the RAW file
            rawFile = RawFileReaderAdapter.FileFactory(os.path.join(working_dir, files[i]))

            # Check if the RAW file is accessible from RawFileReader library
            if not rawFile.IsOpen or rawFile.IsError:
                print('Unable to access the RAW file using the RawFileReader class!')
                sys.exit()

            # Check for any errors in the RAW file
            if rawFile.IsError:
                print('Error opening ({}) - {}'.format(rawFile.FileError, files[i]))
                sys.exit()

            # Check if the RAW file is being acquired
            if rawFile.InAcquisition:
                print('RAW file still being acquired - {}'.format(files[i]))
                sys.exit()

            #Select mass instrument
            rawFile.SelectInstrument(Device.MS, 1)

            # Get the first and last scan from the RAW file
            startScan = rawFile.RunHeaderEx.FirstSpectrum
            endScan = rawFile.RunHeaderEx.LastSpectrum

            # Define the settings for getting the Base Peak chromatogram
            settings = ChromatogramTraceSettings(TraceType.BasePeak)

            # Get the chromatogram from the RAW file.
            data = rawFile.GetChromatogramData([settings], startScan, endScan)

            # Split the data into the chromatograms
            trace = ChromatogramSignal.FromChromatogramData(data)
            end_preprocessing=time.time()
            #Extract MS data from the raw file
            ret_times_ftms, scans_ftms, masses_ftms, intensities_ftms, charges_ftms = load_file_MS(rawFile)
            end_load_file_MS = time.time()
            #Flag which indicates if ret_times are in the csv_file_mass file
            ret_times_available = True
            
            #Read CSV file path
            df_file_mass = pd.read_csv(csv_file_mass, sep=";")
            if("Mass" in df_file_mass.columns): #Check if there is a column named "Mass"
                #Extract masses to search for in the scans
                ions_masses = np.array(df_file_mass['Mass'].tolist())
                                
                if(not all(isinstance(x, float) for x in ions_masses)):
                    print("Error : some ion masses in " + csv_file_mass + " are not real numbers")
                    exit(1)
                
                delta_C12_C13 = 1.003355
                
                ions_masses_first_isotope = ions_masses + delta_C12_C13
            else:
                print("Error : There is no column 'Mass' in " + csv_file_mass)
                exit(1)
            
            if("RT" in df_file_mass.columns):
                #Extract corresponding RT to search for in the scans
                expected_ret_times = np.array(df_file_mass['RT'].tolist())
                
                if(len(ions_masses)!=len(expected_ret_times)):
                    print("Error : the number of ion masses in "+ csv_file_mass +" does not correspond to the number of RT")
                    exit(1)
                
                if(not all(isinstance(x, float) for x in expected_ret_times)):
                    print("Error : some ion RT in " + csv_file_mass + " are not real numbers")
                    exit(1)

                if(len(ions_masses)!=len(expected_ret_times)): #Check that the arrays have the same size
                    #Display error message if sizes are different
                    print("Error : The number of ion masses to search for is not the same as the number of expected retention times")
                    #Exit to avoid doing anything
                    exit(1)
            else:
                #There is no approximation of the retention time in csv_file_mass file
                ret_times_available = False
                
            #Final results array initialization
            ret_times=np.zeros(len(ions_masses), dtype=np.float32)
            scans=np.zeros(len(ions_masses), dtype=np.int32)
            intensity_mass = np.zeros(len(ions_masses), dtype=np.float32)
            charges=np.zeros(len(ions_masses), dtype=np.int32)
            first_iso_intensity_mass = np.zeros(len(ions_masses), dtype=np.float32)
            first_iso_charges=np.zeros(len(ions_masses), dtype=np.int32)
            first_over_mono_iso_ratio=np.zeros(len(ions_masses), dtype=object)

            #Compute the indexes corresponding to the rt boundaries for all ion masses to search for
            idx_boundaries_expected_rt = []
            
            if(ret_times_available):
                for j in range(len(ions_masses)):
                    bound_inf = np.searchsorted(ret_times_ftms, expected_ret_times[j]-delta_rt) #binary search for sorted arrays
                    bound_sup = np.searchsorted(ret_times_ftms, expected_ret_times[j]+delta_rt)-1 #binary search for sorted arrays
                    idx_boundaries_expected_rt.append([bound_inf, bound_sup])
            else: # No approximative retention times were provided
                for j in range(len(ions_masses)):
                    bound_inf = ret_times_ftms[0] #First available RT in MS
                    bound_sup = ret_times_ftms[-1] #Last available RT in MS
                    idx_boundaries_expected_rt.append([bound_inf, bound_sup])

            for j in range(len(scans_ftms)): # For each MS
                #Summed intensities of the ions for current scan +/- delta
                intensity_mass_temp = np.zeros(len(ions_masses), dtype=np.float32)
                #Intensity of the corresponding first isotopes
                intensity_mass_first_iso_temp = np.zeros(len(ions_masses), dtype=np.float32)
                #Max intensity peak of the ions for current scan +/- delta
                max_intensity_peak = np.zeros(len(ions_masses), dtype=np.float32)
                #Max intensity peak of the first isotope for current scan +/- delta
                max_intensity_first_isotopic_peak = np.zeros(len(ions_masses), dtype=np.float32)
                #Max intensity peak charge for every ion mass to search for
                charge_max_intensity_peak = np.zeros(len(ions_masses), dtype=np.int32)
                #Max intensity peak charge for first isotopic peaks of every ion mass to search for
                charge_max_intensity_first_isotopic_peak = np.zeros(len(ions_masses), dtype=np.int32)
                #First isotope mass intensity backup
                first_iso_intensity_mass_backup = np.zeros(len(ions_masses), dtype=np.float32)
                #First isotope charges backup
                first_iso_charges_backup = np.zeros(len(ions_masses), dtype=np.int32)
                
                #Find the ion masses that are concerned by this MS
                ions_masses_idx_rt = []
                for k in range(len(ions_masses)):
                    if(j>=idx_boundaries_expected_rt[k][0] and j<=idx_boundaries_expected_rt[k][1]):
                        ions_masses_idx_rt.append(k)

                for idx_mass in range(len(masses_ftms[j])): # For each mass in the MS
                    for idx_ion in ions_masses_idx_rt: # For each ion_mass to search for in MS
                        if(abs(ions_masses[idx_ion]-masses_ftms[j][idx_mass])/ions_masses[idx_ion]<(delta_mass_ppm/1e6) and intensities_ftms[j][idx_mass]>detection_threshold): # if it matches
                            intensity_mass_temp[idx_ion]+=intensities_ftms[j][idx_mass] # Update summed intensity_mass
                            if(intensities_ftms[j][idx_mass]>max_intensity_peak[idx_ion]): #If the intensity is higher than maximal known intensity 
                                max_intensity_peak[idx_ion]=intensities_ftms[j][idx_mass] #Update maximal known intensity
                                charge_max_intensity_peak[idx_ion]=charges_ftms[j][idx_mass] #Record maximal intensity known peak charge

                for idx_ion in ions_masses_idx_rt: # For each ion mass to search for in the MS
                    if(max_intensity_peak[idx_ion]>intensity_mass[idx_ion]): # If the summed intensities +/- delta_mass is the highest known
                        intensity_mass[idx_ion]=max_intensity_peak[idx_ion] # Update highest known summed intensities +/- delta_mass
                        ret_times[idx_ion]=rawFile.RetentionTimeFromScanNumber(scans_ftms[j]) # Set the new retention time
                        scans[idx_ion]=scans_ftms[j] # Set the new scan number
                        charges[idx_ion]=charge_max_intensity_peak[idx_ion] # Set the new corresponding charge


            for idx_ion in range(len(ions_masses)): #For each ion mass
                if(scans[idx_ion]!=0): # If the monoisotopic peak was found
                    #Check the presence and intensity of the first isotope
                    ret_times_iso = ret_times[idx_ion]
                    scan_iso = scans[idx_ion]
                    mass_iso = ions_masses_first_isotope[idx_ion] #Start at the same mass
                    idx_scan = scans_ftms.index(scan_iso)
                    idx_mass_iso = np.searchsorted(masses_ftms[idx_scan], ions_masses[idx_ion], side='left')

                    while(idx_mass_iso < len(masses_ftms[scan_iso])-1): #While not at the end of the masses list of the scan
                        if((masses_ftms[scan_iso][idx_mass_iso]-ions_masses_first_isotope[idx_ion])/ions_masses_first_isotope[idx_ion]<(delta_mass_ppm/1e6)): #It is still possible to detect the first isotope
                            if(intensities_ftms[scan_iso][idx_mass_iso]>detection_threshold): # First isotope mass matched
                                intensity_mass_first_iso_temp[idx_ion]+=intensities_ftms[scan_iso][idx_mass_iso] # Update summed intensity for the first isotope
                                if(intensities_ftms[scan_iso][idx_mass_iso]>max_intensity_first_isotopic_peak[idx_ion]): #If the intensity is higher than maximal known first iso intensity
                                    max_intensity_first_isotopic_peak[idx_ion]=intensities_ftms[scan_iso][idx_mass_iso] #Update maximal known intensity
                                    charge_max_intensity_first_isotopic_peak[idx_ion]=charges_ftms[scan_iso][idx_mass_iso] #Update maximal known intensity
                        else:
                            break #Exit loop when its not possible to detect the first isotope anymore
                        idx_mass_iso+=1
                        
                    first_iso_intensity_mass[idx_ion] = max_intensity_first_isotopic_peak[idx_ion]
                    first_iso_charges[idx_ion] = charge_max_intensity_first_isotopic_peak[idx_ion]

            global center
            global intensity
            
            #Compute first isotope / monoisotopic peak ratios with an advanced interpolation / integration technique
            for idx_ion in range(len(ions_masses)):
                integral_mono = 0
                integral_first_iso = 0
                scan_iso = scans[idx_ion]

                if(scans[idx_ion]!=0): # If the monoisotopic peak was found
                    idx_scan = scans_ftms.index(scan_iso)
                    params=[1]
                    #Absolute value of the difference between searched ion mass and existing masses in the spectrum
                    err_ppm = abs(np.array(masses_ftms[idx_scan])-ions_masses[idx_ion])/ions_masses[idx_ion]
                    #Index of the closest mass
                    idx_err_ppm = np.argsort(err_ppm)[0]
                    
                    #Get all masses points which are closer than delta_mass_ppm ppm
                    bound_mass_inf = idx_err_ppm
                    
                    while(err_ppm[bound_mass_inf] > (delta_mass_ppm/1e6)):
                        bound_mass_inf-=1
                        
                    bound_mass_inf-=1

                    bound_mass_sup = idx_err_ppm #np.max(idx_err_ppm)
                        
                    while(err_ppm[bound_mass_sup] > (delta_mass_ppm/1e6)):
                        bound_mass_sup+=1
                    
                    bound_mass_sup+=1
                    
                    #Backup for the masses/intensities which correspond to the searched ion mass +/- delta_mass_ppm
                    mass_back = masses_ftms[idx_scan][bound_mass_inf:bound_mass_sup]                  
                    int_back = intensities_ftms[idx_scan][bound_mass_inf:bound_mass_sup]
                    
                    #Monoisotopic centroid mass
                    center_mono=mass_back[int_back.index(np.max(int_back))]
                    #Set monoisotopic peak mass center for the gaussian fitting
                    center = center_mono
                    #Monoisotopic centroid intensity 
                    intensity_mono=int_back[mass_back.index(center_mono)]
                    #Set monoisotopic peak intensity for the gaussian fitting
                    intensity = intensity_mono
                    #Set the half maximal intensity of the monoisotopic peak
                    half_max_intensity_mono = intensity_mono/2.0
                    
                    #Linear interpolation which doubles the number of points
                    mass_fitting_augmented = np.zeros(2*len(mass_back))
                    int_fitting_augmented = np.zeros(2*len(int_back))
                    
                    for j in range(len(mass_fitting_augmented)-1):
                        if(j%2==0):
                            mass_fitting_augmented[j]=mass_back[int(j/2)]
                        else:
                            mass_fitting_augmented[j]=(mass_back[int((j-1)/2)]+mass_back[int((j+1)/2)])/2.0

                    for j in range(len(int_fitting_augmented)-1):
                        if(j%2==0):
                            int_fitting_augmented[j]=int_back[int(j/2)]
                        else:
                            int_fitting_augmented[j]=(int_back[int((j-1)/2)]+int_back[int((j+1)/2)])/2.0
                            
                    #Fit the corresponding data with a gaussian                    
                    result_mono = leastsq(peak_fitting, params, (mass_back, int_back))
                    yfit_mono = np.array([intensity_mono]*len(mass_fitting_augmented))*np.exp(-pow(np.array(mass_fitting_augmented)-np.array([center_mono]*len(mass_fitting_augmented)),2)/(2*(pow(result_mono[0][0],2) if pow(result_mono[0][0],2) > 0 else 2*sys.float_info.min)))
                    #Select the indexes which correspond to intensities greater than half monoisotopic peak maximum
                    indexes_mono_augmented = np.where(np.array(int_fitting_augmented) > half_max_intensity_mono)[0]
                    
                    mass_fitting_augmented = np.take(mass_fitting_augmented, indexes_mono_augmented)
                    int_fitting_augmented = np.take(int_fitting_augmented, indexes_mono_augmented)
                    
                    #Linear interpolation which doubles the number of points again
                    mass_fitting_augmented_integral = np.zeros(2*len(mass_fitting_augmented))

                    for j in range(len(mass_fitting_augmented_integral)-1):
                        if(j%2==0):
                            mass_fitting_augmented_integral[j]=mass_fitting_augmented[int(j/2)]
                        else:
                            mass_fitting_augmented_integral[j]=(mass_fitting_augmented[int((j-1)/2)]+mass_fitting_augmented[int((j+1)/2)])/2.0
                    
                    #Gaussian fitting of the augmented interpolation to increase accuracy
                    yfit_mono_integral = np.array([intensity_mono]*len(mass_fitting_augmented_integral))*np.exp(-pow(np.array(mass_fitting_augmented_integral)-np.array([center_mono]*len(mass_fitting_augmented_integral)),2)/(2*(pow(result_mono[0][0],2) if pow(result_mono[0][0],2) > 0 else 2*sys.float_info.min)))
                    
                    #Integrate mono over the interval
                    integral_mono = np.trapz(yfit_mono_integral, mass_fitting_augmented_integral)
                    
                    #First isotopic peak absolute values errors along the scan spectrum                    
                    err_ppm_first_iso = abs(np.array(masses_ftms[idx_scan])-ions_masses_first_isotope[idx_ion])/ions_masses_first_isotope[idx_ion]

                    #Closest ion mass peak index
                    idx_err_ppm_first_iso = np.argsort(err_ppm_first_iso)[0]
                    
                    #Select masses which correspond to the masses closer than delta_mass_ppm
                    bound_mass_inf_first_iso = idx_err_ppm_first_iso
                    
                    while(err_ppm_first_iso[bound_mass_inf_first_iso] > (delta_mass_ppm/1e6)):
                        bound_mass_inf_first_iso-=1
                        
                    bound_mass_inf_first_iso-=1
                    
                    bound_mass_sup_first_iso = idx_err_ppm_first_iso
                    
                    while(err_ppm_first_iso[bound_mass_sup_first_iso] > (delta_mass_ppm/1e6)):
                        bound_mass_sup_first_iso+=1
                        
                    bound_mass_sup_first_iso+=1
                    
                    #Masses backup for the first isotopic peaks values closer than delta_mass_ppm
                    mass_back_first_iso = masses_ftms[idx_scan][bound_mass_inf_first_iso:bound_mass_sup_first_iso]
                    #Intensities backup for the first isotopic peaks values closer than delta_mass_ppm
                    int_back_first_iso = intensities_ftms[idx_scan][bound_mass_inf_first_iso:bound_mass_sup_first_iso]
                    #Maximal intensity mass among those closer than delta_mass_ppm
                    center_first_iso=np.max(mass_back_first_iso)
                    #Set first isotopic peak mass center for the gaussian fitting
                    center = center_first_iso
                    #Maximal intensity first isotopic peak value
                    intensity_first_iso=int_back_first_iso[mass_back_first_iso.index(center_first_iso)]
                    #Set maximal intensity first isotopic peak value for the gaussian fitting
                    intensity = intensity_first_iso
                    
                    #Linear interpolation which doubles the number of points
                    mass_back_first_iso_augmented = np.zeros(2*len(mass_back_first_iso))
                    int_back_first_iso_augmented = np.zeros(2*len(int_back_first_iso))
                    
                    for j in range(len(mass_back_first_iso_augmented)-1):
                        if(j%2==0):
                            mass_back_first_iso_augmented[j]=mass_back_first_iso[int(j/2)]
                        else:
                            mass_back_first_iso_augmented[j]=(mass_back_first_iso[int((j-1)/2)]+mass_back_first_iso[int((j+1)/2)])/2.0

                    for j in range(len(int_back_first_iso_augmented)-1):
                        if(j%2==0):
                            int_back_first_iso_augmented[j]=int_back_first_iso[int(j/2)]
                        else:
                            int_back_first_iso_augmented[j]=(int_back_first_iso[int((j-1)/2)]+int_back_first_iso[int((j+1)/2)])/2.0
                    
                    #Gaussian fitting of the first isotopic peak
                    result_first_iso = leastsq(peak_fitting, params, (mass_back_first_iso, int_back_first_iso))
                    #Gaussian approximation within the augmented masses window
                    yfit_first_iso = np.array([intensity_first_iso]*len(mass_back_first_iso_augmented))*np.exp(-pow(np.array(mass_back_first_iso_augmented)-np.array([center_first_iso]*len(mass_back_first_iso_augmented)),2)/(2*(pow(result_first_iso[0][0],2) if pow(result_first_iso[0][0],2) > 0 else 2*sys.float_info.min)))
                    #Compute first isotopic peak half maximal intensity
                    half_max_intensity_first_iso = intensity_first_iso/2.0
                    #Select the indexes of the first isotopic peak intensities greater than the half maximum
                    indexes_first_iso_augmented = np.where(np.array(int_back_first_iso_augmented) > half_max_intensity_first_iso)[0]
                                        
                    mass_back_first_iso_augmented = np.take(mass_back_first_iso_augmented, indexes_first_iso_augmented)
                    int_back_first_iso_augmented = np.take(int_back_first_iso_augmented, indexes_first_iso_augmented)

                    #Linear interpolation again to improve accuracy
                    mass_fitting_first_iso_augmented_integral = np.zeros(2*len(mass_back_first_iso_augmented))
                
                    for j in range(len(mass_fitting_first_iso_augmented_integral)-1):
                        if(j%2==0):
                            mass_fitting_first_iso_augmented_integral[j]=mass_back_first_iso_augmented[int(j/2)]
                        else:
                            mass_fitting_first_iso_augmented_integral[j]=(mass_back_first_iso_augmented[int((j-1)/2)]+mass_back_first_iso_augmented[int((j+1)/2)])/2.0
                            
                    #Compute first isotopic peak integral
                    yfit_first_iso_integral = np.array([intensity_first_iso]*len(mass_fitting_first_iso_augmented_integral))*np.exp(-pow(np.array(mass_fitting_first_iso_augmented_integral)-np.array([center_first_iso]*len(mass_fitting_first_iso_augmented_integral)),2)/(2*(pow(result_first_iso[0][0],2) if pow(result_first_iso[0][0],2) > 0 else 2*sys.float_info.min)))
                    
                    #Integrate first isotopic peak
                    integral_first_iso = np.trapz(yfit_first_iso_integral, mass_fitting_first_iso_augmented_integral)
                
                #Compute first isotopic peak integration over monoisotopic peak integration
                first_over_mono_iso_ratio[idx_ion] = (integral_first_iso / integral_mono if (integral_mono != 0 and integral_first_iso != 0) else "N/A")
                
            print("first_over_mono_iso_ratio=", first_over_mono_iso_ratio)
                
            #Build the MS basepeak array
            basepeak = []
            for j in range(len(ions_masses)): #For each ion
                idx=0
                if(scans[j]>0): #If a scan was really found
                    #Find the basepeak index for this scan
                    while(trace[0].Scans[idx]!=scans[j]):
                        idx+=1
                    #Insert basepeak intensity based on the previsouly computed index
                    basepeak.append(trace[0].Intensities[idx])
                else: # If nothing was found for this ion
                    basepeak.append(0.0) #Put basepeak at zero
            #Create file name list
            file_name = [files[i]]*len(ions_masses)
            
            file_name_list.append(file_name)
            ions_masses_list.append(ions_masses)
            charges_list.append(charges)
            first_iso_charges_list.append(first_iso_charges)
            scans_list.append(scans)
            ret_times_list.append(ret_times)
            intensity_mass_list.append(intensity_mass)
            first_iso_intensity_mass_list.append(first_iso_intensity_mass)
            first_over_mono_iso_ratio_list.append(first_over_mono_iso_ratio)
            basepeak_list.append(basepeak)
            end_data_extraction=time.time()
            #Insert total processing time in the total_file_execution_time array
            total_file_execution_time=np.append(total_file_execution_time, end_data_extraction-begin_preprocessing)
            #Display execution times
            print("Preprocessing=", end_preprocessing-begin_preprocessing)
            print("Load file FTMS=", end_load_file_MS-end_preprocessing)
            print("Data extraction=", end_data_extraction-end_load_file_MS)
            print("Total file processing execution time=", end_data_extraction-begin_preprocessing)
    
    file_name_list = np.array(list(itertools.chain.from_iterable(file_name_list)))
    ions_masses_list = np.array(list(itertools.chain.from_iterable(ions_masses_list)))
    charges_list = np.array(list(itertools.chain.from_iterable(charges_list)))
    first_iso_charges_list = np.array(list(itertools.chain.from_iterable(first_iso_charges_list)))
    scans_list = np.array(list(itertools.chain.from_iterable(scans_list)))
    ret_times_list = np.array(list(itertools.chain.from_iterable(ret_times_list)))
    intensity_mass_list = np.array(list(itertools.chain.from_iterable(intensity_mass_list)))
    first_iso_intensity_mass_list = np.array(list(itertools.chain.from_iterable(first_iso_intensity_mass_list)))
    basepeak_list = np.array(list(itertools.chain.from_iterable(basepeak_list)))
    first_over_mono_iso_ratio_list = np.array(list(itertools.chain.from_iterable(first_over_mono_iso_ratio_list)))
    
    raw_files_list_no_duplicate = np.array(list(dict.fromkeys(file_name_list)))
    ions_masses_list_no_duplicate = np.array(list(dict.fromkeys(ions_masses_list)))
    nb_ions_masses = len(ions_masses_list_no_duplicate)
    nb_raw_files = len(raw_files_list_no_duplicate)
    
    intensities_cols = []

    for i in range(nb_raw_files):
        intensities_cols.append(intensity_mass_list[i*nb_ions_masses:(i+1)*nb_ions_masses])
    
    data_int = {'Mass':ions_masses_list_no_duplicate}
    
    for i in range(nb_raw_files):
        data_int.update({raw_files_list_no_duplicate[i]:intensities_cols[i]})

    cols=list(data_int.keys())
    
    df_intensities_arr = pd.DataFrame(data_int, columns=cols)
    
    df_intensities_arr.to_csv(os.path.join(working_dir, "FTMS_int_search_output_trapz.csv"), index=False, sep=';')
    
    #Create pandas dataframe with the extracted data
    data_FTMS = {'File name':file_name_list, 'Mass':ions_masses_list, 'Charges':charges_list, 'Scans':scans_list, 'Ret.time':ret_times_list, 'Intensity':intensity_mass_list, 'First iso intensity':first_iso_intensity_mass_list, 'First iso charges':first_iso_charges_list, 'First over Mono iso ratio': first_over_mono_iso_ratio_list, 'Basepeak':basepeak_list}
    cols=['File name', 'Mass', 'Charges', 'Scans', 'Ret.time', 'Intensity', 'First iso intensity', 'First iso charges', 'First over Mono iso ratio', 'Basepeak']
    df_final_result = pd.DataFrame(data_FTMS, columns=cols)
    #Save the data
    df_final_result.to_csv(os.path.join(working_dir, "FTMS_search_output_trapz.csv"), index=False, sep=';')

    #Sum all individual files execution time
    all_files_execution_time = total_file_execution_time.sum()
    #Display the total execution time
    print("All files execution time=", all_files_execution_time)

#Verify if .NET SDK v8.x is installed
if(not check_dotnet_version()):
    print(".NET has to be >=8.x")
    print("Exiting")
    exit(1)

#Verify if Pythonnet 2.5.x is installed
if(not check_pythonnet_version()):
    print("Pythonnet has to be >= 2.5.x")
    print("Exiting")
    exit(1)

#Call the main function
scan_files(working_dir)
