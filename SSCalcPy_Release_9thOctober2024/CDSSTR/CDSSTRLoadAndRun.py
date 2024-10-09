# -*- coding: utf-8 -*-
"""
Created on Fri Oct 20 17:14:02 2023

@author: au513
"""
import numpy as np
import sys

from CDSSTRPy_Functions import CDSSTRPyFunc #CDSSTRFunctions.py contains the CDSSTR functions




def LoadAndRunCDSSTR(delimIs, CD_or_IR, AU_or_BBK, RefDataLowWL, RefDataHighWL, MyProteinDir, MyProteinFile, ShouldIPlot, ShouldIAlsoPlotSelcon2, scaleCD = 1):
    
    #Tmp fix
    Use_SS_Classification = "AU-F128_T"

    if (CD_or_IR == 'CD'):
        #RefDataDir = "../CDRefData/" #this is NOT the relative path to the GUI .py file
        RefDataDir = "CDRefDataAU/" #this is the relative path to the GUI .py file
        RefDataDirBBK = "CDRefDataBBK/"
        if (RefDataLowWL == 175 and AU_or_BBK == 'BBK'):
            LimWL = 175
            F = np.loadtxt(RefDataDirBBK + "F71.txt", dtype='f', delimiter='\t')  #SP175      #print(F.shape)
            A = np.loadtxt(RefDataDirBBK + "A71.txt", dtype='f', delimiter='\t')  #SP175      #print(A1.shape)
            
        if (RefDataLowWL == 180 and AU_or_BBK == 'BBK'):
            LimWL = 180
            F = np.loadtxt(RefDataDirBBK + "F128.txt", dtype='f', delimiter='\t') #SMP180
            A = np.loadtxt(RefDataDirBBK + "A128.txt", dtype='f', delimiter='\t') #SMP180
        
        if (RefDataLowWL == 175 and AU_or_BBK == 'AU'):
            LimWL = 175
            F1 =  np.loadtxt(RefDataDir + "AU-F128_T-Nov22.txt", dtype='f', delimiter='\t')  #SP175      #print(F.shape)
            A1 = np.loadtxt(RefDataDir + "AU-A128_PCDDB-Nov22.txt", dtype='f', delimiter='\t')  #SP175      #print(A1.shape)
            F = F1[:,:71] #All 66 WL (240-175), first 71 spectra (SP175)
            A = A1[:,:71] #print('F shape: ', F.shape, '. A shape: ', A.shape)
            
        if (RefDataLowWL == 180 and AU_or_BBK == 'AU'):
            LimWL = 180
            F1 =  np.loadtxt(RefDataDir + "AU-F128_T-Nov22.txt", dtype='f', delimiter='\t') #SMP180
            A1 = np.loadtxt(RefDataDir + "AU-A128_PCDDB-Nov22.txt", dtype='f', delimiter='\t') #SMP180
            F = F1[:61,:] #First 61 WL (240-180), all 128 spectra (SMP180)
            A = A1[:61,:] #print('F shape: ', F.shape, '. A shape: ', A.shape)

    NumRefCDSpectra = A.shape[1]
    ##IR Here
    
    #First load the spectrum
    #Find out of there are header lines in the protein data file
    file1 = open(MyProteinDir + MyProteinFile, 'r')
    Lines = file1.readlines()
    file1.close()
     
    skip = 0; readRows = 0
    SkipChars = [";", "#", "X"] #List of first characters in line to be skipped
    for line in Lines:
        if(line[0] in SkipChars or line[0].isnumeric() == False):
            if (readRows == 0):
                skip += 1
        else:
            readRows += 1
    print('Skipping first', skip, 'rows of protein data file', MyProteinFile)
    print('Number of rows to read: ', readRows) #Read max_rows rows of content after skiprows lines.
    
    q1 = np.loadtxt(MyProteinDir + MyProteinFile, dtype='f', delimiter=delimIs, skiprows=skip, max_rows=readRows)
    
    ##MORE CD
    if (CD_or_IR == 'CD'):
        #Find out which format the protein CD spectrum data has (single or double rows)
        try:
            #If the file loaded intp q1 is just a single row, then if (q1.shape[1]>0): fail and an exception is thrown
            #Otherwise assume the file is two rows with WL in first row and Delta Epsilon in second row
            if (q1.shape[1]>0):
                minWLFloat = np.min(q1[:,0])
                LimWL = max(int(minWLFloat), RefDataLowWL) #minWL
                print('LimWL is determined to be=',LimWL)
                wl_increasing = np.linspace(LimWL, RefDataHighWL, int(1+(RefDataHighWL-LimWL)))
                if (q1[-1,0] < q1[0,0]):
                    q_interp = np.interp(wl_increasing, np.flip(q1[:, 0]), np.flip(q1[:, 1]))
                else:
                    q_interp = np.interp(wl_increasing, q1[:, 0], q1[:, 1])
                #reverse q order to fit order of CD reference data
                q = np.flip(q_interp)
                
                
                # #find which data point is 240 nm/RefDataHighWL: iRefDataHighWL
                # iRefDataHighWL=-1
                # for i in range(q1.shape[0]):
                #     if(q1[i,0]==RefDataHighWL):
                #         iRefDataHighWL=i
                # if iRefDataHighWL>=0:
                #     qLength = q1.shape[0]-iRefDataHighWL
                #     q=np.zeros(qLength)
                #     q=q1[iRefDataHighWL:,1] #second row is De and first row is WL
                #     LimWL=max(int(q1[-1,0]),RefDataLowWL)
                #     print('LimWL is determined to be=',LimWL)
                # else:
                #     print('Failed to find the position of the 240 nm data point')
                #     q=q1
                
        except:
            #Assume that the file is just a single row of delta Epsilon starting at 240 nm
            #and limited to LimWL as set above
            tmpLimWL = LimWL
            LimWL = max(tmpLimWL, RefDataLowWL)
            q=q1
        #Finally scale the spectrum
        q = q * scaleCD
    #END: if (CD_or_IR == 'CD'):

    if (CD_or_IR == 'CD'):
        pass
    
    CDSSTRResultsList = CDSSTRPyFunc(A, F, q, LimWL, RefDataLowWL, NumRefCDSpectra, MyProteinFile, MyProteinDir, AU_or_BBK, CD_or_IR, Use_SS_Classification)
    
    #Write the results to file
    
    # Filenames are
    # Prot_CDSSTR-SP_AU-PCDDB_xxx
    # Prot_CDSSTR-SMP_AU-PCDDB_xxx
    # Prot_CDSSTR-SP_BBK_xxx
    # Prot_CDSSTR-SMP_BBK_xxx
    # (Prot_CDSSTR-IR_xxx)
    
    UseRefSet = '_' + AU_or_BBK  #UseRefSet = AU_or_BBK
    if (CD_or_IR == 'IR'):
        UseRefSet = '-IR'
    
    if (AU_or_BBK == 'AU'):
        UseRefSet += '-PCDDB'
    
    SP_or_SMP = ''
    if (RefDataLowWL == 175):
        SP_or_SMP = '-SP'
    elif (RefDataLowWL == 180):
        SP_or_SMP = '-SMP'
    else:
        pass
    
    if (len(CDSSTRResultsList) > 0):
        #with open(MyProteinDir + MyProteinFile[0:-4]+'_CDSSTR-' + str(RefDataLowWL) + '-' + UseRefSet +'_out.txt', 'w') as f:
        with open(MyProteinDir + MyProteinFile[0:-4]+'_CDSSTR' + SP_or_SMP + '' + UseRefSet +'_out.txt', 'w') as f:
            for item in CDSSTRResultsList:
                f.write("%s\n" % item)
            if (scaleCD != 1.0):
                f.write("Spectrum was scaled %.4f\n" %scaleCD)
    #END: if (len(CDSSTRResultsList) > 0):
     
#end Def ,...

# #Test function:
# CD_or_IR = 'CD'
# AU_or_BBK = 'AU'
# RefDataLowWL = 175
# RefDataHighWL = 240
# MyProteinDir = "../TestProteins/"
# MyProteinFile = "ConA.txt"
# ShouldIPlot = 1
# ShouldIAlsoPlotSelcon2 = 0

# LoadAndRunCDSSTR(CD_or_IR, AU_or_BBK, RefDataLowWL, RefDataHighWL, MyProteinDir, MyProteinFile, ShouldIPlot, ShouldIAlsoPlotSelcon2)

