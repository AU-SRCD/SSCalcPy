# -*- coding: utf-8 -*-
"""
Created on Fri Sep 22 12:30:48 2023

@author: au513
"""

#Selcon3LoadAndRun.py

import numpy as np
#import sys
from SelconsFunction import SelconsPy #SelconsFunction.py contains the SelconPy function


def LoadAndRun(delimIs, CD_or_IR, AU_or_BBK, RefDataLowWL, RefDataHighWL, MyProteinDir, MyProteinFile, ShouldIPlot, ShouldIAlsoPlotSelcon2, scaleCD = 1):
    
    print("loading data and Reference set...")
    
    if (CD_or_IR == 'CD'):
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
        
    #*** SPECIAL code for IR spectra fitting ****
    if (CD_or_IR == 'IR'):
        #Note on reference set: Reference publication
        #Marco Pinto Corujo, Adewale Olamoyesan, Anastasiia Tukova, Dale Ang,Erik Goormaghtigh, Jason Peterson, Victor Sharov, Nikola Chmel and Alison Rodger. Front. Chem. 9:784625 (2022) doi: 10.3389/fchem.2021.784625
        A_IRfile = 'IRRefData' + '/' + 'A50-AR-EG-IR_Sept2023.txt' # 'IR_ref/A50-AR-EG-IR_Sept2023.txt'
        F_IRfile = 'IRRefData' + '/' + 'F50-AR-EG-IR_Sept2023.txt' # 'IR_ref/F50-AR-EG-IR_Sept2023.txt'
        #IR notes: 
        #Ref data are from 1600 to 1800 cm-1 in increasing order
        #Protein IR spectrum should also start from 1600 to 1800. We will interpolate the spectrum to steps of 1 cm-1
        IR_LowKInRefSet = 1600; IR_HighKInRefSet = 1800; IR_RefSetSteps = 1
        IR_RefSetNumDataPoints = 1 + ((IR_HighKInRefSet-IR_LowKInRefSet)/IR_RefSetSteps); print('IR_RefSetNumDataPoints: ', IR_RefSetNumDataPoints)
        
        LimWL = 1800 #is this right?
        RefDataLowWL = 1800 #hack to make it fit 240-174 nm range thinking in the code
        
        # skip = 0 #for now assume the IR data are starting at first line (line 0)
        # q1 = np.loadtxt(MyProteinDir + MyProteinFile, dtype='f', delimiter='\t', skiprows=skip)
        
        #Extract the k and the IR values from the spectrum
        k_values = q1[:, 0]; IR_values =q1[:, 1]
        #.... and interpolate
        k_InRefSet = np.linspace(IR_LowKInRefSet, IR_HighKInRefSet, int(IR_RefSetNumDataPoints))
        q_interp = np.interp(k_InRefSet, k_values, IR_values) #Note: k_InRefSet values out side range of k_values are given as IR_values[0] and IR_values[-1]
        q = q_interp
        q = q_interp / max(q_interp) #Should this be for all files???
    
        #Now laod the reference data
        F =  np.loadtxt(F_IRfile, dtype='f', delimiter='\t')  ; print('F matrix shape: ', F.shape)
        A = np.loadtxt(A_IRfile, dtype='f', delimiter='\t') ; print('A matrix shape: ', A.shape)
        #sys.exit()
        
    #END: if (CD_or_IR == 'IR'):
    #***END: SPECIAL code for IR spectra fitting ****
    
    #*** SPECIAL code for CD-IR spectra fitting ****
    if 'CD-IR' in CD_or_IR: #if (CD_or_IR == 'CD-IR'):
        #The IR data might need scaling. The IR scale is passed via scaleCD
        scaleIR = scaleCD
        
        A_CD_IRfile = 'CD-IRRefDataAU' + '/' + 'AU-A28_CD-IR_Aug24.txt'
        if '2cm' in CD_or_IR:
            A_CD_IRfile = 'CD-IRRefDataAU' + '/' + 'AU-A28_CD-IR2cm-1_Aug24.txt'
        F_CD_IRfile = 'CD-IRRefDataAU' + '/' + 'AU-F28_CD-IR_Aug24.txt'
        LimWL = 175 #for CD data
        #we have RefDataLowWL, RefDataHighWL set to 175, 240 in startSSCalc
        
       
        #R reduce A to only contain the CD data from 240 to actual LimWL
        #--- first load the data and find who many CD datapoints there are.
        
        #Load the CD and IR data 240, 239,...LimWLCD, 1600, 1601, ... 1800
        #... but the CD data could also start at a higher WL, followed by IR data that might also not be in only the range 1600-1800(1720)
        #Read the CD data and interpolate
        #...first find lowest WL datapoint
        minWLFloat = np.min(q1[:,0]) #even for mixed CD and IR data, this is the lowest CD data wl
        LimWL = max(int(minWLFloat), RefDataLowWL) #minWL
        print('LimWL is determined to be=',LimWL)
        CD_NumDataPoints = int(1+(RefDataHighWL-LimWL))
        wl_increasingCD = np.linspace(LimWL, RefDataHighWL, CD_NumDataPoints)
        #q1CD = q1[0:CD_NumDataPoints,:] #! Doesn't work if maxWL>RefDataHighWL. extract the CD data

        #Find the place in the q1[:,0] where there is a big jump in WL (where IR data starts)
        lastCDwl = q1[-1,0] #temp set to last element (although this is most likely an IR k value)
        for count, qvalue in enumerate(q1[1:,0]):
            #note: count is 1 less than the index for the qvalue, i.e. starts at zero. So q1[count,:]=q1[0,:] for the first cycle in the for loop
            if (abs(q1[count,0]-qvalue)>1):
                lastCDwl = q1[count,0]
                break
        CD_NumDataPointsInq1 = int(1+abs(q1[0,0]-lastCDwl))
        q1CD = q1[0:CD_NumDataPointsInq1,:] #extract the CD data
        print("CD part of data: First WL, Last WL ", q1CD[0,0],q1CD[-1,0])
        if (q1CD[-1,0] < q1CD[0,0]):
            q_interp = np.interp(wl_increasingCD, np.flip(q1CD[:, 0]), np.flip(q1CD[:, 1]))
        else:
            q_interp = np.interp(wl_increasingCD, q1CD[:, 0], q1CD[:, 1])
        #reverse q order to fit order of CD reference data
        qCD = np.flip(q_interp)
        #wlCD = np.arange(RefDataHighWL, LimWL-1, -1) #240:-1:LimWL; %Vector of wavelengths
        
        #Next extract the IR data
        pass
        IR_LowKInRefSet = 1600; IR_HighKInRefSet = 1800; IR_RefSetSteps = 1
        if '2cm' in CD_or_IR:
            IR_LowKInRefSet = 1600; IR_HighKInRefSet = 1720; IR_RefSetSteps = 2
        
        IR_RefSetNumDataPoints = 1 + ((IR_HighKInRefSet-IR_LowKInRefSet)/IR_RefSetSteps); print('IR_RefSetNumDataPoints: ', IR_RefSetNumDataPoints)
        #Extract the k and the IR values from the spectrum
        #k_values = q1[CD_NumDataPoints:, 0]; IR_values =q1[CD_NumDataPoints:, 1]
        k_values = q1[CD_NumDataPointsInq1:, 0]; IR_values =q1[CD_NumDataPointsInq1:, 1]
        print("IR part of data: First k, Last k ", k_values[0],k_values[-1])
        #.... and interpolate
        k_InRefSet = np.linspace(IR_LowKInRefSet, IR_HighKInRefSet, int(IR_RefSetNumDataPoints))
        q_interp = np.interp(k_InRefSet, k_values, IR_values) #Note: k_InRefSet values out side range of k_values are given as IR_values[0] and IR_values[-1]
        qIR = q_interp * (scaleIR / max(q_interp)) #note in comb CD-IR the IR might be scaled, passed via scaleCD
        #wlIR=np.arange(IR_LowKInRefSet, IR_HighKInRefSet+1, +1) #from 1600 to 1800
        
        #Put the IR and CD data together
        #wl= wlCD q = qCD wl= wlIR q = qIR
        #wl = np.concatenate([wlCD,wlIR])
        q = np.concatenate([qCD,qIR])
        
        #Now laod the reference data
        F1 =  np.loadtxt(F_CD_IRfile, dtype='f', delimiter='\t')  ; print('F1 matrix shape: ', F1.shape)
        A1 = np.loadtxt(A_CD_IRfile, dtype='f', delimiter='\t') ; print('A1 matrix shape: ', A1.shape)
        
        A_CD = A1[:CD_NumDataPoints,:]
        CD_NumDataPointsInA1 = int(1+(RefDataHighWL-RefDataLowWL))
        A_IR = A1[CD_NumDataPointsInA1:,:] * scaleIR #in combined CD-IR the IR might be scaled. Scale both IR parts of q and A
        
        A = np.concatenate([A_CD,A_IR]); print('A matrix shape: ', A.shape)
        F = F1; print('F matrix shape: ', F.shape)
        
        if (scaleIR == 0):
            #Special case where we shold ignore the IR part of the spectrum
            q = qCD
            A = A_CD; print('A matrix shape (scaleIR=0): ', A.shape)
        
        #Now we got q, A and F
        #How to deal with in LimWL and RefDataLowWL?
        #LimWL is the real low WL for the CD data and RefDataLowWL = 175 in the CD-IR case

        #SelconsPy(A,F,q, ShouldIPlot, ShouldIAlsoPlotSelcon2, LimWL, RefDataLowWL, MyProteinFile, MyProteinDir,AU_or_BBK, CD_or_IR)
        #return
    
    # END: if (CD_or_IR == 'CD-IR'):
    #***END: SPECIAL code for CD-IR spectra fitting ****
    
    if (CD_or_IR == 'CD'):
        #Find out which format the protein CD spectrum data has (single or double rows)
        try:
            #If the file loaded intp q1 is just a single colunm, then the following "if (q1.shape[1]>0):"
            # ...will fail and an exception is thrown. I that case just assume file is just a single column of delta Epsilon starting at 240 nm
            #Otherwise, assume the file is two columns with WL in first colunm and Delta Epsilon in second row
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
                #     print('240 nm data is in line', iRefDataHighWL, ' of the data lines in the file. First line is line = 0')
                #     qLength = q1.shape[0]-iRefDataHighWL
                #     q=np.zeros(qLength)
                #     q=q1[iRefDataHighWL:,1] #second row is De and first row is WL
                #     LimWL=max(int(q1[-1,0]),RefDataLowWL)
                #     print('LimWL is determined to be=',LimWL)
                # else:
                #     print('Failed to find the position of the 240 nm data point')
                #     q=q1
                
        except:
            #Assume that the file is just a single colunm of delta Epsilon starting at 240 nm
            q=q1
        #Finally scale the spectrum
        q = q * scaleCD
    #END: if (CD_or_IR == 'CD'):
    
    
    #Call the function found in SelconsFunction.py
    Selcon3ResultsList = SelconsPy(A,F,q, ShouldIPlot, ShouldIAlsoPlotSelcon2, LimWL, RefDataLowWL, MyProteinFile, MyProteinDir,AU_or_BBK, CD_or_IR)
    
    # Filenames are
    # Prot_SELCON3-SP_AU-PCDDB_xxx
    # Prot_SELCON3-SMP_AU-PCDDB_xxx
    # Prot_SELCON3-SP_BBK_xxx
    # Prot_SELCON3-SMP_BBK_xxx
    # Prot_SELCON3-IR_xxx
    # Prot_SELCON3-CD-IR_xxx
    # Prot_SELCON3-CD-IR2cm_xxx
    
    DataSet = '_' + AU_or_BBK
    if (CD_or_IR == 'IR'):
        DataSet = '-IR'
    if (CD_or_IR == 'CD-IR'):
        DataSet = '-CD-IR'
    if (CD_or_IR == 'CD-IR2cm'):
        DataSet = '-CD-IR2cm'
    
    if (AU_or_BBK == 'AU' and CD_or_IR == 'CD'):
        DataSet += '-PCDDB'
        
    SP_or_SMP = ''
    if (RefDataLowWL == 175):
        SP_or_SMP = '-SP'
    elif (RefDataLowWL == 180):
        SP_or_SMP = '-SMP'
    else:
        pass
    
    #Write the Selcon3 results to a file whihc is the ProtFilemane_out.txt
    #with open(MyProteinDir + MyProteinFile[0:-4]+'_selcon3-'+ str(RefDataLowWL) + '-'+ DataSet + '_out.txt', 'w') as f:
    #fOut = open(MyProteinDir + MyProteinFile[0:-4]+'_SELCON3-'+ str(RefDataLowWL) + '-' + DataSet + '_out.txt', 'w')
    fOut = open(MyProteinDir + MyProteinFile[0:-4]+'_SELCON3'+ SP_or_SMP + '' + DataSet + '_out.txt', 'w')
    for item in Selcon3ResultsList:
        fOut.write("%s\n" % item)
    if (scaleCD != 1.0):
        if 'CD-IR' in CD_or_IR:
            fOut.write("IR Spectrum was scaled %.4f\n" %scaleCD)
        else:
            fOut.write("CD Spectrum was scaled %.4f\n" %scaleCD)
    fOut.close()
    
#END: def LoadAndRun(CD_or_IR, AU_or_BBK, RefDataLowWL, LimWL):