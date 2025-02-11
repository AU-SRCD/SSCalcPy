# -*- coding: utf-8 -*-
"""
Created on Wed Sep 20 12:19:05 2023

@author: au513
"""

import numpy as np
import sys
#Check that the file SelconsFunction.py is in the path given below
sys.path.append("../Selcon3Py")
from SelconsFunction import SelconsPy #SelconsFunction.py contains the SelconPy function

#Remember to give the correct protein file name and path here:
#Path or directory
#MyProteinDir = "../Data/timW/"
#MyProteinDir = "../TestProteins/"
MyProteinDir = "IR_SpectraFromA/"
#MyProteinDir = "tmpMayBeDeleted/"
MyProteinDir = "CD-IR_SpectraFromA/"
MyProteinDir = "SpectraFromA/"

CD_or_IR = 'CD-IR2cm'
scaleIR = 15 #15.0
#CD_or_IR = 'IR'
#CD_or_IR = 'CD'

AU_or_BBK = ''

ShouldIPlot = 1
ShouldIAlsoPlotSelcon2 = 0

if (CD_or_IR == 'IR'):
    #Note on reference set: Reference publication
    #Marco Pinto Corujo, Adewale Olamoyesan, Anastasiia Tukova, Dale Ang,Erik Goormaghtigh, Jason Peterson, Victor Sharov, Nikola Chmel and Alison Rodger. Front. Chem. 9:784625 (2022) doi: 10.3389/fchem.2021.784625
    #A_IRfile = '../IR_ref/A50-AR-EG-IR_Sept2023.txt'
    A_IRfile = '../IRRefData/A50-AR-EG-IR_Sept2023.txt'
    F_IRfile = '../IRRefData/F50-AR-EG-IR_Sept2023.txt'
    #IR notes: 
    #Ref data are from 1600 to 1800 cm-1 in increasing order
    #Protein IR spectrum should also start from 1600 to 1800. We will interpolate the spectrum to steps of 1 cm-1
    IR_LowKInRefSet = 1600; IR_HighKInRefSet = 1800; IR_RefSetSteps = 1
    IR_RefSetNumDataPoints = 1 + ((IR_HighKInRefSet-IR_LowKInRefSet)/IR_RefSetSteps); print('IR_RefSetNumDataPoints: ', IR_RefSetNumDataPoints)
    
    LimWL = 1800 #is this right
    RefDataLowWL = 1800 #hack to make it fit 240-174 nm range thinking in the code
    #Now laod the reference data
    F1 =  np.loadtxt(F_IRfile, dtype='f', delimiter='\t')  ; print('F matrix shape: ', F1.shape)
    A1 = np.loadtxt(A_IRfile, dtype='f', delimiter='\t') ; print('A matrix shape: ', A1.shape)
    #sys.exit()
    
#END: if (CD_or_IR == 'IR'):

if 'CD-IR' in CD_or_IR:
    A_CD_IRfile = '../'+'CD-IRRefDataAU' + '/' + 'AU-A28_CD-IR_Aug24.txt'
    if '2cm' in CD_or_IR:
        A_CD_IRfile = '../'+'CD-IRRefDataAU' + '/' + 'AU-A28_CD-IR2cm-1_Aug24.txt'
    F_CD_IRfile = '../'+'CD-IRRefDataAU' + '/' + 'AU-F28_CD-IR_Aug24.txt'
    LimWL = 175 #for CD data
    RefDataHighWL = 240
    RefDataLowWL = 175
    
    
    F1 = np.loadtxt(F_CD_IRfile, dtype='f', delimiter='\t')  ; print('F matrix shape: ', F1.shape)
    A2 = np.loadtxt(A_CD_IRfile, dtype='f', delimiter='\t') ; print('A matrix shape: ', A2.shape)
    
    #Scale IR part of A and thus extracted q
    IR_LowKInRefSet = 1600; IR_HighKInRefSet = 1800; IR_RefSetSteps = 1
    if '2cm' in CD_or_IR:
        IR_LowKInRefSet = 1600; IR_HighKInRefSet = 1720; IR_RefSetSteps = 2
    
    #The CD data are in pos 0 to 240-175 = 65 (66th element)
    CD_NumDataPointsInA2 = int(1+(RefDataHighWL-LimWL))
    A1_CD = A2[:CD_NumDataPointsInA2,:]
    A1_IR = A2[CD_NumDataPointsInA2:,:] * scaleIR
    A1 = np.concatenate([A1_CD,A1_IR]);
    
    if scaleIR == 0:
        #only use the CD part of the A
        A1 = A1_CD
    
    
    
    
#END: if 'CD-IR' in CD_or_IR:

#loop over all proteins in the refrence data set (A1 / F1)
for RefProteinNo in range(A1.shape[1]):
    #Make new A and F matrix where REfProteinNo is left out.
    #q=np.zeros(A1.shape[0])
    MyProteinFile = "ProteinNo_{:d}.txt".format(RefProteinNo)
    
    A=np.delete(A1,RefProteinNo,1) #1=second axis or in this case the colomn
    F=np.delete(F1,RefProteinNo,1)
    q=A1[:,RefProteinNo]

    #Run the selcon3 function
    #def is              SelconsPy(A1,F,q1, ShouldIPlot, ShouldIAlsoPlotSelcon2, LimWL, RefDataLowWL, MyProteinFile, MyProteinDir,AU_or_BBK,CD_or_IR):
    Selcon3ResultsList = SelconsPy(A ,F, q, ShouldIPlot, ShouldIAlsoPlotSelcon2, LimWL, RefDataLowWL, MyProteinFile, MyProteinDir,AU_or_BBK, CD_or_IR)
    #IR_SpectraFromA_Dir = "IR_SpectraFromA/"; DataSet = 'IR'
    IR_SpectraFromA_Dir = MyProteinDir; DataSet = 'IR'
    if 'CD-IR' in CD_or_IR:
        #IR_SpectraFromA_Dir = "CD-IR_SpectraFromA/"; DataSet = 'CD-IR'
        IR_SpectraFromA_Dir = MyProteinDir; DataSet = 'CD-IR'
        if '2cm' in CD_or_IR:
            DataSet = 'CD-IR2cm'

#    with open(IR_SpectraFromA_Dir + MyProteinFile[0:-4]+'_selcon3-'+ DataSet + '.txt', 'w') as f:
#        for item in q:
#            f.write("%s\n" % item)
    if (RefProteinNo==0):
        with open(IR_SpectraFromA_Dir + 'All_Results'+'_selcon3-'+ DataSet + '_out.txt', 'w') as f:
            item0=Selcon3ResultsList[-2]
            f.write("%s\n" % item0)
    #END: if (RefProteinNo==0):
    with open(IR_SpectraFromA_Dir + 'All_Results'+'_selcon3-'+  DataSet + '_out.txt', 'a') as f:
        item1=Selcon3ResultsList[-1]
        f.write("%s\n" % item1)
   
#END: for RefProteinNo in range(A1.shape[1]):