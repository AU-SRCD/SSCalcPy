# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
scaleIR = 20 #1
skipLines = 1
MyProteinDir = "CD-IR_SpectraFromA/" #"IR_SpectraFromA/"
MyRefirDir = MyProteinDir + "LOOV_All-Results_refit_Scale" + str(scaleIR) + "_20240814/"
#ProteinNo_0_SELCON3-CD-IR2cm_refit.txt
FilePart1 = "ProteinNo_"
FilePart2 = "_SELCON3-CD-IR2cm_refit.txt"

CDstartRow = 1-1
CDEndRow = 66-1
IRStartRow = CDEndRow+1
IREndRow = IRStartRow + 61-1

rmsdAll=[]
rmsdCD=[]
rmsdIR=[]

numProteins = 28 #3


for RefProteinNo in range(0, numProteins, 1):
    File = MyRefirDir + FilePart1 + str(RefProteinNo) + FilePart2
    #print('RefProteinNo',RefProteinNo)
    #q1 = np.loadtxt(MyProteinDir + MyProteinFile, dtype='f', delimiter=delimIs, skiprows=skip, max_rows=readRows)
    q1 = np.loadtxt(File, dtype='f', delimiter='\t', skiprows=skipLines) #; print(q1.shape[0])
    #print("CD start/end WL " +str(q1[CDstartRow,0])+ "/"+ str(q1[CDEndRow,0])+ " IR start/end k "+ str(q1[IRStartRow,0])+"/"+ str(q1[IREndRow,0]))
    if RefProteinNo==0:
        if (q1.shape[0] <= IRStartRow):
            print('q1.shape[0] <= IRStartRow')
            #Fix: just calc CD RMSD for IR
            IRStartRow = CDstartRow
            IREndRow = CDEndRow
        print("File", File)
        print('A matrix shape: ', q1.shape)
        print("CD start/end WL",q1[CDstartRow,0],"/",q1[CDEndRow,0],"IR start/end k",q1[IRStartRow,0],"/",q1[IREndRow,0])
    rmsdAll.append(np.sqrt(((q1[:,1] - q1[:,2]) ** 2).mean()))
    rmsdCD.append(np.sqrt(((q1[CDstartRow:CDEndRow+1,1] - q1[CDstartRow:CDEndRow+1,2]) ** 2).mean()))
    rmsdIR.append(np.sqrt(((q1[IRStartRow:IREndRow+1,1] - q1[IRStartRow:IREndRow+1,2]) ** 2).mean()))
    print('rmsdAll/CD/IR: %.4f\t%.4f\t%.4f' %(rmsdAll[RefProteinNo],rmsdCD[RefProteinNo],rmsdIR[RefProteinNo]))
    #print("rmsdAll/CD/IR", RefProteinNo, ":", rmsdAll[RefProteinNo], rmsdCD[RefProteinNo], rmsdIR[RefProteinNo])

    if (RefProteinNo==0):
        with open(MyRefirDir + 'All_RMSD_Results' + FilePart2[:-4] + '_Scale' + str(scaleIR) + '_out.txt', 'w') as f:
            item0='RMSD\n' + 'All\tCD\tIR'
            f.write("%s\n" % item0)
    #END: if (RefProteinNo==0):
    with open(MyRefirDir + 'All_RMSD_Results' + FilePart2[:-4] + '_Scale' + str(scaleIR) + '_out.txt', 'a') as f:
        item1='%.4f\t%.4f\t%.4f' %(rmsdAll[RefProteinNo],rmsdCD[RefProteinNo],rmsdIR[RefProteinNo])
        f.write("%s\n" % item1)
#END: for RefProteinNo in range(0, numProteins, 1):


#np.savetxt(MyProteinDir + MyProteinFile[0:-4]+'_SELCON3-'+ UseRefSet +'_refit.txt', np.column_stack((wl_k, q, Mean_refit_prot)), fmt='%d \t %0.4f \t %0.4f', header='WL/k \t Spectrum \t Mean refit', comments='') #https://numpy.org/doc/stable/reference/generated/numpy.savetxt.html
    