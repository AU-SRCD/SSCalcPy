# -*- coding: utf-8 -*-
"""
Created on Tue Aug 30 11:19:53 2022

@author: au513
"""

#C    ------------------------------------------------------------------    
#C                      
#C    Program Name: CDsstr.for (Variable Selection)                      
#C    Written by Manavalan Parthasarathy, Arazdordi Toumadje, and W Curtis Johnson                    
#C                      
#C     ------------------------------------------------------------------    
#C
#C       Version 1.8.2 (formerly VSbsCD8.for)
#C
#C  Written for an IBM computer and to be compiled by Microsoft Fortran
#C  Powerstation
#C
#C  This version incorporates random choice into variable selection, uses a minimum
#c     basis set and self consistency, and chooses the combinations that give alpha
#C     helix similar to that predicted by Hennessey and Johnson SVD
#C     with self consistency.
#C
#C
#$NOTSTRICT                     
#C                     
#C                      
#C    THIS PROGRAM COMPUTES THE SECONDARY STRUCTURE OF A PROTEIN            
#C    FROM ITS CD, BY USING THE METHOD OF J. P. HENNESSEY, Jr.             
#C    AND W. C. JOHNSON, Jr. (Biochem. (1981) 20, 1085-1094) AND THE  
#C    SVD ALGORITHM ( L. A. Compton and W. C. Johnson, Jr, Anal.   
#C    Biochemistry (1986) 155, 155-167), FOR COMBINATIONS OF THE PROTEINS IN
#C    A MINIMUM BASIS SET (P. Manavalan  
#C    and W. C. Johnson, Jr., Anal Biochem (1987) 167, 76-85).  IT IS SELF CONSISTENT
#C    AS (N. Sreerama and R. W. Woody, Anal Biochem (1993) 209, 32-44), USES A
#C    MINIMUM BASIS SET AS IS IMPLICIT IN (I. H. M. van Stokkum, H. J. W. Spoelder,
#C    M. Bloemendal, R. van Grondelle, and F. C. A. Groen, Anal Biochem (1990)
#C    191, 110-118), AND CHOOSES THIS MINIMUM BASIS SET RANDOMLY AS SUGGESTED BY
#C    (B. Dalmas and W. H. Bannister, Anal Biochem (1995) 225, 39-48).                
#C                      
#C    THIS PROGRAM CALLS THE SUBROUTINE SVD1 (SINGULAR VALUE DECOMP).       
#C                      
#C    This version analyzes the CD of a set of up to test 100 proteins in file proCD.dta
#C    using basis proteins in file basCD.dta and their corresponding secondary stucture
#C    in file secstr.dta.  It randomly chooses NbsCDf basis proteins from the NbsCDi
#C    proteins in the basis set to analyze the test protein. It is self-consistent with
#C    the test protein in the SVD matrix. The first ncomb calculated sets of secondary
#C    structures with 0.952<=Tot=>1.05 and RMSE<=RmseMx are searched for alpha helix
#C    similar to the alpha helix predicted by H and J SVD.  The average of the choosen
#C    combinations are stored in file anal.out.  The average of the reconstructed CD
#C    for the choosen combinations are stored in reconCD.out. 
#C                      
#C        INPUT DATA VARIABLES                      
#C        ---------------------                      
#C                     
#C  BasCD:    CD spectra for the proteins in the basis set                   
#C            to be used in the analysis
#C  FssMin:   The greatest negative value allowed for a sec structure
#C  Fsstr:    Fractions of sec structure for each protein in basis set
#c  H:        Fraction of alpha helix from H and J SVD       
#C  icombf:   Program stops after icombf combinations if the number of successful
#C            combinations does not reach ncomb           
#C  icombi:   Initial combination number for this analysis, here set = 1
#C  Matxu and Matxv =1 to calc U and V in subroutine SVD                     
#C  NbasCD:   Number of proteins in basis set                      
#C  NbsCDv:   Number of basis CD vectors.  ONLY CD BASIS VECTORS WITH        
#C            SIGNIFICANT SINGULAR VALUES SHOULD BE USED IN AN ANALYSIS      
#C  NbsCD1:   Dimensioned number of CD basis proteins                      
#C  NbsCDf:   Final number basis proteins in set after some eliminated       
#C            for this analysis
#C  NCD:      Number of CD spectra undergoing SVD.  With self-consistency
#C            this is NbsCDf+1 (including the test protein)
#C  ncomb:    The number of successful combinations to be averaged
#C  Npro:     The number of test proteins in proCD.dta                               
#C  Nsstr:    Number of different secondary structures(eg. H,G,E,T,P,O =6)
#C  Nsstr1:   Dimensioned number of secondary structures     
#C  Numbas:   Numbers corresponding to the basis proteins being used in SVD                    
#C  Nwave:    Number of CD data points (number of wavelengths)              
#C  Nwave1:   Dimensioned number of CD data points                      
#C  Proname:  Name of each protein in input data (max is 40 characters)      
#C  ProCD:    CD data for protein to be analysed                      
#C  ProFss:   Predicted fraction of sec struc for protein analyzed; the average
#C            of all the successful combinations
#C  Psstr:    Predicted fractions stored for each successful combination
#C  RecnCD:   Average of the reconstructed CD for all the successful combinations
#C  RecnCD1:  Reconstructed CD stored for each successful combination          
#C  RMSE:     Average of Root Mean Square of Error at each wavelength        
#C            when comparing measured ProCD to CD reconstructed        
#C            from analysis for sec struc of one combination
#C  RmseMx:   The maximum value of RMSE for which a calculated set of
#C            secondary structures will be considered
#C
#C            *** The program will accept spectral data with the wavelengths
#C            in any direction, but all data must be consistent.
#C            The basis CD spectra supplied with this new version have the 
#C            INITIAL WAVELENGTH as the LONGEST WAVELENGTH, so your protein
#C            CD must begin with the longest wavelength if you are using the
#C            new basis CD spectra.
#C      
#C            SECONDARY STRUCTURES ARE:                      
#C            H:alpha-Helix; G:3/10-Helix; E:Extended-beta-strand;T:beta-Turns;
#C P:Polyproline-like 3/1-helix; O:OTHERS.  TOT is the total of all secondary
#C structures with any negative values excluded.        
#C 

import numpy as np
import sys

from CDSSTRPy_Functions import BVPRO4, CDSSTRPyFunc #CDSSTRFunctions.py contains the CDSSTR functions

#Define a functions just like
#def SelconsPy(
#   A1,F,q1, ShouldIPlot, ShouldIAlsoPlotSelcon2, 
#   LimWL, RefDataLowWL, MyProteinFile, MyProteinDir,
#   AU_or_BBK,CD_or_IR):    
#What is needed here:
#   LimWL, RefDataLowWL
#   AU_or_BBK  
#   MyProteinFile, MyProteinDir  

# **** USER INPUT HERE: ******

#**** Select protein CD file to be anayzed ****
#Data directory location with reference to py files
MyProteinDir = "../Data/timW/"
MyProteinDir = "TestProteins/"
MyProteinDir = "../TestProteins/"

#File name of spectrum to be analysed
MyProteinFile = "PeaGlobulin_1pt25.txt"
MyProteinFile = "ConA.txt"
MyProteinFile = "Myo.txt"
#MyProteinFile = "q.txt"

#**** Select Reference set. Uncomment only the set you use ****
#i.e. Select one of the four "UseRefSet = " lines below to be uncommented (no hash before the line)
#NCJ added Aug2023
#Enter 1, 2, 3, or 4 for the relevant number for the reference set
# 1 :SP175
# 2 :SMP180
# 3 :SP_AU-PCDDB
# 4 :SMP_AU-PCDDB
UseRefSetNum = 3 
AU_or_BBK =""
if (UseRefSetNum == 1):
    UseRefSet = "SP175" #SP175
    AU_or_BBK ="BBK"
if (UseRefSetNum == 2):
    UseRefSet = "SMP180" #SMP180
    AU_or_BBK ="BBK"
if (UseRefSetNum == 3):
    UseRefSet = "SP_AU-PCDDB"  #A_PCDDB, AU-F128_T/AU-F128_TS2  #Uses only the first 71 spectra in the reference file i.e. those in the SP175 set
    AU_or_BBK ="AU"
if (UseRefSetNum == 4):
    UseRefSet = "SMP_AU-PCDDB"      #A_PCDDB, AU-F128_T/AU-F128_TS2  #Uses all spectra in the reference file i.e. those in the SMP180 set      
    AU_or_BBK ="AU"
#End:NCJ added Aug2023

#Original data sets from B.A. Wallace /  A.M. Miles, Birckbeck
#UseRefSet = "SP175"
#UseRefSet = "SMP180"

#Updated data sets made by AU-SRCD from PCDDB and DSSP server 
#UseRefSet = "SP_AU-PCDDB"      #A_PCDDB, AU-F128_T/AU-F128_TS2  #Uses only the first 71 spectra in the reference file i.e. those in the SP175 set
#UseRefSet = "SMP_AU-PCDDB"      #A_PCDDB, AU-F128_T/AU-F128_TS2  #Uses all spectra in the reference file i.e. those in the SMP180 set

#Secondary structure classification:
Use_SS_Classification = "AU-F128_T"    # Turns: All T's in DSSP - Same definition as used in Dichroweb
#Use_SS_Classification = "AU-F128_TS2"  # Turns: 2 or more consecutive T and S in DSSP

#AU files for A and F matrix. Update when there are new files generated by AU-SRCD
# OLD name: A_AU_PCCDB_File = "A_AU-PCCDB-Nov22.txt" A_AU_PCCDB_File = "AU-A128_PCDDB-Nov22.txt"
CDRefDataPath = "../CDRefData/"
AU_A128_PCDDB_File = CDRefDataPath+"AU-A128_PCDDB-Nov22.txt"
AU_F128_T_File = CDRefDataPath+"AU-F128_T-Nov22.txt"
AU_F128_TS2_File = CDRefDataPath+"AU-F128_TS2-Nov22.txt"


#**** Lower limit of WL in the protein spectrum to be analyzed
#If you supply a two colomn dataset with WL in the 1st colomn, just leave it as 175 nm
LimWL = 175
#Note: A single colomn data set assumes the data starts at 240 nm and ends at LimWL.

# **** END OF USER INPUT ******

#Special output for Nyk: =1 or =0
specialOutputNyk = 1 #0
CDSSTRResultsList=[]

#For debugging
printMore = False #True

RefDataHighWL = 240
#Default Ref set is SP175
RefDataLowWL = 175 #SP175, this is the default
NumRefCDSpectra = 71 #SP175, this is the default
if (UseRefSet == "SMP180"):
    RefDataLowWL = 180 #SMP180
    NumRefCDSpectra = 128
    LimWL = 180
if (UseRefSet == "Another"):
    RefDataLowWL = 175 #Another
    NumRefCDSpectra = 71
if (UseRefSet == "SP_AU-PCDDB"):
    RefDataLowWL = 175 #SP_AU-PCDDB
    NumRefCDSpectra = 71
if (UseRefSet == "SMP_AU-PCDDB"):
    RefDataLowWL = 180 #SMP_AU-PCDDB
    NumRefCDSpectra = 128
    LimWL = 180


#The lower WL limt of the protein CD data to be analyzed
#LimWL = 175 #check what the limit actually is when reading in the query protein file.

if (UseRefSet =="SP175"):
    #Note: when loading witn np-loadtxt the type becomes a <class 'numpy.ndarray'>
    #Load the array with all the 71 secondary structures in SP175
    F =  np.loadtxt(CDRefDataPath+"F71.txt", dtype='f', delimiter='\t')  #SP175      #print(F.shape)
    #Load the array with all the 71 delta epsilon CD spectra in SP175
    A = np.loadtxt(CDRefDataPath+"A71.txt", dtype='f', delimiter='\t')  #SP175      #print(A1.shape)
elif (UseRefSet =="SMP180"):
    F =  np.loadtxt(CDRefDataPath+"F128.txt", dtype='f', delimiter='\t') #SMP180
    A = np.loadtxt(CDRefDataPath+"A128.txt", dtype='f', delimiter='\t') #SMP180
elif (UseRefSet == "SP_AU-PCDDB"):
    A2 = np.loadtxt(AU_A128_PCDDB_File, dtype='f', delimiter='\t')  # A2 = np.loadtxt(A_AU_PCCDB_File, dtype='f', delimiter='\t') #PCDDB smo data of SP175 + SMP180 protein spectra
    A=A2[:66, :71] #Rows 0->65 are 66 WLs 240 -> 175, Columns 0 -> 70 are 71 first SP175 proteins
    #update when we have AU-F128_T / AU-F128_TS2
    #F2 =  np.loadtxt("F71.txt", dtype='f', delimiter='\t') #AU-F128_T or AU-F128_TS2
    #F = F2
    if (Use_SS_Classification == "AU-F128_T"):
        F2 =  np.loadtxt(AU_F128_T_File, dtype='f', delimiter='\t') #AU-F128_T or AU-F128_TS2
        F=F2[:, :71]
    if (Use_SS_Classification == "AU-F128_TS2"):
        F2 =  np.loadtxt(AU_F128_TS2_File, dtype='f', delimiter='\t') #AU-F128_T or AU-F128_TS2
        F=F2[:, :71]
elif (UseRefSet == "SMP_AU-PCDDB"):
    A2 = np.loadtxt(AU_A128_PCDDB_File, dtype='f', delimiter='\t')  #A2 = np.loadtxt(A_AU_PCCDB_File, dtype='f', delimiter='\t') #PCDDB smo data of SP175 + SMP180 protein spectra
    A=A2[:61, :128] #Rows 0->60 are 61 WLs 240 -> 180, Columns 0 -> 127 are 128 SMP180 proteins
    #update when we have AU-F128_T / AU-F128_TS2
    #F2 =  np.loadtxt("F128.txt", dtype='f', delimiter='\t') #AU-F128_T or AU-F128_TS2
    #F = F2
    if (Use_SS_Classification == "AU-F128_T"):
        F2 =  np.loadtxt(AU_F128_T_File, dtype='f', delimiter='\t') #AU-F128_T or AU-F128_TS2
        F=F2[:, :128]
    if (Use_SS_Classification == "AU-F128_TS2"):
        F2 =  np.loadtxt(AU_F128_TS2_File, dtype='f', delimiter='\t') #AU-F128_T or AU-F128_TS2
        F=F2[:, :128]
else:
    print("No dataset selected. Stopping script")
    sys.exit()


#Find out of there are header lines in the protein data file
file1 = open(MyProteinDir + MyProteinFile, 'r')
Lines = file1.readlines()
 
skip = 0
# Strips the newline character
for line in Lines:
    if(line[0] == ";"):
        skip += 1
    else:
        pass 
print('Skipping first', skip, 'rows of protein data file', MyProteinDir + MyProteinFile)

#load in the delta Epsilon data of the protein you want to analyze
q1 = np.loadtxt(MyProteinDir + MyProteinFile, dtype='f', delimiter='\t', skiprows=skip)

#Find out which format the protein CD spectrum data has (single or double rows)
try:
    #If the file loaded intp q1 is just a single row, then if (q1.shape[1]>0): fail and an exception is thrown
    #Otherwise assume the file is two rows with WL in first row and Delta Epsilon in second row
    if (q1.shape[1]>0):
        #find which data point is 240 nm/RefDataHighWL: iRefDataHighWL
        iRefDataHighWL=-1
        for i in range(q1.shape[0]):
            if(q1[i,0]==RefDataHighWL):
                iRefDataHighWL=i
        if iRefDataHighWL>=0:
            qLength = q1.shape[0]-iRefDataHighWL
            q=np.zeros(qLength)
            q=q1[iRefDataHighWL:,1] #second row is De and first row is WL
            LimWL=max(int(q1[-1,0]),RefDataLowWL)
            print('LimWL is determined to be=',LimWL)
        else:
            print('Failed to find the position of the 240 nm data point')
            q=q1
        
except:
    #Assume that the file is just a single row of delta Epsilon starting at 240 nm
    #and limited to LimWL as set above
    tmpLimWL = LimWL
    LimWL = max(tmpLimWL, RefDataLowWL)
    q=q1

#Call the function and pray
CDSSTRResultsList = CDSSTRPyFunc(A, F, q, LimWL, RefDataLowWL, NumRefCDSpectra, MyProteinFile, MyProteinDir, AU_or_BBK, Use_SS_Classification)

#Write the results to file
if (len(CDSSTRResultsList) > 0):
    with open(MyProteinDir + MyProteinFile[0:-4]+'_CDSSTR-'+ UseRefSet +'_out.txt', 'w') as f:
        for item in CDSSTRResultsList:
            f.write("%s\n" % item)
     

#At this point we have 
#A, F, q, limWl, RefDataLowWL
#MyProteinFile is needed for plotting
#Also use NumRefCDSpectra, Use_SS_Classification, UseRefSet

#Old code commented out from here:
    #(and to end)
#def CDSSTRPyFunc(A, F, q, limWl, RefDataLowWL, MyProteinFile, MyProteinDir, AU_or_BBK):
#    pass

# #Parameter (Nwave1=100, NbsCD1=35, Nsstr1=6, Matxu=1, Matxv=1)
# Nwave1=240-LimWL+1 #C  Nwave1:   Dimensioned number of CD data points 
# NbsCD1=NumRefCDSpectra #71 #C  NbsCD1:   Dimensioned number of CD basis proteins 
# Nsstr1=6 #C  Nsstr1:   Dimensioned number of secondary structures 
# #Needed??? Matxu=1; Matxv=1

# #Parameter (NbsCDv=5,Nsstr=6,combi=1)
# NbsCDv=5 #C  NbsCDv:   Number of basis CD vectors.  ONLY CD BASIS VECTORS WITH SIGNIFICANT SINGULAR VALUES SHOULD BE USED IN AN ANALYSIS
# Nsstr=6 #C  Nsstr:    Number of different secondary structures(eg. H,G,E,T,P,O =6)
# combi=1 #?? icombi:   Initial combination number for this analysis, here set = 1
# #Parameter (RmseMx=0.25,FssMin=-0.03,NbsCDf=8)
# RmseMx=0.25 #RmseMx:   The maximum value of RMSE for which a calculated set of secondary structures will be considered
# FssMin=-0.025 #-0.03 #C  FssMin:   The greatest negative value allowed for a sec structure
# SumRule = 0.05 #Sum of SS should be within 5%
# NbsCDf=8 #C  NbsCDf:   Final number basis proteins in set after some eliminated for this analysis

# #Data which are input in CDsstr.for program:
# NbasCD = NumRefCDSpectra #71 #C  NbasCD:   Number of proteins in basis set
# Nwave = Nwave1 #the number of wavelengths (Nwave) in basCD.dta
# Npro = 1 # the number of protein CD spectra in proCD.dta (number of prots to be analyzed, we do one at a time!!!)
# ncomb = 400 #the number of successful combinations to be considered (ncomb),usually 100, but for all alpha proteins, use 400
# icombf = 100000 #the maximum number of trial combinations allowed', before the program stops trying to analyze a protein, usually 100000

# #Should include the option of limiting the lower WL
# #F1 = F
# A1 = A
# q2 = q
# NumWL = A1.shape[0] #size(A1,1); %Should be 66
# #NumWL = NumWL-(LimWL-175) #; %Should be 66 or lower
# NumWL = NumWL-(LimWL-RefDataLowWL) #; %Should be 66 or lower
# print('NumWL = ', NumWL)
# #Notice: A should have indices 0 to 66-1=65. The slice :NumWl (= :66) takes all indicies 0 to NumWL-1 which is correct.
# A=A1[:NumWL,:] #; %New ref prot matrix limited to LimWL
# q=q2[:NumWL] #; %Query prot matrix limited to LimWL
   
# #******************HJ******************
# #carry out svd on the data
# u, s, vh = np.linalg.svd(A, full_matrices=True)
# #       print ('u-shape: ',u.shape,'s-shape: ', s.shape,'vh-shape: ', vh.shape)
# # This gives u, s and vH such that A = u*S*vh
# # In this case U is 66x66 and vh is 71x71
# # However, s is a verctor with the number of data points (66) singular values. 
# # This has to be made into a matrix, in this case a number of data points x number of ref spectra (66x71)
# S = np.diag(s)
# S = np.zeros((u.shape[0],vh.shape[0])) #66x71
# S[:u.shape[0],:u.shape[0]] = np.diag(s) #S[:66,:66] = ...
# #       print('S-shape: ', S.shape)

# # For A = u*S*vh then U*S the basis spectra of the reference dataset 
# # ...and vh contains the least squares coefficients that refit the basis CD spectra to the experimental CD spectra

# #******************HJ5******************
# #Do the matrix multiplication, using the first 5 eigen vectors
# X = F @ np.transpose(vh)[:,:5] @ np.linalg.pinv(S[:5,:5]) @ np.transpose(u[:,:5])
# hj5 = np.matmul(X,q)

# #print('hj5 result')
# ##SStruct = ['Alpha-r','Alpha-d', 'Beta-r', 'Beta-d', 'Turns', 'Unord']
# ##format_list = [SStruct[hj5.tolist().index(item)]+'\t'+'{:.4f}' for item in hj5]
# ##strFormat = '\t'.join(format_list)
# ##print ('hj5 result\t\t  ', strFormat.format(*hj5), 'sum = ', '{:.4f}'.format(hj5.sum()))

# H = hj5[0] #H=ProFss(1)
# #Use both reg and distorted Helix
# H = hj5[0]+hj5[1] #H=ProFss(1)

# print(' Helix from H and J SVD with all basis proteins is ', '{:.4f}'.format(H)) #     Write(7,2060) H 2060 Format(' Helix from H and J SVD with all basis proteins is',F6.3) 


# ##Dimension ProCD(Nwave1),BasCD(Nwave1,NbsCD1),Fsstr(NbsCD1,Nsstr1)
# ##Dimension Psstr(Nsstr1,400),ProFss(Nsstr1)
# ##Dimension RecnCD(Nwave1),RecnCD1(Nwave1,400)
# ##Character*80 Proname 

# #w, h = 400, Nsstr1
# Psstr = np.zeros((Nsstr1, 400))
# if (printMore): print('Psstr-shape: ', Psstr.shape)
# RecnCD = np.zeros(Nwave1)
# if (printMore): print('RecnCD-shape: ', RecnCD.shape) #[0 for x in range(Nwave1)]
# RecnCD1 = np.zeros((Nwave1,400))
# if (printMore): print('RecnCD1-shape: ', RecnCD1.shape)
# ChosenRefsInValidSol=np.zeros((NbsCDf,400),  dtype=int)



# #C     Set up seed for generation random numbers
# #
# IproCD = 1 #Only one protein to be analyzed!!! = Npro, IproCD is a counter for Protein number IproCD
# NbsCDi=NbasCD
# icnt=0 #number of successful solutions
# BasCD = A; Fsstr = F; ProCD = q
# #print('Searching random comibations of ',NbsCDf,'spectra from the CD Reference set. Max iterations: ', icombf)
# #print('Looking for up to ',ncomb,' solutions')
# writeStr = 'Searching random comibations of ' +str(NbsCDf) +' spectra from the CD Reference set. Max iterations: ' +str(icombf) +'\n'
# writeStr += 'Looking for up to '+ str(ncomb) + ' solutions'
# print(writeStr)
# CDSSTRResultsList.append(writeStr)


# for icomb in range(1,icombf+1): #DO 1200 icomb=1,icombf
#     # Call BVPRO4(Nwave,BasCD,Fsstr,Nsstr,Psstr,RecnCD1,IproCD,RmseMx,FssMin,NbsCDv,ProCD,NbsCDi,NbsCDf,icnt,icomb,icombf)
#     #BasCD = A; Fsstr = F; ProCD = q
 
#     #Input  : Nwave,BasCD,Fsstr,Nsstr, IproCD, RmseMx, FssMin, NbsCDv, ProCD, NbsCDi, NbsCDf, icnt, icomb, icomb
#     #Return : Psstr, RecnCD1, icnt, selfConTrials
#     Psstr, RecnCD1, icnt, selfConTrials, ChosenRefsInValidSol = BVPRO4(Nwave,BasCD,Fsstr,Nsstr,Psstr,RecnCD1,IproCD,RmseMx,FssMin,SumRule,NbsCDv,ProCD,NbsCDi,NbsCDf,icnt,icomb,icombf,ChosenRefsInValidSol) # Call BVPRO4(Nwave,BasCD,Fsstr,Nsstr,Psstr,RecnCD1,IproCD,RmseMx,FssMin,NbsCDv,ProCD,NbsCDi,NbsCDf,icnt,icomb,icombf)
#     if (icomb % 1000 == 0):
#         print('Random Combination number:', icomb, '. Successful solutions so far: ', icnt, '. Selfcon itereations',selfConTrials)
#         #print('Refs used',ChosenRefsInValidSol[:,icnt-1])
#     if (icnt >= ncomb):  #If(icnt.ge.ncomb) GoTo 1205
#         break
# #END: for icomb in .... # 1200 Continue

# #print('Number of sum/fract rule Solutions is',icnt, 'found after', icomb,'combinations')
# writeStr = 'Number of Sum('+str(SumRule)+')/Fract('+str(FssMin)+')/RMSD('+str(RmseMx)+') rule solutions is ' + str(icnt) + ' found after '+ str(icomb) + ' combinations'
# print(writeStr)
# CDSSTRResultsList.append(writeStr)


# #How many solutions icntf where finally found
# if (icomb >= icombf):
#     icntf=icnt
# else:
#     icntf=ncomb #i.e. all desired 400 soliutions found
# #Endif

# if (icntf<5):
#     print('The number of acceptable solutions is less than five')
#     print(' There is a problem. Are you using delta epsilon units?')
#     sys.exit() #GoTo 1201
# #Endif
    

# # C
# # C  The program now searches the ncomb successful combinations for ones that have
# # C  a fraction of alpha helix similar to that predicted by Hennessey and Johnson
# # C  SVD (H).  For most proteins with H less than 0.15 it doesn't matter, but a few
# # C  have alpha helix closer to the minimum in the successful combinations while
# # C  others have alpha closer to H.  We average these two possibilities to get the
# # C  basis for the search.  Proteins with alpha greater than 0.65 need all the alpha
# # C  they can get, and we search the successful combinations for the maximum alpha.
# # C  For most proteins with H between 0.25 and 0.65 it doesn't matter, but a few
# # C  have alpha closer to the maximum in the successful combinations while
# # C  others have alpha closer to H.  We average these two possibilities to get the
# # C  basis for this search.  The search of the combinations for proteins with H
# # C  between 0.15 and 0.25 simply uses H as the basis.
# # C




# min=icntf/20
# if (min < 4): min=4
# print('Will search for at least ', min, ' solutions with the highest H(r)+H(d), unless H<0.25')

# Hmin = 0.0
# Hmax = 1.0
# hel = Psstr[0,:icntf]
# hel = Psstr[0,:icntf]  + Psstr[1,:icntf] #both dist and reg helix
# #First for H<0.15
# if (H < 0.15):
#     #Find the min Helix in Psstr
#     writeStr = 'H = %8.3f' % H + ' is H <= 0.15'
#     Hmin = np.min(hel)
#     if (Hmin < 0.0): Hmin=0.0
#     Havg=(Hmin+H)/2.0
#     Hmax=Havg+0.03
#     Hmin=Havg-0.03
# elif (H > 0.25):
#     #print('H = %8.3f' % H,' is > 0.25')
#     Hmax = np.max(hel) #Do 50 I=1,icntf 50     If(Psstr(1,I).gt.Hmax) Hmax=Psstr(1,I)
    
#     #How many are >= Hmax?
#     icnt=0 #51
#     while (icnt < min):
#         for I in range(icntf): #Do 52 I=1,icntf
#             if (hel[I] >= Hmax): icnt=icnt+1 #IF(Psstr(1,I).ge.Hmax) icnt=icnt+1
#         #Next I  /  52   Continue
#         if (icnt < min): #Do it again with a smaller Hmax. Note: first time it will only find one, icnt=1
#             Hmax=Hmax-0.001
#             icnt = 0
#     #End while / GoTo 51
#     #Now the Hmax is set such that at least 20 solutions have hel higher than Hmax
    
#     if (H < 0.65): #53     If(H.lt.0.65) then
#         #print('H = %8.3f' % H,' is 0.25 < H < 0.65')
#         writeStr ='H = %8.3f' % H + ' is 0.25 < H < 0.65'
#         Hmax=(Hmax+H)/2.0
#         Hmin=Hmax-0.03
#         Hmax=Hmax+0.03
#     else:
#         #print('H = %8.3f' % H,' is >= 0.65')
#         writeStr = 'H = %8.3f' % H + ' is >= 0.65'
#         Hmin=Hmax-0.03
#         Hmax=1.0
#     #Endif
# else:
#     #In this case the H must be between 0.15 and 0.25
#     #print('H =', H,' is 0.15 <= H <= 0.25')
#     writeStr = 'H = %8.3f' % H + ' is 0.15 <= H <= 0.25'
#     Hmin=H-0.03
#     Hmax=H+0.03
# #Endif
# #print('Use Hmin, Hmax = %8.3f %8.3f ' % (Hmin,Hmax))
# writeStr += '\n' + 'Use Hmin, Hmax = %8.3f %8.3f ' % (Hmin,Hmax)
# print(writeStr)
# CDSSTRResultsList.append(writeStr)

    
# #Now Hmin and Hmax if found for the different H cases: Find the solutions
# ProFss = np.zeros(6)
# CDSSTR3_solutions=[] #All solutions of sec struct which adhere to the helix rules (like SelCon3)
# CDSSTR3_RecnCD=[]
# CDSSTR3_RefsUsed=[] #basenumbers on the choosen ref with valid solution: ChosenRefsInValidSol
# isstr=0 #numnber of good solutions for Helix test
# for I in range(icntf): #Do 45 I=1,icntf
#     if ((hel[I] >= Hmin) and (hel[I] <= Hmax)): #If((Psstr(1,I).ge.Hmin).and.(Psstr(1,I).le.Hmax)) then
#         ProFss = ProFss + Psstr[:,I] # Do 46 K=1,Nsstr  46     ProFss(K)=Psstr(K,I)+ProFss(K)
#         RecnCD = RecnCD + RecnCD1[:,I] # Do 141 J=1,Nwave   141     RecnCD(J)=RecnCD(J)+RecnCD1(J,I)
#         CDSSTR3_solutions.append(Psstr[:,I])
#         CDSSTR3_RecnCD.append(RecnCD1[:,I])
#         CDSSTR3_RefsUsed.append(ChosenRefsInValidSol[:,I])
#         isstr=isstr+1
#     #Endif
# #Next I / 45   Continue
   
# if (isstr == 0):
#     print('Random selection failed to give good helix, try again with larger ncomb')
#     sys.exit()    #Goto 1201
# #Endif
    
# ProFss = ProFss / isstr #Do 47 K=1,Nsstr 47   ProFss(K)=ProFss(K)/isstr
# RecnCD = RecnCD / isstr #Do 142 J=1,Nwave 142   RecnCD(J)=RecnCD(J)/isstr

# #Write(9,1010) (RecnCD(J),J=1,Nwave) #Write(7,2011) (ProFss(K),K=1,Nsstr),isstr #C
    
# #Are there valid solutions using the helix rule?
# if (len(CDSSTR3_solutions) == 0):
#     print('No valid CDSSTR solutions based on the helix rule')
# else:
#     CDSSTR3_solutionsMat = np.transpose(np.array(CDSSTR3_solutions))
#     CDSSTR3_RecnCDMat = np.transpose(np.array(CDSSTR3_RecnCD))
         
#     #print('Number of solutions in the Hmin, Hmax range ', len(CDSSTR3_solutions),'(',isstr,')')
#     writeStr = 'Number of solutions in the Hmin, Hmax range is ' + str(len(CDSSTR3_solutions))# +' ('+str(isstr)+') '
#     CDSSTRResultsList.append(writeStr)
#     print(writeStr)
    
#     #NCJ 26/1/2023
#     TS_Classification = "All T's in DSSP"
#     if  (UseRefSet == "SM_AU-PCDDB") or (UseRefSet == "SMP_AU-PCDDB"):
#           if (Use_SS_Classification == "AU-F128_T"):
#               TS_Classification = "All T's in DSSP"
#           else:
#               TS_Classification = "2 or more consecutive T and S in DSSP"
            
#     writeStr = 'Reference set: ' + UseRefSet + "  Classification: " + TS_Classification
#     CDSSTRResultsList.append(writeStr)
#     print(writeStr)
    
#     #Show the solution including standart deviation
#     stDev = np.std(CDSSTR3_solutionsMat, axis=1)
#     StructureName  = ['Alpha(r)', 'Alpha(d)', 'Beta(r) ', 'Beta(d) ', 'Turns   ', 'Other   ']
#     for i in range(len(ProFss)):
#         writeStr = StructureName[i] + ' %4.1f ' % (100.0*ProFss[i]) + '%' + ' +/- %4.1f' % (100.0*stDev[i])
#         print(writeStr)
#         CDSSTRResultsList.append(writeStr)
#         #print(StructureName[i], ' %4.1f ' % (100.0*ProFss[i]) , '%' ,' +/- %4.1f' % (100.0*stDev[i]))
#     writeStr = 'Sum = %4.1f ' % (100.0*ProFss.sum())
#     print(writeStr)
#     CDSSTRResultsList.append(writeStr)
#     #print('Sum = %4.1f ' % (100.0*ProFss.sum()))

#    #RMS:
#     rmsdRefit = np.sqrt(((RecnCD - ProCD) ** 2).mean()) #rmsdRefit=rmsd(Mean_refit_prot,q);
#     writeStr = 'RMSD between Query and Refitted spectrum ' + '{:.4f}'.format(rmsdRefit)
#     print(writeStr)
#     CDSSTRResultsList.append(writeStr)
    
#    #NRMSD:
#     #Wallace et.al. (2003) Prot. Sci., Mao and Wallace (1984)
#     #SVH: Wallace Prot. Sci. above must be "Analyses of circular dichroism spectra of membrane proteins" Prot. Sci. 2009 !!!
#     #NRMSD = Sum((Exp-Cal)**2) / Sum(Exp**2) SVH: and the square root of that.
#     NRMSD = ((RecnCD - ProCD) ** 2).sum() 
#     NRMSD = np.sqrt(NRMSD / ((ProCD ** 2).sum()))
#     writeStr = 'NRMSD between Query and Refitted spectrum ' + '{:.4f}'.format(NRMSD)
#     print(writeStr)
#     CDSSTRResultsList.append(writeStr)
    
#     #average number of helices and strands per 100 residues and average helix and strand lengths
#     AvgHelicesPer100Res = 100.0*ProFss[1]/4 #Distored is index=1 (index=0 is regular)
#     AvgStrandsPer100Res = 100.0*ProFss[3]/2
#     AvgHelLength = 100*(ProFss[0] + ProFss[1])/AvgHelicesPer100Res
#     AvgStrandLength = 100*(ProFss[2] + ProFss[3])/AvgStrandsPer100Res
#     #Write the results
#     writeStr1 = 'Average Number of Helices per 100 residues: ' + '%5.2f' % AvgHelicesPer100Res
#     writeStr2 = 'Average Number of Strands per 100 residues: ' + '%5.2f' % AvgStrandsPer100Res
#     writeStr3 = 'Average length of Helices: ' + '%5.2f' % AvgHelLength
#     writeStr4 = 'Average length of Strands: ' + '%5.2f' % AvgStrandLength
#     print(writeStr1); print(writeStr2); print(writeStr3); print(writeStr4)
#     CDSSTRResultsList.append(writeStr1); CDSSTRResultsList.append(writeStr2); CDSSTRResultsList.append(writeStr3); CDSSTRResultsList.append(writeStr4)
    
#     #For Alison: Which proteins where used in solution: CDSSTR3_RefsUsed from ChosenRefsInValidSol
#     CDSSTR3_RefsUsedList = []
#     CDSSTR3_RefsUsedList.extend(CDSSTR3_RefsUsed[0]) #add the first solutions ref spectra numbers. Not as an element (using append) but as a list (using Extend)
#     for I in range(1,len(CDSSTR3_RefsUsed)): #CDSSTR3_RefsUsed[up to icntf] contains [8] size lists
#         CDSSTR3_RefsUsedList.extend(CDSSTR3_RefsUsed[I])
#     CDSSTR3_RefsUsedListUnique = np.unique(np.array(CDSSTR3_RefsUsedList)) #CDSSTR3_RefsUsedListTmpArr)
#     #Write the results
#     writeStr1 = 'List of reference proteins numbers used in solutions:'
#     writeStr2 = CDSSTR3_RefsUsedListUnique
#     writeStr3 = 'Length of list: ' + '{}'.format(len(CDSSTR3_RefsUsedListUnique))
#     print(writeStr1); print(writeStr2); print(writeStr3); CDSSTRResultsList.append(writeStr1); CDSSTRResultsList.append(writeStr2); CDSSTRResultsList.append(writeStr3)
    
#     #Finally make the format for Nyk
#     if (specialOutputNyk == 1):
#         formatSpec = '{:4.1f}'
#         strFormat = formatSpec + '\t'
#         writeStrNyk = ""
#         #All the secondary structure solutions
#         for i in range(len(ProFss)):
#             writeStrNyk += strFormat.format(100*ProFss[i])
#         #The sum of solutions
#         writeStrNyk += strFormat.format(100*np.sum(ProFss))
#         #All the stDev of the solutions
#         for i in range(len(ProFss)):
#             writeStrNyk += strFormat.format(100*stDev[i])
#         #Which selcon
#         writeStrNyk += 'CDSSTR3'+'\t'
#         #number of Selcon3 solutions
#         writeStrNyk += str(len(CDSSTR3_solutions))+'\t'
#         #Add the RMSD for Nyk
#         writeStrNyk += '{:.4f} \t'.format(rmsdRefit)
#         writeStrNyk += '{:.4f}'.format(NRMSD)+'\t'
#         writeStrNyk += str(icomb) 
    
#     if (specialOutputNyk == 1):
        
#         #The long data string for Nyk is the last added to the return list.
#         CDSSTRResultsList.append(' ')
#         CDSSTRResultsList.append('Output for Nyk:')
#         listStr = 'AlphaR'+'\t'+'AlphaD'+'\t'+'BetaR'+'\t'+'BetaD'+'\t'+'Turns'+'\t'+'Unord'+'\t'+'Sum'+'\t'
#         listStr += 'AR Dev'+'\t'+'AD Dev'+'\t'+'BR Dev'+'\t'+'BD Dev'+'\t'+'T Dev'+'\t'+'Un Dev'+'\t'
#         listStr += 'CDSSTR?'+'\t'+'#solu'+'\t'+'RMSD'+'\t'+'NRMSD'+'\t'+'Combinations'
#         CDSSTRResultsList.append(listStr)
#         CDSSTRResultsList.append(writeStrNyk)
    
# #end of else in if (len(CDSSTR3_solutions) == 0):

# #Write the results to file
# if (len(CDSSTRResultsList) > 0):
#     with open(MyProteinDir + MyProteinFile[0:-4]+'_CDSSTR-'+ UseRefSet +'_out.txt', 'w') as f:
# #    with open(MyProteinDir + MyProteinFile[0:-4]+'_out.txt', 'w') as f:
#         for item in CDSSTRResultsList:
#             f.write("%s\n" % item)
    
    
# #Plotting. Only if valid Helix validated solutions.
# if (len(CDSSTR3_solutions) > 0):
#     import matplotlib.pyplot as plt
#     import matplotlib as mpl        
    
#     #Makesure that the plots are in separate windows and not inline in IPython. For inline, use gui='inline'
#     #(Command directly in consle is either %matplotlib qt or %matplotlib inline)
#     try:
#         import IPython
#         shell = IPython.get_ipython()
#         shell.enable_matplotlib(gui='qt')
#     except:
#         pass
    
#     #%Plot solutions
#     plt.close('all')
#     #mpl.rcParams['toolbar'] = 'None'
#     px = 1/mpl.rcParams['figure.dpi']  # pixel in inches
#     figStartX = 100
#     figStartY = 50
#     figWidth = 900
#     figHeight = 600
#     figHSpacer = 5 #25
#     figWSpacer = 50
#     mpl.rcParams['font.size'] = 16
#     mpl.rcParams['legend.fontsize'] = 12
       
#     #Define WL range for plot
#     wl=np.arange(240, LimWL-1, -1) #240:-1:LimWL; %Vector of wavelengths
    
#     #PLOT 1:  %PLOT query prot and all refits
#     fig1 = plt.figure(1, figsize=[figWidth*px,figHeight*px], frameon = False) #figure(1); %Define figure
#     fig1.canvas.manager.window.move(figStartX,figStartY)
#     plt.axhline(linewidth=1, color='grey', dashes=[2,2])
#     plt.plot(wl, ProCD, 'bo', markersize  = 5.0, mfc='none', label = MyProteinFile)
#     #plt.plot(wl, ProCD, 'bo', markersize  = 5.0, mfc='none', label = "Protein")
#     for i in range(CDSSTR3_RecnCDMat.shape[1]):
#         if (i==0):
#             plt.plot(wl, CDSSTR3_RecnCDMat[:,i], '-', linewidth = 1.0, label='Refitted spectra')
#         else:
#             plt.plot(wl, CDSSTR3_RecnCDMat[:,i], '-', linewidth = 1.0)
#     plt.legend()
    
#     #%Define lots of things like linewidths, markers.. and add labels and titles        
#     plt.ylabel(r'$\Delta \epsilon$ [M$^{-1}$ cm$^{-1}$ residue$^{-1}$]') #        ylabel('\Delta \epsilon [mdeg /M /cm /residue]');
#     plt.xlabel('Wavelength [nm]') #        xlabel('Wavelength [nm]');
#     plt.title('Spectra of Query protein and all valid CDSSTR3 refits')
#     fig1.show()
    
#     #PLOT 2: %PLOT query prot and MEAN of refits
#     fig2 = plt.figure(2, figsize=[figWidth*px,figHeight*px], frameon = False)   #figure(1); %Define figure
#     fig2.canvas.manager.window.move(figStartX + figWidth + figHSpacer, figStartY)
#     plt.axhline(linewidth=1, color='grey', dashes=[2,2])
#     plt.plot(wl, ProCD, 'bo', markersize  = 5.0, mfc='none', label = MyProteinFile)
#     #plt.plot(wl, ProCD, 'bo', markersize  = 5.0, mfc='none', label = "Protein")
#     plt.plot(wl, RecnCD, '-', linewidth = 1.0, label="Mean refit")
#     plt.legend()
    
#     #%Define lots of things like linewidths, markers.. and add labels and titles        
#     plt.ylabel(r'$\Delta \epsilon$ [M$^{-1}$ cm$^{-1}$ residue$^{-1}$]') #        ylabel('\Delta \epsilon [mdeg /M /cm /residue]');
#     plt.xlabel('Wavelength [nm]') #        xlabel('Wavelength [nm]');
#     plt.title('Spectra of Query protein and mean of all valid CDSSTR3 refits')
#     fig2.show()
#     plt.savefig(MyProteinDir + MyProteinFile[0:-4]+'_CDSSTR-'+ UseRefSet +'_fig.png')
    
    

