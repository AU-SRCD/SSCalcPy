# -*- coding: utf-8 -*-
"""
Created on Tue Sep 13 15:05:32 2022

@author: au513
"""

import numpy as np
import sys

printMore = True
printMore = False

def CDSSTRPyFunc(A, F, q, LimWL, RefDataLowWL, NumRefCDSpectra, MyProteinFile, MyProteinDir, AU_or_BBK, CD_or_IR, Use_SS_Classification):
    #This also uses NumRefCDSpectra, Use_SS_Classification, UseRefSet
    
    specialOutputNyk = 1 #0  #Special output for Nyk: =1 or =0
    includeProteinsUsedForCalc = 0 #1
    CDSSTRResultsList=[]
    
    # Filenames are
    # Prot_SELCON3-SP_AU-PCDDB_xxx (SP_AU-PCDDB)
    # Prot_SELCON3-SMP_AU-PCDDB_xxx (SMP_AU-PCDDB)
    # Prot_SELCON3-SP_BBK_xxx
    # Prot_SELCON3-SMP_BBK_xxx
    # Prot_SELCON3-IR_xxx
    
    UseRefSet ="SP_AU-PCDDB" #default is AU SP175 set
    if (RefDataLowWL == 180):
        UseRefSet ="SMP_AU-PCDDB"
    if (AU_or_BBK == 'BBK'):
        #UseRefSet = "SP175" #SP175
        UseRefSet = "SP_BBK" #SP175
        if (RefDataLowWL == 180):
            UseRefSet = "SMP_BBK" #SMP180
            #UseRefSet = "SMP180" #SMP180
    if (CD_or_IR == 'IR'):
        UseRefSet ="IR"
    
    
    CD_or_IR = 'CD' #for now only CD
    if (CD_or_IR == 'IR'):
        IR_LowKInRefSet = 1600; IR_HighKInRefSet = 1800; IR_RefSetSteps = 1; IR_RefSetNumDataPoints = 1 + ((IR_HighKInRefSet-IR_LowKInRefSet)/IR_RefSetSteps)
        wl_k = np.linspace(IR_LowKInRefSet, IR_HighKInRefSet, int(IR_RefSetNumDataPoints))
    else:
        RefDataHighWL = 240
        wl_k = np.arange(RefDataHighWL, LimWL-1, -1)
    
    #Parameter (Nwave1=100, NbsCD1=35, Nsstr1=6, Matxu=1, Matxv=1)
    Nwave1=240-LimWL+1 #C  Nwave1:   Dimensioned number of CD data points 
    #Not used: NbsCD1=NumRefCDSpectra #71 #C  NbsCD1:   Dimensioned number of CD basis proteins 
    Nsstr1=6 #C  Nsstr1:   Dimensioned number of secondary structures 
    #Needed??? Matxu=1; Matxv=1
    
    #Parameter (NbsCDv=5,Nsstr=6,combi=1)
    NbsCDv=5 #C  NbsCDv:   Number of basis CD vectors.  ONLY CD BASIS VECTORS WITH SIGNIFICANT SINGULAR VALUES SHOULD BE USED IN AN ANALYSIS
    Nsstr=6 #C  Nsstr:    Number of different secondary structures(eg. H,G,E,T,P,O =6)
    #Not used: combi=1 #?? icombi:   Initial combination number for this analysis, here set = 1
    #Parameter (RmseMx=0.25,FssMin=-0.03,NbsCDf=8)
    RmseMx=0.25 #RmseMx:   The maximum value of RMSE for which a calculated set of secondary structures will be considered
    FssMin=-0.025 #-0.03 #C  FssMin:   The greatest negative value allowed for a sec structure
    SumRule = 0.05 #Sum of SS should be within 5%
    NbsCDf=8 #C  NbsCDf:   Final number basis proteins in set after some eliminated for this analysis
    
    #Data which are input in CDsstr.for program:
    NbasCD = NumRefCDSpectra #71 #C  NbasCD:   Number of proteins in basis set
    Nwave = Nwave1 #the number of wavelengths (Nwave) in basCD.dta
    #Not used: Npro = 1 # the number of protein CD spectra in proCD.dta (number of prots to be analyzed, we do one at a time!!!)
    ncomb = 400 #the number of successful combinations to be considered (ncomb),usually 100, but for all alpha proteins, use 400
    icombf = 100000 #the maximum number of trial combinations allowed', before the program stops trying to analyze a protein, usually 100000
    
    #Should include the option of limiting the lower WL
    #F1 = F
    A1 = A
    q2 = q
    NumWL = A1.shape[0] #size(A1,1); %Should be 66
    #NumWL = NumWL-(LimWL-175) #; %Should be 66 or lower
    NumWL = NumWL-(LimWL-RefDataLowWL) #; %Should be 66 or lower
    print('NumWL = ', NumWL)
    #Notice: A should have indices 0 to 66-1=65. The slice :NumWl (= :66) takes all indicies 0 to NumWL-1 which is correct.
    A=A1[:NumWL,:] #; %New ref prot matrix limited to LimWL
    q=q2[:NumWL] #; %Query prot matrix limited to LimWL
       
    #******************HJ******************
    #carry out svd on the data
    u, s, vh = np.linalg.svd(A, full_matrices=True)
    #       print ('u-shape: ',u.shape,'s-shape: ', s.shape,'vh-shape: ', vh.shape)
    # This gives u, s and vH such that A = u*S*vh
    # In this case U is 66x66 and vh is 71x71
    # However, s is a verctor with the number of data points (66) singular values. 
    # This has to be made into a matrix, in this case a number of data points x number of ref spectra (66x71)
    S = np.diag(s)
    S = np.zeros((u.shape[0],vh.shape[0])) #66x71
    S[:u.shape[0],:u.shape[0]] = np.diag(s) #S[:66,:66] = ...
    #       print('S-shape: ', S.shape)
    
    # For A = u*S*vh then U*S the basis spectra of the reference dataset 
    # ...and vh contains the least squares coefficients that refit the basis CD spectra to the experimental CD spectra
    
    #******************HJ5******************
    #Do the matrix multiplication, using the first 5 eigen vectors
    X = F @ np.transpose(vh)[:,:5] @ np.linalg.pinv(S[:5,:5]) @ np.transpose(u[:,:5])
    hj5 = np.matmul(X,q)
    
    #print('hj5 result')
    ##SStruct = ['Alpha-r','Alpha-d', 'Beta-r', 'Beta-d', 'Turns', 'Unord']
    ##format_list = [SStruct[hj5.tolist().index(item)]+'\t'+'{:.4f}' for item in hj5]
    ##strFormat = '\t'.join(format_list)
    ##print ('hj5 result\t\t  ', strFormat.format(*hj5), 'sum = ', '{:.4f}'.format(hj5.sum()))
    
    H = hj5[0] #H=ProFss(1)
    #Use both reg and distorted Helix
    H = hj5[0]+hj5[1] #H=ProFss(1)
    
    print(' Helix from H and J SVD with all basis proteins is ', '{:.4f}'.format(H)) #     Write(7,2060) H 2060 Format(' Helix from H and J SVD with all basis proteins is',F6.3) 
    
    
    ##Dimension ProCD(Nwave1),BasCD(Nwave1,NbsCD1),Fsstr(NbsCD1,Nsstr1)
    ##Dimension Psstr(Nsstr1,400),ProFss(Nsstr1)
    ##Dimension RecnCD(Nwave1),RecnCD1(Nwave1,400)
    ##Character*80 Proname 
    
    #w, h = 400, Nsstr1
    Psstr = np.zeros((Nsstr1, 400))
    if (printMore): print('Psstr-shape: ', Psstr.shape)
    RecnCD = np.zeros(Nwave1)
    if (printMore): print('RecnCD-shape: ', RecnCD.shape) #[0 for x in range(Nwave1)]
    RecnCD1 = np.zeros((Nwave1,400))
    if (printMore): print('RecnCD1-shape: ', RecnCD1.shape)
    ChosenRefsInValidSol=np.zeros((NbsCDf,400),  dtype=int)
    
    
    
    #C     Set up seed for generation random numbers
    #
    IproCD = 1 #Only one protein to be analyzed!!! = Npro, IproCD is a counter for Protein number IproCD
    NbsCDi=NbasCD
    icnt=0 #number of successful solutions
    BasCD = A; Fsstr = F; ProCD = q
    #print('Searching random comibations of ',NbsCDf,'spectra from the CD Reference set. Max iterations: ', icombf)
    #print('Looking for up to ',ncomb,' solutions')
    writeStr = 'Searching random comibations of ' +str(NbsCDf) +' spectra from the CD Reference set. Max iterations: ' +str(icombf) +'\n'
    writeStr += 'Looking for up to '+ str(ncomb) + ' solutions'
    print(writeStr)
    CDSSTRResultsList.append(writeStr)
    
    
    for icomb in range(1,icombf+1): #DO 1200 icomb=1,icombf
        # Call BVPRO4(Nwave,BasCD,Fsstr,Nsstr,Psstr,RecnCD1,IproCD,RmseMx,FssMin,NbsCDv,ProCD,NbsCDi,NbsCDf,icnt,icomb,icombf)
        #BasCD = A; Fsstr = F; ProCD = q
     
        #Input  : Nwave,BasCD,Fsstr,Nsstr, IproCD, RmseMx, FssMin, NbsCDv, ProCD, NbsCDi, NbsCDf, icnt, icomb, icomb
        #Return : Psstr, RecnCD1, icnt, selfConTrials
        Psstr, RecnCD1, icnt, selfConTrials, ChosenRefsInValidSol = BVPRO4(Nwave,BasCD,Fsstr,Nsstr,Psstr,RecnCD1,IproCD,RmseMx,FssMin,SumRule,NbsCDv,ProCD,NbsCDi,NbsCDf,icnt,icomb,icombf,ChosenRefsInValidSol) # Call BVPRO4(Nwave,BasCD,Fsstr,Nsstr,Psstr,RecnCD1,IproCD,RmseMx,FssMin,NbsCDv,ProCD,NbsCDi,NbsCDf,icnt,icomb,icombf)
        if (icomb % 1000 == 0):
            print('Random Combination number:', icomb, '. Successful solutions so far: ', icnt, '. Selfcon itereations',selfConTrials)
            #print('Refs used',ChosenRefsInValidSol[:,icnt-1])
        if (icnt >= ncomb):  #If(icnt.ge.ncomb) GoTo 1205
            break
    #END: for icomb in .... # 1200 Continue
    
    #print('Number of sum/fract rule Solutions is',icnt, 'found after', icomb,'combinations')
    writeStr = 'Number of Sum('+str(SumRule)+')/Fract('+str(FssMin)+')/RMSD('+str(RmseMx)+') rule solutions is ' + str(icnt) + ' found after '+ str(icomb) + ' combinations'
    print(writeStr)
    CDSSTRResultsList.append(writeStr)
    
    
    #How many solutions icntf where finally found
    if (icomb >= icombf):
        icntf=icnt
    else:
        icntf=ncomb #i.e. all desired 400 soliutions found
    #Endif
    
    if (icntf<5):
        print('The number of acceptable solutions is less than five')
        print(' There is a problem. Are you using delta epsilon units?')
        sys.exit() #GoTo 1201
    #Endif
        
    
    # C
    # C  The program now searches the ncomb successful combinations for ones that have
    # C  a fraction of alpha helix similar to that predicted by Hennessey and Johnson
    # C  SVD (H).  For most proteins with H less than 0.15 it doesn't matter, but a few
    # C  have alpha helix closer to the minimum in the successful combinations while
    # C  others have alpha closer to H.  We average these two possibilities to get the
    # C  basis for the search.  Proteins with alpha greater than 0.65 need all the alpha
    # C  they can get, and we search the successful combinations for the maximum alpha.
    # C  For most proteins with H between 0.25 and 0.65 it doesn't matter, but a few
    # C  have alpha closer to the maximum in the successful combinations while
    # C  others have alpha closer to H.  We average these two possibilities to get the
    # C  basis for this search.  The search of the combinations for proteins with H
    # C  between 0.15 and 0.25 simply uses H as the basis.
    # C
    
    
    
    
    min=icntf/20
    if (min < 4): min=4
    print('Will search for at least ', min, ' solutions with the highest H(r)+H(d), unless H<0.25')
    
    Hmin = 0.0
    Hmax = 1.0
    hel = Psstr[0,:icntf]
    hel = Psstr[0,:icntf]  + Psstr[1,:icntf] #both dist and reg helix
    #First for H<0.15
    if (H < 0.15):
        #Find the min Helix in Psstr
        writeStr = 'H = %8.3f' % H + ' is H <= 0.15'
        Hmin = np.min(hel)
        if (Hmin < 0.0): Hmin=0.0
        Havg=(Hmin+H)/2.0
        Hmax=Havg+0.03
        Hmin=Havg-0.03
    elif (H > 0.25):
        #print('H = %8.3f' % H,' is > 0.25')
        Hmax = np.max(hel) #Do 50 I=1,icntf 50     If(Psstr(1,I).gt.Hmax) Hmax=Psstr(1,I)
        
        #How many are >= Hmax?
        icnt=0 #51
        while (icnt < min):
            for I in range(icntf): #Do 52 I=1,icntf
                if (hel[I] >= Hmax): icnt=icnt+1 #IF(Psstr(1,I).ge.Hmax) icnt=icnt+1
            #Next I  /  52   Continue
            if (icnt < min): #Do it again with a smaller Hmax. Note: first time it will only find one, icnt=1
                Hmax=Hmax-0.001
                icnt = 0
        #End while / GoTo 51
        #Now the Hmax is set such that at least 20 solutions have hel higher than Hmax
        
        if (H < 0.65): #53     If(H.lt.0.65) then
            #print('H = %8.3f' % H,' is 0.25 < H < 0.65')
            writeStr ='H = %8.3f' % H + ' is 0.25 < H < 0.65'
            Hmax=(Hmax+H)/2.0
            Hmin=Hmax-0.03
            Hmax=Hmax+0.03
        else:
            #print('H = %8.3f' % H,' is >= 0.65')
            writeStr = 'H = %8.3f' % H + ' is >= 0.65'
            Hmin=Hmax-0.03
            Hmax=1.0
        #Endif
    else:
        #In this case the H must be between 0.15 and 0.25
        #print('H =', H,' is 0.15 <= H <= 0.25')
        writeStr = 'H = %8.3f' % H + ' is 0.15 <= H <= 0.25'
        Hmin=H-0.03
        Hmax=H+0.03
    #Endif
    #print('Use Hmin, Hmax = %8.3f %8.3f ' % (Hmin,Hmax))
    writeStr += '\n' + 'Use Hmin, Hmax = %8.3f %8.3f ' % (Hmin,Hmax)
    print(writeStr)
    CDSSTRResultsList.append(writeStr)
    
        
    #Now Hmin and Hmax if found for the different H cases: Find the solutions
    ProFss = np.zeros(6)
    CDSSTR3_solutions=[] #All solutions of sec struct which adhere to the helix rules (like SelCon3)
    CDSSTR3_RecnCD=[]
    CDSSTR3_RefsUsed=[] #basenumbers on the choosen ref with valid solution: ChosenRefsInValidSol
    isstr=0 #numnber of good solutions for Helix test
    for I in range(icntf): #Do 45 I=1,icntf
        if ((hel[I] >= Hmin) and (hel[I] <= Hmax)): #If((Psstr(1,I).ge.Hmin).and.(Psstr(1,I).le.Hmax)) then
            ProFss = ProFss + Psstr[:,I] # Do 46 K=1,Nsstr  46     ProFss(K)=Psstr(K,I)+ProFss(K)
            RecnCD = RecnCD + RecnCD1[:,I] # Do 141 J=1,Nwave   141     RecnCD(J)=RecnCD(J)+RecnCD1(J,I)
            CDSSTR3_solutions.append(Psstr[:,I])
            CDSSTR3_RecnCD.append(RecnCD1[:,I])
            CDSSTR3_RefsUsed.append(ChosenRefsInValidSol[:,I])
            isstr=isstr+1
        #Endif
    #Next I / 45   Continue
       
    if (isstr == 0):
        print('Random selection failed to give good helix, try again with larger ncomb')
        sys.exit()    #Goto 1201
    #Endif
        
    ProFss = ProFss / isstr #Do 47 K=1,Nsstr 47   ProFss(K)=ProFss(K)/isstr
    RecnCD = RecnCD / isstr #Do 142 J=1,Nwave 142   RecnCD(J)=RecnCD(J)/isstr
    
    #Write(9,1010) (RecnCD(J),J=1,Nwave) #Write(7,2011) (ProFss(K),K=1,Nsstr),isstr #C
        
    #Are there valid solutions using the helix rule?
    if (len(CDSSTR3_solutions) == 0):
        print('No valid CDSSTR solutions based on the helix rule')
    else:
        CDSSTR3_solutionsMat = np.transpose(np.array(CDSSTR3_solutions))
        CDSSTR3_RecnCDMat = np.transpose(np.array(CDSSTR3_RecnCD))
             
        #print('Number of solutions in the Hmin, Hmax range ', len(CDSSTR3_solutions),'(',isstr,')')
        writeStr = 'Number of solutions in the Hmin, Hmax range is ' + str(len(CDSSTR3_solutions))# +' ('+str(isstr)+') '
        CDSSTRResultsList.append(writeStr)
        print(writeStr)
        
        #NCJ 26/1/2023
        TS_Classification = "All T's in DSSP"
        if  (UseRefSet == "SP_AU-PCDDB") or (UseRefSet == "SMP_AU-PCDDB"):
              if (Use_SS_Classification == "AU-F128_T"):
                  TS_Classification = "All T's in DSSP"
              else:
                  TS_Classification = "2 or more consecutive T and S in DSSP"
                
        writeStr = 'Reference set: ' + UseRefSet + "  Classification: " + TS_Classification
        CDSSTRResultsList.append(writeStr)
        print(writeStr)
        
        #Show the solution including standart deviation
        stDev = np.std(CDSSTR3_solutionsMat, axis=1)
        #StructureName  = ['Alpha(r)', 'Alpha(d)', 'Beta(r) ', 'Beta(d) ', 'Turns   ', 'Other   ']
        StructureName  = ['AlphaR', 'AlphaD', 'BetaR ', 'BetaD ', 'Turns   ', 'Other   ']
        for i in range(len(ProFss)):
            writeStr = StructureName[i] + ' %4.1f ' % (100.0*ProFss[i]) + '%' + ' +/- %4.1f' % (100.0*stDev[i])
            print(writeStr)
            CDSSTRResultsList.append(writeStr)
            #print(StructureName[i], ' %4.1f ' % (100.0*ProFss[i]) , '%' ,' +/- %4.1f' % (100.0*stDev[i]))
        writeStr = 'Sum = %4.1f ' % (100.0*ProFss.sum())
        print(writeStr)
        CDSSTRResultsList.append(writeStr)
        #print('Sum = %4.1f ' % (100.0*ProFss.sum()))
    
       #RMS:
        rmsdRefit = np.sqrt(((RecnCD - ProCD) ** 2).mean()) #rmsdRefit=rmsd(Mean_refit_prot,q);
        writeStr = 'RMSD between Query and Refitted spectrum ' + '{:.4f}'.format(rmsdRefit)
        print(writeStr)
        CDSSTRResultsList.append(writeStr)
        
       #NRMSD:
        #Wallace et.al. (2003) Prot. Sci., Mao and Wallace (1984)
        #SVH: Wallace Prot. Sci. above must be "Analyses of circular dichroism spectra of membrane proteins" Prot. Sci. 2009 !!!
        #NRMSD = Sum((Exp-Cal)**2) / Sum(Exp**2) SVH: and the square root of that.
        NRMSD = ((RecnCD - ProCD) ** 2).sum() 
        NRMSD = np.sqrt(NRMSD / ((ProCD ** 2).sum()))
        writeStr = 'NRMSD between Query and Refitted spectrum ' + '{:.4f}'.format(NRMSD)
        print(writeStr)
        CDSSTRResultsList.append(writeStr)
        
        #average number of helices and strands per 100 residues and average helix and strand lengths
        AvgHelicesPer100Res = 100.0*ProFss[1]/4 #Distored is index=1 (index=0 is regular)
        AvgStrandsPer100Res = 100.0*ProFss[3]/2
        AvgHelLength = 100*(ProFss[0] + ProFss[1])/AvgHelicesPer100Res
        AvgStrandLength = 100*(ProFss[2] + ProFss[3])/AvgStrandsPer100Res
        #Write the results
        writeStr1 = 'Average Number of Helices per 100 residues: ' + '%5.2f' % AvgHelicesPer100Res
        writeStr2 = 'Average Number of Strands per 100 residues: ' + '%5.2f' % AvgStrandsPer100Res
        writeStr3 = 'Average length of Helices: ' + '%5.2f' % AvgHelLength
        writeStr4 = 'Average length of Strands: ' + '%5.2f' % AvgStrandLength
        print(writeStr1); print(writeStr2); print(writeStr3); print(writeStr4)
        CDSSTRResultsList.append(writeStr1); CDSSTRResultsList.append(writeStr2); CDSSTRResultsList.append(writeStr3); CDSSTRResultsList.append(writeStr4)
        
        #For Alison: Which proteins where used in solution: CDSSTR3_RefsUsed from ChosenRefsInValidSol
        if (includeProteinsUsedForCalc == 1):
            CDSSTR3_RefsUsedList = []
            CDSSTR3_RefsUsedList.extend(CDSSTR3_RefsUsed[0]) #add the first solutions ref spectra numbers. Not as an element (using append) but as a list (using Extend)
            for I in range(1,len(CDSSTR3_RefsUsed)): #CDSSTR3_RefsUsed[up to icntf] contains [8] size lists
                CDSSTR3_RefsUsedList.extend(CDSSTR3_RefsUsed[I])
            CDSSTR3_RefsUsedListUnique = np.unique(np.array(CDSSTR3_RefsUsedList)) #CDSSTR3_RefsUsedListTmpArr)
            #Write the results
            writeStr1 = 'List of reference proteins numbers used in solutions:'
            writeStr2 = CDSSTR3_RefsUsedListUnique
            writeStr3 = 'Length of list: ' + '{}'.format(len(CDSSTR3_RefsUsedListUnique))
            print(writeStr1); print(writeStr2); print(writeStr3); CDSSTRResultsList.append(writeStr1); CDSSTRResultsList.append(writeStr2); CDSSTRResultsList.append(writeStr3)
        
        #Finally make the format for Nyk
        if (specialOutputNyk == 1):
            formatSpec = '{:4.1f}'
            strFormat = formatSpec + '\t'
            writeStrNyk = ""
            #All the secondary structure solutions
            for i in range(len(ProFss)):
                writeStrNyk += strFormat.format(100*ProFss[i])
            #The sum of solutions
            writeStrNyk += strFormat.format(100*np.sum(ProFss))
            #All the stDev of the solutions
            for i in range(len(ProFss)):
                writeStrNyk += strFormat.format(100*stDev[i])
            #Which selcon
            writeStrNyk += 'CDSSTR3'+'\t'
            #number of Selcon3 solutions
            writeStrNyk += str(len(CDSSTR3_solutions))+'\t'
            #Add the RMSD for Nyk
            writeStrNyk += '{:.4f} \t'.format(rmsdRefit)
            writeStrNyk += '{:.4f}'.format(NRMSD)+'\t'
            writeStrNyk += str(icomb)+'\t' 
            writeStrNyk +=  UseRefSet+'\t' #added 10/1/2024
            writeStrNyk +=  MyProteinFile #added 10/1/2024
           
        
        if (specialOutputNyk == 1):
            
            #The long data string for Nyk is the last added to the return list.
            CDSSTRResultsList.append(' ')
            CDSSTRResultsList.append('Output for Nyk:')
            #listStr = 'AlphaR'+'\t'+'AlphaD'+'\t'+'BetaR'+'\t'+'BetaD'+'\t'+'Turns'+'\t'+'Unord'+'\t'+'Sum'+'\t'
            listStr = 'AlphaR'+'\t'+'AlphaD'+'\t'+'BetaR'+'\t'+'BetaD'+'\t'+'Turns'+'\t'+'Other'+'\t'+'Sum'+'\t'
            #listStr += 'AR Dev'+'\t'+'AD Dev'+'\t'+'BR Dev'+'\t'+'BD Dev'+'\t'+'T Dev'+'\t'+'Un Dev'+'\t'
            listStr += 'AR Dev'+'\t'+'AD Dev'+'\t'+'BR Dev'+'\t'+'BD Dev'+'\t'+'T Dev'+'\t'+'Ot Dev'+'\t'
            listStr += 'CDSSTR?'+'\t'+'#solu'+'\t'+'RMSD'+'\t'+'NRMSD'+'\t'+'Combinations'
            CDSSTRResultsList.append(listStr)
            CDSSTRResultsList.append(writeStrNyk)
        
    #end of else in if (len(CDSSTR3_solutions) == 0):
    
        
    #Plotting. Only if valid Helix validated solutions.
    if (len(CDSSTR3_solutions) > 0):
        import matplotlib.pyplot as plt
        import matplotlib as mpl
        #Makesure that the plots are in separate windows and not inline in IPython
        mpl.use('qtagg') #matplotlib.use('qtagg')
        
        #%Plot solutions
        plt.close('all')
        #mpl.rcParams['toolbar'] = 'None'
        px = 1/mpl.rcParams['figure.dpi']  # pixel in inches
        figStartX = 400 #100
        figStartY = 50
        figWidth = 900
        figHeight = 600
        figHSpacer = 5 #25
        figWSpacer = 50
        mpl.rcParams['font.size'] = 16
        mpl.rcParams['legend.fontsize'] = 12
           
        #Define WL range for plot
        wl=np.arange(240, LimWL-1, -1) #240:-1:LimWL; %Vector of wavelengths
        
        #PLOT 1:  %PLOT query prot and all refits
        fig1 = plt.figure(1, figsize=[figWidth*px,figHeight*px])#, frameon = False) #figure(1); %Define figure
        fig1.canvas.manager.window.move(figStartX,figStartY)
        plt.axhline(linewidth=1, color='grey', dashes=[2,2])
        plt.plot(wl, ProCD, 'bo', markersize  = 5.0, mfc='none', label = MyProteinFile)
        #plt.plot(wl, ProCD, 'bo', markersize  = 5.0, mfc='none', label = "Protein")
        for i in range(CDSSTR3_RecnCDMat.shape[1]):
            if (i==0):
                plt.plot(wl, CDSSTR3_RecnCDMat[:,i], '-', linewidth = 1.0, label='Refitted spectra')
            else:
                plt.plot(wl, CDSSTR3_RecnCDMat[:,i], '-', linewidth = 1.0)
        plt.legend()
        
        #%Define lots of things like linewidths, markers.. and add labels and titles        
        plt.ylabel(r'$\Delta \epsilon$ [M$^{-1}$ cm$^{-1}$ residue$^{-1}$]') #        ylabel('\Delta \epsilon [mdeg /M /cm /residue]');
        plt.xlabel('Wavelength [nm]') #        xlabel('Wavelength [nm]');
        plt.title('Spectra of Query protein and all valid CDSSTR3 refits')
        fig1.show()
        
        #PLOT 2: %PLOT query prot and MEAN of refits
        fig2 = plt.figure(2, figsize=[figWidth*px,figHeight*px])#, frameon = False)   #figure(1); %Define figure
        fig2.canvas.manager.window.move(figStartX + figWidth + figHSpacer, figStartY)
        plt.axhline(linewidth=1, color='grey', dashes=[2,2])
        plt.plot(wl, ProCD, 'bo', markersize  = 5.0, mfc='none', label = MyProteinFile)
        plt.plot(wl, RecnCD, '-', linewidth = 1.0, label="Mean refit")
        plt.legend()
        
        #%Define lots of things like linewidths, markers.. and add labels and titles        
        plt.ylabel(r'$\Delta \epsilon$ [M$^{-1}$ cm$^{-1}$ residue$^{-1}$]') #        ylabel('\Delta \epsilon [mdeg /M /cm /residue]');
        plt.xlabel('Wavelength [nm]') #        xlabel('Wavelength [nm]');
        plt.title('Spectra of Query protein and mean of all valid CDSSTR3 refits')
        fig2.show()
        plt.savefig(MyProteinDir + MyProteinFile[0:-4]+'_CDSSTR-'+ UseRefSet +'_fig.png')

    #END: if (len(CDSSTR3_solutions) > 0):

    #Before ending, save the RecnCD (mean reconstructed) solution together with the org CD spectrum ProCD (= q)
    np.savetxt(MyProteinDir + MyProteinFile[0:-4]+'_CDSSTR-'+ UseRefSet +'_refit.txt', np.column_stack((wl_k, q, RecnCD)), fmt='%d \t %0.4f \t %0.4f', header='WL/k \t Spectrum \t Mean refit', comments='') #https://numpy.org/doc/stable/reference/generated/numpy.savetxt.html
    
    return CDSSTRResultsList
# END: def CDSSTRPyFunc():

#c********************************************************************
#c Subroutine to Calculate the Number of Segments
#	Subroutine Segment(Str,MxSS,AHel,BStr,HelLn,BetLn)
#	dimension Str(MxSS)
#	AHel = Str(2) * 100.0 / 4.0
#	BStr = Str(4) * 100.0 / 2.0
#	If(AHel.LT.0.0) AHel = 0.0
#	If(BStr.LT.0.0) BStr = 0.0
#	HelLn = 0.0
#	If(AHel.NE.0.0)HelLn = (Str(1) + Str(2)) * 100.0 /  AHel
#	BetLn = 0.0
#	If(BStr.NE.0.0)BetLn = (Str(3) + Str(4)) * 100.0 /  BStr    
#	End
#c ******************************************************************
def Segment(Str):
    pass


#def Rndm(NbsCDi,NbsCDf,Nwave,Nsstr,BasCD,Fsstr,Numbas,CD,SS):
def Rndm(NbsCDi,NbsCDf,Nwave,Nsstr,BasCD,Fsstr): #,Numbas,CD,SS):
    #    Parameter (Nwave1=100, NbsCD1=35, Nsstr1=6)
    #Nwave1=100; NbsCD1=71; Nsstr1=6
    #      Dimension Numbas(NbsCD1),BasCD(Nwave1,NbsCD1),Fsstr(NbsCD1,Nsstr1),CD(Nwave1,NbsCD1),SS(NbsCD1,Nsstr1)
    Numbas = np.zeros(NbsCDf,  dtype=int) #np.zeros(NbsCD1)
    #CD = np.zeros((Nwave1,NbsCD1)); SS=np.zeros((Nsstr1,NbsCD1)) #SS=np.zeros((NbsCD1,Nsstr1))
    CD = np.zeros((Nwave,NbsCDi)); SS=np.zeros((Nsstr,NbsCDi))
    ChosenRefs = np.zeros(NbsCDf) #Which are the randomly chosen reference spectra
    
    import random
    for I in range(0,NbsCDf): #range(1,NbsCDf+1):     #      DO 100 I=1,NbsCDf
        #I1=0 #   10 I1=0
        #      Call Random(ranval)
        #scrap=ranval*NbsCDi #      scrap=ranval*Real(NbsCDi)
        getNewRandomBase = True
        while getNewRandomBase:
            Numbas[I]= random.randint(0, NbsCDi-1) #Int(scrap)+1 #      Numbas(I)=Int(scrap)+1
            na = Numbas[I]
            #print('na = Numbas[',I,'] = ',na)
            #is na in Numbas[:I]?
            if na in Numbas[:I]:
                getNewRandomBase = True
            else:
                getNewRandomBase = False
        #end while
        ChosenRefs[I] = na
        CD[:,I] = BasCD[:,na]
        SS[:,I] = Fsstr[:,na]
    #end for
    return Numbas, CD, SS, ChosenRefs
        
        #I1=I1+1     #   20 I1=I1+1
        #      IF(I1.ge.I) Goto 30
        #      IF(Numbas(I).eq.Numbas(I1))Goto 10      
        #      Goto 20
        #   30 na=Numbas(I)
        #      DO 3 J=1,Nwave
        #    3 CD(J,I)=BasCD(J,na)
        #      DO 4 K=1,Nsstr
        #    4 SS(I,K)=Fsstr(na,K)
    #  100 Continue
    #C
    #C      Write(8,1000) icnt,(Numbas(I),I=1,NbsCDf)
    # 1000 Format('icnt=',I4,2x,30I3)
#      Return
#      End


def BVPRO4(Nwave,BasCD,Fsstr,Nsstr,Psstr,RecnCD1,IproCD,RmseMx,FssMin,SumRule,NbsCDv,ProCD,NbsCDi,NbsCDf,icnt,icomb,icombf,ChosenRefsInValidSol):
    #Input Nwave,BasCD,Fsstr,Nsstr, RmseMx, FssMin, NbsCDv, ProCD, NbsCDi, NbsCDf, icnt, icomb, icomb, IproCD=protein number (=1)
    #return Psstr, RecnCD1
    
    #    Parameter (Nwave1=100, NbsCD1=35, Nsstr1=6, Matxu=1, Matxv=1)
    #Nwave1=240-175+1 #C  Nwave1:   Dimensioned number of CD data points 
    #NbsCD1=71 #C  NbsCD1:   Dimensioned number of CD basis proteins 
    #Nsstr1=6 #C  Nsstr1:   Dimensioned number of secondary structures 
    #    Dimension CD(Nwave1,NbsCD1), U(Nwave1,NbsCD1), V(Nwave1,NbsCD1)       
    #    Dimension S(NbsCD1), Work(NbsCD1),RecnCD1(Nwave1,400)                     
    #    Dimension BasCD(Nwave1,NbsCD1), Fsstr(NbsCD1,Nsstr1)       
    #    Dimension Numbas(NbsCD1), SS(NbsCD1,Nsstr1), Psstr(Nsstr1,400)              
    #    Dimension C(NbsCD1), ProFss(Nsstr1), ProCD(Nwave1),RecnCD(Nwave1)      
    #    LOGICAL Matu, Matv   
    #CD = np.zeros((Nwave1,NbsCD1))
    #Numbas = np.zeros(NbsCD1)
    #SS = np.zeros((NbsCD1,Nsstr1))
    #Psstr = np.zeros((Nsstr,400))
    #RecnCD1 = np.zeros((Nwave,400))
    
    #??IF(icnt.gt.icombf) STOP
    Numbas, CD, SS, ChosenRefs = Rndm(NbsCDi,NbsCDf,Nwave,Nsstr,BasCD,Fsstr) #,Numbas,CD,SS, ChosenRefs)
    NCD=NbsCDf+1
    CD[:,NCD-1] = ProCD[:] #; print('CD-shape: ', CD.shape) # DO 15 J=1,Nwave,    15 CD(J,NCD)=ProCD(J)
    #Initial "guess". Why these values??? Can we do better
    
    SS[0,NCD-1] = 0.36 #(NCD,1)=0.36
    SS[1,NCD-1] = 0.06 #SS(NCD,2)=0.06
    SS[2,NCD-1] = 0.18 #SS(NCD,3)=0.18
    SS[3,NCD-1] = 0.12 #SS(NCD,4)=0.12
    SS[4,NCD-1] = 0.07 #SS(NCD,5)=0.07
    SS[5,NCD-1] = 0.21 #SS(NCD,6)=0.21 
    #Note according to CDsstr_unix, this should be replaced by SS[k,NCD-1]=0.17 i.e. same for all structures
    #print('SS-shape: ', SS.shape); print('SS[:,:NCD]-shape: ', SS[:,:NCD].shape)
    
    #Call SVD(Nwave,NCD,CD,S,Matu,U,Matv,V,ierr,Work)
    #truncate the matrix to be inverted to only have the NCD first elements 0, 1, ..., NCD-1
    u, s, vh = np.linalg.svd(CD[:,:NCD], full_matrices=False) #u, s, vh = np.linalg.svd(CD[:,:NCD], full_matrices=True)
    #print('u-shape: ', u.shape); print('s-shape: ', s.shape); print('vh-shape: ', vh.shape)
    #construct the S matrix
    
    #S = np.zeros((u.shape[0],vh.shape[0])) #66x71
    #S[:s.shape[0],:s.shape[0]] = np.diag(s)
    S=np.diag(s)
    
    #The new sec struct ProFss = X * ProCD where X = SS * v * (s^-1) * uT
    
    #Pre calc some of X, Xsub = v*s^-1*uT
    #Xsub = np.transpose(vh)[:,:NbsCDv] @ np.linalg.pinv(S[:NbsCDv,:NbsCDv]) @ np.transpose(u[:,:NbsCDv]) #; print('X-shape: ', X.shape)
    
    #No need to use the Penrose Inverse pinv as S is already a square matrix
    #Xsub = np.transpose(vh)[:,:NbsCDv] @ np.linalg.inv(S[:NbsCDv,:NbsCDv]) @ np.transpose(u[:,:NbsCDv]) #; print('X-shape: ', X.shape)
    
    #do not transpose the entier vH and then select the NbsCDv colomn, but select first, then transpose
    Xsub = np.transpose(vh[:NbsCDv,:]) @ np.linalg.inv(S[:NbsCDv,:NbsCDv]) @ np.transpose(u[:,:NbsCDv]) #; print('X-shape: ', X.shape)
    
    #Another way of multiplying 1/s onto vHT is to np.multiply (or *) the vHT with the vector 1/s = (1/s1, 1/s2, ...)
    #Xsub1=np.transpose(vh[:NbsCDv,:]) * (1/(s[:NbsCDv]))
    #Xsub = Xsub1 @ np.transpose(u[:,:NbsCDv])
    
    #Repeat up to 20 times for selfconsistency
    selfConTrials = 0
    for trial in range(20):
        #predict the structure based on
        #ProFss = Fsstr * V * S(inverse) * U(transpose) * ProCD
        #Do the matrix multiplication, using the all tjhe NbsCDv eigen vectors

        #S[:u.shape[0],:u.shape[0]] = np.diag(s) #S[:66,:66] = ...
        #X = SS[:,:NCD] @ np.transpose(vh)[:,:NbsCDv] @ np.linalg.pinv(S[:NbsCDv,:NbsCDv]) @ np.transpose(u[:,:NbsCDv]) #; print('X-shape: ', X.shape)
        X = SS[:,:NCD] @ Xsub #Note that @ is shorthand for np.matmul
        ProFss = np.matmul(X,ProCD)
        #How close is the new SS compated to the old? Are we self con?
        rmdsSelfCon = np.sqrt(((ProFss - SS[:,NCD-1]) ** 2).mean())
        #Update SS to new value calc value for structure of the protein
        SS[:,NCD-1] = ProFss #;print('trial ', trial, 'has ProFss =', ProFss)
        
        #%check to see if selfconsistent
        if rmdsSelfCon < 0.0025: #        if((rmsd(F(:,1),mean(np_solutions,2)))<0.0025) %If the new solution is close to the previous, we have selfconsistency, SVH 7/12-17
            #if (icomb % 1000 == 0):
            #    print('Random Combination number:', icomb, '. Iterations for Self const: ', trial, 'rmdsSelfCon = ', rmdsSelfCon)
            selfConTrials = trial
            break
        
        #Now repeat with the new values for the guess of the proteins Sec Struc
    #Next trial
    
    #print('After 20 selfcon,  ProFss =', ProFss)
    #calculate the sum of all fractions. Should not include if a fraction is <0!
    Tot = np.sum(SS[:,NCD-1]) #; print('After 20 selfcon,  Tot =', Tot)
    
    #C     Reconstructed CD spectrum corresponding to predicted sec struc        
    
    #tmp: Is tot within range
    sum_rule = SumRule #0.05
    fract_rule = FssMin #-0.025
    testMinFrac = SS[:,NCD-1].min() 
    if ((abs(1-Tot) <= sum_rule) and (testMinFrac >=fract_rule)): #If((Tot.ge.0.952).and.(Tot.le.1.05)
        #print('New ok sum/fract rule solution ',icnt) #, 'ProFss=',ProFss)
        #Reconstruct the spectrum
        
        
        # SSE=0.0
        # for J in range(1,Nwave+1): 
        #     RecnCD[J-1]=0.0
        #     for JA in range(1,Nwave+1):
        #         for K in range(1,NbsCDv+1):
        #             RecnCD[J-1]=RecnCD[J-1]+u[J-1,K-1]*u[JA-1,K-1]*ProCD[JA-1] # U(J,K)*U(JA,K)*ProCD(JA)
        #     se = RecnCD[J-1]-ProCD[JA-1]
        #     SSE = SSE + se*se
        # #continue J
        # rmse=np.sqrt(SSE/Nwave)
            
        # DO 90 J=1,Nwave
        #     RecnCD(J)=0.0
        #     DO 88 JA=1,Nwave                      
        #     DO 88 K=1,NbsCDv                      
        #     88 RecnCD(J)=RecnCD(J)+U(J,K)*U(JA,K)*ProCD(JA)    
            
        #     se=RecnCD(J)-ProCD(J)                                           
        #     SSE=SSE+se*se                      
        # 90 Continue                      
        # RMSE=SQRT(SSE/FLOAT(Nwave))
        
        
        
        
        refitCD = u[:,:NbsCDv] @ S[:NbsCDv,:NbsCDv] @ vh[:NbsCDv,:] ##refit = U(:,1:basis)*S(1:basis,1:basis)*Vt(1:basis,:);        
        rmse = np.sqrt(((refitCD[:,-1] - ProCD) ** 2).mean()) #refits(i)=rmsd(refit(:,1),q); %SVH 22/9-2006 
        
        if printMore: print('New ok sum/fract rule solution with rmse', rmse) #,' (',icnt,')') #, 'ProFss=',ProFss)
        
#        refitCD_Sreerama_uuT = u[:,:NbsCDv] @ np.transpose(u[:,:NbsCDv])
#        refitCD_Sreerama = np.matmul(refitCD_Sreerama_uuT,ProCD)
#        rmse_Sreerama = np.sqrt(((refitCD_Sreerama - ProCD) ** 2).mean())
#        
#        print('New ok sum/fract rule solution with rmse_Sreerama', rmse_Sreerama)
        
        # with open('out_'+str(icomb)+'.txt', 'w') as f:
        #     for t in range(refitCD.shape[0]):
        #         f.write("{:.4f}".format(refitCD[t,0])+'\n')
        if rmse <= RmseMx:
            icnt = icnt+1
            Psstr[:,icnt-1]= ProFss
            RecnCD1[:,icnt-1]=refitCD[:,-1] #RecnCD1(J,icnt)=RecnCD(J)
            ChosenRefsInValidSol[:,icnt-1] = ChosenRefs
            if printMore: print('New ok sum/fract/RmseMx rule solution ',icnt) #, 'ProFss=',ProFss)
#            with open("OutFiles/"+'outRmseLow'+str(icnt)+'.txt', 'w') as f:
#                for t in range(refitCD.shape[0]): #for t in range(RecnCD.shape[0]): 
#                    f.write("{:.4f}".format(refitCD[t,-1])+'\n')
        
    
    
    return Psstr, RecnCD1, icnt, selfConTrials, ChosenRefsInValidSol
    
 