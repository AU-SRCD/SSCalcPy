# -*- coding: utf-8 -*-
"""
Created on Fri Aug 13 14:18:52 2021

@author: au513

The code for this function was originally made in Matlab in 2005 by the research group of
Prof. B.A.Wallace, Birckbeck Collage, London

The MatLab code was updated in 2006-2007 to allow plotting and calculation of
the mean refitted spectrum to the query protein by S.V. Hoffmann, Aarhus University, Denmark

This code was adapted to Python by S.V. Hoffmann, Aarhus University, Denmark in 2021

"""
import numpy as np
#def SelconsPy(A1,F,q1, ShouldIPlot, ShouldIAlsoPlotSelcon2, LimWL, RefDataLowWL):
#def SelconsPy(A1,F,q1, ShouldIPlot, ShouldIAlsoPlotSelcon2, LimWL, RefDataLowWL, MyProteinFile):
def SelconsPy(A1,F,q1, ShouldIPlot, ShouldIAlsoPlotSelcon2, LimWL, RefDataLowWL, MyProteinFile, MyProteinDir,AU_or_BBK,CD_or_IR):

    #Selcon2 solutions depend on the RMSD between the protein spectrum and the reconstructed spectrum
    #This is traditionally 0.25 delta epsilon units for CD
    #For IR we need a new rule
    refit_rule=0.25
    if (CD_or_IR == 'IR'): refit_rule=0.025 #0.25 / 10
    #The rule for CD_or_IR == 'CD-IR' is the same as CD for now (IR scaled)
    
    #Special for CD-IR: LimWL, RefDataLowWL = 175,176,... (depend on input data) / 175
    #      ... and tHe CD parts of A and q passed are already truncated at CD LimWL
    #Special for IR: LimWL, RefDataLowWL = 1800 / 1800
    #For IR and CD-IR make the definition of the k values here for use later
    IR_LowKInRefSet = 1600; IR_HighKInRefSet = 1800; IR_RefSetSteps = 1
    if 'CD-IR' in CD_or_IR: #if (CD_or_IR == 'CD-IR'):
        IR_LowKInRefSet = 1600; IR_HighKInRefSet = 1800; IR_RefSetSteps = 1
        if '2cm' in CD_or_IR:
            IR_LowKInRefSet = 1600; IR_HighKInRefSet = 1720; IR_RefSetSteps = 2
    IR_RefSetNumDataPoints = 1 + ((IR_HighKInRefSet-IR_LowKInRefSet)/IR_RefSetSteps)
    
    
    # Filenames are
    # Prot_SELCON3-SP_AU-PCDDB_xxx
    # Prot_SELCON3-SMP_AU-PCDDB_xxx
    # Prot_SELCON3-SP_BBK_xxx
    # Prot_SELCON3-SMP_BBK_xxx
    # Prot_SELCON3-IR_xxx
    # Prot_SELCON3-CD-IR_xxx
    
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
    if (CD_or_IR == 'CD-IR'):
        UseRefSet ="CD-IR"
    if (CD_or_IR == 'CD-IR2cm'):
        UseRefSet = "CD-IR2cm"
     
    # DataSet = AU_or_BBK
    # if AU_or_BBK == 'AU':
    #     DataSet +='-PCDDB'
    # if (CD_or_IR == 'IR'):
    #     DataSet = 'IR'
    
    #Special output for Nyk: =1 or =0
    specialOutputNyk = 1
    includeProteinsUsedForCalc = 0 #1 #= 1: will add the Protein list to the return string AND if =2 also print to console
    #returns HJ, selcon, selcon2, and selcon3 solutions 
    #A data matrix of reference protein CD spectra
    #F structure matrix
    #q query protein spectrum
    #LimWL is the lower limit WL. Should be 175 and up to say 200 nm
    
    if (CD_or_IR == 'CD'):
        #Should include the option of limiting the lower WL (CD ONLY)
        NumWL = A1.shape[0] #size(A1,1); %Should be 66
        #NumWL = NumWL-(LimWL-175) #; %Should be 66 or lower
        NumWL = NumWL-(LimWL-RefDataLowWL) #; %Should be 66 or lower
        #Notice: A should have indices 0 to 66-1=65. The slice :NumWl (= :66) takes all indicies 0 to NumWL-1 which is correct.
        A=A1[:NumWL,:] #; %New ref prot matrix limited to LimWL
        q=q1[:NumWL] #; %Query prot matrix limited to LimWL
    else:
        A=A1
        q=q1
        
    #Also make a wl/k vector for saving later on
    if (CD_or_IR == 'IR'):
        #IR_LowKInRefSet = 1600; IR_HighKInRefSet = 1800; IR_RefSetSteps = 1; IR_RefSetNumDataPoints = 1 + ((IR_HighKInRefSet-IR_LowKInRefSet)/IR_RefSetSteps)
        k = np.linspace(IR_LowKInRefSet, IR_HighKInRefSet, int(IR_RefSetNumDataPoints))
        wl_k = k
    else:
        RefDataHighWL = 240
        wl = np.arange(RefDataHighWL, LimWL-1, -1)
        wl_k = wl
        #print(wl_k)
    
    if 'CD-IR' in CD_or_IR: #if (CD_or_IR == 'CD-IR'):
        k = np.linspace(IR_LowKInRefSet, IR_HighKInRefSet, int(IR_RefSetNumDataPoints))
        #wl is already calculated above
        wl_k = np.concatenate([wl,k])
    
    if 'CD-IR' in CD_or_IR:
        #Detect if the IR part has been scaled to zero
        if (len(q)<len(wl_k)):
            wl_k = wl
    
        
        
        
    
   
    #******************HJ******************
    #carry out svd on the data
    u, s, vh = np.linalg.svd(A, full_matrices=True)
    #   print ('u-shape: ',u.shape,'s-shape: ', s.shape,'vh-shape: ', vh.shape)
    # This gives u, s and vH such that A = u*S*vh
    # In this case U is 66x66 and vh is 71x71
    # However, s is a verctor with the number of data points (66) singular values. 
    # This has to be made into a matrix, in this case a number of data points x number of ref spectra (66x71)
    S = np.diag(s)
    S = np.zeros((u.shape[0],vh.shape[0])) #66x71 (IR: 200x50)
    s_size = min(u.shape[0],vh.shape[0])
    S[:s_size,:s_size] = np.diag(s) #S[:66,:66] = ... #S[:u.shape[0],:u.shape[0]] = np.diag(s) #S[:66,:66] = ...
    #       print('S-shape: ', S.shape)
    
    # For A = u*S*vh then U*S the basis spectra of the reference dataset 
    # ...and vh contains the least squares coefficients that refit the basis CD spectra to the experimental CD spectra

    #******************HJ5******************
    #Do the matrix multiplication, using the first 5 eigen vectors
    X = F @ np.transpose(vh)[:,:5] @ np.linalg.pinv(S[:5,:5]) @ np.transpose(u[:,:5])
    hj5 = np.matmul(X,q)
    
    #print('hj5 result')
    #SStruct = ['Alpha-r','Alpha-d', 'Beta-r', 'Beta-d', 'Turns', 'Unord']
    SStruct = ['AlphaR','AlphaD', 'BetaR', 'BetaD', 'Turns', 'Other']
    #SSDev = ['AR Dev', 'AD Dev', 'BR Dev', 'BD Dev', 'T Dev', 'Un Dev']
    SSDev = ['AR Dev', 'AD Dev', 'BR Dev', 'BD Dev', 'T Dev', 'Ot Dev']
    if 'IR' in CD_or_IR: #  if (CD_or_IR == 'IR' or CD_or_IR == 'CD-IR'):
        SStruct = ['Helix','Sheet', 'Turns', 'Rand C']
        SSDev = ['H Dev', 'S Dev', 'T Dev', 'RC Dev']
    format_list = [SStruct[hj5.tolist().index(item)]+'\t'+'{:.4f}' for item in hj5]
    strFormat = '\t'.join(format_list)
    print ('hj5 result\t\t  ', strFormat.format(*hj5), 'sum = ', '{:.4f}'.format(hj5.sum()))

    
    #******************Arrange A (and F) according to similarity to q******************
    #%This is a 1x71 (1x#ref proteins) matrix
    rmsds=np.zeros(A.shape[1]) #rmsds=zeros(1,size(A,2));
    #%calculate index of rmsd values of query to data values
    for i in range(A.shape[1]):  #for i= 1:size(A,2); 
        rmsds[i] = np.sqrt(((q - A[:,i]) ** 2).mean())  #rmsds(i)=rmsd(q,A(:,i)); %SVH 20-9-2006
    
    # Get the indexes of the rmsds sorted so that the smallest rmsds comes first
    ix = np.argsort(rmsds) #[sorted, ix] = sort(rmsds);
    Asort = np.zeros((A.shape[0],A.shape[1])) #Asort = zeros(size(A,1),size(A,2));
    Fsort = np.zeros((F.shape[0],F.shape[1])) #Fsort = zeros(size(F,1),size(F,2));
    #ProtsSort=prots;%  SVH 25/9-06
    
    # Arrange A and F according the the rmsd 
    for i in range(A.shape[1]): #for i= 1:size(A,2);
        Asort[:,i]=A[:,ix[i]] #Asort(:,i)=A(:,ix(i)); 
        Fsort[:,i]=F[:,ix[i]]  #Fsort(:,i)=F(:,ix(i)); 
        #ProtsSort(i)=prots(ix(i));
    #end
    
    #%concatanate q to A
    A=np.c_[q, Asort] #A=[q Asort;];
    #%make an intial guess about ss
    F=np.c_[Fsort[:,0], Fsort] #F=[Fsort(:,1) Fsort;];
    
    #pred.Asort=A;
    #pred.guess = Fsort(:,1);
    
    #******************SELCON1******************
    #%now do selcon1 
    
    #%for %min number of proteins to max  number of proteins
    #%for min number of basis spectra to max number of basis spectra
    #%calculate solution and store it in a matrix to be filtered by sum and
    #%fract rule
    
    selfcon =0
    attempts =0
    while selfcon==0:
        
        #%empty all of the solutions from the variable selection
        solutionsList=[]
        
        #%for various number of cd spectra
        for i in range(2, A.shape[1]): #this is 0 to (71+1)-1    #for i=3:size(A,2)
            #%make the truncated matrix B
            u, s, vh = np.linalg.svd(A[:,:i+1], full_matrices=True)  #[U,S,V] = svd(A(:,1:i),0); #print ('u-shape: ',u.shape,'s-shape: ', s.shape,'vh-shape: ', vh.shape)
            S = np.zeros((u.shape[0],vh.shape[0])) #66x71 #print('S-shape: ', S.shape)
            S[:s.shape[0],:s.shape[0]] = np.diag(s) #S[:66,:66] = ... Note :66 is 0 to 65
            #%calculate the maximum number of basis spectra to use
            maxbasis = 7
            if (i+1 < 7):
                maxbasis = i+1
            #%for various numbers of basis spectra
            for j in range(maxbasis): #this goes from 0 to maxbasis-1 #for j=1:maxbasis
                newSolution = F[:,:i+1] @ np.transpose(vh)[:,:j+1] @ np.linalg.pinv(S[:j+1,:j+1]) @ np.transpose(u[:,:j+1]) @ q #((F(:,1:i)*V(:,1:j) * pinv(S(1:j,1:j)) * U(:,1:j)') *q)
                solutionsList.append(newSolution) #solutions =[solutions ((F(:,1:i)*V(:,1:j) * pinv(S(1:j,1:j)) * U(:,1:j)') *q);];
            #end
        #end
        
        #convert to a matrix (#print('Solutions list')     #print(solutionsList))
        solutions = np.transpose(np.array(solutionsList))
        
        #%now apply the sum and fract rules to the solutions
        #%relaxing as we go to make sure we get a minimum of 1 solution
        
        #%keep track of how far weve iterated along all solutions
        solution_number =1
        
        #%track all valid solutions and parameters not useful for selcon but might
        #%give better results so try in the future
        valid_solutions = []
        #valid_params=[]
        
        #%store the parameters that have obtained the best solutions for each np
        np_solutions =[]
        np_params=[]
        #%init the sum and fract rule that may need to be relaxed later
        sum_rule = 0.05
        fract_rule = -0.025
        
        #%bool to see if weve found a minimum of one solution after filtering with
        #%the sum and fract rules
        found_one = 0
        while found_one == 0:
            solution_number = 1;    
            #%for each value of Np
            for i in range (2, A.shape[1]): #for i=3:size(A,2)
                #% a way of checking if we have the minimum sum error
                sum_best =10000
                #% store the best parameter of the basis spectra
                best_basis=0
                #% store the best parameters np and basis for the best solutions
                #best_param=[]
                #% store the best solution for this particular number of np
                best_solution=[]
                #% determine the maximum number of basis spectra for this np
                maxbasis =7
                if (i+1 <7):
                    maxbasis=i+1    
                #%for each number of basis spectra
                for j in range(maxbasis): #for j=1:maxbasis
                    #%if the solution satisfies the sum and fraction rule
                    testSum = solutions[:,solution_number-1].sum()
                    testMinFrac = solutions[:,solution_number-1].min()   #print('i,j', i,j,'testSum: ',testSum,' testMinFrac ',testMinFrac)
                    if ((abs(1-testSum) <= sum_rule) and (testMinFrac >=fract_rule)): #if (((abs(1-(sum(solutions(:,solution_number)))))<=sum_rule) && min(solutions(:,solution_number))>=fract_rule);
                        #print('testSum', testSum, 'abs(1-testSum)', abs(1-testSum),'testMinFrac', testMinFrac)
                        #%store the solutions and valid parameters
                        valid_solutions.append(solutions[:,solution_number-1])  ##valid_solutions=[valid_solutions solutions(:,solution_number);];
                        #valid_params.append(np.array([i, j])) #!!not sure this is used later!!  : valid_params.append([i;j;]) ##valid_params = [valid_params [i;j;];];    
                       
                        #%now keep track of the best solution and its basis number
                        #%for this particular np
                        if abs(1-testSum)<= sum_best: #if (abs(1-(sum(solutions(:,solution_number)))) < sum_best)
                            best_basis = j #best_basis =j;
                            best_solution = solutions[:,solution_number-1] #best_solution = solutions(:,solution_number);
                            sum_best = abs(1-testSum) #sum_best = abs(1-(sum(solutions(:,solution_number))));
                        #end
                    #end
                    #%increment the solution number
                    solution_number = solution_number +1;

                #end %for nj (num of basis spectra, SVH 22/9-2006)
                #%now go through and only include the closest solution to 1 from each value of Np. This is the best_solution
                #% if we have a valid solution for this np then record basis and solution
                if best_basis > 0: #if (best_basis >0)
                    np_params.append(np.array([i, best_basis])) #    np_params= [np_params [i;best_basis;];]; 
                    np_solutions.append(best_solution) #    np_solutions= [np_solutions best_solution;];
                #end
            
            #end%for np (#ref prot used)
            
            if len(valid_solutions) > 0: #if(size(valid_solutions,2)  > 0)
                found_one = 1        #print('used fract_rule', fract_rule, 'sum_rule= ', sum_rule)
            #end
            #%relax
            fract_rule= fract_rule - 0.005
            sum_rule= sum_rule + 0.01
        #end %while(found_one==0), SVH 22/9-2006

        

        #pred.np = np_solutions;
        #pred.np_params = np_params;
        
        
        #np_solutions is a list of solutions[:,solution_number-1] arrays so create this as a matrix
        np_solutionsMat = np.transpose(np.array(np_solutions))      #print('np_solutionsMat shape: ', np_solutionsMat.shape)
        sel1 = np.mean(np_solutionsMat,axis = 1) #pred.sel1 = mean(np_solutions,2);
        #pred.jon1 = mean(valid_solutions,2);
        
        #%pred.valid = valid_solutions;
        #%pred.valid_params = valid_params;
        
        #disp(['For selfconsist attempt ' num2str(attempts) ' rmsd is ' num2str(rmsd(F(:,1),mean(np_solutions,2)))]) %SVH 22-9-2006
        format_list = [SStruct[sel1.tolist().index(item)]+'\t'+'{:.4f}' for item in sel1]
        strFormat = '\t'.join(format_list)
        print ('For selfconsist attempt ', attempts, strFormat.format(*sel1), 'sum = ', '{:.4f}'.format(sel1.sum()))
        
        #%check to see if selfconsistent
        if np.sqrt(((F[:,0] - sel1) ** 2).mean()) <0.0025: #        if((rmsd(F(:,1),mean(np_solutions,2)))<0.0025) %If the new solution is close to the previous, we have selfconsistency, SVH 7/12-17
            selfcon=1;
        #end;
        
        #%can change this bit just to check if using all valid solutions is better
        F[:,0] = sel1 #mean(np_solutions,2); %Make a new guess for the solution, and run through the solution space again, SVH 7/12-17
        attempts =attempts +1;
        if (attempts >=50):
            print('bail*********************************************************')
            selfcon=1
        #end
        
        #selfcon = 1 #dummy end while loop
    #end  %while(selfcon==0)   %while not selfcon
    
    
    #******************SELCON2******************
    #%now we can do selcon 2
    #%filer out solutions using the refi errors of the successfuk solutuins
    sel2_solutions = []
    np_param_Sel2=[]  #; %Keep track on valid #Np and #Basis for valid Sel2 solutions, SVH 25/9-06
    #SVH: Note that np_params is a list of 1x2 arrays each array is np and basis
    refits = np.zeros(len(np_params))  #;refits = zeros(1,size(np_params,2));
    
    for i in range(len(np_params)): #for i=1:size(np_params,2)
        np_=np_params[i][0] #np=np_params(1,i);
        basis=np_params[i][1] #basis=np_params(2,i);
        #Note: np is in the range 0 to 71-1, and basis is in the range 0 to 7-1.
        #But when slicing e.g. :3 means 0 to 2, so to include the i and j use +1
        #So use them as indices +1 (old wrong comment: as they are (not +1))
        u, s, vh = np.linalg.svd(A[:,:np_+1], full_matrices=True) #[U,S,V] = svd(A(:,1:np),0);         #Vt=V';
        S = np.zeros((u.shape[0],vh.shape[0])) 
        S[:s.shape[0],:s.shape[0]] = np.diag(s)
        refit = u[:,:basis+1] @ S[:basis+1,:basis+1] @ vh[:basis+1,:] ##refit = U(:,1:basis)*S(1:basis,1:basis)*Vt(1:basis,:);        
        refits[i] = np.sqrt(((refit[:,0] - q) ** 2).mean()) #refits(i)=rmsd(refit(:,1),q); %SVH 22/9-2006        
    #end  

    
    got_refit=0
    #refit_rule is defined in top of the function
    #           refit_rule=0.25  #if (CD_or_IR == 'IR'): refit_rule=0.025
    while got_refit==0:
        for i in range(len(np_params)): #for i=1:size(np_params,2) %size(np_params,2) is the number of selcon1 solutions, SVH 7/12-2017
            if (refits[i] <= refit_rule):
                sel2_solutions.append(np_solutions[i]) #sel2_solutions=[sel2_solutions np_solutions(:,i);];
                np_param_Sel2.append(np_params[i]) #np_param_Sel2=[np_param_Sel2 np_params(:,i);]; %SVH 25/9-06
            #end
        #end; 
        if (len(sel2_solutions) > 0): #if (size(sel2_solutions,2) > 0)
            got_refit=1
        else:
            refit_rule = refit_rule + 0.01 #; %Relax refit rule
        #end
    #end %while(got_refit==0)
    sel2_solutionsMat = np.transpose(np.array(sel2_solutions))
    np_param_Sel2Mat = np.transpose(np.array(np_param_Sel2))
    
    sel2 = np.mean(sel2_solutionsMat,axis = 1)  #pred.sel2 = mean(sel2_solutions,2);
    
    format_list = [SStruct[sel2.tolist().index(item)]+'\t'+'{:.4f}' for item in sel2]
    strFormat = '\t'.join(format_list)
    print ('Sel2 solution:\t\t  ', strFormat.format(*sel2), 'sum = ', '{:.4f}'.format(sel2.sum()))
    print ('Sel2 solutions finally used RMSD rule:', refit_rule)
    
    #pred.np_param_Sel2=np_param_Sel2;
    #pred.sel2Sol = sel2_solutions; %SVH 7/12-2017
    
    #******************SELCON3******************
    #%and finally we carry out application of the helix rule (selcon3)
    sel3_solutions=[] #; 
    np_param_Sel3=[] #; %SVH 25/9-06
    
    #%************** change these two lines****
    #%not disordered helix
    if 'IR' in CD_or_IR:
        hjh = hj5[0] #hjh = hj5(1);
        hel = sel2_solutionsMat[0,:] #hel = sel2_solutions(1,:);
    else:
        #%these two lines need to be uncommented if ad (Aplha Disorted) and (Alpha Regular) assignment is used
        hjh = hj5[0] + hj5[1] #hjh = pred.hj5(1) + pred.hj5(2);
        hel = sel2_solutionsMat[1,:] + sel2_solutionsMat[0,:] #hel = sel2_solutions(2,:) + sel2_solutions(1,:);
    
    #%**************************************************************************
    
    hel_max = np.max(hel)  #hel_max = max(hel);
    hel_min = np.min(hel) #hel_min = min(hel);
    hel_ave = np.mean(hel) #hel_ave = mean(hel);
    
    if (hjh >0.65):
        for i in range(len(sel2_solutions)): #for i=1:size(sel2_solutions,2)
            helix = hel[i]
            if (helix > 0.65):
                sel3_solutions.append(sel2_solutions[i]) #sel3_solutions= [sel3_solutions sel2_solutions(:,i);];
                np_param_Sel3.append(np_param_Sel2[i]) #np_param_Sel3=[np_param_Sel3 np_param_Sel2(:,i);]; %SVH 25/9-06
            #end
        #end
    #end
    
    
    if (hjh <= 0.65 and hjh >= 0.25):
        for i in range(len(sel2_solutions)): #for i=1:size(sel2_solutions,2)
            helix = hel[i]
            if (helix <= (((hjh + hel_max)/2)+0.03) and helix >= (((hjh + hel_max)/2)-0.03)):
                sel3_solutions.append(sel2_solutions[i]) #sel3_solutions= [sel3_solutions sel2_solutions(:,i);];
                np_param_Sel3.append(np_param_Sel2[i]) #np_param_Sel3=[np_param_Sel3 np_param_Sel2(:,i);]; %SVH 25/9-06
            #end
        #end
    #end
    
    if (hjh < 0.25 and hjh >= 0.15):
        for i in range(len(sel2_solutions)): #for i=1:size(sel2_solutions,2)
            helix = hel[i]
            if (helix <= (((hjh + hel_ave)/2)+0.03) and helix >= (((hjh + hel_ave)/2)-0.03)):
                sel3_solutions.append(sel2_solutions[i]) #sel3_solutions= [sel3_solutions sel2_solutions(:,i);];
                np_param_Sel3.append(np_param_Sel2[i]) #np_param_Sel3=[np_param_Sel3 np_param_Sel2(:,i);]; %SVH 25/9-06
            #end
        #end
    #end
    
    if (hjh < 0.15):
        for i in range(len(sel2_solutions)): #for i=1:size(sel2_solutions,2)
            helix = hel[i]
            if (helix <= (((hjh + hel_min)/2)+0.03) and helix >= (((hjh + hel_min)/2)-0.03)):
                sel3_solutions.append(sel2_solutions[i]) #sel3_solutions= [sel3_solutions sel2_solutions(:,i);];
                np_param_Sel3.append(np_param_Sel2[i]) #np_param_Sel3=[np_param_Sel3 np_param_Sel2(:,i);]; %SVH 25/9-06
            #end
        #end
    #end
    
    #print('sel3_solutions ',sel3_solutions)
    #print('np_param_Sel3 ',np_param_Sel3)
    
    
    refit_np_param=[] #; %SVH 25/9-06
    #plusMinus6 = np.matrix(" +/- ; +/- ; +/- ; +/- ; +/- ; +/- ")  #[' +/- '; ' +/- '; ' +/- '; ' +/- '; ' +/- '; ' +/- ';]
    #plusMinus6List = [' +/- ', ' +/- ', ' +/- ', ' +/- ', ' +/- ', ' +/- ']  #plusMinus6=[' +/- '; ' +/- '; ' +/- '; ' +/- '; ' +/- '; ' +/- ';]
    #percentList = ['%', '%', '%', '%', '%', '%']# ['%'; '%'; '%'; '%'; '%'; '%';] #; %The ; makes it into a vertical vector
    #print(plusMinus6List, percentList)
    
    return_List=[]
    
    if (len(sel3_solutions) > 0): #if(size(sel3_solutions,2) > 0)
        sel3_solutionsMat = np.transpose(np.array(sel3_solutions))
        sel3 = np.mean(sel3_solutionsMat,axis = 1) #  pred.sel3 = mean(sel3_solutions,2);
        
        refit_np_param = np_param_Sel3 #; %SVH 25/9-06
        stDev = np.std(sel3_solutionsMat, axis=1) #std(sel3_solutions,1,2) #;%SVH 8/12/2017 flag=1: divide by n (flag=0 divide by n-1). 2= along second dimension
        #print(stDev)
        #pred.sel3rms=stDev;%SVH 8/12/2017
        #pred.sel3sol=sel3_solutions;
        writeStr = 'There are ' + str(len(sel3_solutions)) + ' Selcon3 solutions'
        return_List.append(writeStr)
        print(writeStr)
        #print('There are ', len(sel3_solutions), ' Selcon3 solutions') #disp(['There are ' num2str(size(sel3_solutions,2)) ' Selcon3 solutions']);
        writeStr = '\t Mean     Std'
        return_List.append(writeStr)
        print(writeStr)
        #print('\t Mean     Std') #disp ([' Mean'  '     Std']);
        formatSpec = '{:4.1f}'
        for i in range(len(sel3)): #disp ([num2str(100*mean(sel3_solutions,2),formatSpec) percent plusMinus6 num2str(100*stDev,formatSpec) percent;]);
            strFormat = SStruct[i]+'\t'+ formatSpec+' +/-' + formatSpec
            writeStr = strFormat.format(100*sel3[i], 100*stDev[i])
            return_List.append(writeStr)
            print(writeStr)
            #print(strFormat.format(100*sel3[i], 100*stDev[i])) 
        writeStr = 'sum =\t' + formatSpec.format(100*np.sum(sel3)) + '%'
        return_List.append(writeStr)
        print(writeStr)
        #print('sum =\t', formatSpec.format(100*np.sum(sel3)), '%')   #disp(['sum = ' num2str(sum(100*mean(sel3_solutions,2)),formatSpec) '%']);
        
        #Finally make the format for Nyk
        if (specialOutputNyk == 1):
            formatSpec = '{:4.1f}'
            strFormat = formatSpec + '\t'
            writeStrNyk = ""
            #All the secondary structure solutions
            for i in range(len(sel3)):
                writeStrNyk += strFormat.format(100*sel3[i])
            #The sum of solutions
            #writeStrNyk += strFormat.format(100*np.sum(sel2))
            writeStrNyk += strFormat.format(100*np.sum(sel3))
            #All the stDev of the solutions
            for i in range(len(sel3)):
                writeStrNyk += strFormat.format(100*stDev[i])
            #Which selcon
            writeStrNyk += 'Selcon3'+'\t'
            #number of Selcon3 solutions
            writeStrNyk += str(len(sel3_solutions))+'\t'
        
        
    else:
        #sel2_solutionsMat and np_param_Sel2Mat already exists, as well as sel2
        sel3 = np.mean(sel2_solutionsMat,axis = 1) # pred.sel3= mean(sel2_solutions,2);
        
        refit_np_param=np_param_Sel2 #; %SVH 25/9-06
        stDev = np.std(sel2_solutionsMat, axis=1) #stDev=std(sel2_solutions,1,2);%SVH 8/12/2017
        #pred.sel2rms=stDev;%SVH 8/12/2017
        #pred.sel2sol=sel2_solutions;
        #fprintf ( 2, 'No Selcon3 solutions!\n' )#;

        print('No Selcon3 solutions!\n', 'Using the Selcon2 solutions as result') #disp('Using the Selcon2 solutions as result');
        
        writeStr = 'There are ' + str(len(sel2_solutions)) + ' Selcon2 solutions'
        return_List.append(writeStr)
        print(writeStr)
        #print('There are ', len(sel2_solutions), ' Selcon2 solutions') #disp(['There are ' num2str(size(sel2_solutions,2)) ' Selcon2 solutions']);
        
        writeStr = '\t Mean     Std'
        return_List.append(writeStr)
        print(writeStr)
        #print('\t Mean     Std') #disp ([' Mean'  '     Std']);
        
        formatSpec = '{:4.1f}' #formatSpec = '%4.1f';
        for i in range(len(sel2)): #disp ([num2str(100*mean(sel2_solutions,2),formatSpec) percent plusMinus6 num2str(100*stDev,formatSpec) percent;]);
            strFormat = SStruct[i]+'\t'+ formatSpec+' +/-' + formatSpec
            writeStr = strFormat.format(100*sel2[i], 100*stDev[i])
            return_List.append(writeStr)
            print(writeStr)
            #print(strFormat.format(100*sel2[i], 100*stDev[i])) 
        writeStr = 'sum =\t' + formatSpec.format(100*np.sum(sel2)) + '%'
        return_List.append(writeStr)
        print(writeStr)
        #print('sum =\t', formatSpec.format(100*np.sum(sel2)), '%')   #disp(['sum = ' num2str(sum(100*mean(sel2_solutions,2)),formatSpec) '%']);
        
        #Finally make the format for Nyk
        if (specialOutputNyk == 1):
            formatSpec = '{:4.1f}'
            strFormat = formatSpec + '\t'
            writeStrNyk = ""
            #All the secondary structure solutions
            for i in range(len(sel2)):
                writeStrNyk += strFormat.format(100*sel2[i])
            #The sum of solutions
            writeStrNyk += strFormat.format(100*np.sum(sel2))
            #All the stDev of the solutions
            for i in range(len(sel2)):
                writeStrNyk += strFormat.format(100*stDev[i])
            #Which selcon
            writeStrNyk += 'Selcon2'+'\t'
            #number of Selcon2 solutions
            writeStrNyk += str(len(sel2_solutions))+'\t'
      #      writeStrNyk += UseRefSet+'\t' #added 10/1-2024
      #      writeStrNyk += MyProteinFile+'\t' #added 10/1-2024
        
    #end
    
    #%Now refit Selcon solutions, SVH 25/9-06
    
    
    #print('refit_np_param ', refit_np_param)
    refit_prot = np.zeros((A.shape[0],len(refit_np_param))) #number of datapoints (66) x number of np for refit
    for i in range(len(refit_np_param)): #for(i=1:size(refit_np_param,2))
        #Note the values stores in np_params are the i and j indices
        np_=refit_np_param[i][0] #np_=refit_np_param(1,i);
        basis=refit_np_param[i][1] #basis=refit_np_param(2,i);
        #Note: np is in the range 0 to 71-1, and basis is in the range 0 to 7-1. So use them as indices as they are (not +1)
        u, s, vh = np.linalg.svd(A[:,:np_+1], full_matrices=True) #[U,S,V] = svd(A(:,1:np),0);         #Vt=V';
        S = np.zeros((u.shape[0],vh.shape[0])) 
        S[:s.shape[0],:s.shape[0]] = np.diag(s)
        refit = u[:,:basis+1] @ S[:basis+1,:basis+1] @ vh[:basis+1,:] #        refit = U(:,1:basis)*S(1:basis,1:basis)*Vt(1:basis,:);
        refit_prot[:,i]=refit[:,0] #        refit_prot(:,i)=refit(:,1); %This is the refit prot spec for given #Np and #Basis           

    #end
    #print(refit_prot) 

    Mean_refit_prot = np.mean(refit_prot,axis=1) #mean(refit_prot,2); %Vector of mean refitted prot spectrum
    #pred.MeanRefitProt=Mean_refit_prot;
    #pred.RefitProt=refit_prot; %Array of all refited spectra for all good Sel3 #Np and #Basis
    
    #%Calc. RMSD value
    rmsdRefit = np.sqrt(((Mean_refit_prot - q) ** 2).mean()) #rmsdRefit=rmsd(Mean_refit_prot,q);
    writeStr = 'RMSD between Query and Refitted spectra ' + '{:.4f}'.format(rmsdRefit)
    return_List.append(writeStr)
    print(writeStr)
    #print('RMSD between Query and Refitted spectra ', '{:.4f}'.format(rmsdRefit)) #disp(['RMSD between Query and Refitted spectra ' num2str(rmsdRefit)]);
    
    #Add the RMSD for Nyk
    if (specialOutputNyk == 1):
        writeStrNyk += '{:.4f}'.format(rmsdRefit)+'\t'
        writeStrNyk += UseRefSet+'\t' #added 10/1-2024
        writeStrNyk += MyProteinFile+'\t' #added 10/1-2024
        #The long data string for Nyk is the last added to the return list.
        return_List.append(' ')
        return_List.append('Output for Nyk:')
        listStr = ''
        for item in SStruct:
            listStr += item + '\t'
        listStr += 'Sum'+'\t'
        for item in SSDev:
            listStr += item + '\t'
        listStr += 'selcon?'+'\t'+'#solu'+'\t'+'RMSD'
        return_List.append(listStr)
        return_List.append(writeStrNyk)
   
    np.savetxt(MyProteinDir + MyProteinFile[0:-4]+'_SELCON3-'+ UseRefSet +'_refit.txt', np.column_stack((wl_k, q, Mean_refit_prot)), fmt='%d \t %0.4f \t %0.4f', header='WL/k \t Spectrum \t Mean refit', comments='') #https://numpy.org/doc/stable/reference/generated/numpy.savetxt.html
    
    if (ShouldIPlot==1):
        import matplotlib.pyplot as plt
        import matplotlib as mpl   
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
        figWSpacer = 5
            
        mpl.rcParams['font.size'] = 16
        mpl.rcParams['legend.fontsize'] = 12
        NpMax = 0
        
        #Define WL range for plot
        # wl=np.arange(240, LimWL-1, -1) #240:-1:LimWL; %Vector of wavelengths
        # if (CD_or_IR == 'IR'):
        #     wl=np.arange(1600, LimWL+1, +1) #from 1600 to 1800
        #USE: wl_k calculated in start of function
        
        #%Define lots of things like linewidths, markers.. and add labels and titles
        MyyLabel = r'$\Delta \epsilon$ [M$^{-1}$ cm$^{-1}$ residue$^{-1}$]'
        MyxLabel = 'Wavelength [nm]'
        if (CD_or_IR == 'IR'):
            MyxLabel = 'Wavenumbers [cm-1]'
            MyyLabel = 'IR absorbance'
        if 'CD-IR' in CD_or_IR: #if (CD_or_IR == 'CD-IR'):
            MyxLabelIR = 'Wavenumbers [cm-1]'
            MyyLabelIR = 'IR absorbance'
        
        #Setting that change between CD/IR and CD-IR plots
        wlCD = wl_k
        qCD = q
        refit_protCD = refit_prot
        Mean_refit_protCD = Mean_refit_prot
        numCols = 1
        widthMult = 1.0
        if 'CD-IR' in CD_or_IR: #if (CD_or_IR == 'CD-IR'):
            #RefDataLowWL = 175 in the CD-IR case
            #minWLFloat = np.min(wl) #even for mixed CD and IR data, this is the lowest CD data wl
            #Already defined from function call: LimWL = max(int(minWLFloat), RefDataLowWL) #minWL
            CD_NumDataPoints = int(1+(RefDataHighWL-LimWL))
            #CD data
            wlCD = wl_k[:CD_NumDataPoints]; qCD = q[:CD_NumDataPoints]
            refit_protCD = refit_prot[:CD_NumDataPoints,:]; Mean_refit_protCD = Mean_refit_prot[:CD_NumDataPoints]
            #IR data
            wlIR = wl_k[CD_NumDataPoints:]; qIR = q[CD_NumDataPoints:]
            refit_protIR = refit_prot[CD_NumDataPoints:,:]
            Mean_refit_protIR = Mean_refit_prot[CD_NumDataPoints:]
            
            numCols = 2; widthMult = 1.5
        
        #PLOT 1:  %PLOT query prot and all refits
        fig1 = plt.figure(figsize=[widthMult*figWidth*px,figHeight*px])
        fig1.canvas.manager.window.move(figStartX,figStartY)
        gs = fig1.add_gridspec(1, numCols, wspace=0) #Create subplot 1 row and 2 cols with no spacing
        axs = gs.subplots(sharey= 'row') #sharey = True)
        if 'CD-IR' in CD_or_IR: #if (CD_or_IR == 'CD-IR'):
            axs0=axs[0]
        else:
            axs0=axs
       
        axs0.plot(wlCD, qCD, 'bo', markersize  = 5.0, mfc='none', label = MyProteinFile)
        for i in range(refit_prot.shape[1]):
            #plt.plot(wl, refit_prot[:,i], '-', markersize  = 2.0, label='Refir '+str(i)) #tempStr=['Refitted Np=' num2str(refit_np_param(1,i-1)) ', #Basis=' num2str(refit_np_param(2,i-1))];
            axs0.plot(wlCD, refit_protCD[:,i], '-', linewidth = 1.0, label='Refitted Np='+str(refit_np_param[i][0]+1)+' #Basis='+str(refit_np_param[i][1]+1)) #tempStr=['Refitted Np=' num2str(refit_np_param(1,i-1)) ', #Basis=' num2str(refit_np_param(2,i-1))];
            if (refit_np_param[i][0]+1 > NpMax):
                NpMax = refit_np_param[i][0]+1
        
        axs0.legend()
        axs0.set(ylabel =MyyLabel, xlabel = MyxLabel)
        axs0.axhline(linewidth=1, color='grey', dashes=[2,2])
        fig1.suptitle('Spectra of Query protein and all valid selcon3 refits', y = 0.94)
        
        if 'CD-IR' in CD_or_IR: #if (CD_or_IR == 'CD-IR'):
            axs[1].plot(wlIR, qIR, 'bo', markersize  = 5.0, mfc='none', label = MyProteinFile)
            for i in range(refit_prot.shape[1]):
                axs[1].plot(wlIR, refit_protIR[:,i], '-', linewidth = 1.0, label='Refitted Np='+str(refit_np_param[i][0]+1)+' #Basis='+str(refit_np_param[i][1]+1)) #tempStr=['Refitted Np=' num2str(refit_np_param(1,i-1)) ', #Basis=' num2str(refit_np_param(2,i-1))];

            axs[1].set(ylabel =MyyLabelIR, xlabel = MyxLabelIR)
            axs[1].yaxis.set_label_position("right")
            axs[0].spines['right'].set_visible(False)
            axs[1].spines['left'].set_visible(False)
            axs[1].yaxis.tick_right()
            axs[1].axhline(linewidth=1, color='grey', dashes=[2,2])
            axs[1].set_ylim(axs[0].get_ylim())
        #plt.tigth_layout() #Doesnøt work
        #gs.tigth_layout()
        #plt.rcParams['axes.titley'] = 1.0    # y is in axes-relative coordinates.
        #plt.rcParams['axes.titlepad'] = -14  # pad is in points..
        fig1.show()
        
        
        # fig1 = plt.figure(1, figsize=[figWidth*px,figHeight*px])#, frameon = False) #figure(1); %Define figure
        # fig1.canvas.manager.window.move(figStartX,figStartY)
        # plt.axhline(linewidth=1, color='grey', dashes=[2,2])
        # #plt.plot(wl, q, 'bo', markersize  = 5.0, mfc='none', label = "Protein")
        # plt.plot(wl_k, q, 'bo', markersize  = 5.0, mfc='none', label = MyProteinFile)
        # for i in range(refit_prot.shape[1]):
        #     #plt.plot(wl, refit_prot[:,i], '-', markersize  = 2.0, label='Refir '+str(i)) #tempStr=['Refitted Np=' num2str(refit_np_param(1,i-1)) ', #Basis=' num2str(refit_np_param(2,i-1))];
        #     plt.plot(wl_k, refit_prot[:,i], '-', linewidth = 1.0, label='Refitted Np='+str(refit_np_param[i][0]+1)+' #Basis='+str(refit_np_param[i][1]+1)) #tempStr=['Refitted Np=' num2str(refit_np_param(1,i-1)) ', #Basis=' num2str(refit_np_param(2,i-1))];
        #     if (refit_np_param[i][0]+1 > NpMax):
        #         NpMax = refit_np_param[i][0]+1
        # plt.legend()

        # plt.ylabel(MyyLabel) #plt.ylabel(r'$\Delta \epsilon$ [M$^{-1}$ cm$^{-1}$ residue$^{-1}$]') #        ylabel('\Delta \epsilon [mdeg /M /cm /residue]');
        # plt.xlabel(MyxLabel) #plt.xlabel('Wavelength [nm]') #        xlabel('Wavelength [nm]');
        # plt.title('Spectra of Query protein and all valid selcon3 refits')
        # fig1.show()
        
        #PLOT 2: %PLOT query prot and MEAN of refits
        fig2 = plt.figure(figsize=[widthMult*figWidth*px,figHeight*px])
        fig2.canvas.manager.window.move(figStartX + int(widthMult*figWidth) + figHSpacer, figStartY)
        gs = fig2.add_gridspec(1, numCols, wspace=0) #Create subplot 1 row and 2 cols with no spacing
        axs = gs.subplots(sharey='row') #sharey = True)
        if 'CD-IR' in CD_or_IR: #if (CD_or_IR == 'CD-IR'):
            axs0=axs[0]
        else:
            axs0=axs
       
        axs0.plot(wlCD, qCD, 'bo', markersize  = 5.0, mfc='none', label = MyProteinFile)
        axs0.plot(wlCD, Mean_refit_protCD, '-', linewidth = 1.0, label="Mean refit")
        axs0.legend()
        axs0.set(ylabel =MyyLabel, xlabel = MyxLabel)
        axs0.axhline(linewidth=1, color='grey', dashes=[2,2])
        fig2.suptitle('Spectra of Query protein and mean of all valid selcon3 refits', y = 0.94)
        
        if 'CD-IR' in CD_or_IR: #if (CD_or_IR == 'CD-IR'):
            axs[1].plot(wlIR, qIR, 'bo', markersize  = 5.0, mfc='none', label = MyProteinFile)
            axs[1].plot(wlIR, Mean_refit_protIR, '-', linewidth = 1.0, label="Mean refit")
            axs[1].set(ylabel =MyyLabelIR, xlabel = MyxLabelIR)
            axs[1].yaxis.set_label_position("right")
            axs[0].spines['right'].set_visible(False)
            axs[1].spines['left'].set_visible(False)
            axs[1].yaxis.tick_right()
            axs[1].axhline(linewidth=1, color='grey', dashes=[2,2])
            axs[1].set_ylim(axs[0].get_ylim())
        fig2.show()
        
        # fig2 = plt.figure(2, figsize=[figWidth*px,figHeight*px])#, frameon = False)   #figure(1); %Define figure
        # fig2.canvas.manager.window.move(figStartX + figWidth + figHSpacer, figStartY)
        # plt.axhline(linewidth=1, color='grey', dashes=[2,2])
        # plt.plot(wl_k, q, 'bo', markersize  = 5.0, mfc='none', label = MyProteinFile)
        # plt.plot(wl_k, Mean_refit_prot, '-', linewidth = 1.0, label="Mean refit")
        # plt.legend()
        
        # #%Define lots of things like linewidths, markers.. and add labels and titles        
        # plt.ylabel(MyyLabel) #plt.ylabel(r'$\Delta \epsilon$ [M$^{-1}$ cm$^{-1}$ residue$^{-1}$]') #        ylabel('\Delta \epsilon [mdeg /M /cm /residue]');
        # plt.xlabel(MyxLabel) #plt.xlabel('Wavelength [nm]') #        xlabel('Wavelength [nm]');
        # plt.title('Spectra of Query protein and mean of all valid selcon3 refits')
        
        #fig2.show()

        plt.savefig(MyProteinDir + MyProteinFile[0:-4]+'_SELCON3-'+ UseRefSet +'_fig.png')
  
        #ShouldIAlsoPlotSelcon2 is not implemented for combined CD-IR
        if(ShouldIAlsoPlotSelcon2 == 1 and 'CD-IR' not in CD_or_IR): #CD_or_IR != 'CD-IR' and CD_or_IR != 'CD-IR2m'):
            #%Calc. selcon2 solutions
            refit_prot_Sel2 = np.zeros((A.shape[0],len(np_param_Sel2))) #number of datapoints (66) x number of np for refit
            for i in range(len(np_param_Sel2)): #for(i=1:size(refit_np_param,2))
                np_=np_param_Sel2[i][0] #np=np_param_Sel2(1,i);
                basis=np_param_Sel2[i][1] #basis=np_param_Sel2(2,i);
                #Note: np is in the range 0 to 71-1, and basis is in the range 0 to 7-1. So use them as indices as they are (not +1)
                u, s, vh = np.linalg.svd(A[:,:np_+1], full_matrices=True) #[U,S,V] = svd(A(:,1:np),0);         #Vt=V';
                S = np.zeros((u.shape[0],vh.shape[0])) 
                S[:s.shape[0],:s.shape[0]] = np.diag(s)
                refit = u[:,:basis+1] @ S[:basis+1,:basis+1] @ vh[:basis+1,:] #        refit = U(:,1:basis)*S(1:basis,1:basis)*Vt(1:basis,:);
                refit_prot_Sel2[:,i]=refit[:,0] #        refit_prot(:,i)=refit(:,1); %This is the refit prot spec for given #Np and #Basis  
            
            #PLOT 3: %Plot Selcon2 solutions
            fig3 = plt.figure(3, figsize=[figWidth*px,figHeight*px])#, frameon = False) #figure(1); %Define figure
            fig3.canvas.manager.window.move(figStartX, figStartY + figHeight + figWSpacer)
            #plt.plot(wl, q, 'bo', markersize  = 5.0, mfc='none', label = "Protein")
            plt.plot(wl_k, q, 'bo', markersize  = 5.0, mfc='none', label = MyProteinFile)
            for i in range(refit_prot_Sel2.shape[1]):
                #labelTxt='Refitted Np='+str(np_param_Sel2[i][0]+1)+' #Basis='+str(np_param_Sel2[i][1]+1)
                labelTxt = ""
                plt.plot(wl_k, refit_prot_Sel2[:,i], '-', linewidth = 1.0, label=labelTxt) #plt.plot(wl, refit_prot_Sel2[:,i], '-', markersize  = 2.0, label='refit '+str(i))
            plt.legend()
            
            #%Define lots of things like linewidths, markers.. and add labels and titles  
            plt.ylabel(MyyLabel) #plt.ylabel(r'$\Delta \epsilon$ [M$^{-1}$ cm$^{-1}$ residue$^{-1}$]') #        ylabel('\Delta \epsilon [mdeg /M /cm /residue]');
            plt.xlabel(MyxLabel) #plt.xlabel('Wavelength [nm]') #        xlabel('Wavelength [nm]');
            plt.title('Spectra of Query protein and all valid selcon2 refits')
            fig3.show()
        

    #end %if(ShouldIPlot==1) %So plot solutions
    
    if (includeProteinsUsedForCalc>=1):
        #Which proteins where included in the solutions
        writeStr ='\n'
        writeStr +='****** Which proteins used for the calculation? ******'+'\n'
        writeStr += 'Max Np ' + '{:d}'.format(NpMax)+'\n'
        
        #Default (CD_or_IR)
        try:
            #Lbl =  np.loadtxt("CDRefData/Labels-SMP180_Andy-A71Order.txt", dtype='str', delimiter='\t')
            if (AU_or_BBK == 'AU'):
                Lbl =  np.loadtxt("CDRefDataAU/Labels-SMP180_PCDDBOrder.txt", dtype='str', delimiter='\t')  
            else:
                #For BBK ref set
                Lbl =  np.loadtxt("CDRefDataBBK/Labels-SMP180_Andy-A71Order.txt", dtype='str', delimiter='\t')
                #    Lbl =  np.loadtxt("Labels-SMP180_Andy-A71Order.txt", dtype='str', delimiter='\t')  
            if (CD_or_IR == 'IR'):
                Lbl =  np.loadtxt("IRRefData/Lbl_A50-AR-EG-IR_Sept2023.txt", dtype='str', delimiter='\t')  
            if 'CD-IR' in CD_or_IR:
                Lbl =  np.loadtxt("CD-IRRefDataAU/Lbl_A28-CD-IR_Aug2024.txt", dtype='str', delimiter='\t') 
            writeStr += 'Proteins used in solution in RMSD order:\n'
            for i in range(min(NpMax,len(ix))):
                writeStr += '{:d}: '.format(ix[i])
                writeStr += Lbl[ix[i]]+'\n'
        except:
            writeStr +="Could not find the Label file"
            print(writeStr)
        #END: Try/Except block

    
    if (includeProteinsUsedForCalc>=1):
        return_List.append(writeStr)
        if (includeProteinsUsedForCalc==2):
            print(writeStr)
        else:
            print('List of Proteins used for the calculation is included in the OUT file')    

   
    return return_List
        