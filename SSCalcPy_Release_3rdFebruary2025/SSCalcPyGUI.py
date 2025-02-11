# -*- coding: utf-8 -*-
"""
Created on Thu Sep 21 16:54:58 2023

@author: au513
"""
#Selcon3GUI.py

#Global variables
# figStartX = 100
# figStartY = 50
# figWidth = 900
# figHeight = 600
# figHSpacer = 5
# figWSpacer = 5

import sys
sys.path.append("Selcon3Py")
sys.path.append("CDSSTR")
sys.path.append("C:\ProgramData\Anaconda3")


from PyQt5 import QtWidgets, uic #, QtGui
#from PyQt5 import QtCore
from Selcon3LoadAndRun import LoadAndRun
from CDSSTRLoadAndRun import LoadAndRunCDSSTR
import numpy as np

import matplotlib.pyplot as plt
import matplotlib as mpl 
mpl.use('qtagg') #matplotlib.use('qtagg')

#gui file from Anaconda Promt as designer.exe
#see also https://realpython.com/python-pyqt-gui-calculator/
#qtcreator_file  = "Selcon3GUI.ui" # "tax_calc.ui" # Enter file here.
qtcreator_file  = "SSCalcPyGUI.ui"
Ui_MainWindow, QtBaseClass = uic.loadUiType(qtcreator_file)


class MyApp(QtWidgets.QMainWindow, Ui_MainWindow):
    
    def __init__(self): #, CD_OR_IR):
        #self.CD_OR_IR = CD_OR_IR
        
        QtWidgets.QMainWindow.__init__(self)
        Ui_MainWindow.__init__(self)
        
        self.setupUi(self)

        figStartX = 25
        figStartY = 50
        self.move(figStartX,figStartY)
        
        #Get the window to be in front when launching the GUI for the first time
        #all the activatewindow and setFocus doesn't work. Use the hack to minimize and showNormal
        #self.activateWindow()
        #self.setFocus()
        #QtWidgets.QMainWindow.activateWindow
        #self.Ui_MainWindow.setFocus(Qt.NoFocusReason)    ## self.setFocus(QtWidgets.NoFocusReason) #  Qt.NoFocusReason)
        self.showMinimized()
        self.showNormal()
        
        
        self.labelFileName.setText("")
        
        self.pushButtonStartCalc.clicked.connect(self.startSSCalc)
        self.pushButtonSelectFile.clicked.connect(self.selectFile)
        
        #The BBK reference set is unfortunately not an option for redistributed versions
        self.radioButtonRefSetBBK.setHidden(True)
        #self.radioButtonRefSetBBK.setHidden(False) #NCJ Dec. 2023 change to be able to run with BBK ref sets.
    
    def selectFile(self):
        fileName = ""
        CurrentPath = self.labelMyProteinDir.text()
        if (CurrentPath == ""):
            CurrentPath = "./."
        fileTypes = "Ascii files (*.txt *.dat)"
        if self.radioButtonCSVXY.isChecked():
            fileTypes = "CSV XY files (*.csv)"
        
        #Clear the filenames
        self.labelFileName.setText(""); self.labelMyProteinFile.setText(""); self.labelMyProteinDir.setText("")
        #Open dialog file
        dialog = QtWidgets.QFileDialog()
        fileName, _ = dialog.getOpenFileName(self, "Select Protein Spectrum File", CurrentPath, fileTypes)
        self.labelFileName.setText(fileName)
        
        
        #Things that need to be set
        MyProteinFile = fileName.split('/')[-1]
        self.labelMyProteinFile.setText(MyProteinFile)
        MyProteinDir = fileName[:(len(fileName) - len(MyProteinFile))]
        print('MyProteinDir: ' + MyProteinDir)
        self.labelMyProteinDir.setText(MyProteinDir)
        if (MyProteinFile != ""):
            RefDataHighWL = 240
            RefDataLowWL = 175
            LimWL = 175
            if self.radioButtonRefSetSMP180.isChecked():
                LimWL = 180
                RefDataLowWL = 180
            if self.radioButtonIR.isChecked():
                LimWL = 1800
                RefDataLowWL = 1800
            if self.radioButtonCDIR.isChecked():
                LimWL = 175
                #How to deal with IR ref set. Hard code to 1800 cm-1???
            
            
            isOK = self.plotData(MyProteinDir, MyProteinFile, LimWL, RefDataHighWL, RefDataLowWL)
            if isOK:
                info = "Protein spectrum is selected \n ... and plotted"
                self.labelMessage.setText(info)
                self.labelMessage.show()
            else:
                pass  #Do not show, error in plotData
                self.labelFileName.setText('')
                self.labelMyProteinDir.setText('')
                self.labelMyProteinFile.setText('')
            
    def loadData(self, MyProteinDir, MyProteinFile):
        print("loading data...")
        #These are the CD reference set defaults. Set otherwise if IR
        RefDataHighWL = 240
        RefDataLowWL = 175
        LimWL = 175
        if self.radioButtonRefSetSMP180.isChecked():
            LimWL = 180
            RefDataLowWL = 180
        if self.radioButtonCDIR.isChecked():
            LimWL = 175
            RefDataLowWL = 175
        
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
        
        #load in the delta Epsilon data of the protein you want to analyze
        delimIs ='\t' # radioButtonASCIIXY / radioButtonCSVXY
        if self.radioButtonCSVXY.isChecked():
            delimIs = ','

        q1 = np.loadtxt(MyProteinDir + MyProteinFile, dtype='f', delimiter=delimIs, skiprows=skip, max_rows=readRows)    #q1 = np.loadtxt(MyProteinDir + MyProteinFile, dtype='f', delimiter='\t', skiprows=skip)
        
        if self.radioButtonCD.isChecked():
            #Find out which format the protein CD spectrum data has (single or double rows)
            try:
                #If the file loaded intp q1 is just a single row, then if (q1.shape[1]>0): fail and an exception is thrown
                #Otherwise assume the file is two rows with WL in first row and Delta Epsilon in second row
                if (q1.shape[1]>0):
                    #Interpolate, first find lowest WL datapoint
                    minWLFloat = np.min(q1[:,0])
                    LimWL = max(int(minWLFloat), RefDataLowWL) #minWL
                    print('LimWL is determined to be=',LimWL)
                    wl_increasing = np.linspace(LimWL, RefDataHighWL, int(1+(RefDataHighWL-LimWL)))
                    #print(wl_increasing)
                    if (q1[-1,0] < q1[0,0]):
                        q_interp = np.interp(wl_increasing, np.flip(q1[:, 0]), np.flip(q1[:, 1]))
                    else:
                        q_interp = np.interp(wl_increasing, q1[:, 0], q1[:, 1])
                    #reverse q order to fit order of CD reference data
                    q = np.flip(q_interp)
            except:
                #Assume that the file is just a single row of delta Epsilon starting at 240 nm
                q=q1
            wl = np.arange(RefDataHighWL, LimWL-1, -1) #240:-1:LimWL; %Vector of wavelengths
        #END: if self.radioButtonCD.isChecked():

        if self.radioButtonIR.isChecked():
            IR_LowKInRefSet = 1600; IR_HighKInRefSet = 1800; IR_RefSetSteps = 1
            #if self.checkBoxUseCDIR2cm.isChecked():
            #    IR_LowKInRefSet = 1600; IR_HighKInRefSet = 1720; IR_RefSetSteps = 2
            IR_RefSetNumDataPoints = 1 + ((IR_HighKInRefSet-IR_LowKInRefSet)/IR_RefSetSteps); print('IR_RefSetNumDataPoints: ', IR_RefSetNumDataPoints)
            
            LimWL = 1800 #is this right
            RefDataLowWL = 1800 #hack to make it fit 240-174 nm range thinking in the code
             
            #Extract the k and the IR values from the spectrum
            k_values = q1[:, 0]; IR_values =q1[:, 1]
            #.... and interpolate
            k_InRefSet = np.linspace(IR_LowKInRefSet, IR_HighKInRefSet, int(IR_RefSetNumDataPoints))
            q_interp = np.interp(k_InRefSet, k_values, IR_values) #Note: k_InRefSet values out side range of k_values are given as IR_values[0] and IR_values[-1]
            q = q_interp / max(q_interp)
            wl=np.arange(1600, LimWL+1, +1) #from 1600 to 1800
        #END: if self.radioButtonIR.isChecked():
        
        if self.radioButtonCDIR.isChecked():
            #First read the CD data
            #Interpolate, first find lowest WL datapoint
            minWLFloat = np.min(q1[:,0]) #even for mixed CD and IR data, this is the lowest CD data wl
            LimWL = max(int(minWLFloat), RefDataLowWL) #minWL
            print('LimWL is determined to be=',LimWL)
            CD_NumDataPoints = int(1+(RefDataHighWL-LimWL))
            #if reading a pure IR dataset (mistake) CD_NumDataPoints is most likely negative
            if (CD_NumDataPoints <0):
                print("Error in defining CD part of CD-IR spectrum")
                #Return dummy wl, q
                return [0],[0]
            wl_increasingCD = np.linspace(LimWL, RefDataHighWL, CD_NumDataPoints) #wl_increasingCD = np.linspace(LimWL, RefDataHighWL, int(1+(RefDataHighWL-LimWL)))
            #q1CD = q1[0:CD_NumDataPoints,:] #Doesn't work if maxWL>RefDataHighWL. extract the CD data
            
            #Next, find the place in the q1[:,0] where there is a big jump in WL (where IR data starts)
            lastCDwl = q1[-1,0] #temp set to last element (although this is most likely an IR k value)
            for count, qvalue in enumerate(q1[1:,0]):
                #note: count is 1 less than the index for the qvalue, i.e. starts at zero. So q1[count,:]=q1[0,:] for the first cycle in the for loop
                if (abs(q1[count,0]-qvalue)>1):
                    lastCDwl = q1[count,0]
                    break
            CD_NumDataPointsInq1 = int(1+abs(q1[0,0]-lastCDwl))
            q1CD = q1[0:CD_NumDataPointsInq1,:] #extract the CD data
                
            print("CD part of data: First WL, Last WL ", q1CD[0,0],q1CD[-1,0])
            if (q1CD[1,0] < q1CD[0,0]): #Is the second wl smaller than the first, then in decending order. if (q1CD[-1,0] < q1CD[0,0]):
                q_interp = np.interp(wl_increasingCD, np.flip(q1CD[:, 0]), np.flip(q1CD[:, 1]))
            else:
                q_interp = np.interp(wl_increasingCD, q1CD[:, 0], q1CD[:, 1])
            #reverse q order to fit order of CD reference data
            qCD = np.flip(q_interp)
            wlCD = np.arange(RefDataHighWL, LimWL-1, -1) #240:-1:LimWL; %Vector of wavelengths
            
            #Next extract the IR data
            IR_LowKInRefSet = 1600; IR_HighKInRefSet = 1800; IR_RefSetSteps = 1
            if self.checkBoxUseCDIR2cm.isChecked():
                IR_LowKInRefSet = 1600; IR_HighKInRefSet = 1720; IR_RefSetSteps = 2
            IR_RefSetNumDataPoints = 1 + ((IR_HighKInRefSet-IR_LowKInRefSet)/IR_RefSetSteps); print('IR_RefSetNumDataPoints: ', IR_RefSetNumDataPoints)
            #Extract the k and the IR values from the spectrum
            #k_values = q1[CD_NumDataPoints:, 0]; IR_values =q1[CD_NumDataPoints:, 1]
            k_values = q1[CD_NumDataPointsInq1:, 0]; IR_values =q1[CD_NumDataPointsInq1:, 1]
            print("IR part of data: First k, Last k ", k_values[0],k_values[-1])
            #.... and interpolate
            k_InRefSet = np.linspace(IR_LowKInRefSet, IR_HighKInRefSet, int(IR_RefSetNumDataPoints))
            q_interp = np.interp(k_InRefSet, k_values, IR_values) #Note: k_InRefSet values out side range of k_values are given as IR_values[0] and IR_values[-1]
            
            
            #Get the scale factor for IR in combined CD and IR
            scaleIRText = self.lineEditIRAbsMax.text()
            self.lineEditIRAbsMax.setText(scaleIRText.replace(",", "."))
            scaleIR = float(self.lineEditIRAbsMax.text())
            #End: get the scaling factor for IR in combined IR and CD
            qIR = q_interp * (scaleIR / max(q_interp))
            
            wlIR=np.arange(IR_LowKInRefSet, IR_HighKInRefSet+1, IR_RefSetSteps) #+1) #from 1600 to 1800 or 1600 to 1720 in 2cm-1 steps
            
            #Put the IR and CD data together
            #wl= wlCD q = qCD wl= wlIR q = qIR
            wl = np.concatenate([wlCD,wlIR])
            q = np.concatenate([qCD,qIR])
        #END: if self.radioButtonCDIR.isChecked():
            
        print("RefDataHighWL, RefDataLowWL, LimWL: ", RefDataHighWL, RefDataLowWL, LimWL)
        return wl, q
    
    def plotData(self, MyProteinDir, MyProteinFile, LimWL, RefDataHighWL, RefDataLowWL):
        
        wl, q = self.loadData(MyProteinDir, MyProteinFile)
        #print("shapes of wl and q", len(wl), ", ", len(q), ". Max/Min of q: ", np.max(q),"/",np.min(q))
        if (len(wl) != len(q) or abs(np.max(q) - np.min(q)) < 0.0001):
            Message = "ERROR: The protein spectrum couldn't be plotted !"
            Message += "\n" + "Did you choose the correct type of file (CD/IR) for the spectrum you selected?"
            plt.close('all')
            self.labelMessage.setText(Message)
            self.labelMessage.show()
            return False #None #exit functions
        
        
        
        #%Plot solutions
        plt.close('all')
        #mpl.rcParams['toolbar'] = 'None'
        px = 1/mpl.rcParams['figure.dpi']  # pixel in inches
        figStartX = 25
        figStartY = 50
        figWidth = 900
        figHeight = 600
        #figHSpacer = 5 #25
        figWSpacer = 5 #50
        mpl.rcParams['font.size'] = 16
        mpl.rcParams['legend.fontsize'] = 12
        #mpl.rcParams['figure.facecolor'] = 'yellow'
        
        #%Define lots of things like linewidths, markers.. and add labels and titles
        MyyLabel = r'$\Delta \epsilon$ [M$^{-1}$ cm$^{-1}$ residue$^{-1}$]'
        MyxLabel = 'Wavelength [nm]'
        #if (CD_or_IR == 'IR'):
        if self.radioButtonIR.isChecked():
            MyxLabel = 'Wavenumbers [cm-1]'
            MyyLabel = 'IR absorbance'
        if self.radioButtonCDIR.isChecked():
            MyxLabelIR = 'Wavenumbers [cm-1]'
            MyyLabelIR = 'IR absorbance'
            
        #General plot, both for CD, IR alone and combined CD-IR
        wlCD = wl
        qCD = q
        numCols = 1
        widthMult = 1.0
        if self.radioButtonCDIR.isChecked():
            #RefDataLowWL = 175 in the CD-IR case
            minWLFloat = np.min(wl) #even for mixed CD and IR data, this is the lowest CD data wl
            LimWL = max(int(minWLFloat), RefDataLowWL) #minWL
            CD_NumDataPoints = int(1+(RefDataHighWL-LimWL))
            wlCD = wl[:CD_NumDataPoints]; qCD = q[:CD_NumDataPoints]
            wlIR = wl[CD_NumDataPoints:]; qIR = q[CD_NumDataPoints:]
            numCols = 2; widthMult = 1.5
        
        fig6 = plt.figure(figsize=[widthMult*figWidth*px,figHeight*px])
        fig6.canvas.manager.window.move(figStartX, figStartY + figHeight + figWSpacer)
        gs = fig6.add_gridspec(1, numCols, wspace=0) #Create subplot 1 row and 2 cols with no spacing
        axs = gs.subplots(sharey='row') #sharey = True)
        if self.radioButtonCDIR.isChecked():
            axs0=axs[0]
        else:
            axs0=axs
       
        axs0.plot(wlCD, qCD, 'bo', markersize  = 5.0, mfc='none', label = MyProteinFile)
        axs0.legend()
        axs0.set(ylabel =MyyLabel, xlabel = MyxLabel)
        axs0.axhline(linewidth=1, color='grey', dashes=[2,2])
        
        if self.radioButtonCDIR.isChecked():
            axs[1].plot(wlIR, qIR, 'bo', markersize  = 5.0, mfc='none', label = MyProteinFile)
            axs[1].set(ylabel =MyyLabelIR, xlabel = MyxLabelIR)
            axs[1].yaxis.set_label_position("right")
            axs[0].spines['right'].set_visible(False)
            axs[1].spines['left'].set_visible(False)
            axs[1].yaxis.tick_right()
            axs[1].axhline(linewidth=1, color='grey', dashes=[2,2])
            #axs[1].set_ylim(axs[0].get_ylim()) #ymin, ymax = axs[0].get_ylim()  #axs[1].set_ylim(ymin,ymax)
        fig6.show()
        
        #Tmp save plot
        #plt.savefig(MyProteinFile[0:-4]+'_Loaded' + '_fig.png') #, facecolor='yellow')
        
        return True
    
    def startSSCalc(self):
        self.labelMessage.setText('')
        #self.show()
        self.show(); QtWidgets.QApplication.processEvents()
        AU_or_BBK = ''
        CD_or_IR = ''
        if self.radioButtonRefSetBBK.isChecked():
            AU_or_BBK = "BBK"
        elif self.radioButtonRefSetAU.isChecked():
            AU_or_BBK = "AU"

        if self.radioButtonCD.isChecked():
            CD_or_IR ='CD'
        elif self.radioButtonIR.isChecked():
            CD_or_IR ='IR'
        elif self.radioButtonCDIR.isChecked():
            CD_or_IR ='CD-IR'
            if self.checkBoxUseCDIR2cm.isChecked():
                CD_or_IR ='CD-IR2cm'
        
        delimIs ='\t' # radioButtonASCIIXY / radioButtonCSVXY
        if self.radioButtonCSVXY.isChecked():
            delimIs = ','
        info = "" #self.labelMessage.text()
        newInfo = "Which reference set: " + AU_or_BBK #newInfo = info +'\n' + "AU_or_BBK is: " + AU_or_BBK
        newInfo += '\n' + "Spectrum type is: " + CD_or_IR
        self.labelMessage.setText(newInfo)
        self.show(); QtWidgets.QApplication.processEvents()
        
        fileName = self.labelFileName.text()
        MyProteinFile = fileName.split('/')[-1]
        MyProteinDir = fileName[:(len(fileName) - len(MyProteinFile))]
        print('MyProteinDir' + MyProteinDir)
        if (fileName == ""):
            info = self.labelMessage.text()
            newInfo = info + '\n' + 'NO FILE LOADED. CANT RUN Selcon3/CDSSTR'
            self.labelMessage.setText(newInfo)
            return None

        info = self.labelMessage.text()
        newInfo = info + '\n' + 'File is loaded. Running Selcon3'
        if self.radioButtonCDSSTR.isChecked():
            newInfo = info + '\n' + 'File is loaded. Running CDSSTR'
        self.labelMessage.setText(newInfo)
        QtWidgets.QApplication.processEvents()
        #Defaults
        RefDataHighWL = 240
        RefDataLowWL = 175 #SP175
        ShouldIPlot = 0; ShouldIAlsoPlotSelcon2 = 0
        if self.radioButtonRefSetSMP180.isChecked():
            RefDataLowWL = 180 #SMP180
        if self.checkBoxShouldIPlot.isChecked():
            ShouldIPlot = 1
        if self.checkBoxShouldIAlsoPlotSelcon2.isChecked():
            ShouldIAlsoPlotSelcon2 = 1
        
        okToRunAnalysis = True
        #Scaling
        scaleCD = 1
        if self.radioButtonCD.isChecked():
            if self.checkBoxScaleCD.isChecked():
                #First check for comma and replace with punktum
                scaleCDText = self.lineEditScaleFactor.text()
                self.lineEditScaleFactor.setText(scaleCDText.replace(",", "."))
                scaleCD = float(self.lineEditScaleFactor.text())
                if (scaleCD < 0.9 or scaleCD > 1.1): 
                    okToRunAnalysis = False
                    self.lineEditScaleFactor.setText('1.0')
                    self.checkBoxScaleCD.setChecked(False)
                    #info = self.labelMessage.text()
                    newInfo = info + '\n \n' + 'Error: CD scale factor outside limit.' '\n Analysis aborted'
                    #self.labelMessage.setText(newInfo)
                    #QtWidgets.QApplication.processEvents()
            pass
        else:
            #Ignore CD scaling for IR and CD-IR
            self.lineEditScaleFactor.setText('1.0')
            self.checkBoxScaleCD.setChecked(False)
        #END Scaling: if self.radioButtonCD.isChecked():
            
        #scaling IR in combined IR and CD
        if self.radioButtonCDIR.isChecked():
            scaleIRText = self.lineEditIRAbsMax.text()
            self.lineEditIRAbsMax.setText(scaleIRText.replace(",", "."))
            scaleIR = float(self.lineEditIRAbsMax.text())
            if (scaleIR < 1.0 and scaleIR !=0):
                okToRunAnalysis = False
                self.lineEditIRAbsMax.setText('1.0')
                newInfo = info + '\n \n' + 'Error: IR scale factor has to be >=1.' '\n Analysis aborted'
            #For now pass the scaling factor to SelconPy via the CD scaling fator
            #I.e. it is NOT possible to scale the CD in combined CD and IR mode i.e. if(CD_or_IR == 'CD-IR'):
            if (scaleIR == 0):
                info = self.labelMessage.text()
                newInfo = info + '\n' + 'IR spectrum scale = 0. Ignoring IR part of spectrum'
                self.labelMessage.setText(newInfo)
                QtWidgets.QApplication.processEvents()
            scaleCD = scaleIR
        #End: scaling IR in combined IR and CD
        
        
        if okToRunAnalysis:
            if self.radioButtonCDSSTR.isChecked():
                if (self.radioButtonIR.isChecked() or self.radioButtonCDIR.isChecked()):
                    info = self.labelMessage.text()
                    newInfo = info + '\n \n' + 'Error: Analysis of IR spectra with CDSSTR is not implemented'
                    self.labelMessage.setText(newInfo)
                    newInfo = ''
                else:
                    LoadAndRunCDSSTR(delimIs, CD_or_IR, AU_or_BBK, RefDataLowWL, RefDataHighWL, MyProteinDir, MyProteinFile, ShouldIPlot, ShouldIAlsoPlotSelcon2, scaleCD)
                    newInfo = '\n \n' + 'Analysis is completed'
            elif self.radioButtonSelcon3.isChecked():
                #def LoadAndRun(CD_or_IR, AU_or_BBK, RefDataLowWL, RefDataHighWL, MyProteinDir, MyProteinFile, ShouldIPlot, ShouldIAlsoPlotSelcon2):
                LoadAndRun(delimIs, CD_or_IR, AU_or_BBK, RefDataLowWL, RefDataHighWL, MyProteinDir, MyProteinFile, ShouldIPlot, ShouldIAlsoPlotSelcon2, scaleCD)
                newInfo = '\n \n' + 'Analysis is completed'
            else:
                #Really no other options
                pass
        #END: if okToRunAnalysis:
        
        #Info that it is done.
        info = self.labelMessage.text()
        #newInfo = info + '\n \n' + 'Analysis is completed'
        self.labelMessage.setText(info + newInfo)
        self.show()
        #But did it go ok?
        print("End of calculations")
        
    #END: def startSSCalc(self):    
        
    def closeEvent(self, event):
        plt.close('all')
        #sys.exit(0) #Uncomment if exit doesn't work
        

if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    #app.setWindowIcon(QtGui.QIcon('Logo_Icons/SS-Icon.png')) #self.setWindowIcon(QtGui.QIcon('Logo_Icons/SS-Icon.png'))
    #app.setQuitOnLastWindowClosed(True)
    window = MyApp()
    
    window.show()
    #window.setWindowState(window.windowState() & ~QtCore.Qt.WindowMinimized | QtCore.Qt.WindowActive)
    window.raise_()
    window.activateWindow()
    sys.exit(app.exec())