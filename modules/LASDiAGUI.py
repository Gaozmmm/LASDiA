# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file '.\modules\LASDiAGUI.ui'
#
# Created by: PyQt5 UI code generator 5.7
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_LASDiAGUI(object):
    def setupUi(self, LASDiAGUI):
        LASDiAGUI.setObjectName("LASDiAGUI")
        LASDiAGUI.resize(1487, 931)
        LASDiAGUI.setStyleSheet("background-color:rgb(255, 255, 255)")
        self.experimentPanel = QtWidgets.QTabWidget(LASDiAGUI)
        self.experimentPanel.setGeometry(QtCore.QRect(760, 10, 701, 381))
        self.experimentPanel.setAutoFillBackground(False)
        self.experimentPanel.setStyleSheet("background-color:rgb(97, 126, 255)")
        self.experimentPanel.setObjectName("experimentPanel")
        self.fileTab = QtWidgets.QWidget()
        self.fileTab.setObjectName("fileTab")
        self.dataPathName = QtWidgets.QTextEdit(self.fileTab)
        self.dataPathName.setGeometry(QtCore.QRect(140, 10, 381, 31))
        self.dataPathName.setStyleSheet("background-color:rgb(255, 255, 255);color:rgb(0, 0, 0)")
        self.dataPathName.setObjectName("dataPathName")
        self.dataPath = QtWidgets.QPushButton(self.fileTab)
        self.dataPath.setGeometry(QtCore.QRect(10, 10, 115, 32))
        self.dataPath.setStyleSheet("color:rgb(255, 255, 255);border: 2px solid rgb(28, 31, 255)")
        self.dataPath.setObjectName("dataPath")
        self.dataBox = QtWidgets.QGroupBox(self.fileTab)
        self.dataBox.setGeometry(QtCore.QRect(10, 70, 331, 101))
        self.dataBox.setAutoFillBackground(False)
        self.dataBox.setStyleSheet("background-color:rgb(97, 126, 255);color:rgb(255, 255, 255)")
        self.dataBox.setObjectName("dataBox")
        self.importData = QtWidgets.QPushButton(self.dataBox)
        self.importData.setGeometry(QtCore.QRect(240, 20, 71, 31))
        self.importData.setAutoFillBackground(False)
        self.importData.setStyleSheet("background-color:rgb(97, 126, 255);border: 2px solid rgb(28, 31, 255)")
        self.importData.setObjectName("importData")
        self.dataFileName = QtWidgets.QTextEdit(self.dataBox)
        self.dataFileName.setGeometry(QtCore.QRect(111, 20, 111, 31))
        self.dataFileName.setStyleSheet("background-color:rgb(255, 255, 255);color:rgb(0, 0, 0)")
        self.dataFileName.setObjectName("dataFileName")
        self.removePeaksData = QtWidgets.QPushButton(self.dataBox)
        self.removePeaksData.setGeometry(QtCore.QRect(10, 60, 115, 32))
        self.removePeaksData.setStyleSheet("background-color:rgb(97, 126, 255);border: 2px solid rgb(28, 31, 255)")
        self.removePeaksData.setObjectName("removePeaksData")
        self.dataFileName_2 = QtWidgets.QLabel(self.dataBox)
        self.dataFileName_2.setGeometry(QtCore.QRect(20, 30, 46, 16))
        self.dataFileName_2.setObjectName("dataFileName_2")
        self.bkgBox = QtWidgets.QGroupBox(self.fileTab)
        self.bkgBox.setGeometry(QtCore.QRect(350, 70, 331, 101))
        self.bkgBox.setAutoFillBackground(False)
        self.bkgBox.setStyleSheet("color:rgb(255, 255, 255)")
        self.bkgBox.setObjectName("bkgBox")
        self.importBkg = QtWidgets.QPushButton(self.bkgBox)
        self.importBkg.setGeometry(QtCore.QRect(240, 20, 71, 31))
        self.importBkg.setStyleSheet("background-color:rgb(97, 126, 255);border: 2px solid rgb(28, 31, 255)")
        self.importBkg.setObjectName("importBkg")
        self.removePeaksBkg = QtWidgets.QPushButton(self.bkgBox)
        self.removePeaksBkg.setGeometry(QtCore.QRect(20, 60, 115, 32))
        self.removePeaksBkg.setStyleSheet("background-color:rgb(97, 126, 255);border: 2px solid rgb(28, 31, 255)")
        self.removePeaksBkg.setObjectName("removePeaksBkg")
        self.bkgFileName_2 = QtWidgets.QLabel(self.bkgBox)
        self.bkgFileName_2.setGeometry(QtCore.QRect(20, 30, 61, 16))
        self.bkgFileName_2.setObjectName("bkgFileName_2")
        self.bkgFileName = QtWidgets.QTextEdit(self.bkgBox)
        self.bkgFileName.setGeometry(QtCore.QRect(110, 20, 111, 31))
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.bkgFileName.sizePolicy().hasHeightForWidth())
        self.bkgFileName.setSizePolicy(sizePolicy)
        self.bkgFileName.setStyleSheet("background-color:rgb(255, 255, 255);color:rgb(0, 0, 0)")
        self.bkgFileName.setObjectName("bkgFileName")
        self.interpolationBox = QtWidgets.QGroupBox(self.fileTab)
        self.interpolationBox.setGeometry(QtCore.QRect(10, 200, 671, 141))
        self.interpolationBox.setAutoFillBackground(False)
        self.interpolationBox.setStyleSheet("background-color:rgb(97, 126, 255);color:rgb(255, 255, 255)")
        self.interpolationBox.setObjectName("interpolationBox")
        self.interpolationPoints = QtWidgets.QSpinBox(self.interpolationBox)
        self.interpolationPoints.setGeometry(QtCore.QRect(100, 40, 111, 22))
        self.interpolationPoints.setStyleSheet("background-color:rgb(255, 255, 255);color:rgb(0, 0, 0)")
        self.interpolationPoints.setMaximum(10000)
        self.interpolationPoints.setProperty("value", 2048)
        self.interpolationPoints.setObjectName("interpolationPoints")
        self.interpolationPoints_2 = QtWidgets.QLabel(self.interpolationBox)
        self.interpolationPoints_2.setGeometry(QtCore.QRect(10, 40, 91, 16))
        self.interpolationPoints_2.setObjectName("interpolationPoints_2")
        self.QrangeBox = QtWidgets.QGroupBox(self.interpolationBox)
        self.QrangeBox.setGeometry(QtCore.QRect(260, 20, 171, 111))
        self.QrangeBox.setStyleSheet("")
        self.QrangeBox.setObjectName("QrangeBox")
        self.verticalLayout = QtWidgets.QVBoxLayout(self.QrangeBox)
        self.verticalLayout.setObjectName("verticalLayout")
        self.gridLayout = QtWidgets.QGridLayout()
        self.gridLayout.setObjectName("gridLayout")
        self.minQ_2 = QtWidgets.QLabel(self.QrangeBox)
        self.minQ_2.setObjectName("minQ_2")
        self.gridLayout.addWidget(self.minQ_2, 0, 0, 1, 1)
        self.minQ = QtWidgets.QDoubleSpinBox(self.QrangeBox)
        self.minQ.setStyleSheet("background-color:rgb(255, 255, 255);color:rgb(0, 0, 0)")
        self.minQ.setMaximum(1000.0)
        self.minQ.setProperty("value", 3.0)
        self.minQ.setObjectName("minQ")
        self.gridLayout.addWidget(self.minQ, 0, 1, 1, 1)
        self.maxQ_2 = QtWidgets.QLabel(self.QrangeBox)
        self.maxQ_2.setObjectName("maxQ_2")
        self.gridLayout.addWidget(self.maxQ_2, 1, 0, 1, 1)
        self.maxQ = QtWidgets.QDoubleSpinBox(self.QrangeBox)
        self.maxQ.setStyleSheet("background-color:rgb(255, 255, 255);color:rgb(0, 0, 0)")
        self.maxQ.setMaximum(1000.0)
        self.maxQ.setProperty("value", 109.0)
        self.maxQ.setObjectName("maxQ")
        self.gridLayout.addWidget(self.maxQ, 1, 1, 1, 1)
        self.verticalLayout.addLayout(self.gridLayout)
        self.dataInterpolation = QtWidgets.QPushButton(self.interpolationBox)
        self.dataInterpolation.setGeometry(QtCore.QRect(500, 50, 115, 32))
        self.dataInterpolation.setStyleSheet("background-color:rgb(97, 126, 255);border: 2px solid rgb(28, 31, 255)")
        self.dataInterpolation.setObjectName("dataInterpolation")
        self.experimentPanel.addTab(self.fileTab, "")
        self.geometryTab = QtWidgets.QWidget()
        self.geometryTab.setObjectName("geometryTab")
        self.xrayEnergy_2 = QtWidgets.QLabel(self.geometryTab)
        self.xrayEnergy_2.setGeometry(QtCore.QRect(30, 30, 121, 16))
        self.xrayEnergy_2.setStyleSheet("color:rgb(255, 255, 255)")
        self.xrayEnergy_2.setObjectName("xrayEnergy_2")
        self.xrayEnergy = QtWidgets.QDoubleSpinBox(self.geometryTab)
        self.xrayEnergy.setGeometry(QtCore.QRect(160, 30, 66, 24))
        self.xrayEnergy.setStyleSheet("background-color:rgb(255, 255, 255)")
        self.xrayEnergy.setObjectName("xrayEnergy")
        self.dacCorrection = QtWidgets.QGroupBox(self.geometryTab)
        self.dacCorrection.setGeometry(QtCore.QRect(20, 90, 461, 51))
        self.dacCorrection.setStyleSheet("background-color:rgb(97, 126, 255);color:rgb(255, 255, 255)")
        self.dacCorrection.setObjectName("dacCorrection")
        self.dacThickness = QtWidgets.QDoubleSpinBox(self.dacCorrection)
        self.dacThickness.setGeometry(QtCore.QRect(140, 20, 66, 24))
        self.dacThickness.setStyleSheet("background-color:rgb(255, 255, 255);color:rgb(0, 0, 0)")
        self.dacThickness.setDecimals(4)
        self.dacThickness.setMaximum(99.9999)
        self.dacThickness.setObjectName("dacThickness")
        self.dacThickness_2 = QtWidgets.QLabel(self.dacCorrection)
        self.dacThickness_2.setGeometry(QtCore.QRect(10, 20, 131, 16))
        self.dacThickness_2.setObjectName("dacThickness_2")
        self.dacAngle = QtWidgets.QDoubleSpinBox(self.dacCorrection)
        self.dacAngle.setGeometry(QtCore.QRect(370, 20, 66, 24))
        self.dacAngle.setStyleSheet("background-color:rgb(255, 255, 255);color:rgb(0, 0, 0)")
        self.dacAngle.setObjectName("dacAngle")
        self.dacAngle_2 = QtWidgets.QLabel(self.dacCorrection)
        self.dacAngle_2.setGeometry(QtCore.QRect(250, 20, 111, 16))
        self.dacAngle_2.setObjectName("dacAngle_2")
        self.MCCCorrection = QtWidgets.QGroupBox(self.geometryTab)
        self.MCCCorrection.setGeometry(QtCore.QRect(20, 180, 461, 91))
        self.MCCCorrection.setStyleSheet("background-color:rgb(97, 126, 255);color:rgb(255, 255, 255)")
        self.MCCCorrection.setObjectName("MCCCorrection")
        self.ws1_2 = QtWidgets.QLabel(self.MCCCorrection)
        self.ws1_2.setGeometry(QtCore.QRect(10, 20, 61, 16))
        self.ws1_2.setObjectName("ws1_2")
        self.ws1 = QtWidgets.QDoubleSpinBox(self.MCCCorrection)
        self.ws1.setGeometry(QtCore.QRect(80, 20, 66, 24))
        self.ws1.setStyleSheet("background-color:rgb(255, 255, 255);color:rgb(0, 0, 0)")
        self.ws1.setDecimals(4)
        self.ws1.setMaximum(99.9999)
        self.ws1.setObjectName("ws1")
        self.r1 = QtWidgets.QDoubleSpinBox(self.MCCCorrection)
        self.r1.setGeometry(QtCore.QRect(80, 50, 66, 24))
        self.r1.setStyleSheet("background-color:rgb(255, 255, 255);color:rgb(0, 0, 0)")
        self.r1.setDecimals(4)
        self.r1.setMaximum(99.9999)
        self.r1.setObjectName("r1")
        self.r2 = QtWidgets.QDoubleSpinBox(self.MCCCorrection)
        self.r2.setGeometry(QtCore.QRect(230, 50, 66, 24))
        self.r2.setStyleSheet("background-color:rgb(255, 255, 255);color:rgb(0, 0, 0)")
        self.r2.setDecimals(4)
        self.r2.setMaximum(99.9999)
        self.r2.setObjectName("r2")
        self.d_2 = QtWidgets.QLabel(self.MCCCorrection)
        self.d_2.setGeometry(QtCore.QRect(320, 20, 41, 16))
        self.d_2.setObjectName("d_2")
        self.ws2_2 = QtWidgets.QLabel(self.MCCCorrection)
        self.ws2_2.setGeometry(QtCore.QRect(160, 20, 61, 16))
        self.ws2_2.setObjectName("ws2_2")
        self.ws2 = QtWidgets.QDoubleSpinBox(self.MCCCorrection)
        self.ws2.setGeometry(QtCore.QRect(230, 20, 66, 24))
        self.ws2.setStyleSheet("background-color:rgb(255, 255, 255);color:rgb(0, 0, 0)")
        self.ws2.setDecimals(4)
        self.ws2.setMaximum(99.9999)
        self.ws2.setObjectName("ws2")
        self.r1_2 = QtWidgets.QLabel(self.MCCCorrection)
        self.r1_2.setGeometry(QtCore.QRect(10, 50, 51, 16))
        self.r1_2.setObjectName("r1_2")
        self.d = QtWidgets.QDoubleSpinBox(self.MCCCorrection)
        self.d.setGeometry(QtCore.QRect(370, 20, 66, 24))
        self.d.setStyleSheet("background-color:rgb(255, 255, 255);color:rgb(0, 0, 0)")
        self.d.setDecimals(4)
        self.d.setMaximum(99.9999)
        self.d.setObjectName("d")
        self.r2_2 = QtWidgets.QLabel(self.MCCCorrection)
        self.r2_2.setGeometry(QtCore.QRect(160, 50, 51, 16))
        self.r2_2.setObjectName("r2_2")
        self.experimentPanel.addTab(self.geometryTab, "")
        self.sampleTab = QtWidgets.QWidget()
        self.sampleTab.setObjectName("sampleTab")
        self.importXYZFile = QtWidgets.QPushButton(self.sampleTab)
        self.importXYZFile.setGeometry(QtCore.QRect(10, 260, 115, 32))
        self.importXYZFile.setStyleSheet("color:rgb(255, 255, 255);border: 2px solid rgb(28, 31, 255)")
        self.importXYZFile.setObjectName("importXYZFile")
        self.sampleComposition = QtWidgets.QPushButton(self.sampleTab)
        self.sampleComposition.setGeometry(QtCore.QRect(10, 20, 75, 23))
        self.sampleComposition.setStyleSheet("color:rgb(255, 255, 255);border: 2px solid rgb(28, 31, 255)")
        self.sampleComposition.setObjectName("sampleComposition")
        self.xyzFileName = QtWidgets.QTextEdit(self.sampleTab)
        self.xyzFileName.setGeometry(QtCore.QRect(160, 260, 331, 31))
        self.xyzFileName.setStyleSheet("background-color:rgb(255, 255, 255);color:rgb(0, 0, 0)")
        self.xyzFileName.setObjectName("xyzFileName")
        self.composition = QtWidgets.QTextEdit(self.sampleTab)
        self.composition.setGeometry(QtCore.QRect(160, 20, 331, 31))
        self.composition.setStyleSheet("background-color:rgb(255, 255, 255);color:rgb(0, 0, 0)")
        self.composition.setObjectName("composition")
        self.experimentPanel.addTab(self.sampleTab, "")
        self.displayTab = QtWidgets.QWidget()
        self.displayTab.setObjectName("displayTab")
        self.middlePlotBox = QtWidgets.QGroupBox(self.displayTab)
        self.middlePlotBox.setGeometry(QtCore.QRect(20, 30, 291, 111))
        self.middlePlotBox.setStyleSheet("background-color:rgb(97, 126, 255);color:rgb(255, 255, 255)")
        self.middlePlotBox.setObjectName("middlePlotBox")
        self.formFactorCheck = QtWidgets.QCheckBox(self.middlePlotBox)
        self.formFactorCheck.setGeometry(QtCore.QRect(20, 30, 81, 17))
        self.formFactorCheck.setObjectName("formFactorCheck")
        self.SQCheck = QtWidgets.QCheckBox(self.middlePlotBox)
        self.SQCheck.setGeometry(QtCore.QRect(20, 70, 70, 17))
        self.SQCheck.setObjectName("SQCheck")
        self.incohCheck = QtWidgets.QCheckBox(self.middlePlotBox)
        self.incohCheck.setGeometry(QtCore.QRect(160, 30, 121, 17))
        self.incohCheck.setObjectName("incohCheck")
        self.QiQCheck = QtWidgets.QCheckBox(self.middlePlotBox)
        self.QiQCheck.setGeometry(QtCore.QRect(160, 70, 70, 17))
        self.QiQCheck.setObjectName("QiQCheck")
        self.bottomPlotBox = QtWidgets.QGroupBox(self.displayTab)
        self.bottomPlotBox.setGeometry(QtCore.QRect(20, 200, 111, 111))
        self.bottomPlotBox.setStyleSheet("background-color:rgb(97, 126, 255);color:rgb(255, 255, 255)")
        self.bottomPlotBox.setObjectName("bottomPlotBox")
        self.FrCheck = QtWidgets.QCheckBox(self.bottomPlotBox)
        self.FrCheck.setGeometry(QtCore.QRect(30, 30, 41, 17))
        self.FrCheck.setObjectName("FrCheck")
        self.grCheck = QtWidgets.QCheckBox(self.bottomPlotBox)
        self.grCheck.setGeometry(QtCore.QRect(30, 70, 41, 17))
        self.grCheck.setObjectName("grCheck")
        self.calcSQ_4 = QtWidgets.QPushButton(self.displayTab)
        self.calcSQ_4.setGeometry(QtCore.QRect(550, 100, 111, 23))
        self.calcSQ_4.setStyleSheet("background-color:rgb(97, 126, 255);border: 2px solid rgb(28, 31, 255)")
        self.calcSQ_4.setObjectName("calcSQ_4")
        self.calcSQ = QtWidgets.QPushButton(self.displayTab)
        self.calcSQ.setGeometry(QtCore.QRect(440, 140, 75, 23))
        self.calcSQ.setStyleSheet("background-color:rgb(97, 126, 255);border: 2px solid rgb(28, 31, 255)")
        self.calcSQ.setObjectName("calcSQ")
        self.calcSQ_2 = QtWidgets.QPushButton(self.displayTab)
        self.calcSQ_2.setGeometry(QtCore.QRect(550, 140, 75, 23))
        self.calcSQ_2.setObjectName("calcSQ_2")
        self.calcSQ_3 = QtWidgets.QPushButton(self.displayTab)
        self.calcSQ_3.setGeometry(QtCore.QRect(440, 100, 75, 23))
        self.calcSQ_3.setStyleSheet("background-color:rgb(97, 126, 255);border: 2px solid rgb(28, 31, 255)")
        self.calcSQ_3.setObjectName("calcSQ_3")
        self.calcFr = QtWidgets.QPushButton(self.displayTab)
        self.calcFr.setGeometry(QtCore.QRect(540, 230, 75, 23))
        self.calcFr.setObjectName("calcFr")
        self.calcFr_2 = QtWidgets.QPushButton(self.displayTab)
        self.calcFr_2.setGeometry(QtCore.QRect(430, 230, 75, 23))
        self.calcFr_2.setObjectName("calcFr_2")
        self.experimentPanel.addTab(self.displayTab, "")
        self.overlayTab = QtWidgets.QWidget()
        self.overlayTab.setObjectName("overlayTab")
        self.experimentPanel.addTab(self.overlayTab, "")
        self.analysisPanel = QtWidgets.QGroupBox(LASDiAGUI)
        self.analysisPanel.setGeometry(QtCore.QRect(760, 410, 701, 511))
        self.analysisPanel.setStyleSheet("background-color:rgb(97, 126, 255);color:rgb(255, 255, 255)")
        self.analysisPanel.setObjectName("analysisPanel")
        self.chi2_plot = mplwidget.MplWidget(self.analysisPanel)
        self.chi2_plot.setGeometry(QtCore.QRect(300, 130, 391, 361))
        self.chi2_plot.setObjectName("chi2_plot")
        self.refThickness_2 = QtWidgets.QLabel(self.analysisPanel)
        self.refThickness_2.setGeometry(QtCore.QRect(15, 400, 101, 16))
        self.refThickness_2.setObjectName("refThickness_2")
        self.refThickness = QtWidgets.QDoubleSpinBox(self.analysisPanel)
        self.refThickness.setGeometry(QtCore.QRect(140, 400, 91, 24))
        self.refThickness.setStyleSheet("background-color:rgb(255, 255, 255);color:rgb(0, 0, 0)")
        self.refThickness.setDecimals(4)
        self.refThickness.setMaximum(1000.0)
        self.refThickness.setObjectName("refThickness")
        self.sampThickness_2 = QtWidgets.QLabel(self.analysisPanel)
        self.sampThickness_2.setGeometry(QtCore.QRect(15, 450, 121, 16))
        self.sampThickness_2.setObjectName("sampThickness_2")
        self.sampThickness = QtWidgets.QDoubleSpinBox(self.analysisPanel)
        self.sampThickness.setGeometry(QtCore.QRect(140, 450, 91, 24))
        self.sampThickness.setStyleSheet("background-color:rgb(255, 255, 255);color:rgb(0, 0, 0)")
        self.sampThickness.setDecimals(4)
        self.sampThickness.setMaximum(1000.0)
        self.sampThickness.setObjectName("sampThickness")
        self.optimize = QtWidgets.QPushButton(self.analysisPanel)
        self.optimize.setGeometry(QtCore.QRect(300, 80, 101, 31))
        self.optimize.setStyleSheet("background-color:rgb(97, 126, 255);border: 2px solid rgb(28, 31, 255)")
        self.optimize.setObjectName("optimize")
        self.QmaxIntegrate = QtWidgets.QDoubleSpinBox(self.analysisPanel)
        self.QmaxIntegrate.setGeometry(QtCore.QRect(610, 30, 81, 24))
        self.QmaxIntegrate.setStyleSheet("background-color:rgb(255, 255, 255);color:rgb(0, 0, 0)")
        self.QmaxIntegrate.setMaximum(1000.0)
        self.QmaxIntegrate.setProperty("value", 98.0)
        self.QmaxIntegrate.setObjectName("QmaxIntegrate")
        self.smoothingFactor = QtWidgets.QDoubleSpinBox(self.analysisPanel)
        self.smoothingFactor.setGeometry(QtCore.QRect(380, 30, 91, 24))
        self.smoothingFactor.setStyleSheet("background-color:rgb(255, 255, 255);color:rgb(0, 0, 0)")
        self.smoothingFactor.setMaximum(1000.0)
        self.smoothingFactor.setProperty("value", 0.25)
        self.smoothingFactor.setObjectName("smoothingFactor")
        self.QmaxIntegrate_2 = QtWidgets.QLabel(self.analysisPanel)
        self.QmaxIntegrate_2.setGeometry(QtCore.QRect(500, 30, 101, 31))
        self.QmaxIntegrate_2.setObjectName("QmaxIntegrate_2")
        self.smoothingFactor_2 = QtWidgets.QLabel(self.analysisPanel)
        self.smoothingFactor_2.setGeometry(QtCore.QRect(280, 30, 91, 16))
        self.smoothingFactor_2.setObjectName("smoothingFactor_2")
        self.refThicknessCheck = QtWidgets.QCheckBox(self.analysisPanel)
        self.refThicknessCheck.setGeometry(QtCore.QRect(250, 400, 41, 17))
        self.refThicknessCheck.setText("")
        self.refThicknessCheck.setChecked(True)
        self.refThicknessCheck.setObjectName("refThicknessCheck")
        self.sampThicknessCheck = QtWidgets.QCheckBox(self.analysisPanel)
        self.sampThicknessCheck.setGeometry(QtCore.QRect(250, 450, 31, 17))
        self.sampThicknessCheck.setText("")
        self.sampThicknessCheck.setChecked(False)
        self.sampThicknessCheck.setObjectName("sampThicknessCheck")
        self.dampingFunctBox = QtWidgets.QGroupBox(self.analysisPanel)
        self.dampingFunctBox.setGeometry(QtCore.QRect(10, 30, 241, 91))
        self.dampingFunctBox.setObjectName("dampingFunctBox")
        self.dampingFactor = QtWidgets.QDoubleSpinBox(self.dampingFunctBox)
        self.dampingFactor.setGeometry(QtCore.QRect(130, 20, 91, 24))
        self.dampingFactor.setStyleSheet("background-color:rgb(255, 255, 255);color:rgb(0, 0, 0)")
        self.dampingFactor.setMaximum(1000.0)
        self.dampingFactor.setProperty("value", 1.0)
        self.dampingFactor.setObjectName("dampingFactor")
        self.dampingFactor_2 = QtWidgets.QLabel(self.dampingFunctBox)
        self.dampingFactor_2.setGeometry(QtCore.QRect(10, 20, 81, 16))
        self.dampingFactor_2.setObjectName("dampingFactor_2")
        self.dampingFunction = QtWidgets.QComboBox(self.dampingFunctBox)
        self.dampingFunction.setGeometry(QtCore.QRect(130, 60, 91, 22))
        self.dampingFunction.setStyleSheet("background-color:rgb(255, 255, 255);color:rgb(0, 0, 0)")
        self.dampingFunction.setObjectName("dampingFunction")
        self.dampingFunction.addItem("")
        self.dampingFunction.addItem("")
        self.FrOptBox = QtWidgets.QGroupBox(self.analysisPanel)
        self.FrOptBox.setGeometry(QtCore.QRect(10, 250, 241, 111))
        self.FrOptBox.setObjectName("FrOptBox")
        self.rmin = QtWidgets.QDoubleSpinBox(self.FrOptBox)
        self.rmin.setGeometry(QtCore.QRect(130, 70, 91, 24))
        self.rmin.setStyleSheet("background-color:rgb(255, 255, 255);color:rgb(0, 0, 0)")
        self.rmin.setDecimals(3)
        self.rmin.setMaximum(1000.0)
        self.rmin.setProperty("value", 0.24)
        self.rmin.setObjectName("rmin")
        self.rmin_2 = QtWidgets.QLabel(self.FrOptBox)
        self.rmin_2.setGeometry(QtCore.QRect(10, 70, 71, 16))
        self.rmin_2.setObjectName("rmin_2")
        self.iterations = QtWidgets.QSpinBox(self.FrOptBox)
        self.iterations.setGeometry(QtCore.QRect(130, 20, 66, 24))
        self.iterations.setStyleSheet("background-color:rgb(255, 255, 255);color:rgb(0, 0, 0)")
        self.iterations.setMaximum(99)
        self.iterations.setProperty("value", 2)
        self.iterations.setObjectName("iterations")
        self.iterations_2 = QtWidgets.QLabel(self.FrOptBox)
        self.iterations_2.setGeometry(QtCore.QRect(10, 20, 61, 31))
        self.iterations_2.setObjectName("iterations_2")
        self.initialValueBox = QtWidgets.QGroupBox(self.analysisPanel)
        self.initialValueBox.setGeometry(QtCore.QRect(10, 130, 241, 101))
        self.initialValueBox.setObjectName("initialValueBox")
        self.sfValue_2 = QtWidgets.QLabel(self.initialValueBox)
        self.sfValue_2.setGeometry(QtCore.QRect(10, 60, 59, 16))
        self.sfValue_2.setObjectName("sfValue_2")
        self.rho0Value = QtWidgets.QDoubleSpinBox(self.initialValueBox)
        self.rho0Value.setGeometry(QtCore.QRect(140, 20, 91, 24))
        self.rho0Value.setStyleSheet("background-color:rgb(255, 255, 255);color:rgb(0, 0, 0)")
        self.rho0Value.setDecimals(5)
        self.rho0Value.setMaximum(1000.0)
        self.rho0Value.setProperty("value", 25.0584)
        self.rho0Value.setObjectName("rho0Value")
        self.sfValue = QtWidgets.QDoubleSpinBox(self.initialValueBox)
        self.sfValue.setGeometry(QtCore.QRect(140, 60, 91, 24))
        self.sfValue.setStyleSheet("background-color:rgb(255, 255, 255);color:rgb(0, 0, 0)")
        self.sfValue.setDecimals(5)
        self.sfValue.setMaximum(1000.0)
        self.sfValue.setProperty("value", 0.58536)
        self.sfValue.setObjectName("sfValue")
        self.rho0Value_2 = QtWidgets.QLabel(self.initialValueBox)
        self.rho0Value_2.setGeometry(QtCore.QRect(10, 20, 111, 16))
        self.rho0Value_2.setObjectName("rho0Value_2")
        self.graphicPanel = QtWidgets.QGroupBox(LASDiAGUI)
        self.graphicPanel.setGeometry(QtCore.QRect(0, 0, 751, 921))
        self.graphicPanel.setAutoFillBackground(False)
        self.graphicPanel.setStyleSheet("background-color:rgb(97, 126, 255);color:rgb(255, 255, 255)")
        self.graphicPanel.setObjectName("graphicPanel")
        self.factorPlot = mplwidget.MplWidget(self.graphicPanel)
        self.factorPlot.setGeometry(QtCore.QRect(10, 320, 731, 291))
        self.factorPlot.setObjectName("factorPlot")
        self.rawDataPlot = mplwidget.MplWidget(self.graphicPanel)
        self.rawDataPlot.setGeometry(QtCore.QRect(10, 20, 731, 301))
        self.rawDataPlot.setStyleSheet("")
        self.rawDataPlot.setObjectName("rawDataPlot")
        self.distfuncPlot = mplwidget.MplWidget(self.graphicPanel)
        self.distfuncPlot.setGeometry(QtCore.QRect(10, 610, 731, 291))
        self.distfuncPlot.setObjectName("distfuncPlot")
        self.analysisPanel.raise_()
        self.graphicPanel.raise_()
        self.experimentPanel.raise_()

        self.retranslateUi(LASDiAGUI)
        self.experimentPanel.setCurrentIndex(0)
        self.dampingFunction.setCurrentIndex(0)
        QtCore.QMetaObject.connectSlotsByName(LASDiAGUI)

    def retranslateUi(self, LASDiAGUI):
        _translate = QtCore.QCoreApplication.translate
        LASDiAGUI.setWindowTitle(_translate("LASDiAGUI", "LASDiA"))
        self.dataPath.setText(_translate("LASDiAGUI", "Data Path"))
        self.dataBox.setTitle(_translate("LASDiAGUI", "Sample"))
        self.importData.setText(_translate("LASDiAGUI", "Load"))
        self.removePeaksData.setText(_translate("LASDiAGUI", "Remove Peaks"))
        self.dataFileName_2.setText(_translate("LASDiAGUI", "File Name"))
        self.bkgBox.setTitle(_translate("LASDiAGUI", "Background"))
        self.importBkg.setText(_translate("LASDiAGUI", "Load"))
        self.removePeaksBkg.setText(_translate("LASDiAGUI", "Remove Peaks"))
        self.bkgFileName_2.setText(_translate("LASDiAGUI", "File Name"))
        self.interpolationBox.setTitle(_translate("LASDiAGUI", "Interpolation"))
        self.interpolationPoints_2.setText(_translate("LASDiAGUI", "Num of Points"))
        self.QrangeBox.setTitle(_translate("LASDiAGUI", "Q Range (nm-1)"))
        self.minQ_2.setText(_translate("LASDiAGUI", "<html><head/><body><p>Minimum</p></body></html>"))
        self.maxQ_2.setText(_translate("LASDiAGUI", "<html><head/><body><p>Maximim</p></body></html>"))
        self.dataInterpolation.setText(_translate("LASDiAGUI", "Interpolate Data"))
        self.experimentPanel.setTabText(self.experimentPanel.indexOf(self.fileTab), _translate("LASDiAGUI", "File"))
        self.xrayEnergy_2.setText(_translate("LASDiAGUI", "X-Ray Energy (keV)"))
        self.dacCorrection.setTitle(_translate("LASDiAGUI", "DAC Correction"))
        self.dacThickness_2.setText(_translate("LASDiAGUI", "Anvil thickness (cm)"))
        self.dacAngle_2.setText(_translate("LASDiAGUI", "DAC Angle (deg)"))
        self.MCCCorrection.setTitle(_translate("LASDiAGUI", "MCC Correction"))
        self.ws1_2.setText(_translate("LASDiAGUI", "WS1 (cm)"))
        self.d_2.setText(_translate("LASDiAGUI", "D (cm)"))
        self.ws2_2.setText(_translate("LASDiAGUI", "WS2 (cm)"))
        self.r1_2.setText(_translate("LASDiAGUI", "R1 (cm)"))
        self.r2_2.setText(_translate("LASDiAGUI", "R2 (cm)"))
        self.experimentPanel.setTabText(self.experimentPanel.indexOf(self.geometryTab), _translate("LASDiAGUI", "Geometry"))
        self.importXYZFile.setText(_translate("LASDiAGUI", "Import XYZ file"))
        self.sampleComposition.setText(_translate("LASDiAGUI", "Sample"))
        self.experimentPanel.setTabText(self.experimentPanel.indexOf(self.sampleTab), _translate("LASDiAGUI", "Sample"))
        self.middlePlotBox.setTitle(_translate("LASDiAGUI", "Middle Plot"))
        self.formFactorCheck.setText(_translate("LASDiAGUI", "Form Factor"))
        self.SQCheck.setText(_translate("LASDiAGUI", "S(Q)"))
        self.incohCheck.setText(_translate("LASDiAGUI", "Incoherent Intensity"))
        self.QiQCheck.setText(_translate("LASDiAGUI", "Qi(Q)"))
        self.bottomPlotBox.setTitle(_translate("LASDiAGUI", "Bottom Plot"))
        self.FrCheck.setText(_translate("LASDiAGUI", "F(r)"))
        self.grCheck.setText(_translate("LASDiAGUI", "g(r)"))
        self.calcSQ_4.setText(_translate("LASDiAGUI", "Incoherent Intensity"))
        self.calcSQ.setText(_translate("LASDiAGUI", "S(Q)"))
        self.calcSQ_2.setText(_translate("LASDiAGUI", "Qi(Q)"))
        self.calcSQ_3.setText(_translate("LASDiAGUI", "Form Factor"))
        self.calcFr.setText(_translate("LASDiAGUI", "F(r)"))
        self.calcFr_2.setText(_translate("LASDiAGUI", "g(r)"))
        self.experimentPanel.setTabText(self.experimentPanel.indexOf(self.displayTab), _translate("LASDiAGUI", "Display"))
        self.experimentPanel.setTabText(self.experimentPanel.indexOf(self.overlayTab), _translate("LASDiAGUI", "Overlay"))
        self.analysisPanel.setTitle(_translate("LASDiAGUI", "Analysis Panel"))
        self.refThickness_2.setText(_translate("LASDiAGUI", "Ref Thickness (cm)"))
        self.sampThickness_2.setText(_translate("LASDiAGUI", "Samp Thinckness (cm)"))
        self.optimize.setText(_translate("LASDiAGUI", "Optimize"))
        self.QmaxIntegrate_2.setText(_translate("LASDiAGUI", "<html><head/><body><p>Max  Q value <br> for integration (nm<span style=\" vertical-align:super;\">-1</span>)</p></body></html>\n"
""))
        self.smoothingFactor_2.setText(_translate("LASDiAGUI", "Smoothing Factor"))
        self.dampingFunctBox.setTitle(_translate("LASDiAGUI", "Damping Function"))
        self.dampingFactor_2.setText(_translate("LASDiAGUI", "Damping Factor"))
        self.dampingFunction.setItemText(0, _translate("LASDiAGUI", "Exponential"))
        self.dampingFunction.setItemText(1, _translate("LASDiAGUI", "Lorch Function"))
        self.FrOptBox.setTitle(_translate("LASDiAGUI", "F(r) Optimization"))
        self.rmin_2.setText(_translate("LASDiAGUI", "r cutoff (nm)"))
        self.iterations_2.setText(_translate("LASDiAGUI", "Iterations"))
        self.initialValueBox.setTitle(_translate("LASDiAGUI", "Initial values"))
        self.sfValue_2.setText(_translate("LASDiAGUI", "Scale Factor"))
        self.rho0Value_2.setText(_translate("LASDiAGUI", "Density (atm/nm<sup>3</sup>)"))
        self.graphicPanel.setTitle(_translate("LASDiAGUI", "Graphic Panel"))

from modules import mplwidget
