# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'gui.ui'
#
# Created: Sun Apr  3 19:17:53 2016
#      by: pyside-uic 0.2.15 running on PySide 1.2.4
#
# WARNING! All changes made in this file will be lost!

from PySide import QtCore, QtGui

class Ui_LASDiAGui(object):
    def setupUi(self, LASDiAGui):
        LASDiAGui.setObjectName("LASDiAGui")
        LASDiAGui.resize(961, 459)
        self.layoutWidget = QtGui.QWidget(LASDiAGui)
        self.layoutWidget.setGeometry(QtCore.QRect(20, 10, 204, 218))
        self.layoutWidget.setObjectName("layoutWidget")
        self.verticalLayout = QtGui.QVBoxLayout(self.layoutWidget)
        self.verticalLayout.setContentsMargins(0, 0, 0, 0)
        self.verticalLayout.setObjectName("verticalLayout")
        self.LoadData = QtGui.QPushButton(self.layoutWidget)
        self.LoadData.setObjectName("LoadData")
        self.verticalLayout.addWidget(self.LoadData)
        self.CalcSQ = QtGui.QPushButton(self.layoutWidget)
        self.CalcSQ.setObjectName("CalcSQ")
        self.verticalLayout.addWidget(self.CalcSQ)
        self.CalcFr = QtGui.QPushButton(self.layoutWidget)
        self.CalcFr.setObjectName("CalcFr")
        self.verticalLayout.addWidget(self.CalcFr)
        self.groupBox_2 = QtGui.QGroupBox(LASDiAGui)
        self.groupBox_2.setGeometry(QtCore.QRect(330, 280, 201, 161))
        self.groupBox_2.setObjectName("groupBox_2")
        self.Optimization = QtGui.QPushButton(self.groupBox_2)
        self.Optimization.setGeometry(QtCore.QRect(0, 120, 144, 32))
        self.Optimization.setObjectName("Optimization")
        self.label_6 = QtGui.QLabel(self.groupBox_2)
        self.label_6.setGeometry(QtCore.QRect(14, 34, 51, 16))
        self.label_6.setObjectName("label_6")
        self.label_7 = QtGui.QLabel(self.groupBox_2)
        self.label_7.setGeometry(QtCore.QRect(14, 90, 50, 16))
        self.label_7.setObjectName("label_7")
        self.Iteration = QtGui.QSpinBox(self.groupBox_2)
        self.Iteration.setGeometry(QtCore.QRect(77, 34, 47, 24))
        self.Iteration.setObjectName("Iteration")
        self.rcutoff = QtGui.QDoubleSpinBox(self.groupBox_2)
        self.rcutoff.setGeometry(QtCore.QRect(77, 90, 66, 24))
        self.rcutoff.setLocale(QtCore.QLocale(QtCore.QLocale.C, QtCore.QLocale.AnyCountry))
        self.rcutoff.setObjectName("rcutoff")
        self.groupBox = QtGui.QGroupBox(LASDiAGui)
        self.groupBox.setGeometry(QtCore.QRect(700, 160, 251, 281))
        self.groupBox.setObjectName("groupBox")
        self.Minimization = QtGui.QPushButton(self.groupBox)
        self.Minimization.setGeometry(QtCore.QRect(60, 250, 115, 32))
        self.Minimization.setObjectName("Minimization")
        self.layoutWidget1 = QtGui.QWidget(self.groupBox)
        self.layoutWidget1.setGeometry(QtCore.QRect(20, 32, 179, 198))
        self.layoutWidget1.setObjectName("layoutWidget1")
        self.horizontalLayout = QtGui.QHBoxLayout(self.layoutWidget1)
        self.horizontalLayout.setContentsMargins(0, 0, 0, 0)
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.verticalLayout_3 = QtGui.QVBoxLayout()
        self.verticalLayout_3.setObjectName("verticalLayout_3")
        self.label_11 = QtGui.QLabel(self.layoutWidget1)
        self.label_11.setObjectName("label_11")
        self.verticalLayout_3.addWidget(self.label_11)
        self.label_12 = QtGui.QLabel(self.layoutWidget1)
        self.label_12.setObjectName("label_12")
        self.verticalLayout_3.addWidget(self.label_12)
        self.label_14 = QtGui.QLabel(self.layoutWidget1)
        self.label_14.setObjectName("label_14")
        self.verticalLayout_3.addWidget(self.label_14)
        self.label_13 = QtGui.QLabel(self.layoutWidget1)
        self.label_13.setObjectName("label_13")
        self.verticalLayout_3.addWidget(self.label_13)
        self.horizontalLayout.addLayout(self.verticalLayout_3)
        self.verticalLayout_2 = QtGui.QVBoxLayout()
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        self.rho0LowerLim = QtGui.QDoubleSpinBox(self.layoutWidget1)
        self.rho0LowerLim.setLocale(QtCore.QLocale(QtCore.QLocale.C, QtCore.QLocale.AnyCountry))
        self.rho0LowerLim.setObjectName("rho0LowerLim")
        self.verticalLayout_2.addWidget(self.rho0LowerLim)
        self.rho0UpperLim = QtGui.QDoubleSpinBox(self.layoutWidget1)
        self.rho0UpperLim.setLocale(QtCore.QLocale(QtCore.QLocale.C, QtCore.QLocale.AnyCountry))
        self.rho0UpperLim.setObjectName("rho0UpperLim")
        self.verticalLayout_2.addWidget(self.rho0UpperLim)
        self.sLowerLim = QtGui.QDoubleSpinBox(self.layoutWidget1)
        self.sLowerLim.setLocale(QtCore.QLocale(QtCore.QLocale.C, QtCore.QLocale.AnyCountry))
        self.sLowerLim.setObjectName("sLowerLim")
        self.verticalLayout_2.addWidget(self.sLowerLim)
        self.sUpperLim = QtGui.QDoubleSpinBox(self.layoutWidget1)
        self.sUpperLim.setLocale(QtCore.QLocale(QtCore.QLocale.C, QtCore.QLocale.AnyCountry))
        self.sUpperLim.setObjectName("sUpperLim")
        self.verticalLayout_2.addWidget(self.sUpperLim)
        self.horizontalLayout.addLayout(self.verticalLayout_2)
        self.groupBox_5 = QtGui.QGroupBox(LASDiAGui)
        self.groupBox_5.setGeometry(QtCore.QRect(330, 10, 281, 91))
        self.groupBox_5.setObjectName("groupBox_5")
        self.verticalLayoutWidget = QtGui.QWidget(self.groupBox_5)
        self.verticalLayoutWidget.setGeometry(QtCore.QRect(160, 30, 114, 50))
        self.verticalLayoutWidget.setObjectName("verticalLayoutWidget")
        self.verticalLayout_10 = QtGui.QVBoxLayout(self.verticalLayoutWidget)
        self.verticalLayout_10.setContentsMargins(0, 0, 0, 0)
        self.verticalLayout_10.setObjectName("verticalLayout_10")
        self.label = QtGui.QLabel(self.verticalLayoutWidget)
        self.label.setObjectName("label")
        self.verticalLayout_10.addWidget(self.label)
        self.molDistance = QtGui.QDoubleSpinBox(self.verticalLayoutWidget)
        self.molDistance.setLocale(QtCore.QLocale(QtCore.QLocale.C, QtCore.QLocale.AnyCountry))
        self.molDistance.setDecimals(6)
        self.molDistance.setMaximum(99.999999)
        self.molDistance.setObjectName("molDistance")
        self.verticalLayout_10.addWidget(self.molDistance)
        self.verticalLayoutWidget_2 = QtGui.QWidget(self.groupBox_5)
        self.verticalLayoutWidget_2.setGeometry(QtCore.QRect(20, 30, 91, 49))
        self.verticalLayoutWidget_2.setObjectName("verticalLayoutWidget_2")
        self.verticalLayout_11 = QtGui.QVBoxLayout(self.verticalLayoutWidget_2)
        self.verticalLayout_11.setContentsMargins(0, 0, 0, 0)
        self.verticalLayout_11.setObjectName("verticalLayout_11")
        self.label_17 = QtGui.QLabel(self.verticalLayoutWidget_2)
        self.label_17.setObjectName("label_17")
        self.verticalLayout_11.addWidget(self.label_17)
        self.Molecule = QtGui.QComboBox(self.verticalLayoutWidget_2)
        self.Molecule.setObjectName("Molecule")
        self.Molecule.addItem("")
        self.Molecule.addItem("")
        self.Molecule.addItem("")
        self.verticalLayout_11.addWidget(self.Molecule)
        self.groupBox_6 = QtGui.QGroupBox(LASDiAGui)
        self.groupBox_6.setGeometry(QtCore.QRect(710, 10, 241, 131))
        self.groupBox_6.setObjectName("groupBox_6")
        self.layoutWidget2 = QtGui.QWidget(self.groupBox_6)
        self.layoutWidget2.setGeometry(QtCore.QRect(10, 30, 211, 96))
        self.layoutWidget2.setObjectName("layoutWidget2")
        self.horizontalLayout_4 = QtGui.QHBoxLayout(self.layoutWidget2)
        self.horizontalLayout_4.setContentsMargins(0, 0, 0, 0)
        self.horizontalLayout_4.setObjectName("horizontalLayout_4")
        self.verticalLayout_7 = QtGui.QVBoxLayout()
        self.verticalLayout_7.setObjectName("verticalLayout_7")
        self.label_2 = QtGui.QLabel(self.layoutWidget2)
        self.label_2.setObjectName("label_2")
        self.verticalLayout_7.addWidget(self.label_2)
        self.label_3 = QtGui.QLabel(self.layoutWidget2)
        self.label_3.setObjectName("label_3")
        self.verticalLayout_7.addWidget(self.label_3)
        self.label_4 = QtGui.QLabel(self.layoutWidget2)
        self.label_4.setObjectName("label_4")
        self.verticalLayout_7.addWidget(self.label_4)
        self.horizontalLayout_4.addLayout(self.verticalLayout_7)
        self.verticalLayout_8 = QtGui.QVBoxLayout()
        self.verticalLayout_8.setObjectName("verticalLayout_8")
        self.minQ = QtGui.QDoubleSpinBox(self.layoutWidget2)
        self.minQ.setLocale(QtCore.QLocale(QtCore.QLocale.C, QtCore.QLocale.AnyCountry))
        self.minQ.setDecimals(6)
        self.minQ.setMaximum(999.999999)
        self.minQ.setObjectName("minQ")
        self.verticalLayout_8.addWidget(self.minQ)
        self.maxQ = QtGui.QDoubleSpinBox(self.layoutWidget2)
        self.maxQ.setLocale(QtCore.QLocale(QtCore.QLocale.C, QtCore.QLocale.AnyCountry))
        self.maxQ.setDecimals(6)
        self.maxQ.setMaximum(999.999999)
        self.maxQ.setObjectName("maxQ")
        self.verticalLayout_8.addWidget(self.maxQ)
        self.QmaxIntegrate = QtGui.QDoubleSpinBox(self.layoutWidget2)
        self.QmaxIntegrate.setLocale(QtCore.QLocale(QtCore.QLocale.C, QtCore.QLocale.AnyCountry))
        self.QmaxIntegrate.setDecimals(6)
        self.QmaxIntegrate.setMaximum(999.999999)
        self.QmaxIntegrate.setObjectName("QmaxIntegrate")
        self.verticalLayout_8.addWidget(self.QmaxIntegrate)
        self.horizontalLayout_4.addLayout(self.verticalLayout_8)
        self.groupBox_3 = QtGui.QGroupBox(LASDiAGui)
        self.groupBox_3.setGeometry(QtCore.QRect(330, 130, 191, 131))
        self.groupBox_3.setObjectName("groupBox_3")
        self.layoutWidget3 = QtGui.QWidget(self.groupBox_3)
        self.layoutWidget3.setGeometry(QtCore.QRect(0, 30, 177, 81))
        self.layoutWidget3.setObjectName("layoutWidget3")
        self.horizontalLayout_2 = QtGui.QHBoxLayout(self.layoutWidget3)
        self.horizontalLayout_2.setContentsMargins(0, 0, 0, 0)
        self.horizontalLayout_2.setObjectName("horizontalLayout_2")
        self.verticalLayout_4 = QtGui.QVBoxLayout()
        self.verticalLayout_4.setObjectName("verticalLayout_4")
        self.label_15 = QtGui.QLabel(self.layoutWidget3)
        self.label_15.setObjectName("label_15")
        self.verticalLayout_4.addWidget(self.label_15)
        self.label_16 = QtGui.QLabel(self.layoutWidget3)
        self.label_16.setObjectName("label_16")
        self.verticalLayout_4.addWidget(self.label_16)
        self.horizontalLayout_2.addLayout(self.verticalLayout_4)
        self.verticalLayout_5 = QtGui.QVBoxLayout()
        self.verticalLayout_5.setObjectName("verticalLayout_5")
        self.sValue = QtGui.QDoubleSpinBox(self.layoutWidget3)
        self.sValue.setLocale(QtCore.QLocale(QtCore.QLocale.C, QtCore.QLocale.AnyCountry))
        self.sValue.setDecimals(6)
        self.sValue.setMaximum(99.999999)
        self.sValue.setObjectName("sValue")
        self.verticalLayout_5.addWidget(self.sValue)
        self.rho0Value = QtGui.QDoubleSpinBox(self.layoutWidget3)
        self.rho0Value.setLocale(QtCore.QLocale(QtCore.QLocale.C, QtCore.QLocale.AnyCountry))
        self.rho0Value.setDecimals(6)
        self.rho0Value.setObjectName("rho0Value")
        self.verticalLayout_5.addWidget(self.rho0Value)
        self.horizontalLayout_2.addLayout(self.verticalLayout_5)
        self.groupBox_4 = QtGui.QGroupBox(LASDiAGui)
        self.groupBox_4.setGeometry(QtCore.QRect(40, 260, 161, 91))
        self.groupBox_4.setObjectName("groupBox_4")
        self.layoutWidget4 = QtGui.QWidget(self.groupBox_4)
        self.layoutWidget4.setGeometry(QtCore.QRect(10, 20, 138, 62))
        self.layoutWidget4.setObjectName("layoutWidget4")
        self.horizontalLayout_3 = QtGui.QHBoxLayout(self.layoutWidget4)
        self.horizontalLayout_3.setContentsMargins(0, 0, 0, 0)
        self.horizontalLayout_3.setObjectName("horizontalLayout_3")
        self.verticalLayout_9 = QtGui.QVBoxLayout()
        self.verticalLayout_9.setObjectName("verticalLayout_9")
        self.label_5 = QtGui.QLabel(self.layoutWidget4)
        self.label_5.setObjectName("label_5")
        self.verticalLayout_9.addWidget(self.label_5)
        self.label_8 = QtGui.QLabel(self.layoutWidget4)
        self.label_8.setObjectName("label_8")
        self.verticalLayout_9.addWidget(self.label_8)
        self.horizontalLayout_3.addLayout(self.verticalLayout_9)
        self.verticalLayout_6 = QtGui.QVBoxLayout()
        self.verticalLayout_6.setObjectName("verticalLayout_6")
        self.smoothFactor = QtGui.QDoubleSpinBox(self.layoutWidget4)
        self.smoothFactor.setLocale(QtCore.QLocale(QtCore.QLocale.C, QtCore.QLocale.AnyCountry))
        self.smoothFactor.setDecimals(2)
        self.smoothFactor.setMaximum(99.99)
        self.smoothFactor.setObjectName("smoothFactor")
        self.verticalLayout_6.addWidget(self.smoothFactor)
        self.dampingFactor = QtGui.QDoubleSpinBox(self.layoutWidget4)
        self.dampingFactor.setLocale(QtCore.QLocale(QtCore.QLocale.C, QtCore.QLocale.AnyCountry))
        self.dampingFactor.setDecimals(2)
        self.dampingFactor.setMaximum(99.99)
        self.dampingFactor.setObjectName("dampingFactor")
        self.verticalLayout_6.addWidget(self.dampingFactor)
        self.horizontalLayout_3.addLayout(self.verticalLayout_6)

        self.retranslateUi(LASDiAGui)
        QtCore.QMetaObject.connectSlotsByName(LASDiAGui)

    def retranslateUi(self, LASDiAGui):
        LASDiAGui.setWindowTitle(QtGui.QApplication.translate("LASDiAGui", "LASDiAGui", None, QtGui.QApplication.UnicodeUTF8))
        self.LoadData.setText(QtGui.QApplication.translate("LASDiAGui", "Load Data", None, QtGui.QApplication.UnicodeUTF8))
        self.CalcSQ.setText(QtGui.QApplication.translate("LASDiAGui", "S(Q)", None, QtGui.QApplication.UnicodeUTF8))
        self.CalcFr.setText(QtGui.QApplication.translate("LASDiAGui", "F(r)", None, QtGui.QApplication.UnicodeUTF8))
        self.groupBox_2.setTitle(QtGui.QApplication.translate("LASDiAGui", "F(r) Optimization Parameters", None, QtGui.QApplication.UnicodeUTF8))
        self.Optimization.setText(QtGui.QApplication.translate("LASDiAGui", "Optimization", None, QtGui.QApplication.UnicodeUTF8))
        self.label_6.setText(QtGui.QApplication.translate("LASDiAGui", "Iteration", None, QtGui.QApplication.UnicodeUTF8))
        self.label_7.setText(QtGui.QApplication.translate("LASDiAGui", "r cut-off", None, QtGui.QApplication.UnicodeUTF8))
        self.groupBox.setTitle(QtGui.QApplication.translate("LASDiAGui", "Parameters for minimization", None, QtGui.QApplication.UnicodeUTF8))
        self.Minimization.setText(QtGui.QApplication.translate("LASDiAGui", "Minimization", None, QtGui.QApplication.UnicodeUTF8))
        self.label_11.setText(QtGui.QApplication.translate("LASDiAGui", "rho0 lower limit", None, QtGui.QApplication.UnicodeUTF8))
        self.label_12.setText(QtGui.QApplication.translate("LASDiAGui", "rho0 upper limit", None, QtGui.QApplication.UnicodeUTF8))
        self.label_14.setText(QtGui.QApplication.translate("LASDiAGui", "s lower limit", None, QtGui.QApplication.UnicodeUTF8))
        self.label_13.setText(QtGui.QApplication.translate("LASDiAGui", "s upper limit", None, QtGui.QApplication.UnicodeUTF8))
        self.groupBox_5.setTitle(QtGui.QApplication.translate("LASDiAGui", "Sample Composition", None, QtGui.QApplication.UnicodeUTF8))
        self.label.setText(QtGui.QApplication.translate("LASDiAGui", "Molecule distance", None, QtGui.QApplication.UnicodeUTF8))
        self.label_17.setText(QtGui.QApplication.translate("LASDiAGui", "Molecule", None, QtGui.QApplication.UnicodeUTF8))
        self.Molecule.setItemText(0, QtGui.QApplication.translate("LASDiAGui", "Ar", None, QtGui.QApplication.UnicodeUTF8))
        self.Molecule.setItemText(1, QtGui.QApplication.translate("LASDiAGui", "CO2", None, QtGui.QApplication.UnicodeUTF8))
        self.Molecule.setItemText(2, QtGui.QApplication.translate("LASDiAGui", "N2", None, QtGui.QApplication.UnicodeUTF8))
        self.groupBox_6.setTitle(QtGui.QApplication.translate("LASDiAGui", "Limits", None, QtGui.QApplication.UnicodeUTF8))
        self.label_2.setText(QtGui.QApplication.translate("LASDiAGui", "minQ", None, QtGui.QApplication.UnicodeUTF8))
        self.label_3.setText(QtGui.QApplication.translate("LASDiAGui", "maxQ", None, QtGui.QApplication.UnicodeUTF8))
        self.label_4.setText(QtGui.QApplication.translate("LASDiAGui", "QmaxIntegrate", None, QtGui.QApplication.UnicodeUTF8))
        self.groupBox_3.setTitle(QtGui.QApplication.translate("LASDiAGui", "Initial values", None, QtGui.QApplication.UnicodeUTF8))
        self.label_15.setText(QtGui.QApplication.translate("LASDiAGui", "s value", None, QtGui.QApplication.UnicodeUTF8))
        self.label_16.setText(QtGui.QApplication.translate("LASDiAGui", "rho0 value", None, QtGui.QApplication.UnicodeUTF8))
        self.groupBox_4.setTitle(QtGui.QApplication.translate("LASDiAGui", "Factors", None, QtGui.QApplication.UnicodeUTF8))
        self.label_5.setText(QtGui.QApplication.translate("LASDiAGui", "Smooth", None, QtGui.QApplication.UnicodeUTF8))
        self.label_8.setText(QtGui.QApplication.translate("LASDiAGui", "Damping", None, QtGui.QApplication.UnicodeUTF8))
