from PyQt5 import QtCore, QtWidgets
from PyQt5.QtWidgets import QApplication, QMainWindow
from PyQt5.QtCore import pyqtSlot

from myMassEl_v16 import *
import Utility

# from modules.myMassEl_v16 import *
# from modules import Utility

def reset_(self):
    self.rbr = 0
    self.n = 0
    self.num = 0
    self.ui.label_el_1.setText("-")
    self.ui.label_el_2.setText("-")
    self.ui.label_el_3.setText("-")
    self.ui.label_el_4.setText("-")
    self.ui.spinBox_el_1.setValue(0)
    self.ui.spinBox_el_2.setValue(0)
    self.ui.spinBox_el_3.setValue(0)
    self.ui.spinBox_el_4.setValue(0)
    self.mr_el_1 = 0
    self.mr_el_2 = 0
    self.mr_el_3 = 0
    self.mr_el_4 = 0
    self.ui.checkBox_Ru.setChecked(False)
    self.ui.checkBox_Re.setChecked(False)
    self.ui.checkBox_Rf.setChecked(False)
    self.ui.checkBox_Rg.setChecked(False)
    self.ui.checkBox_Ra.setChecked(False)
    self.ui.checkBox_Rb.setChecked(False)
    self.ui.checkBox_Rn.setChecked(False)
    self.ui.checkBox_Rh.setChecked(False)
    self.ui.checkBox_Be.setChecked(False)
    self.ui.checkBox_Ba.setChecked(False)
    self.ui.checkBox_Bh.setChecked(False)
    self.ui.checkBox_Bi.setChecked(False)
    self.ui.checkBox_Bk.setChecked(False)
    self.ui.checkBox_Br.setChecked(False)
    self.ui.checkBox_H.setChecked(False)
    self.ui.checkBox_P.setChecked(False)
    self.ui.checkBox_Os.setChecked(False)
    self.ui.checkBox_Es.setChecked(False)
    self.ui.checkBox_Hg.setChecked(False)
    self.ui.checkBox_Ge.setChecked(False)
    self.ui.checkBox_Gd.setChecked(False)
    self.ui.checkBox_Ga.setChecked(False)
    self.ui.checkBox_Pr.setChecked(False)
    self.ui.checkBox_Pt.setChecked(False)
    self.ui.checkBox_Pu.setChecked(False)
    self.ui.checkBox_C.setChecked(False)
    self.ui.checkBox_Pb.setChecked(False)
    self.ui.checkBox_Pa.setChecked(False)
    self.ui.checkBox_Pd.setChecked(False)
    self.ui.checkBox_Cd.setChecked(False)
    self.ui.checkBox_Po.setChecked(False)
    self.ui.checkBox_Pm.setChecked(False)
    self.ui.checkBox_Hs.setChecked(False)
    self.ui.checkBox_Uup.setChecked(False)
    self.ui.checkBox_Uus.setChecked(False)
    self.ui.checkBox_Uuo.setChecked(False)
    self.ui.checkBox_Ho.setChecked(False)
    self.ui.checkBox_Hf.setChecked(False)
    self.ui.checkBox_K.setChecked(False)
    self.ui.checkBox_He.setChecked(False)
    self.ui.checkBox_Md.setChecked(False)
    self.ui.checkBox_Mg.setChecked(False)
    self.ui.checkBox_Mo.setChecked(False)
    self.ui.checkBox_Mn.setChecked(False)
    self.ui.checkBox_O.setChecked(False)
    self.ui.checkBox_Mt.setChecked(False)
    self.ui.checkBox_S.setChecked(False)
    self.ui.checkBox_W.setChecked(False)
    self.ui.checkBox_Zn.setChecked(False)
    self.ui.checkBox_Eu.setChecked(False)
    self.ui.checkBox_Zr.setChecked(False)
    self.ui.checkBox_Er.setChecked(False)
    self.ui.checkBox_Ni.setChecked(False)
    self.ui.checkBox_No.setChecked(False)
    self.ui.checkBox_Na.setChecked(False)
    self.ui.checkBox_Nb.setChecked(False)
    self.ui.checkBox_Nd.setChecked(False)
    self.ui.checkBox_Ne.setChecked(False)
    self.ui.checkBox_Np.setChecked(False)
    self.ui.checkBox_Fr.setChecked(False)
    self.ui.checkBox_Fe.setChecked(False)
    self.ui.checkBox_Fl.setChecked(False)
    self.ui.checkBox_Fm.setChecked(False)
    self.ui.checkBox_B.setChecked(False)
    self.ui.checkBox_F.setChecked(False)
    self.ui.checkBox_Sr.setChecked(False)
    self.ui.checkBox_N.setChecked(False)
    self.ui.checkBox_Kr.setChecked(False)
    self.ui.checkBox_Si.setChecked(False)
    self.ui.checkBox_Sn.setChecked(False)
    self.ui.checkBox_Sm.setChecked(False)
    self.ui.checkBox_V.setChecked(False)
    self.ui.checkBox_Sc.setChecked(False)
    self.ui.checkBox_Sb.setChecked(False)
    self.ui.checkBox_Sg.setChecked(False)
    self.ui.checkBox_Se.setChecked(False)
    self.ui.checkBox_Co.setChecked(False)
    self.ui.checkBox_Cn.setChecked(False)
    self.ui.checkBox_Cm.setChecked(False)
    self.ui.checkBox_Cl.setChecked(False)
    self.ui.checkBox_Ca.setChecked(False)
    self.ui.checkBox_Cf.setChecked(False)
    self.ui.checkBox_Ce.setChecked(False)
    self.ui.checkBox_Xe.setChecked(False)
    self.ui.checkBox_Lu.setChecked(False)
    self.ui.checkBox_Cs.setChecked(False)
    self.ui.checkBox_Cr.setChecked(False)
    self.ui.checkBox_Cu.setChecked(False)
    self.ui.checkBox_La.setChecked(False)
    self.ui.checkBox_Li.setChecked(False)
    self.ui.checkBox_Lv.setChecked(False)
    self.ui.checkBox_Tl.setChecked(False)
    self.ui.checkBox_Tm.setChecked(False)
    self.ui.checkBox_Lr.setChecked(False)
    self.ui.checkBox_Th.setChecked(False)
    self.ui.checkBox_Ti.setChecked(False)
    self.ui.checkBox_Te.setChecked(False)
    self.ui.checkBox_Tb.setChecked(False)
    self.ui.checkBox_Tc.setChecked(False)
    self.ui.checkBox_Ta.setChecked(False)
    self.ui.checkBox_Yb.setChecked(False)
    self.ui.checkBox_Db.setChecked(False)
    self.ui.checkBox_Dy.setChecked(False)
    self.ui.checkBox_Ds.setChecked(False)
    self.ui.checkBox_I.setChecked(False)
    self.ui.checkBox_U.setChecked(False)
    self.ui.checkBox_Y.setChecked(False)
    self.ui.checkBox_Ac.setChecked(False)
    self.ui.checkBox_Ag.setChecked(False)
    self.ui.checkBox_Uut.setChecked(False)
    self.ui.checkBox_Ir.setChecked(False)
    self.ui.checkBox_Am.setChecked(False)
    self.ui.checkBox_Al.setChecked(False)
    self.ui.checkBox_As.setChecked(False)
    self.ui.checkBox_Ar.setChecked(False)
    self.ui.checkBox_Au.setChecked(False)
    self.ui.checkBox_At.setChecked(False)
    self.ui.checkBox_In.setChecked(False)
    
def populate(self, elem, mr_el):
    if self.n == 0:
        self.n = 1
        self.ui.label_el_1.setText(elem)
    elif  self.n == 1:
        self.n = 2
        self.ui.label_el_2.setText(elem)
    elif  self.n == 2:
        self.n = 3
        self.ui.label_el_3.setText(elem)
    elif  self.n == 3:
        self.n = 4
        self.ui.label_el_4.setText(elem)
    else:
        reset_(self)     

class Calculate(QMainWindow): 
    def __init__(self, parent = None):
        super(Calculate, self).__init__(parent)
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)
        reset_(self)
        self.elems__ = {"H":1.008, "He":4.0026, "Li":6.94, "Be":9.0122, "B":10.81, "C":12.011, "N":14.007, "O":15.999, "F":18.998,
                        "Ne":20.180, "Na":22.990, "Mg":24.305, "Al":26.982, "Si":28.085, "P":30.974, "S":32.06, "Cl":35.45,
                        "Ar":39.948, "K":39.098, "Ca":40.078, "Sc":44.956, "Ti":47.867, "V":50.942, "Cr":51.996, "Mn":54.938, "Fe":55.845,
                        "Co":58.933, "Ni":58.693, "Cu":63.546, "Zn":65.38, "Ga":69.723, "Ge":72.63, "As":74.922, "Se":78.96, "Br":79.904,
                        "Kr":83.798, "Rb":85.468, "Sr":87.62, "Y":88.906, "Zr":91.224, "Nb":92.906, "Mo":95.96, "Tc":97.91, "Ru":101.07, "Rh":102.91,
                        "Pd":106.42, "Ag":107.87, "Cd":112.41, "In":114.82, "Sn":118.71, "Sb":121.76, "Te":127.60, "I":126.90, "Xe":131.29,
                        "Cs":132.91, "Ba":137.33, "La":138.91, "Ce":140.12, "Pr":140.91, "Nd":144.24, "Pm":144.91, "Sm":150.36, "Eu":151.96,
                        "Gd":157.25, "Tb":158.93, "Dy":162.50, "Ho":164.93, "Er":167.26, "Tm":168.93, "Yb":173.05, "Lu":174.97, "Hf":178.49,
                        "Ta":180.95, "W":183.84, "Re":186.21, "Os":190.23, "Ir":192.22, "Pt":195.08, "Au":196.97, "Hg":200.59, "Tl":204.38,
                        "Pb":207.2, "Bi":208.98, "Po":208.98, "At":209.99, "Rn":222.02, "Fr":223.02, "Ra":226.03, "Ac":227.03, "Th":232.04, "Pa":231.04,
                        "U":238.03, "Np":237.05, "Pu":244.06, "Am":243.06, "Cm":247.07, "Bk":247.07, "Cf":251.08, "Es":252.08, "Fm":257.10, "Md":258.10,
                        "No":259.10, "Lr":262.11, "Rf":265.12, "Db":268.13, "Sg":271.13, "Bh":270, "Hs":277.15, "Mt":276.15, "Ds":281.16, "Rg":280.16,
                        "Cn":285.17, "Uut":284.18, "Fl":289.19, "Uup":288.19,"Lv":293, "Uus":294, "Uuo":294}

        
        # self.connect(self, Qt.SIGNAL('triggered()'), self.closeEvent
    
    
    @pyqtSlot()
    def on_pushButton_close_clicked(self):
        molecule = ""
        if self.ui.label_el_1.text() != "-":
            molecule += self.ui.label_el_1.text() + str(self.ui.spinBox_el_1.value())
        if self.ui.label_el_2.text() != "-":
            molecule += self.ui.label_el_2.text() + str(self.ui.spinBox_el_2.value())
        if self.ui.label_el_3.text() != "-":
            molecule += self.ui.label_el_3.text() + str(self.ui.spinBox_el_3.value())
        if self.ui.label_el_4.text() != "-":
            molecule += self.ui.label_el_4.text() + str(self.ui.spinBox_el_4.value())
        
        # elementList = Utility.molToElemList(molecule)
        print(molecule)
        self.close()
        
    @pyqtSlot()
    def on_pushButton_reset_clicked(self):
        reset_(self)

    @pyqtSlot()
    def on_checkBox_H_clicked(self):
        if self.ui.checkBox_H.isChecked():
            populate(self, "H", self.elems__["H"])    
        if not self.ui.checkBox_H.isChecked():
            reset_(self)
      
    @pyqtSlot()
    def on_checkBox_He_clicked(self):
        if self.ui.checkBox_He.isChecked():
            populate(self, "He", self.elems__["He"])
        if not self.ui.checkBox_He.isChecked():
            reset_(self)

    @pyqtSlot()
    def on_checkBox_Li_clicked(self):
        if self.ui.checkBox_Li.isChecked():
            populate(self, "Li", self.elems__["Li"])
        if not self.ui.checkBox_Li.isChecked():
            reset_(self)

    @pyqtSlot()
    def on_checkBox_Be_clicked(self):
        if self.ui.checkBox_Be.isChecked():
            populate(self, "Be", self.elems__["Be"])
        if not self.ui.checkBox_Be.isChecked():
            reset_(self)

    @pyqtSlot()
    def on_checkBox_B_clicked(self):
        if self.ui.checkBox_B.isChecked():
            populate(self, "B", self.elems__["B"])
        if not self.ui.checkBox_B.isChecked():
            reset_(self)

    @pyqtSlot()
    def on_checkBox_C_clicked(self):
        if self.ui.checkBox_C.isChecked():
            populate(self, "C", self.elems__["C"])    
        if not self.ui.checkBox_C.isChecked():
            reset_(self)
      
    @pyqtSlot()
    def on_checkBox_N_clicked(self):
        if self.ui.checkBox_N.isChecked():
            populate(self, "N", self.elems__["N"])
        if not self.ui.checkBox_N.isChecked():
            reset_(self)

    @pyqtSlot()
    def on_checkBox_O_clicked(self):
        if self.ui.checkBox_O.isChecked():
            populate(self, "O", self.elems__["O"])
        if not self.ui.checkBox_O.isChecked():
            reset_(self)

    @pyqtSlot()
    def on_checkBox_F_clicked(self):
        if self.ui.checkBox_F.isChecked():
            populate(self, "F", self.elems__["F"])
        if not self.ui.checkBox_F.isChecked():
            reset_(self)

    @pyqtSlot()
    def on_checkBox_Ne_clicked(self):
        if self.ui.checkBox_Ne.isChecked():
            populate(self, "Ne", self.elems__["Ne"])
        if not self.ui.checkBox_Ne.isChecked():
            reset_(self)

    @pyqtSlot()
    def on_checkBox_Na_clicked(self):
        if self.ui.checkBox_Na.isChecked():
            populate(self, "Na", self.elems__["Na"])    
        if not self.ui.checkBox_Na.isChecked():
            reset_(self)
      
    @pyqtSlot()
    def on_checkBox_Mg_clicked(self):
        if self.ui.checkBox_Mg.isChecked():
            populate(self, "Mg", self.elems__["Mg"])
        if not self.ui.checkBox_Mg.isChecked():
            reset_(self)

    @pyqtSlot()
    def on_checkBox_Al_clicked(self):
        if self.ui.checkBox_Al.isChecked():
            populate(self, "Al", self.elems__["Al"])
        if not self.ui.checkBox_Al.isChecked():
            reset_(self)

    @pyqtSlot()
    def on_checkBox_Si_clicked(self):
        if self.ui.checkBox_Si.isChecked():
            populate(self, "Si", self.elems__["Si"])
        if not self.ui.checkBox_Si.isChecked():
            reset_(self)

    @pyqtSlot()
    def on_checkBox_P_clicked(self):
        if self.ui.checkBox_P.isChecked():
            populate(self, "P", self.elems__["P"])
        if not self.ui.checkBox_P.isChecked():
            reset_(self)

    @pyqtSlot()
    def on_checkBox_S_clicked(self):
        if self.ui.checkBox_S.isChecked():
            populate(self, "S", self.elems__["S"])    
        if not self.ui.checkBox_S.isChecked():
            reset_(self)
      
    @pyqtSlot()
    def on_checkBox_Cl_clicked(self):
        if self.ui.checkBox_Cl.isChecked():
            populate(self, "Cl", self.elems__["Cl"])
        if not self.ui.checkBox_Cl.isChecked():
            reset_(self)

    @pyqtSlot()
    def on_checkBox_Ar_clicked(self):
        if self.ui.checkBox_Ar.isChecked():
            populate(self, "Ar", self.elems__["Ar"])
        if not self.ui.checkBox_Ar.isChecked():
            reset_(self)

    @pyqtSlot()
    def on_checkBox_K_clicked(self):
        if self.ui.checkBox_K.isChecked():
            populate(self, "K", self.elems__["K"])
        if not self.ui.checkBox_K.isChecked():
            reset_(self)

    @pyqtSlot()
    def on_checkBox_Ca_clicked(self):
        if self.ui.checkBox_Ca.isChecked():
            populate(self, "Ca", self.elems__["Ca"])
        if not self.ui.checkBox_Ca.isChecked():
            reset_(self)

    @pyqtSlot()
    def on_checkBox_Sc_clicked(self):
        if self.ui.checkBox_Sc.isChecked():
            populate(self, "Sc", self.elems__["Sc"])    
        if not self.ui.checkBox_Sc.isChecked():
            reset_(self)
      
    @pyqtSlot()
    def on_checkBox_Ti_clicked(self):
        if self.ui.checkBox_Ti.isChecked():
            populate(self, "Ti", self.elems__["Ti"])
        if not self.ui.checkBox_Ti.isChecked():
            reset_(self)

    @pyqtSlot()
    def on_checkBox_V_clicked(self):
        if self.ui.checkBox_V.isChecked():
            populate(self, "V", self.elems__["V"])
        if not self.ui.checkBox_V.isChecked():
            reset_(self)

    @pyqtSlot()
    def on_checkBox_Cr_clicked(self):
        if self.ui.checkBox_Cr.isChecked():
            populate(self, "Cr", self.elems__["Cr"])
        if not self.ui.checkBox_Cr.isChecked():
            reset_(self)

    @pyqtSlot()
    def on_checkBox_Mn_clicked(self):
        if self.ui.checkBox_Mn.isChecked():
            populate(self, "Mn", self.elems__["Mn"])
        if not self.ui.checkBox_Mn.isChecked():
            reset_(self)

    @pyqtSlot()
    def on_checkBox_Fe_clicked(self):
        if self.ui.checkBox_Fe.isChecked():
            populate(self, "Fe", self.elems__["Fe"])    
        if not self.ui.checkBox_Fe.isChecked():
            reset_(self)
      
    @pyqtSlot()
    def on_checkBox_Co_clicked(self):
        if self.ui.checkBox_Co.isChecked():
            populate(self, "Co", self.elems__["Co"])
        if not self.ui.checkBox_Co.isChecked():
            reset_(self)

    @pyqtSlot()
    def on_checkBox_Ni_clicked(self):
        if self.ui.checkBox_Ni.isChecked():
            populate(self, "Ni", self.elems__["Ni"])
        if not self.ui.checkBox_Ni.isChecked():
            reset_(self)

    @pyqtSlot()
    def on_checkBox_Cu_clicked(self):
        if self.ui.checkBox_Cu.isChecked():
            populate(self, "Cu", self.elems__["Cu"])
        if not self.ui.checkBox_Cu.isChecked():
            reset_(self)

    @pyqtSlot()
    def on_checkBox_Zn_clicked(self):
        if self.ui.checkBox_Zn.isChecked():
            populate(self, "Zn", self.elems__["Zn"])
        if not self.ui.checkBox_Zn.isChecked():
            reset_(self)

    @pyqtSlot()
    def on_checkBox_Ga_clicked(self):
        if self.ui.checkBox_Ga.isChecked():
            populate(self, "Ga", self.elems__["Ga"])    
        if not self.ui.checkBox_Ga.isChecked():
            reset_(self)
      
    @pyqtSlot()
    def on_checkBox_Ge_clicked(self):
        if self.ui.checkBox_Ge.isChecked():
            populate(self, "Ge", self.elems__["Ge"])
        if not self.ui.checkBox_Ge.isChecked():
            reset_(self)

    @pyqtSlot()
    def on_checkBox_As_clicked(self):
        if self.ui.checkBox_As.isChecked():
            populate(self, "As", self.elems__["As"])
        if not self.ui.checkBox_As.isChecked():
            reset_(self)

    @pyqtSlot()
    def on_checkBox_Se_clicked(self):
        if self.ui.checkBox_Se.isChecked():
            populate(self, "Se", self.elems__["Se"])
        if not self.ui.checkBox_Se.isChecked():
            reset_(self)

    @pyqtSlot()
    def on_checkBox_Br_clicked(self):
        if self.ui.checkBox_Br.isChecked():
            populate(self, "Br", self.elems__["Br"])
        if not self.ui.checkBox_Br.isChecked():
            reset_(self)

    @pyqtSlot()
    def on_checkBox_Kr_clicked(self):
        if self.ui.checkBox_Kr.isChecked():
            populate(self, "Kr", self.elems__["Kr"])    
        if not self.ui.checkBox_Kr.isChecked():
            reset_(self)
      
    @pyqtSlot()
    def on_checkBox_Rb_clicked(self):
        if self.ui.checkBox_Rb.isChecked():
            populate(self, "Rb", self.elems__["Rb"])
        if not self.ui.checkBox_Rb.isChecked():
            reset_(self)

    @pyqtSlot()
    def on_checkBox_Sr_clicked(self):
        if self.ui.checkBox_Sr.isChecked():
            populate(self, "Sr", self.elems__["Sr"])
        if not self.ui.checkBox_Sr.isChecked():
            reset_(self)

    @pyqtSlot()
    def on_checkBox_Y_clicked(self):
        if self.ui.checkBox_Y.isChecked():
            populate(self, "Y", self.elems__["Y"])
        if not self.ui.checkBox_Y.isChecked():
            reset_(self)

    @pyqtSlot()
    def on_checkBox_Zr_clicked(self):
        if self.ui.checkBox_Zr.isChecked():
            populate(self, "Zr", self.elems__["Zr"])
        if not self.ui.checkBox_Zr.isChecked():
            reset_(self)

    @pyqtSlot()
    def on_checkBox_Nb_clicked(self):
        if self.ui.checkBox_Nb.isChecked():
            populate(self, "Nb", self.elems__["Nb"])    
        if not self.ui.checkBox_Nb.isChecked():
            reset_(self)
      
    @pyqtSlot()
    def on_checkBox_Mo_clicked(self):
        if self.ui.checkBox_Mo.isChecked():
            populate(self, "Mo", self.elems__["Mo"])
        if not self.ui.checkBox_Mo.isChecked():
            reset_(self)

    @pyqtSlot()
    def on_checkBox_Tc_clicked(self):
        if self.ui.checkBox_Tc.isChecked():
            populate(self, "Tc", self.elems__["Tc"])
        if not self.ui.checkBox_Tc.isChecked():
            reset_(self)

    @pyqtSlot()
    def on_checkBox_Ru_clicked(self):
        if self.ui.checkBox_Ru.isChecked():
            populate(self, "Ru", self.elems__["Ru"])
        if not self.ui.checkBox_Ru.isChecked():
            reset_(self)

    @pyqtSlot()
    def on_checkBox_Rh_clicked(self):
        if self.ui.checkBox_Rh.isChecked():
            populate(self, "Rh", self.elems__["Rh"])
        if not self.ui.checkBox_Rh.isChecked():
            reset_(self)

    @pyqtSlot()
    def on_checkBox_Pd_clicked(self):
        if self.ui.checkBox_Pd.isChecked():
            populate(self, "Pd", self.elems__["Pd"])    
        if not self.ui.checkBox_Pd.isChecked():
            reset_(self)
      
    @pyqtSlot()
    def on_checkBox_Ag_clicked(self):
        if self.ui.checkBox_Ag.isChecked():
            populate(self, "Ag", self.elems__["Ag"])
        if not self.ui.checkBox_Ag.isChecked():
            reset_(self)

    @pyqtSlot()
    def on_checkBox_Cd_clicked(self):
        if self.ui.checkBox_Cd.isChecked():
            populate(self, "Cd", self.elems__["Cd"])
        if not self.ui.checkBox_Cd.isChecked():
            reset_(self)

    @pyqtSlot()
    def on_checkBox_In_clicked(self):
        if self.ui.checkBox_In.isChecked():
            populate(self, "In", self.elems__["In"])
        if not self.ui.checkBox_In.isChecked():
            reset_(self)

    @pyqtSlot()
    def on_checkBox_Sn_clicked(self):
        if self.ui.checkBox_Sn.isChecked():
            populate(self, "Sn", self.elems__["Sn"])
        if not self.ui.checkBox_Sn.isChecked():
            reset_(self)

    @pyqtSlot()
    def on_checkBox_Sb_clicked(self):
        if self.ui.checkBox_Sb.isChecked():
            populate(self, "Sb", self.elems__["Sb"])    
        if not self.ui.checkBox_Sb.isChecked():
            reset_(self)
      
    @pyqtSlot()
    def on_checkBox_Te_clicked(self):
        if self.ui.checkBox_Te.isChecked():
            populate(self, "Te", self.elems__["Te"])
        if not self.ui.checkBox_Te.isChecked():
            reset_(self)

    @pyqtSlot()
    def on_checkBox_I_clicked(self):
        if self.ui.checkBox_I.isChecked():
            populate(self, "I", self.elems__["I"])
        if not self.ui.checkBox_I.isChecked():
            reset_(self)

    @pyqtSlot()
    def on_checkBox_Xe_clicked(self):
        if self.ui.checkBox_Xe.isChecked():
            populate(self, "Xe", self.elems__["Xe"])
        if not self.ui.checkBox_Xe.isChecked():
            reset_(self)

    @pyqtSlot()
    def on_checkBox_Cs_clicked(self):
        if self.ui.checkBox_Cs.isChecked():
            populate(self, "Cs", self.elems__["Cs"])
        if not self.ui.checkBox_Cs.isChecked():
            reset_(self)

    @pyqtSlot()
    def on_checkBox_Ba_clicked(self):
        if self.ui.checkBox_Ba.isChecked():
            populate(self, "Ba", self.elems__["Ba"])    
        if not self.ui.checkBox_Ba.isChecked():
            reset_(self)
      
    @pyqtSlot()
    def on_checkBox_La_clicked(self):
        if self.ui.checkBox_La.isChecked():
            populate(self, "La", self.elems__["La"])
        if not self.ui.checkBox_La.isChecked():
            reset_(self)

    @pyqtSlot()
    def on_checkBox_Ce_clicked(self):
        if self.ui.checkBox_Ce.isChecked():
            populate(self, "Ce", self.elems__["Ce"])
        if not self.ui.checkBox_Ce.isChecked():
            reset_(self)

    @pyqtSlot()
    def on_checkBox_Pr_clicked(self):
        if self.ui.checkBox_Pr.isChecked():
            populate(self, "Pr", self.elems__["Pr"])
        if not self.ui.checkBox_Pr.isChecked():
            reset_(self)

    @pyqtSlot()
    def on_checkBox_Nd_clicked(self):
        if self.ui.checkBox_Nd.isChecked():
            populate(self, "Nd", self.elems__["Nd"])
        if not self.ui.checkBox_Nd.isChecked():
            reset_(self)

    @pyqtSlot()
    def on_checkBox_Pm_clicked(self):
        if self.ui.checkBox_Pm.isChecked():
            populate(self, "Pm", self.elems__["Pm"])    
        if not self.ui.checkBox_Pm.isChecked():
            reset_(self)
      
    @pyqtSlot()
    def on_checkBox_Sm_clicked(self):
        if self.ui.checkBox_Sm.isChecked():
            populate(self, "Sm", self.elems__["Sm"])
        if not self.ui.checkBox_Sm.isChecked():
            reset_(self)

    @pyqtSlot()
    def on_checkBox_Eu_clicked(self):
        if self.ui.checkBox_Eu.isChecked():
            populate(self, "Eu", self.elems__["Eu"])
        if not self.ui.checkBox_Eu.isChecked():
            reset_(self)

    @pyqtSlot()
    def on_checkBox_Gd_clicked(self):
        if self.ui.checkBox_Gd.isChecked():
            populate(self, "Gd", self.elems__["Gd"])
        if not self.ui.checkBox_Gd.isChecked():
            reset_(self)

    @pyqtSlot()
    def on_checkBox_Tb_clicked(self):
        if self.ui.checkBox_Tb.isChecked():
            populate(self, "Tb", self.elems__["Tb"])
        if not self.ui.checkBox_Tb.isChecked():
            reset_(self)

    @pyqtSlot()
    def on_checkBox_Dy_clicked(self):
        if self.ui.checkBox_Dy.isChecked():
            populate(self, "Dy", self.elems__["Dy"])    
        if not self.ui.checkBox_Dy.isChecked():
            reset_(self)
      
    @pyqtSlot()
    def on_checkBox_Ho_clicked(self):
        if self.ui.checkBox_Ho.isChecked():
            populate(self, "Ho", self.elems__["Ho"])
        if not self.ui.checkBox_Ho.isChecked():
            reset_(self)

    @pyqtSlot()
    def on_checkBox_Er_clicked(self):
        if self.ui.checkBox_Er.isChecked():
            populate(self, "Er", self.elems__["Er"])
        if not self.ui.checkBox_Er.isChecked():
            reset_(self)

    @pyqtSlot()
    def on_checkBox_Tm_clicked(self):
        if self.ui.checkBox_Tm.isChecked():
            populate(self, "Tm", self.elems__["Tm"])
        if not self.ui.checkBox_Tm.isChecked():
            reset_(self)

    @pyqtSlot()
    def on_checkBox_Yb_clicked(self):
        if self.ui.checkBox_Yb.isChecked():
            populate(self, "Yb", self.elems__["Yb"])
        if not self.ui.checkBox_Yb.isChecked():
            reset_(self)

    @pyqtSlot()
    def on_checkBox_Lu_clicked(self):
        if self.ui.checkBox_Lu.isChecked():
            populate(self, "Lu", self.elems__["Lu"])    
        if not self.ui.checkBox_Lu.isChecked():
            reset_(self)
      
    @pyqtSlot()
    def on_checkBox_Hf_clicked(self):
        if self.ui.checkBox_Hf.isChecked():
            populate(self, "Hf", self.elems__["Hf"])
        if not self.ui.checkBox_Hf.isChecked():
            reset_(self)

    @pyqtSlot()
    def on_checkBox_Ta_clicked(self):
        if self.ui.checkBox_Ta.isChecked():
            populate(self, "Ta", self.elems__["Ta"])
        if not self.ui.checkBox_Ta.isChecked():
            reset_(self)

    @pyqtSlot()
    def on_checkBox_W_clicked(self):
        if self.ui.checkBox_W.isChecked():
            populate(self, "W", self.elems__["W"])
        if not self.ui.checkBox_W.isChecked():
            reset_(self)

    @pyqtSlot()
    def on_checkBox_Re_clicked(self):
        if self.ui.checkBox_Re.isChecked():
            populate(self, "Re", self.elems__["Re"])
        if not self.ui.checkBox_Re.isChecked():
            reset_(self)

    @pyqtSlot()
    def on_checkBox_Os_clicked(self):
        if self.ui.checkBox_Os.isChecked():
            populate(self, "Os", self.elems__["Os"])    
        if not self.ui.checkBox_Os.isChecked():
            reset_(self)
      
    @pyqtSlot()
    def on_checkBox_Ir_clicked(self):
        if self.ui.checkBox_Ir.isChecked():
            populate(self, "Ir", self.elems__["Ir"])
        if not self.ui.checkBox_Ir.isChecked():
            reset_(self)

    @pyqtSlot()
    def on_checkBox_Pt_clicked(self):
        if self.ui.checkBox_Pt.isChecked():
            populate(self, "Pt", self.elems__["Pt"])
        if not self.ui.checkBox_Pt.isChecked():
            reset_(self)

    @pyqtSlot()
    def on_checkBox_Au_clicked(self):
        if self.ui.checkBox_Au.isChecked():
            populate(self, "Au", self.elems__["Au"])
        if not self.ui.checkBox_Au.isChecked():
            reset_(self)

    @pyqtSlot()
    def on_checkBox_Hg_clicked(self):
        if self.ui.checkBox_Hg.isChecked():
            populate(self, "Hg", self.elems__["Hg"])
        if not self.ui.checkBox_Hg.isChecked():
            reset_(self)

    @pyqtSlot()
    def on_checkBox_Tl_clicked(self):
        if self.ui.checkBox_Tl.isChecked():
            populate(self, "Tl", self.elems__["Tl"])
        if not self.ui.checkBox_Tl.isChecked():
            reset_(self)

    @pyqtSlot()
    def on_checkBox_Pb_clicked(self):
        if self.ui.checkBox_Pb.isChecked():
            populate(self, "Pb", self.elems__["Pb"])
        if not self.ui.checkBox_Pb.isChecked():
            reset_(self)

    @pyqtSlot()
    def on_checkBox_Bi_clicked(self):
        if self.ui.checkBox_Bi.isChecked():
            populate(self, "Bi", self.elems__["Bi"])
        if not self.ui.checkBox_Bi.isChecked():
            reset_(self)

    @pyqtSlot()
    def on_checkBox_Po_clicked(self):
        if self.ui.checkBox_Po.isChecked():
            populate(self, "Po", self.elems__["Po"])    
        if not self.ui.checkBox_Po.isChecked():
            reset_(self)
      
    @pyqtSlot()
    def on_checkBox_At_clicked(self):
        if self.ui.checkBox_At.isChecked():
            populate(self, "At", self.elems__["At"])
        if not self.ui.checkBox_At.isChecked():
            reset_(self)

    @pyqtSlot()
    def on_checkBox_Rn_clicked(self):
        if self.ui.checkBox_Rn.isChecked():
            populate(self, "Rn", self.elems__["Rn"])
        if not self.ui.checkBox_Rn.isChecked():
            reset_(self)

    @pyqtSlot()
    def on_checkBox_Fr_clicked(self):
        if self.ui.checkBox_Fr.isChecked():
            populate(self, "Fr", self.elems__["Fr"])
        if not self.ui.checkBox_Fr.isChecked():
            reset_(self)

    @pyqtSlot()
    def on_checkBox_Ra_clicked(self):
        if self.ui.checkBox_Ra.isChecked():
            populate(self, "Ra", self.elems__["Ra"])
        if not self.ui.checkBox_Ra.isChecked():
            reset_(self)

    @pyqtSlot()
    def on_checkBox_Ac_clicked(self):
        if self.ui.checkBox_Ac.isChecked():
            populate(self, "Ac", self.elems__["Ac"])    
        if not self.ui.checkBox_Ac.isChecked():
            reset_(self)
      
    @pyqtSlot()
    def on_checkBox_Th_clicked(self):
        if self.ui.checkBox_Th.isChecked():
            populate(self, "Th", self.elems__["Th"])
        if not self.ui.checkBox_Th.isChecked():
            reset_(self)

    @pyqtSlot()
    def on_checkBox_Pa_clicked(self):
        if self.ui.checkBox_Pa.isChecked():
            populate(self, "Pa", self.elems__["Pa"])
        if not self.ui.checkBox_Pa.isChecked():
            reset_(self)

    @pyqtSlot()
    def on_checkBox_U_clicked(self):
        if self.ui.checkBox_U.isChecked():
            populate(self, "U", self.elems__["U"])
        if not self.ui.checkBox_U.isChecked():
            reset_(self)

    @pyqtSlot()
    def on_checkBox_Np_clicked(self):
        if self.ui.checkBox_Np.isChecked():
            populate(self, "Np", self.elems__["Np"])
        if not self.ui.checkBox_Np.isChecked():
            reset_(self)

    @pyqtSlot()
    def on_checkBox_Pu_clicked(self):
        if self.ui.checkBox_Pu.isChecked():
            populate(self, "Pu", self.elems__["Pu"])    
        if not self.ui.checkBox_Pu.isChecked():
            reset_(self)
      
    @pyqtSlot()
    def on_checkBox_Am_clicked(self):
        if self.ui.checkBox_Am.isChecked():
            populate(self, "Am", self.elems__["Am"])
        if not self.ui.checkBox_Am.isChecked():
            reset_(self)

    @pyqtSlot()
    def on_checkBox_Cm_clicked(self):
        if self.ui.checkBox_Cm.isChecked():
            populate(self, "Cm", self.elems__["Cm"])
        if not self.ui.checkBox_Cm.isChecked():
            reset_(self)

    @pyqtSlot()
    def on_checkBox_Bk_clicked(self):
        if self.ui.checkBox_Bk.isChecked():
            populate(self, "Bk", self.elems__["Bk"])
        if not self.ui.checkBox_Bk.isChecked():
            reset_(self)

    @pyqtSlot()
    def on_checkBox_Cf_clicked(self):
        if self.ui.checkBox_Cf.isChecked():
            populate(self, "Cf", self.elems__["Cf"])
        if not self.ui.checkBox_Cf.isChecked():
            reset_(self)

    @pyqtSlot()
    def on_checkBox_Es_clicked(self):
        if self.ui.checkBox_Es.isChecked():
            populate(self, "Es", self.elems__["Es"])    
        if not self.ui.checkBox_Es.isChecked():
            reset_(self)
      
    @pyqtSlot()
    def on_checkBox_Fm_clicked(self):
        if self.ui.checkBox_Fm.isChecked():
            populate(self, "Fm", self.elems__["Fm"])
        if not self.ui.checkBox_Fm.isChecked():
            reset_(self)

    @pyqtSlot()
    def on_checkBox_Md_clicked(self):
        if self.ui.checkBox_Md.isChecked():
            populate(self, "Md", self.elems__["Md"])
        if not self.ui.checkBox_Md.isChecked():
            reset_(self)

    @pyqtSlot()
    def on_checkBox_No_clicked(self):
        if self.ui.checkBox_No.isChecked():
            populate(self, "No", self.elems__["No"])
        if not self.ui.checkBox_No.isChecked():
            reset_(self)

    @pyqtSlot()
    def on_checkBox_Lr_clicked(self):
        if self.ui.checkBox_Lr.isChecked():
            populate(self, "Lr", self.elems__["Lr"])
        if not self.ui.checkBox_Lr.isChecked():
            reset_(self)

    @pyqtSlot()
    def on_checkBox_Rf_clicked(self):
        if self.ui.checkBox_Rf.isChecked():
            populate(self, "Rf", self.elems__["Rf"])    
        if not self.ui.checkBox_Rf.isChecked():
            reset_(self)
      
    @pyqtSlot()
    def on_checkBox_Db_clicked(self):
        if self.ui.checkBox_Db.isChecked():
            populate(self, "Db", self.elems__["Db"])
        if not self.ui.checkBox_Db.isChecked():
            reset_(self)

    @pyqtSlot()
    def on_checkBox_Sg_clicked(self):
        if self.ui.checkBox_Sg.isChecked():
            populate(self, "Sg", self.elems__["Sg"])
        if not self.ui.checkBox_Sg.isChecked():
            reset_(self)

    @pyqtSlot()
    def on_checkBox_Bh_clicked(self):
        if self.ui.checkBox_Bh.isChecked():
            populate(self, "Bh", self.elems__["Bh"])
        if not self.ui.checkBox_Bh.isChecked():
            reset_(self)

    @pyqtSlot()
    def on_checkBox_Hs_clicked(self):
        if self.ui.checkBox_Hs.isChecked():
            populate(self, "Hs", self.elems__["Hs"])
        if not self.ui.checkBox_Hs.isChecked():
            reset_(self)

    @pyqtSlot()
    def on_checkBox_Mt_clicked(self):
        if self.ui.checkBox_Mt.isChecked():
            populate(self, "Mt", self.elems__["Mt"])    
        if not self.ui.checkBox_Mt.isChecked():
            reset_(self)
      
    @pyqtSlot()
    def on_checkBox_Ds_clicked(self):
        if self.ui.checkBox_Ds.isChecked():
            populate(self, "Ds", self.elems__["Ds"])
        if not self.ui.checkBox_Ds.isChecked():
            reset_(self)

    @pyqtSlot()
    def on_checkBox_Rg_clicked(self):
        if self.ui.checkBox_Rg.isChecked():
            populate(self, "Rg", self.elems__["Rg"])
        if not self.ui.checkBox_Rg.isChecked():
            reset_(self)

    @pyqtSlot()
    def on_checkBox_Cn_clicked(self):
        if self.ui.checkBox_Cn.isChecked():
            populate(self, "Cn", self.elems__["Cn"])
        if not self.ui.checkBox_Cn.isChecked():
            reset_(self)

    @pyqtSlot()
    def on_checkBox_Uut_clicked(self):
        if self.ui.checkBox_Uut.isChecked():
            populate(self, "Uut", self.elems__["Uut"])
        if not self.ui.checkBox_Uut.isChecked():
            reset_(self)

    @pyqtSlot()
    def on_checkBox_Fl_clicked(self):
        if self.ui.checkBox_Fl.isChecked():
            populate(self, "Fl", self.elems__["Fl"])    
        if not self.ui.checkBox_Fl.isChecked():
            reset_(self)
      
    @pyqtSlot()
    def on_checkBox_Uup_clicked(self):
        if self.ui.checkBox_Uup.isChecked():
            populate(self, "Uup", self.elems__["Uup"])
        if not self.ui.checkBox_Uup.isChecked():
            reset_(self)

    @pyqtSlot()
    def on_checkBox_Lv_clicked(self):
        if self.ui.checkBox_Lv.isChecked():
            populate(self, "Lv", self.elems__["Lv"])
        if not self.ui.checkBox_Lv.isChecked():
            reset_(self)

    @pyqtSlot()
    def on_checkBox_Uus_clicked(self):
        if self.ui.checkBox_Uus.isChecked():
            populate(self, "Uus", self.elems__["Uus"])
        if not self.ui.checkBox_Uus.isChecked():
            reset_(self)

    @pyqtSlot()
    def on_checkBox_Uuo_clicked(self):
        if self.ui.checkBox_Uuo.isChecked():
            populate(self, "Uuo", self.elems__["Uuo"])
        if not self.ui.checkBox_Uuo.isChecked():
            reset_(self)
   
    # @pyqtSlot()
    # def on_pushButton_save_clicked(self):
        # f_name = self.ui.lineEdit_html.text() + ".html"
        # file = open(f_name,"a")
        # file.write(self.ui.textBrowser.toHtml())    
 

if __name__ == "__main__":
    import sys
    app = QApplication(sys.argv)
    form = Calculate()
    form.show()
    app.exec_()
