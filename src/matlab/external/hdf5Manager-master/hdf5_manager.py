# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'hdf5_manager.ui'
#
# Created by: PyQt4 UI code generator 4.11.4
#
# WARNING! All changes made in this file will be lost!

from PyQt4 import QtCore, QtGui

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    def _fromUtf8(s):
        return s

try:
    _encoding = QtGui.QApplication.UnicodeUTF8
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig, _encoding)
except AttributeError:
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig)

class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName(_fromUtf8("MainWindow"))
        MainWindow.resize(1091, 866)
        MainWindow.setMinimumSize(QtCore.QSize(500, 700))
        MainWindow.setMaximumSize(QtCore.QSize(16777215, 16777215))
        MainWindow.setStyleSheet(_fromUtf8("background-color: rgb(226, 226, 226);"))
        self.centralwidget = QtGui.QWidget(MainWindow)
        self.centralwidget.setObjectName(_fromUtf8("centralwidget"))
        self.gridLayout_3 = QtGui.QGridLayout(self.centralwidget)
        self.gridLayout_3.setObjectName(_fromUtf8("gridLayout_3"))
        self.gridFrame = QtGui.QFrame(self.centralwidget)
        self.gridFrame.setMinimumSize(QtCore.QSize(0, 10))
        self.gridFrame.setFrameShape(QtGui.QFrame.StyledPanel)
        self.gridFrame.setObjectName(_fromUtf8("gridFrame"))
        self.gridLayout = QtGui.QGridLayout(self.gridFrame)
        self.gridLayout.setMargin(6)
        self.gridLayout.setObjectName(_fromUtf8("gridLayout"))
        self.horizontalWidget = QtGui.QWidget(self.gridFrame)
        self.horizontalWidget.setObjectName(_fromUtf8("horizontalWidget"))
        self.horizontalLayout_2 = QtGui.QHBoxLayout(self.horizontalWidget)
        self.horizontalLayout_2.setMargin(0)
        self.horizontalLayout_2.setObjectName(_fromUtf8("horizontalLayout_2"))
        self.loadDirectoryBtn = QtGui.QPushButton(self.horizontalWidget)
        self.loadDirectoryBtn.setStyleSheet(_fromUtf8("background-color: rgb(240, 240, 240);"))
        self.loadDirectoryBtn.setObjectName(_fromUtf8("loadDirectoryBtn"))
        self.horizontalLayout_2.addWidget(self.loadDirectoryBtn)
        self.reloadDirectoryBtn = QtGui.QPushButton(self.horizontalWidget)
        self.reloadDirectoryBtn.setStyleSheet(_fromUtf8("background-color: rgb(240, 240, 240);"))
        self.reloadDirectoryBtn.setObjectName(_fromUtf8("reloadDirectoryBtn"))
        self.horizontalLayout_2.addWidget(self.reloadDirectoryBtn)
        self.workingDirectory = QtGui.QLineEdit(self.horizontalWidget)
        self.workingDirectory.setStyleSheet(_fromUtf8("background-color: rgb(255, 255, 255);"))
        self.workingDirectory.setObjectName(_fromUtf8("workingDirectory"))
        self.horizontalLayout_2.addWidget(self.workingDirectory)
        self.gridLayout.addWidget(self.horizontalWidget, 0, 0, 1, 1)
        self.gridLayout_3.addWidget(self.gridFrame, 0, 0, 1, 1)
        self.gridFrame1 = QtGui.QFrame(self.centralwidget)
        self.gridFrame1.setMinimumSize(QtCore.QSize(0, 20))
        self.gridFrame1.setFrameShape(QtGui.QFrame.NoFrame)
        self.gridFrame1.setObjectName(_fromUtf8("gridFrame1"))
        self.gridLayout_7 = QtGui.QGridLayout(self.gridFrame1)
        self.gridLayout_7.setMargin(0)
        self.gridLayout_7.setObjectName(_fromUtf8("gridLayout_7"))
        self.horizontalWidget1 = QtGui.QWidget(self.gridFrame1)
        self.horizontalWidget1.setObjectName(_fromUtf8("horizontalWidget1"))
        self.horizontalLayout_10 = QtGui.QHBoxLayout(self.horizontalWidget1)
        self.horizontalLayout_10.setSizeConstraint(QtGui.QLayout.SetNoConstraint)
        self.horizontalLayout_10.setMargin(0)
        self.horizontalLayout_10.setObjectName(_fromUtf8("horizontalLayout_10"))
        self.gridFrame2 = QtGui.QFrame(self.horizontalWidget1)
        self.gridFrame2.setBaseSize(QtCore.QSize(0, 0))
        self.gridFrame2.setFrameShape(QtGui.QFrame.StyledPanel)
        self.gridFrame2.setLineWidth(1)
        self.gridFrame2.setObjectName(_fromUtf8("gridFrame2"))
        self.gridLayout_6 = QtGui.QGridLayout(self.gridFrame2)
        self.gridLayout_6.setMargin(6)
        self.gridLayout_6.setSpacing(2)
        self.gridLayout_6.setObjectName(_fromUtf8("gridLayout_6"))
        self.label_5 = QtGui.QLabel(self.gridFrame2)
        self.label_5.setMaximumSize(QtCore.QSize(16777215, 20))
        font = QtGui.QFont()
        font.setBold(True)
        font.setWeight(75)
        self.label_5.setFont(font)
        self.label_5.setObjectName(_fromUtf8("label_5"))
        self.gridLayout_6.addWidget(self.label_5, 0, 0, 1, 1)
        self.horizontalWidget2 = QtGui.QWidget(self.gridFrame2)
        self.horizontalWidget2.setObjectName(_fromUtf8("horizontalWidget2"))
        self.horizontalLayout_5 = QtGui.QHBoxLayout(self.horizontalWidget2)
        self.horizontalLayout_5.setContentsMargins(5, 3, 5, 3)
        self.horizontalLayout_5.setSpacing(5)
        self.horizontalLayout_5.setObjectName(_fromUtf8("horizontalLayout_5"))
        self.genericRadioBtn = QtGui.QRadioButton(self.horizontalWidget2)
        self.genericRadioBtn.setObjectName(_fromUtf8("genericRadioBtn"))
        self.plottingChoice = QtGui.QButtonGroup(MainWindow)
        self.plottingChoice.setObjectName(_fromUtf8("plottingChoice"))
        self.plottingChoice.addButton(self.genericRadioBtn)
        self.horizontalLayout_5.addWidget(self.genericRadioBtn)
        self.timeSeriesRadioBtn = QtGui.QRadioButton(self.horizontalWidget2)
        self.timeSeriesRadioBtn.setObjectName(_fromUtf8("timeSeriesRadioBtn"))
        self.plottingChoice.addButton(self.timeSeriesRadioBtn)
        self.horizontalLayout_5.addWidget(self.timeSeriesRadioBtn)
        self.spikesRadioBtn = QtGui.QRadioButton(self.horizontalWidget2)
        self.spikesRadioBtn.setObjectName(_fromUtf8("spikesRadioBtn"))
        self.plottingChoice.addButton(self.spikesRadioBtn)
        self.horizontalLayout_5.addWidget(self.spikesRadioBtn)
        self.threeDRadioBtn = QtGui.QRadioButton(self.horizontalWidget2)
        self.threeDRadioBtn.setObjectName(_fromUtf8("threeDRadioBtn"))
        self.plottingChoice.addButton(self.threeDRadioBtn)
        self.horizontalLayout_5.addWidget(self.threeDRadioBtn)
        self.ImageStackRadioBtn = QtGui.QRadioButton(self.horizontalWidget2)
        self.ImageStackRadioBtn.setObjectName(_fromUtf8("ImageStackRadioBtn"))
        self.plottingChoice.addButton(self.ImageStackRadioBtn)
        self.horizontalLayout_5.addWidget(self.ImageStackRadioBtn)
        self.gridLayout_6.addWidget(self.horizontalWidget2, 1, 0, 1, 1)
        self.horizontalWidget3 = QtGui.QWidget(self.gridFrame2)
        self.horizontalWidget3.setObjectName(_fromUtf8("horizontalWidget3"))
        self.horizontalLayout_3 = QtGui.QHBoxLayout(self.horizontalWidget3)
        self.horizontalLayout_3.setMargin(0)
        self.horizontalLayout_3.setObjectName(_fromUtf8("horizontalLayout_3"))
        self.dataSetSelection = QtGui.QLineEdit(self.horizontalWidget3)
        self.dataSetSelection.setMaximumSize(QtCore.QSize(200, 16777215))
        self.dataSetSelection.setStyleSheet(_fromUtf8("background-color: rgb(255, 255, 255);"))
        self.dataSetSelection.setObjectName(_fromUtf8("dataSetSelection"))
        self.horizontalLayout_3.addWidget(self.dataSetSelection)
        self.plotDataBtn = QtGui.QPushButton(self.horizontalWidget3)
        self.plotDataBtn.setStyleSheet(_fromUtf8("background-color: rgb(240, 240, 240);"))
        self.plotDataBtn.setObjectName(_fromUtf8("plotDataBtn"))
        self.horizontalLayout_3.addWidget(self.plotDataBtn)
        self.addToPlotBtn = QtGui.QPushButton(self.horizontalWidget3)
        self.addToPlotBtn.setStyleSheet(_fromUtf8("background-color: rgb(240, 240, 240);"))
        self.addToPlotBtn.setObjectName(_fromUtf8("addToPlotBtn"))
        self.horizontalLayout_3.addWidget(self.addToPlotBtn)
        self.gridLayout_6.addWidget(self.horizontalWidget3, 3, 0, 1, 1)
        self.horizontalLayout_10.addWidget(self.gridFrame2)
        self.gridFrame3 = QtGui.QFrame(self.horizontalWidget1)
        self.gridFrame3.setFrameShape(QtGui.QFrame.StyledPanel)
        self.gridFrame3.setObjectName(_fromUtf8("gridFrame3"))
        self.gridLayout_5 = QtGui.QGridLayout(self.gridFrame3)
        self.gridLayout_5.setMargin(6)
        self.gridLayout_5.setObjectName(_fromUtf8("gridLayout_5"))
        self.horizontalWidget4 = QtGui.QWidget(self.gridFrame3)
        self.horizontalWidget4.setObjectName(_fromUtf8("horizontalWidget4"))
        self.horizontalLayout_6 = QtGui.QHBoxLayout(self.horizontalWidget4)
        self.horizontalLayout_6.setMargin(0)
        self.horizontalLayout_6.setObjectName(_fromUtf8("horizontalLayout_6"))
        self.launchiPythonBtn = QtGui.QPushButton(self.horizontalWidget4)
        self.launchiPythonBtn.setStyleSheet(_fromUtf8("background-color: rgb(240, 240, 240);"))
        self.launchiPythonBtn.setObjectName(_fromUtf8("launchiPythonBtn"))
        self.horizontalLayout_6.addWidget(self.launchiPythonBtn)
        self.sendToConsoleBtn = QtGui.QPushButton(self.horizontalWidget4)
        self.sendToConsoleBtn.setStyleSheet(_fromUtf8("background-color: rgb(240, 240, 240);"))
        self.sendToConsoleBtn.setObjectName(_fromUtf8("sendToConsoleBtn"))
        self.horizontalLayout_6.addWidget(self.sendToConsoleBtn)
        self.plotNameSpaceBtn = QtGui.QPushButton(self.horizontalWidget4)
        self.plotNameSpaceBtn.setStyleSheet(_fromUtf8("background-color: rgb(240, 240, 240);"))
        self.plotNameSpaceBtn.setObjectName(_fromUtf8("plotNameSpaceBtn"))
        self.horizontalLayout_6.addWidget(self.plotNameSpaceBtn)
        self.gridLayout_5.addWidget(self.horizontalWidget4, 2, 0, 1, 1)
        self.label_6 = QtGui.QLabel(self.gridFrame3)
        font = QtGui.QFont()
        font.setBold(True)
        font.setWeight(75)
        self.label_6.setFont(font)
        self.label_6.setObjectName(_fromUtf8("label_6"))
        self.gridLayout_5.addWidget(self.label_6, 0, 0, 1, 1)
        spacerItem = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.gridLayout_5.addItem(spacerItem, 1, 0, 1, 1)
        self.horizontalLayout_10.addWidget(self.gridFrame3)
        self.gridLayout_7.addWidget(self.horizontalWidget1, 0, 0, 1, 1)
        self.gridLayout_3.addWidget(self.gridFrame1, 3, 0, 1, 2)
        self.gridFrame4 = QtGui.QFrame(self.centralwidget)
        self.gridFrame4.setMinimumSize(QtCore.QSize(0, 10))
        self.gridFrame4.setFrameShape(QtGui.QFrame.StyledPanel)
        self.gridFrame4.setObjectName(_fromUtf8("gridFrame4"))
        self.gridLayout_8 = QtGui.QGridLayout(self.gridFrame4)
        self.gridLayout_8.setContentsMargins(9, 6, 6, 6)
        self.gridLayout_8.setObjectName(_fromUtf8("gridLayout_8"))
        self.horizontalWidget5 = QtGui.QWidget(self.gridFrame4)
        self.horizontalWidget5.setObjectName(_fromUtf8("horizontalWidget5"))
        self.horizontalLayout = QtGui.QHBoxLayout(self.horizontalWidget5)
        self.horizontalLayout.setMargin(0)
        self.horizontalLayout.setSpacing(10)
        self.horizontalLayout.setObjectName(_fromUtf8("horizontalLayout"))
        self.dataSetTree = QtGui.QTreeWidget(self.horizontalWidget5)
        self.dataSetTree.setStyleSheet(_fromUtf8("background-color: rgb(255, 255, 255);"))
        self.dataSetTree.setFrameShadow(QtGui.QFrame.Sunken)
        self.dataSetTree.setMidLineWidth(0)
        self.dataSetTree.setAnimated(True)
        self.dataSetTree.setObjectName(_fromUtf8("dataSetTree"))
        self.dataSetTree.headerItem().setTextAlignment(1, QtCore.Qt.AlignHCenter|QtCore.Qt.AlignVCenter|QtCore.Qt.AlignCenter)
        self.dataSetTree.headerItem().setTextAlignment(2, QtCore.Qt.AlignHCenter|QtCore.Qt.AlignVCenter|QtCore.Qt.AlignCenter)
        self.dataSetTree.header().setVisible(True)
        self.dataSetTree.header().setCascadingSectionResizes(False)
        self.dataSetTree.header().setDefaultSectionSize(200)
        self.dataSetTree.header().setMinimumSectionSize(20)
        self.dataSetTree.header().setStretchLastSection(False)
        self.horizontalLayout.addWidget(self.dataSetTree)
        self.attributesTree = QtGui.QTreeWidget(self.horizontalWidget5)
        self.attributesTree.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        self.attributesTree.setStyleSheet(_fromUtf8("background-color: rgb(255, 255, 255);"))
        self.attributesTree.setFrameShadow(QtGui.QFrame.Sunken)
        self.attributesTree.setEditTriggers(QtGui.QAbstractItemView.AllEditTriggers)
        self.attributesTree.setProperty("showDropIndicator", True)
        self.attributesTree.setDragEnabled(False)
        self.attributesTree.setDragDropOverwriteMode(False)
        self.attributesTree.setDragDropMode(QtGui.QAbstractItemView.NoDragDrop)
        self.attributesTree.setDefaultDropAction(QtCore.Qt.CopyAction)
        self.attributesTree.setSelectionMode(QtGui.QAbstractItemView.ExtendedSelection)
        self.attributesTree.setSelectionBehavior(QtGui.QAbstractItemView.SelectItems)
        self.attributesTree.setAutoExpandDelay(-1)
        self.attributesTree.setIndentation(20)
        self.attributesTree.setItemsExpandable(True)
        self.attributesTree.setAnimated(True)
        self.attributesTree.setColumnCount(3)
        self.attributesTree.setObjectName(_fromUtf8("attributesTree"))
        self.attributesTree.headerItem().setTextAlignment(1, QtCore.Qt.AlignHCenter|QtCore.Qt.AlignVCenter|QtCore.Qt.AlignCenter)
        self.attributesTree.headerItem().setTextAlignment(2, QtCore.Qt.AlignHCenter|QtCore.Qt.AlignVCenter|QtCore.Qt.AlignCenter)
        self.attributesTree.header().setDefaultSectionSize(200)
        self.attributesTree.header().setMinimumSectionSize(20)
        self.attributesTree.header().setSortIndicatorShown(True)
        self.attributesTree.header().setStretchLastSection(False)
        self.horizontalLayout.addWidget(self.attributesTree)
        self.gridLayout_8.addWidget(self.horizontalWidget5, 0, 0, 1, 1)
        self.gridLayout_3.addWidget(self.gridFrame4, 1, 0, 1, 2)
        self.gridFrame5 = QtGui.QFrame(self.centralwidget)
        self.gridFrame5.setMinimumSize(QtCore.QSize(30, 0))
        self.gridFrame5.setFrameShape(QtGui.QFrame.StyledPanel)
        self.gridFrame5.setObjectName(_fromUtf8("gridFrame5"))
        self.gridLayout_2 = QtGui.QGridLayout(self.gridFrame5)
        self.gridLayout_2.setMargin(6)
        self.gridLayout_2.setObjectName(_fromUtf8("gridLayout_2"))
        self.horizontalWidget6 = QtGui.QWidget(self.gridFrame5)
        self.horizontalWidget6.setObjectName(_fromUtf8("horizontalWidget6"))
        self.horizontalLayout_4 = QtGui.QHBoxLayout(self.horizontalWidget6)
        self.horizontalLayout_4.setMargin(0)
        self.horizontalLayout_4.setSpacing(6)
        self.horizontalLayout_4.setObjectName(_fromUtf8("horizontalLayout_4"))
        self.label = QtGui.QLabel(self.horizontalWidget6)
        font = QtGui.QFont()
        font.setBold(True)
        font.setWeight(75)
        self.label.setFont(font)
        self.label.setObjectName(_fromUtf8("label"))
        self.horizontalLayout_4.addWidget(self.label)
        self.currentSelectionValue = QtGui.QLineEdit(self.horizontalWidget6)
        self.currentSelectionValue.setStyleSheet(_fromUtf8("background-color: rgb(255, 255, 255);"))
        self.currentSelectionValue.setObjectName(_fromUtf8("currentSelectionValue"))
        self.horizontalLayout_4.addWidget(self.currentSelectionValue)
        self.gridLayout_2.addWidget(self.horizontalWidget6, 0, 0, 1, 1)
        self.gridLayout_3.addWidget(self.gridFrame5, 2, 0, 1, 2)
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtGui.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 1091, 25))
        self.menubar.setObjectName(_fromUtf8("menubar"))
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtGui.QStatusBar(MainWindow)
        self.statusbar.setObjectName(_fromUtf8("statusbar"))
        MainWindow.setStatusBar(self.statusbar)

        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        MainWindow.setWindowTitle(_translate("MainWindow", "MainWindow", None))
        self.loadDirectoryBtn.setText(_translate("MainWindow", "Load directory", None))
        self.reloadDirectoryBtn.setText(_translate("MainWindow", "Reload", None))
        self.workingDirectory.setPlaceholderText(_translate("MainWindow", "Directory of hdf5 files", None))
        self.label_5.setText(_translate("MainWindow", "Plotting", None))
        self.genericRadioBtn.setText(_translate("MainWindow", "Generic", None))
        self.timeSeriesRadioBtn.setText(_translate("MainWindow", "TimeSeries", None))
        self.spikesRadioBtn.setText(_translate("MainWindow", "Spikes", None))
        self.threeDRadioBtn.setText(_translate("MainWindow", "3D", None))
        self.ImageStackRadioBtn.setText(_translate("MainWindow", "ImageStack", None))
        self.dataSetSelection.setPlaceholderText(_translate("MainWindow", "Subselection e.g. 0,2,4-6", None))
        self.plotDataBtn.setText(_translate("MainWindow", "Plot data", None))
        self.addToPlotBtn.setText(_translate("MainWindow", "Add to plot", None))
        self.launchiPythonBtn.setText(_translate("MainWindow", "Launch iPython", None))
        self.sendToConsoleBtn.setText(_translate("MainWindow", "Send to concole", None))
        self.plotNameSpaceBtn.setText(_translate("MainWindow", "Print NameSpace", None))
        self.label_6.setText(_translate("MainWindow", "Analysis", None))
        self.dataSetTree.setSortingEnabled(True)
        self.dataSetTree.headerItem().setText(0, _translate("MainWindow", "DataSets", None))
        self.dataSetTree.headerItem().setText(1, _translate("MainWindow", "Size", None))
        self.dataSetTree.headerItem().setText(2, _translate("MainWindow", "Type", None))
        self.attributesTree.setSortingEnabled(True)
        self.attributesTree.headerItem().setText(0, _translate("MainWindow", "Attributes", None))
        self.attributesTree.headerItem().setText(1, _translate("MainWindow", "Value", None))
        self.attributesTree.headerItem().setText(2, _translate("MainWindow", "Type", None))
        self.label.setText(_translate("MainWindow", "Selection", None))

