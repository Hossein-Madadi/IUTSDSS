from PyQt4 import QtGui
from PyQt4.QtGui import *
from PyQt4.QtCore import *
from qgis.core import *
from qgis.utils import *
from osgeo import gdal, ogr, gdalconst
import os, platform, tempfile, pickle, processing, itertools, osr
from ui.ui_Footprint_Dialog_base import Ui_FootPrintClassDialogBase
from qgis.analysis import QgsRasterCalculator, QgsRasterCalculatorEntry
import numpy as np
from qgis.analysis import QgsGeometryAnalyzer


# Create the dialog.
class FootprintClassDialog(QtGui.QDialog):
    def __init__(self, iface):
        QtGui.QDialog.__init__(self)
        self.iface = iface

        # Set up the user interface from Designer.
        self.ui = Ui_FootPrintClassDialogBase()
        self.ui.setupUi(self)
        self.ui.btnProcess.setDisabled(True)

        # Set up the signals.
        self.connectSignals()

    # Connections.
    def connectSignals(self):
        self.ui.cBox_line.currentIndexChanged.connect(self.chek_fields)
        self.ui.linOutput.textChanged.connect(self.chek_fields)
        self.ui.pButton_Output.clicked.connect(self.select_output_file)
        self.ui.btnClose.clicked.connect(self.close)
        self.ui.btnProcess.clicked.connect(self.btnProcessClicked)

    # Get point layer into a list.
    def getLines(self):
        layers = QgsMapLayerRegistry.instance().mapLayers().values()
        road_list = []
        for layer in layers:
            if layer.type() == QgsMapLayer.VectorLayer and (
                            layer.wkbType() == QGis.WKBLineString or layer.wkbType() == QGis.WKBMultiLineString):
                road_list.append(layer)
        return road_list

    # Add Line layer fields to combo box.
    def loadRoadField(self):
        layers = self.getLines()
        self.ui.cBox_Width.clear()
        self.ui.cBox_Length.clear()
        self.ui.cBox_Vehicle.clear()
        self.ui.cBox_Fuel.clear()
        selectedLayerIndex = self.ui.cBox_line.currentIndex()
        selectedLayer = layers[selectedLayerIndex]
        source_layer_fields = [self.tr("choose a field")]
        for field in selectedLayer.pendingFields():
            if field.type() == QVariant.Int or field.type() == QVariant.Double or \
                            field.type() == QVariant.LongLong or field.type() == QVariant.ULongLong:
                source_layer_fields.append(unicode(field.name()))

        for f_label in source_layer_fields:
            self.ui.cBox_Width.addItem(f_label)
            self.ui.cBox_Length.addItem(f_label)
            self.ui.cBox_Vehicle.addItem(f_label)
            self.ui.cBox_Fuel.addItem(f_label)

    # Set number type
    def is_number(self, s):
        try:
            float(s)
            return True
        except ValueError:
            return False

    # Warning box.
    def userWarning(self, text, details):
        msg = QMessageBox()
        msg.setIcon(QMessageBox.Warning)
        msg.setText(text)
        msg.setWindowTitle("Warning")
        msg.setDetailedText(details)
        msg.exec_()

    # Update combo box fields for new run.
    def update_fields(self):
        line_list = self.getLines()

        # clear the combos.
        self.ui.cBox_line.clear()
        self.ui.cBox_Width.clear()
        self.ui.cBox_Length.clear()
        self.ui.cBox_Vehicle.clear()
        self.ui.cBox_Fuel.clear()
        # self.ui.Energy.clear()

        # Clear fields.
        self.ui.linOutput.clear()

        # Add layers to combos.
        for layer in line_list:
            self.ui.cBox_line.addItem(layer.name(), layer)
            self.linLayer = self.ui.cBox_line.itemData(self.ui.cBox_line.currentIndex())
            if self.linLayer is not None:
                self.loadRoadField()

    # Check input fields.
    def chek_fields(self):
        # Defining inputs
        self.layer = self.ui.cBox_line.itemData(self.ui.cBox_line.currentIndex())
        self.field1 = self.ui.cBox_Width.itemData(self.ui.cBox_Width.currentIndex())
        self.field2 = self.ui.cBox_Length.itemData(self.ui.cBox_Length.currentIndex())
        self.field3 = self.ui.cBox_Vehicle.itemData(self.ui.cBox_Vehicle.currentIndex())
        self.field4 = self.ui.cBox_Fuel.itemData(self.ui.cBox_Fuel.currentIndex())
        self.constant = self.ui.Energy.text()
        self.outPutName = self.ui.linOutput.text()

        f1 = str(self.ui.cBox_Width.currentText())
        f2 = str(self.ui.cBox_Length.currentText())
        f3 = str(self.ui.cBox_Vehicle.currentText())
        f4 = str(self.ui.cBox_Fuel.currentText())
        '''
        if not self.outPutName:
            if not (self.layer and self.field1 and self.field2 and self.field3 and self.field4):
                self.ui.btnProcess.setEnabled(True)
            elif (self.layer and self.field1 and self.field2 and self.field3 and self.field4) is not None:
                self.ui.btnProcess.setEnabled(True)
        elif self.outPutName is not None:
            if not (self.layer and self.field1 and self.field2 and self.field3 and self.field4):
                self.ui.btnProcess.setEnabled(True)
            elif (self.layer and self.field1 and self.field2 and self.field3 and self.field4) is not None:
                if not self.is_number(self.constant):
                    self.iface.messageBar().pushMessage("Error", "Constant value not numeric",
                                                        level=QgsMessageBar.CRITICAL,
                                                        duration=5)
                else:
                    self.ui.btnProcess.setEnabled(True)
        '''

        if not self.outPutName:
            if not self.layer:
                if (f1 and f2 and f3 and f4) == "choose a field":
                    self.ui.btnProcess.setEnabled(False)
            elif self.layer is not None:
                if (f1 and f2 and f3 and f4) == "choose a field":
                    self.ui.btnProcess.setEnabled(False)
        elif self.outPutName is not None:
            if not self.layer:
                if (f1 and f2 and f3 and f4) == "choose a field":
                    self.ui.btnProcess.setEnabled(False)
            elif self.layer is not None:
                if (f1 and f2 and f3 and f4) == "choose a field":
                    self.ui.btnProcess.setEnabled(False)
                elif (f1 and f2 and f3 and f4) != "choose a field":
                    if not self.is_number(self.constant):
                        self.iface.messageBar().pushMessage("Error", "Constant value not numeric",
                                                            level=QgsMessageBar.CRITICAL,
                                                            duration=5)
                    else:
                        self.ui.btnProcess.setEnabled(True)

    # Get Source Layer to action.
    def road_layer(self):
        for layer in self.getLines():
            if layer.name() == self.ui.cBox_line.currentText():
                self.vlayer = layer
                return self.vlayer

    # Check geometries of input layers.
    # If input data are valid then run model.
    def checkGeodata(self):
        player = self.point_layer()
        rlayer = self.grid_layer()
        lstLyr = []
        lstLyr.extend([player, rlayer])

        # Get crs of input raster layers.
        registry = QgsMapLayerRegistry.instance()
        crs_Lst = []
        for layer in lstLyr:
            lyr = QgsMapLayer.crs(layer)  # registry.mapLayersByName(layer)
            # xLyr = lyr[0].crs().authid()
            crs_Lst.append(lyr)
        if not self.all_same(crs_Lst):
            self.userWarning("Error", "Coordinate Reference System (CRS) "
                                      "must be the same for all input layers."
                                      "CRS of all input Layers is : '%s'."
                             % ((str(crs_Lst))))
        else:
            self.mean_distance()
            self.population_mapping()

    # Check for repeated item in a list.
    def all_same(self, items):
        return all(x == items[0] for x in items)

    # Define output path and name for vulnerability layer.
    def select_output_file(self):
        fileTypes = 'Shapefile (*.shp)'
        filename = QFileDialog.getSaveFileName(self, "Select output file ", "", fileTypes)
        if filename:
            # If the output file has no extension it considered as TIFF file.
            if not os.path.splitext(filename)[1]:
                outputFile = os.path.splitext(filename)[0] + '.shp'
            elif os.path.splitext(filename)[1] != '.tif':
                outputFile = os.path.splitext(filename)[0] + '.shp'
            else:
                outputFile = filename
        self.ui.linOutput.setText(outputFile)

    # Progress bar
    def progdialog(self, progress, txt):
        dialog = QProgressDialog()
        dialog.setCancelButton(None)
        dialog.setWindowTitle("Processing...")
        dialog.setLabelText(txt)
        bar = QProgressBar(dialog)
        bar.setTextVisible(True)
        bar.setValue(progress)
        dialog.setBar(bar)
        dialog.setMinimumWidth(300)
        dialog.show()
        return dialog, bar

    # Called when "Process" button pressed.  ########################################################################
    def btnProcessClicked(self):
        self.buffering()

    # Get temporary folder for intermediate layers.
    def gettempdir(self):
        platform_type = platform.system()
        tempplatform = tempfile.gettempdir()
        tempdir = ''
        if platform_type == 'Linux':
            tempdir = str(tempplatform) + '/tmpSDSSRoad/population'
        elif platform_type == 'Windows':
            tempdir = str(tempplatform) + '/tmpSDSSRoad/population'
        elif platform_type == 'Darwin':
            tempdir = str(tempplatform) + '/tmpSDSSRoad/population'
        return tempdir

    def buffer_size(self):
        layer = self.road_layer()

        # Get input values.
        road_width = self.ui.cBox_Width.currentText()
        road_length = self.ui.cBox_Length.currentText()
        vehicle = self.ui.cBox_Vehicle.currentText()
        fuel = self.ui.cBox_Fuel.currentText()
        constant = self.ui.Energy.text()

        #### Add new field to source layer, "buffer_w" field name
        new_Field = "buffer_m"
        provider = layer.dataProvider()
        caps = provider.capabilities()

        # Check if attribute is already there, return "-1" if not
        idx = provider.fieldNameIndex(new_Field)
        try:
            if idx == -1:
                if caps & QgsVectorDataProvider.AddAttributes:
                    res = provider.addAttributes([QgsField(new_Field, QVariant.Double, "double", 20, 2)])
        except:
            return False

        # Source table start editing.
        layer.updateFields()
        layer.startEditing()
        new_Field_idx = layer.fieldNameIndex(new_Field)

        # Component of formula.
        physical_footprint = "(\"" + str(road_width)+ "\"" + '*' + "\"" + str(road_length) + "\")"
        energy_footprint = ("(\"" + str(vehicle)+ "\"" + '*' + "\"" + str(fuel) + "\"" + '*' + constant + ")")

        # Calculate buffer size in METER unit.
        formula = QgsExpression("(((" + physical_footprint + " + " + energy_footprint + ")" + " / " + "\"" + str(road_length) + "\"" + "))" + " * " + '1000')
        formula.prepare(layer.pendingFields())
        for f in layer.getFeatures():
            value = formula.evaluate(f)
            f[new_Field_idx] = value
            layer.updateFeature(f)

        layer.commitChanges()
        return new_Field_idx

    # Create Buffer map in vector as Footprint.
    def buffering(self):
        layer = self.road_layer()
        idx = self.buffer_size()
        outputfile = self.ui.linOutput.text()

        QgsGeometryAnalyzer().buffer(layer,
                                     outputfile,
                                     -1,
                                     False,
                                     False,
                                     idx)
