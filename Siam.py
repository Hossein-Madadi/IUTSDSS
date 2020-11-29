from PyQt4.QtCore import QFileInfo
from PyQt4 import QtGui
from PyQt4.QtGui import *
from qgis.core import *
from qgis.analysis import QgsRasterCalculator, QgsRasterCalculatorEntry
from qgis.utils import *
from ui.ui_SIAM_Dialog_base import Ui_SIAMClassDialogBase
import collections
from osgeo import gdal, osr
import platform, tempfile
import processing
import math
import struct


# create the dialog.
class SiamClassDialog(QtGui.QDialog):
    def __init__(self, iface):
        QtGui.QDialog.__init__(self)
        self.iface = iface

        # Set up the user interface from Designer.
        self.ui = Ui_SIAMClassDialogBase()
        self.ui.setupUi(self)
        self.ui.btnProcess.setDisabled(True)

        # Set up the signals.
        self.connectSignals()

    # Connections.
    def connectSignals(self):
        self.ui.cBox_Study.currentIndexChanged.connect(self.chek_fields)
        self.ui.cBox_LuLc.currentIndexChanged.connect(self.chek_fields)
        self.ui.cBox_Pop.currentIndexChanged.connect(self.chek_fields)
        self.ui.cBox_Impact_1.currentIndexChanged.connect(self.chek_fields)
        self.ui.cBox_Impact_2.currentIndexChanged.connect(self.chek_fields)
        self.ui.linOutput.textChanged.connect(self.chek_fields)
        self.ui.pButton_Output.clicked.connect(self.select_output_file)
        self.ui.btnClose.clicked.connect(self.close)
        self.ui.btnProcess.clicked.connect(self.btnProcessClicked)

    # Get raster layers into a list.
    def getRasters(self):
        layers = QgsMapLayerRegistry.instance().mapLayers().values()
        raster_list = []
        for layer in layers:
            if layer.type() == QgsMapLayer.RasterLayer:
                raster_list.append(layer)
        return raster_list

    # Update combo box fields for new run.
    def update_fields(self):
        raster_list = self.getRasters()

        # clear the combos.
        self.ui.cBox_Study.clear()
        self.ui.cBox_LuLc.clear()
        self.ui.cBox_Pop.clear()
        self.ui.cBox_Impact_1.clear()
        self.ui.cBox_Impact_2.clear()

        # Clear fields.
        self.ui.linOutput.clear()

        # Add layers to combos.
        for layer in raster_list:
            self.ui.cBox_Study.addItem(layer.name(), layer)
            self.ui.cBox_LuLc.addItem(layer.name(), layer)
            self.ui.cBox_Pop.addItem(layer.name(), layer)
            self.ui.cBox_Impact_1.addItem(layer.name(), layer)
            self.ui.cBox_Impact_2.addItem(layer.name(), layer)

    # Check input fields.
    def chek_fields(self):
        # Defining layer input
        self.layer0 = self.ui.cBox_Study.itemData(self.ui.cBox_Study.currentIndex())
        self.layer1 = self.ui.cBox_LuLc.itemData(self.ui.cBox_LuLc.currentIndex())
        self.layer2 = self.ui.cBox_Pop.itemData(self.ui.cBox_Pop.currentIndex())
        self.layer3 = self.ui.cBox_Impact_1.itemData(self.ui.cBox_Impact_1.currentIndex())
        self.layer4 = self.ui.cBox_Impact_2.itemData(self.ui.cBox_Impact_2.currentIndex())

        # Make a list from raster layers in combo boxes.
        lst = []
        lst.extend([self.layer0, self.layer1, self.layer2, self.layer3, self.layer4])

        # Check layers in combo boxes are individual or not.
        lst2 = [item for item, count in collections.Counter(lst).items() if count > 1]

        # Disable process button if layers iterated.
        self.outPutName = self.ui.linOutput.text()
        if not self.outPutName:
            if not lst2:
                self.ui.btnProcess.setEnabled(False)
            elif lst2 is not None:
                self.ui.btnProcess.setEnabled(False)
        elif self.outPutName is not None:
            if not lst2:
                self.ui.btnProcess.setEnabled(True)
            elif lst2 is not None:
                self.ui.btnProcess.setEnabled(False)

        return lst

    # Add input layers name in a list.
    def readRasters(self):
        # Defining layer input.
        self.layer0 = self.ui.cBox_Study.currentText()
        self.layer1 = self.ui.cBox_LuLc.currentText()
        self.layer2 = self.ui.cBox_Pop.currentText()
        self.layer3 = self.ui.cBox_Impact_1.currentText()
        self.layer4 = self.ui.cBox_Impact_2.currentText()

        # Make a list from raster layers in combo boxes.
        lstLyr = []
        lstLyr.extend(
            [self.layer0, self.layer1, self.layer2, self.layer3, self.layer4])

        return lstLyr

    # Check cell size and geometries of input layers.
    # If input data are valid then run model.
    def checkGeodata(self):
        lstLyr = self.readRasters()

        # Get cell size of input raster layers.
        registry = QgsMapLayerRegistry.instance()
        cellSize_Lst = []

        for layer in lstLyr:
            lyr = registry.mapLayersByName(layer)
            xLyr = lyr[0].rasterUnitsPerPixelX()
            cellSize_Lst.append(xLyr)

        # Get crs of input raster layers.
        crs_Lst = []
        for layer in lstLyr:
            lyr = registry.mapLayersByName(layer)
            xLyr = lyr[0].crs().authid()
            crs_Lst.append(xLyr)

        if not self.all_same(cellSize_Lst):
            if not self.all_same(crs_Lst):
                self.userWarning("Error", "Pixel Size and Coordinate Reference System (CRS) "
                                          "must be the same for all input layers."
                                          "Pixel Size of all input Layers is : '%s' and "
                                          "CRS of all input Layers is : '%s'."
                                 % (str(cellSize_Lst), (str(crs_Lst))))
            elif self.all_same(crs_Lst) is not None:
                self.userWarning("Error", "Pixel Size "
                                          "must be the same for all input layers."
                                          "Pixel Size of all input Layers is : '%s'"
                                 % (str(cellSize_Lst)))
        elif self.all_same(cellSize_Lst) is not None:
            if not self.all_same(crs_Lst):
                self.userWarning("Error", "Coordinate Reference System (CRS)"
                                          " must be the same for all input layers."
                                          "CRS of all input Layers is : '%s'."
                                 % ((str(crs_Lst))))
            elif self.all_same(crs_Lst) is not None:
                self.lulc_unique()

    # Check for repeated item in a list.
    def all_same(self, items):
        return all(x == items[0] for x in items)

    # Define output path and name for vulnerability layer.
    def select_output_file(self):
        fileTypes = 'GeoTIFF Files (*.tif *.tiff)'
        filename = QFileDialog.getSaveFileName(self, "Select output file ", "", fileTypes)
        if filename:
            # If the output file has no extension it considered as TIFF file.
            if not os.path.splitext(filename)[1]:
                outputFile = os.path.splitext(filename)[0] + '.tif'

            elif os.path.splitext(filename)[1] != '.tif':
                outputFile = os.path.splitext(filename)[0] + '.tif'

            else:
                outputFile = filename
        self.ui.linOutput.setText(outputFile)

    # Warning box.
    def userWarning(self, text, details):
        msg = QMessageBox()
        msg.setIcon(QMessageBox.Warning)
        msg.setText(text)
        msg.setWindowTitle("Warning")
        msg.setDetailedText(details)
        msg.exec_()

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

    # Get No Data value
    def get_nodata(self, layer):
        extent = layer.extent()

        provider = layer.dataProvider()

        rows = layer.rasterUnitsPerPixelY()
        cols = layer.rasterUnitsPerPixelX()

        block = provider.block(1, extent, rows, cols)

        return block.noDataValue()

    # Create temporary folder for intermediate layers.
    def gettempdir(self):
        platform_type = platform.system()
        tempplatform = tempfile.gettempdir()
        tempdir = ''
        if platform_type == 'Linux':
            tempdir = str(tempplatform) + '/tmpSDSSRoad/siam'
        elif platform_type == 'Windows':
            tempdir = str(tempplatform) + '/tmpSDSSRoad/siam'
        return tempdir

    # Called when "Process" button pressed.  ###########################################################################
    def btnProcessClicked(self):
        self.siam()

    # Get LULC raster layer to action.
    def lulc_layer(self):
        for layer in self.getRasters():
            if layer.name() == self.ui.cBox_LuLc.currentText():
                self.lulclayer = layer
        return self.lulclayer

    # Get raster map unique codes.
    def unique(self, layer):
        # lulclayer = self.lulc_layer()
        filepath = layer.source()
        # load raster
        gdalData = gdal.Open(filepath)

        # get width and heights of the raster
        xsize = gdalData.RasterXSize
        ysize = gdalData.RasterYSize

        # get number of bands
        bands = gdalData.RasterCount

        # process the raster
        band_i = gdalData.GetRasterBand(1)
        raster = band_i.ReadAsArray()

        # create dictionary for unique values count
        count = {}

        # count unique values for the given band
        for col in range(xsize):
            for row in range(ysize):
                cell_value = raster[row, col]

                # check if cell_value is NaN
                if math.isnan(cell_value):
                    cell_value = 'Null'

                else:
                    cell_value = round(cell_value)

                # add cell_value to dictionary
                try:
                    count[cell_value] += 1
                except:
                    count[cell_value] = 1

        codes = count.items()
        return codes

    # Get Land use/Land Cover types
    def lulc_classess(self):
        lulc_layer = self.lulc_layer()
        codes = self.unique(lulc_layer)
        nodata = self.get_nodata(lulc_layer)
        lst1 = []
        lst2 = []
        lst_nodata = []
        for key, value in codes:
            if 0 <= key < 10:
                lst1.append(int(key))
            elif key == nodata:
                lst_nodata.append(key)
            elif key < 0: # for values such as -99999 or -3.4028231e+38 etc.
                lst_nodata.append(key)
            else: # if not 0 < key < 10:
                lst2.append(key)

        if len(lst2) > 0:
            self.userWarning("Error", "Land Use/Land Cover raster map has "
                                      "incorrect code for types, Please See Help.")
        else:
            pass
        return lst1

    # Get First Impact raster layer to action
    def first_impact_layer(self):
        for layer in self.getRasters():
            if layer.name() == self.ui.cBox_Impact_1.currentText():
                self.firstimpact = layer
        return self.firstimpact

    # Get First impact layer classes
    def first_impact_classes(self):
        first_impact = self.first_impact_layer()
        layer_name = first_impact.name()
        codes = self.unique(first_impact)
        nodata = self.get_nodata(first_impact)
        lst1 = []
        lst2 = []
        lst_nodata = []
        for key, value in codes:
            if 0 <= key <= 10:
                lst1.append(int(key))
            elif key == nodata:
                lst_nodata.append(key)
            else: # if not 0 < key < 10:
                lst2.append(key)

        if len(lst2) > 0:
            self.userWarning("Error", str(layer_name) + " raster map has "
                                      "incorrect code for types, Please See Help.")
        else:
            pass
        return lst1

    # Get Second Impact raster layer to action
    def second_impact_layer(self):
        for layer in self.getRasters():
            if layer.name() == self.ui.cBox_Impact_2.currentText():
                self.secondimpact = layer
        return self.secondimpact

    # Get Second impact layer classes
    def second_impact_classes(self):
        second_impact = self.second_impact_layer()
        layer_name = second_impact.name()
        codes = self.unique(second_impact)
        nodata = self.get_nodata(second_impact)
        lst1 = []
        lst2 = []
        lst_nodata = []
        for key, value in codes:
            if 0 <= key <= 10:
                lst1.append(int(key))
            elif key == nodata:
                lst_nodata.append(key)
            else: # if not 0 < key < 10:
                lst2.append(key)

        if len(lst2) > 0:
            self.userWarning("Error", str(layer_name) + " raster map has "
                                      "incorrect code for types, Please See Help.")
        else:
            pass

        return lst1

    # Get Study area raster layer to action
    def study_layer(self):
        for layer in self.getRasters():
            if layer.name() == self.ui.cBox_Study.currentText():
                self.study = layer
        return self.study

    # Check study area raster layer is in boolean (0, 1).
    def study_classes(self):
        study = self.study_layer()
        layer_name = study.name()
        codes = self.unique(study)
        nodata = self.get_nodata(study)
        lst1 = []
        lst2 = []
        lst_nodata = []
        for key, value in codes:
            if 0 <= key <= 1:
                lst1.append(key)
            elif key == nodata:
                lst_nodata.append(key)
            else:  # if not 0 < key < 10:
                lst2.append(key)

        if len(lst2) > 0:
            self.userWarning("Error", str(layer_name) + " raster map has "
                                                        "incorrect code for types, Please See Help.")
        else:
            pass

    # Convert First Impact raster layer classes to individual maps
    #  (4 class --> 4 raster) and save them in temporary directory.
    def first_impact_rasters(self):
        tempdir = self.gettempdir()
        classes = self.first_impact_classes()
        layer = self.first_impact_layer()

        lyrPath = layer.source()
        lyrName = layer.name()
        lyr2 = QgsRasterLayer(lyrPath, lyrName)

        entries = []
        ras1 = QgsRasterCalculatorEntry()
        ras1.ref = lyrName + '@1'
        ras1.raster = lyr2
        ras1.bandNumber = 1
        entries.append(ras1)

        # List of output rasters for lulc types
        first_list = []

        for i in classes:
            formula = "(" + "\"" + ras1.ref + "\"" + ' = ' + str(i) + " ) " + " * " + str(i)
            output = tempdir + "/" + str(lyrName)+ "_%s.tif" % str(i)
            first_list.append(output)
            calc = QgsRasterCalculator(formula, output, 'GTiff', lyr2.extent(), lyr2.width(), lyr2.height(), entries)
            calc.processCalculation()
        del entries
        return first_list

    # Convert Second Impact raster layer classes to individual maps
    #  (4 class --> 4 raster) and save them in temporary directory.
    def second_impact_rasters(self):
        tempdir = self.gettempdir()
        classes = self.second_impact_classes()
        layer = self.second_impact_layer()

        lyrPath = layer.source()
        lyrName = layer.name()
        lyr2 = QgsRasterLayer(lyrPath, lyrName)

        entries = []
        ras1 = QgsRasterCalculatorEntry()
        ras1.ref = lyrName + '@1'
        ras1.raster = lyr2
        ras1.bandNumber = 1
        entries.append(ras1)

        # List of output rasters for lulc types
        second_list = []

        for i in classes:
            formula = "(" + "\"" + ras1.ref + "\"" + ' = ' + str(i) + " ) " + " * " + str(i)
            output = tempdir + "/" + str(lyrName) + "_%s.tif" % str(i)
            second_list.append(output)
            calc = QgsRasterCalculator(formula, output, 'GTiff', lyr2.extent(), lyr2.width(), lyr2.height(),
                                       entries)
            calc.processCalculation()
        del entries
        return second_list

    #Multiply Study area Layer and LULC layer.
    def lulc_study_area(self):
        tempdir = self.gettempdir()
        lulc_layer = self.lulc_layer()
        study_area = self.study_layer()

        # Get Study area layer info.
        studylyrName = study_area.name()
        entries = []
        studyras = QgsRasterCalculatorEntry()
        studyras.ref = studylyrName + '@1'
        studyras.raster = study_area
        studyras.bandNumber = 1
        entries.append(studyras)

        # Get lulc layer info.
        lulcName = lulc_layer.name()
        lulcras = QgsRasterCalculatorEntry()
        lulcras.ref = lulcName + '@1'
        lulcras.raster = lulc_layer
        lulcras.bandNumber = 1
        entries.append(lulcras)

        formula = "\"" + studyras.ref + "\"" + ' * ' + "\"" + lulcras.ref + "\""
        output_lulc_study = tempdir + "/%s_study.tif" % str(lulcName)
        calc = QgsRasterCalculator(formula, output_lulc_study, 'GTiff', study_area.extent(), study_area.width(),
                                    study_area.height(), entries)
        calc.processCalculation()
        del entries

        return output_lulc_study

    # Create LULC individual layers based on lulc types.
    def lulc_rasters(self):
        tempdir = self.gettempdir()
        classes = self.lulc_classess()
        layerpath = self.lulc_study_area()

        fileInfo = QFileInfo(layerpath)
        baseName = fileInfo.baseName()
        lyr = QgsRasterLayer(layerpath, baseName)

        entries = []
        ras = QgsRasterCalculatorEntry()
        ras.ref = baseName + '@1'
        ras.raster = lyr
        ras.bandNumber = 1
        entries.append(ras)

        # List of output rasters for lulc types
        lulc_list = []

        for i in classes:
            formula = "(" + "\"" + ras.ref + "\"" + ' = ' + str(i) + " ) " + " * " + str(i)
            output = tempdir + "/" + str(baseName) + "_%s.tif" % str(i)
            lulc_list.append(output)
            calc = QgsRasterCalculator(formula, output, 'GTiff', lyr.extent(), lyr.width(), lyr.height(),
                                       entries)
            calc.processCalculation()
        del entries
        return lulc_list

    # Get number of pixels in LULC types of individual raster layer.
    def total_pixel_lulc_classes(self):
        # Get full path of new raster layers extracted from lulc layer.
        lulc_list_fullpath = self.lulc_rasters()

        # Get list of pixel value classes in lulc layer.
        lulc_classes = self.lulc_classess()

        # Dict for lulc classes (as key) and pixel numbers (as value).
        lulc_pixels_dict = {}

        for lyr in lulc_list_fullpath:
            layer = QgsRasterLayer(lyr)
            provider = layer.dataProvider()

            fmttypes = {'Byte': 'B', 'UInt16': 'H', 'Int16': 'h', 'UInt32': 'I', 'Int32': 'i', 'Float32': 'f',
                        'Float64': 'd'}

            lyr_path = provider.dataSourceUri()
            dataset = gdal.Open(lyr_path)
            band = dataset.GetRasterBand(1)
            BandType = gdal.GetDataTypeName(band.DataType)

            # Number of pixels for each class.
            count_value = 0

            for y in range(band.YSize):

                scanline = band.ReadRaster(0, y, band.XSize, 1, band.XSize, 1, band.DataType)
                values = struct.unpack(fmttypes[BandType] * band.XSize, scanline)

                for value in values:
                    for i in lulc_classes:
                        if value == i and value != 0:
                            count_value += 1
            dataset = None

            lulc_pixels_dict[lyr] = count_value
        return lulc_pixels_dict

    # Multiply lulc classes and impact class layers.
    def multiply_lulc_impacts(self):
        tempdir = self.gettempdir()
        first_list = self.first_impact_rasters()
        second_list = self.second_impact_rasters()
        lulc_list = self.lulc_rasters()

        first_impact_lulc_list = []
        second_impact_lulc_list = []

        # Get lulc layers info.
        for lulc_fullpathes in lulc_list:
            fileInfo = QFileInfo(lulc_fullpathes)
            baseName = fileInfo.baseName()
            lulclyr = processing.getObject(str(lulc_fullpathes))

            ras = QgsRasterCalculatorEntry()
            ras.ref = baseName + '@1'
            ras.raster = lulclyr
            ras.bandNumber = 1

            # multiply lulc classes with first impact class layers.
            for fullpathlayer1 in first_list:
                fileInfo1 = QFileInfo(fullpathlayer1)
                baseName1 = fileInfo1.baseName()
                firstlyr = processing.getObject(str(fullpathlayer1))
                ras1 = QgsRasterCalculatorEntry()
                ras1.ref = baseName1 + '@1'
                ras1.raster = firstlyr
                ras1.bandNumber = 1
                entries1 = []
                entries1.append(ras)
                entries1.append(ras1)

                formula1 = "\"" + ras1.ref + "\"" + ' * ' + "\"" + ras.ref + "\""
                output1 = tempdir + "/%s_%s.tif" % (str(baseName1), str(baseName))
                first_impact_lulc_list.append(output1)
                calc1 = QgsRasterCalculator(formula1, output1, 'GTiff', firstlyr.extent(), firstlyr.width(),
                                            firstlyr.height(), entries1)
                calc1.processCalculation()
                del entries1

            # multiply lulc classes with second impact class layers.
            for fullpathlayer2 in second_list:
                fileInfo2 = QFileInfo(fullpathlayer2)
                baseName2 = fileInfo2.baseName()
                secondlyr = processing.getObject(str(fullpathlayer2))
                ras2 = QgsRasterCalculatorEntry()
                ras2.ref = baseName2 + '@1'
                ras2.raster = secondlyr
                ras2.bandNumber = 1
                entries2 = []
                entries2.append(ras)
                entries2.append(ras2)

                formula2 = "\"" + ras2.ref + "\"" + ' * ' + "\"" + ras.ref + "\""
                output2 = tempdir + "/%s_%s.tif" % (str(baseName2), str(baseName))
                second_impact_lulc_list.append(output2)
                calc2 = QgsRasterCalculator(formula2, output2, 'GTiff', secondlyr.extent(), secondlyr.width(),
                                            secondlyr.height(), entries2)
                calc2.processCalculation()
                del entries2

        #print second_impact_lulc_list
        return [first_impact_lulc_list, second_impact_lulc_list]

    # Count Totatl pixel of study area raster layer (as total area).
    def total_study_pixels(self):
        study_area = self.study_layer()
        provider = study_area.dataProvider()

        fmttypes = {'Byte': 'B', 'UInt16': 'H', 'Int16': 'h', 'UInt32': 'I', 'Int32': 'i', 'Float32': 'f',
                    'Float64': 'd'}

        lyr_path = provider.dataSourceUri()
        dataset = gdal.Open(lyr_path)
        band = dataset.GetRasterBand(1)
        BandType = gdal.GetDataTypeName(band.DataType)
        count_value = 0

        for y in range(band.YSize):

            scanline = band.ReadRaster(0, y, band.XSize, 1, band.XSize, 1, band.DataType)
            values = struct.unpack(fmttypes[BandType] * band.XSize, scanline)

            for value in values:
                if value == 1:
                    count_value += 1
        if count_value == 0:
            self.userWarning("Error", "Study Area raster layer must be in boolean (0, 1) format, See Help")
        dataset = None
        return count_value

    # Count total number of pixels for each class of first and second impact layer multiply lulc classes.
    def total_pixel_first_second_impact_classes(self):
        # Get full path of new raster layers extracted from first impact layer classes multiply lulc layer.
        first_impact_multiply_fullpath, second_impact_multiply_fullpath = self.multiply_lulc_impacts()

        # Dict for impact layers (as key) and pixel numbers (as value).
        first_impact_dict = {}
        second_impact_dict = {}

        for lyr1 in first_impact_multiply_fullpath:
            layer1 = QgsRasterLayer(lyr1)
            provider1 = layer1.dataProvider()

            fmttypes = {'Byte': 'B', 'UInt16': 'H', 'Int16': 'h', 'UInt32': 'I', 'Int32': 'i', 'Float32': 'f',
                        'Float64': 'd'}

            lyr_path1 = provider1.dataSourceUri()
            # (root, filename) = path.split(lyr_path)
            dataset1 = gdal.Open(lyr_path1)
            band1 = dataset1.GetRasterBand(1)
            BandType1 = gdal.GetDataTypeName(band1.DataType)

            # Number of pixels for each impact class i.e. (0, 3, 7, 10).
            count_value1 = 0

            for y in range(band1.YSize):

                scanline = band1.ReadRaster(0, y, band1.XSize, 1, band1.XSize, 1, band1.DataType)
                values1 = struct.unpack(fmttypes[BandType1] * band1.XSize, scanline)

                for value1 in values1:
                    #for i in first_impact_classes:
                    if value1 > 0:#if value1 == i and value1 != 0:
                        count_value1 += 1
            dataset = None

            first_impact_dict[lyr1] = count_value1

        for lyr2 in second_impact_multiply_fullpath:
            layer2 = QgsRasterLayer(lyr2)
            provider2 = layer2.dataProvider()

            fmttypes = {'Byte': 'B', 'UInt16': 'H', 'Int16': 'h', 'UInt32': 'I', 'Int32': 'i', 'Float32': 'f',
                        'Float64': 'd'}

            lyr_path2 = provider2.dataSourceUri()
            dataset2 = gdal.Open(lyr_path2)
            band2 = dataset2.GetRasterBand(1)
            BandType2 = gdal.GetDataTypeName(band2.DataType)

            # Number of pixels for each impact class i.e. (0, 3, 7, 10).
            count_value2 = 0

            for y in range(band2.YSize):

                scanline = band2.ReadRaster(0, y, band2.XSize, 1, band2.XSize, 1, band2.DataType)
                values2 = struct.unpack(fmttypes[BandType2] * band2.XSize, scanline)

                for value2 in values2:
                    #for i in second_impact_classes:
                    if value2 > 0:#if value2 == i and value2 != 0:
                        count_value2 += 1
            dataset = None

            second_impact_dict[lyr2] = count_value2
        return [first_impact_dict, second_impact_dict]

    ##################### Population Section
    # Get Study area raster layer to action
    def pop_layer(self):
        for layer in self.getRasters():
            if layer.name() == self.ui.cBox_Pop.currentText():
                self.pop = layer
        return self.pop

    # Extract population in study area.
    def pop_boundary(self):
        tempdir = self.gettempdir()
        study_area = self.study_layer()
        pop_layer = self.pop_layer()
        entries = []

        # Get Study area layer info.
        studylyrName = study_area.name()
        studyras = QgsRasterCalculatorEntry()
        studyras.ref = studylyrName + '@1'
        studyras.raster = study_area
        studyras.bandNumber = 1
        entries.append(studyras)

        # Get Population area layer info.
        poplyrName = pop_layer.name()
        popras = QgsRasterCalculatorEntry()
        popras.ref = poplyrName + '@1'
        popras.raster = pop_layer
        popras.bandNumber = 1
        entries.append(popras)

        formula = "\"" + studyras.ref + "\"" + ' * ' + "\"" + popras.ref + "\""
        #print formula
        pop_output = tempdir + "/%s_study.tif" % str(poplyrName)
        #print pop_output
        calc = QgsRasterCalculator(formula, pop_output, 'GTiff', study_area.extent(), study_area.width(),
                                    study_area.height(), entries)
        calc.processCalculation()
        del entries
        return pop_output

    # Sum all pixel values for a raster layer.
    def sumValue(self, raster):
        provider = raster.dataProvider()
        extent = provider.extent()

        rows = raster.height()
        cols = raster.width()
        block = provider.block(1, extent, cols, rows)
        sum = 0

        for i in range(rows):
            for j in range(cols):
                if block.value(i, j) >= 0:
                    sum += block.value(i, j)
        return sum

    # Sum all population in study area.
    def pop_study(self):
        # Get full path of population map.
        pop_output = self.pop_boundary()
        pop_layer = QgsRasterLayer(pop_output)
        sum_pop = self.sumValue(pop_layer)
        return sum_pop

    # Multiply Population of study area and first impact class layer.
    def multiply_pop_first_impacts(self):
        tempdir = self.gettempdir()
        first_list = self.first_impact_rasters()
        pop_boundary_fullpath = self.pop_boundary()

        # Get Population area layer info.
        fileInfo = QFileInfo(pop_boundary_fullpath)
        poplyrName = fileInfo.baseName()
        poplyr = QgsRasterLayer(pop_boundary_fullpath, poplyrName)

        # poplyrName = pop_layer.name()
        popras = QgsRasterCalculatorEntry()
        popras.ref = poplyrName + '@1'
        popras.raster = poplyr
        popras.bandNumber = 1

        first_pop_list = []
        # multiply lulc classes with first impact class layers.
        for fullpathlayer1 in first_list:
            fileInfo1 = QFileInfo(fullpathlayer1)
            baseName1 = fileInfo1.baseName()
            firstlyr = processing.getObject(str(fullpathlayer1))
            ras1 = QgsRasterCalculatorEntry()
            ras1.ref = baseName1 + '@1'
            ras1.raster = firstlyr
            ras1.bandNumber = 1
            entries1 = []
            entries1.append(popras)
            entries1.append(ras1)

            formula1 = "\"" + ras1.ref + "\"" + ' * ' + "\"" + popras.ref + "\""
            output1 = tempdir + "/%s_%s.tif" % (str(baseName1), str(poplyrName))
            first_pop_list.append(output1)
            calc1 = QgsRasterCalculator(formula1, output1, 'GTiff', firstlyr.extent(), firstlyr.width(),
                                        firstlyr.height(), entries1)
            calc1.processCalculation()
            del entries1
        return first_pop_list

    # Multiply Population of study area and second impact class layer.
    def multiply_pop_second_impacts(self):
        tempdir = self.gettempdir()
        second_list = self.second_impact_rasters()
        pop_boundary_fullpath = self.pop_boundary()

        # Get Population area layer info.
        fileInfo = QFileInfo(pop_boundary_fullpath)
        poplyrName = fileInfo.baseName()
        poplyr = QgsRasterLayer(pop_boundary_fullpath, poplyrName)

        # poplyrName = pop_layer.name()
        popras = QgsRasterCalculatorEntry()
        popras.ref = poplyrName + '@1'
        popras.raster = poplyr
        popras.bandNumber = 1

        second_pop_list = []
        # multiply lulc classes with first impact class layers.
        for fullpathlayer1 in second_list:
            fileInfo1 = QFileInfo(fullpathlayer1)
            baseName1 = fileInfo1.baseName()
            firstlyr = processing.getObject(str(fullpathlayer1))
            ras1 = QgsRasterCalculatorEntry()
            ras1.ref = baseName1 + '@1'
            ras1.raster = firstlyr
            ras1.bandNumber = 1
            entries1 = []
            entries1.append(popras)
            entries1.append(ras1)

            formula1 = "\"" + ras1.ref + "\"" + ' * ' + "\"" + popras.ref + "\""
            output1 = tempdir + "/%s_%s.tif" % (str(baseName1), str(poplyrName))
            second_pop_list.append(output1)
            calc1 = QgsRasterCalculator(formula1, output1, 'GTiff', firstlyr.extent(), firstlyr.width(),
                                        firstlyr.height(), entries1)
            calc1.processCalculation()
            del entries1
        return second_pop_list

    # Count population for each class of first and second impact layer.
    def total_population_first_second_impact_classes(self):
        # Get full path of new raster layers extracted from first impact layer classes multiply study area.
        first_impact_multiply_fullpath = self.multiply_pop_first_impacts()
        second_impact_multiply_fullpath = self.multiply_pop_second_impacts()

        # Dict for impact layers (as key) and pixel numbers (as value).
        pop_first_impact_dict = {}
        pop_second_impact_dict = {}

        for lyr1 in first_impact_multiply_fullpath:
            layer1 = QgsRasterLayer(lyr1)
            sum_pop_first_impact = self.sumValue(layer1)
            pop_first_impact_dict[lyr1] = sum_pop_first_impact

        for lyr2 in second_impact_multiply_fullpath:
            layer2 = QgsRasterLayer(lyr2)
            sum_pop_second_impact = self.sumValue(layer2)
            pop_second_impact_dict[lyr2] = sum_pop_second_impact

        return [pop_first_impact_dict, pop_second_impact_dict]

    # Check results from all functions for input in SIAM model.
    def check_siam_inputs(self):
        # Get weights
        area_weight = self.ui.lineEdit_1.text()
        #print "Area Weight = " + area_weight
        pop_weight = self.ui.lineEdit_2.text()
        #print "Population Weight = " + pop_weight

        # Get full path of new raster layers extracted from first impact layer classes multiply study area.
        first_impact_multiply_fullpath, second_impact_multiply_fullpath = self.multiply_area_impacts()
        #print "Full path First impact layers = " + str(first_impact_multiply_fullpath)
        #print "Full path Second impact layers = " + str(second_impact_multiply_fullpath)

        # Get list of pixel value classes in first impact layer.
        first_impact_classes = self.first_impact_classes()
        #print "First impact classes = " + str(first_impact_classes)

        # Get list of pixel value classes in second impact layer.
        second_impact_classes = self.second_impact_classes()
        #print "Second impact classes = " + str(second_impact_classes)

        # Get total pixel of study area (as total area).
        total_study_pixels = self.total_study_pixels()
        #print "Total Pixel Study Area = " + str(total_study_pixels)

        # Get pixel count of each extracted impact layers (as impacted area in each class).
        first_impact_dict, second_impact_dict = self.total_pixel_first_second_impact_classes()
        #print "Pixel Count First impact classes = " + str(first_impact_dict)
        #print "Pixel Count Second impact classes = " + str(second_impact_dict)

        # Get total Population in study area.
        sum_pop = self.pop_study()
        #print "Totatl Population in study area = " + str(sum_pop)

        # Get Population count of each extracted impact layers (as impacted population in each class).
        pop_first_impact_dict, pop_second_impact_dict = self.total_population_first_second_impact_classes()
        #print "Population count First impact classes = " + str(pop_first_impact_dict)
        second_pop_list = self.multiply_pop_second_impacts()
        #print "Population count Second impact classes = " + str(pop_second_impact_dict)

    # Prepare layers required for Siam calculation.
    def pre_siam(self):
        # Get temporary directory.
        tempdir = self.gettempdir()

        # Get full path of new raster layers extracted from first impact layer classes multiply study area
        # and pixel count of each extracted impact layers (as impacted area in each class).
        first_impact_dict, second_impact_dict = self.total_pixel_first_second_impact_classes()

        # Get total pixel of study area (as total area).
        total_study_pixels = self.total_study_pixels()

        # Get Population count of each extracted impact layers (as impacted population in each class).
        pop_first_impact_dict, pop_second_impact_dict = self.total_population_first_second_impact_classes()

        # Get total Population in study area.
        sum_pop = self.pop_study()

        first_area_list = []
        second_area_list = []
        first_pop_list = []
        second_pop_list = []

        # Refer fullpath as key and pixel count as value. (First impact layers and area)
        entries1 = []
        for fullpath1a, count1a in first_impact_dict.items():
            # Get layer info.
            fileInfo1a = QFileInfo(str(fullpath1a))
            lyr1a_Name = fileInfo1a.baseName()
            impactlyr1a = processing.getObject(str(fullpath1a)) # For prevent QGIS crashes instead of QgsRasterLayer
            firstras_a = QgsRasterCalculatorEntry()
            firstras_a.ref = lyr1a_Name + '@1'
            firstras_a.raster = impactlyr1a
            firstras_a.bandNumber = 1
            entries1.append(firstras_a)
            formula_area1 = "\"" + firstras_a.ref + "\"" + ' * ' + "(" + str(count1a) + '/' + str(total_study_pixels) + ")"
            output_area1 = tempdir + "/" + str(lyr1a_Name) + "_area1.tif"
            first_area_list.append(output_area1)
            calc_area1 = QgsRasterCalculator(formula_area1, output_area1, 'GTiff', impactlyr1a.extent(), impactlyr1a.width(), impactlyr1a.height(), entries1)
            calc_area1.processCalculation()
        del entries1

        # Refer fullpath as key and pixel count as value. (Second impact layers and area)
        entries2 = []
        for fullpath2a, count2a in second_impact_dict.items():
            # Get layer info.
            fileInfo2a = QFileInfo(str(fullpath2a))
            lyr2a_Name = fileInfo2a.baseName()
            impactlyr2a = processing.getObject(str(fullpath2a)) # For prevent QGIS crashes instead of QgsRasterLayer
            secondras_a = QgsRasterCalculatorEntry()
            secondras_a.ref = lyr2a_Name + '@1'
            secondras_a.raster = impactlyr2a
            secondras_a.bandNumber = 1
            entries2.append(secondras_a)
            formula_area2 = "\"" + secondras_a.ref + "\"" + ' * ' + "(" + str(count2a) + '/' + str(total_study_pixels) + ")"
            output_area2 = tempdir + "/" + str(lyr2a_Name) + "_area2.tif"
            second_area_list.append(output_area2)
            calc_area2 = QgsRasterCalculator(formula_area2, output_area2, 'GTiff', impactlyr2a.extent(), impactlyr2a.width(), impactlyr2a.height(), entries2)
            calc_area2.processCalculation()
        del entries2

        # Refer fullpath as key and pixel count as value. (First impact layers and population)
        entries3 = []
        for fullpath1p, count1p in pop_first_impact_dict.items():
            # Get layer info.
            fileInfo1p = QFileInfo(str(fullpath1p))
            lyr1p_Name = fileInfo1p.baseName()
            impactlyr1p = processing.getObject(str(fullpath1p))  # For prevent QGIS crashes instead of QgsRasterLayer
            firstras_p = QgsRasterCalculatorEntry()
            firstras_p.ref = lyr1p_Name + '@1'
            firstras_p.raster = impactlyr1p
            firstras_p.bandNumber = 1
            entries3.append(firstras_p)
            formula_pop1 = "\"" + firstras_p.ref + "\"" + ' * ' + "(" + str(count1p) + '/' + str(
                sum_pop) + ")"
            output_pop1 = tempdir + "/" + str(lyr1p_Name) + "_pop1.tif"
            first_pop_list.append(output_pop1)
            calc_pop1 = QgsRasterCalculator(formula_pop1, output_pop1, 'GTiff', impactlyr1p.extent(),
                                            impactlyr1p.width(), impactlyr1p.height(), entries3)
            calc_pop1.processCalculation()
        del entries3

        # Refer fullpath as key and pixel count as value. (Second impact layers and population)
        entries4 = []
        for fullpath2p, count2p in pop_second_impact_dict.items():
            # Get layer info.
            fileInfo2p = QFileInfo(str(fullpath2p))
            lyr2p_Name = fileInfo2p.baseName()
            impactlyr2p = processing.getObject(str(fullpath2p))  # For prevent QGIS crashes instead of QgsRasterLayer
            secondras_p = QgsRasterCalculatorEntry()
            secondras_p.ref = lyr2p_Name + '@1'
            secondras_p.raster = impactlyr2p
            secondras_p.bandNumber = 1
            entries4.append(secondras_p)
            formula_pop2 = "\"" + secondras_p.ref + "\"" + ' * ' + "(" + str(count2p) + '/' + str(
                sum_pop) + ")"
            output_pop2 = tempdir + "/" + str(lyr2p_Name) + "_pop2.tif"
            second_pop_list.append(output_pop2)
            calc_pop2 = QgsRasterCalculator(formula_pop2, output_pop2, 'GTiff', impactlyr2p.extent(),
                                            impactlyr2p.width(), impactlyr2p.height(), entries4)
            calc_pop2.processCalculation()
        del entries4

        return [first_area_list, second_area_list, first_pop_list, second_pop_list]

    # Sum first impact-area classes.
    def sum_first_impact_area(self):
        # Get temporary directory.
        tempdir = self.gettempdir()

        # Get layers already created in temporary directory.
        first_area_list, second_area_list, first_pop_list, second_pop_list = self.pre_siam()

        entries = []

        for i in range(0, len(first_area_list)):
            raster = first_area_list[i]
            readRst = processing.getObject(str(raster))
            ras1 = QgsRasterCalculatorEntry()
            ras1.raster = readRst
            ras1.ref = "first_impact_area_" + str(i + 1) + "@1"
            ras1.bandNumber = 1
            entries.append(ras1)

        reflist = " + ".join([ent.ref for ent in entries])
        formula = '(' + reflist + ')'

        readRst = QgsRasterLayer(first_area_list[0])
        sum_first_impact_area = tempdir + "/sum_first_impact_area.tif"
        calc = QgsRasterCalculator(formula, sum_first_impact_area, 'GTiff', readRst.extent(), readRst.width(), readRst.height(),
                                   entries)
        calc.processCalculation()
        return sum_first_impact_area

    # Sum second impact-area classes.
    def sum_second_impact_area(self):
        # Get temporary directory.
        tempdir = self.gettempdir()

        # Get layers already created in temporary directory.
        first_area_list, second_area_list, first_pop_list, second_pop_list = self.pre_siam()

        entries = []

        for i in range(0, len(second_area_list)):
            raster = second_area_list[i]
            readRst = processing.getObject(str(raster))
            ras1 = QgsRasterCalculatorEntry()
            ras1.raster = readRst
            ras1.ref = "second_impact_area_" + str(i + 1) + "@1"
            ras1.bandNumber = 1
            entries.append(ras1)

        reflist = " + ".join([ent.ref for ent in entries])
        formula = '(' + reflist + ')'

        readRst = QgsRasterLayer(second_area_list[0])
        sum_second_impact_area = tempdir + "/sum_second_impact_area.tif"
        calc = QgsRasterCalculator(formula, sum_second_impact_area, 'GTiff', readRst.extent(), readRst.width(), readRst.height(),
                                   entries)
        calc.processCalculation()
        return sum_second_impact_area

    # Sum first impact-population classes.
    def sum_first_impact_pop(self):
        # Get temporary directory.
        tempdir = self.gettempdir()

        # Get layers already created in temporary directory.
        first_area_list, second_area_list, first_pop_list, second_pop_list = self.pre_siam()
        entries = []

        for i in range(0, len(first_pop_list)):
            raster = first_pop_list[i]
            readRst = processing.getObject(str(raster))
            ras1 = QgsRasterCalculatorEntry()
            ras1.raster = readRst
            ras1.ref = "first_impact_pop_" + str(i + 1) + "@1"
            ras1.bandNumber = 1
            entries.append(ras1)

        reflist = " + ".join([ent.ref for ent in entries])
        formula = '(' + reflist + ')'

        readRst = QgsRasterLayer(first_pop_list[0])
        sum_first_impact_pop = tempdir + "/sum_first_impact_pop.tif"
        calc = QgsRasterCalculator(formula, sum_first_impact_pop, 'GTiff', readRst.extent(), readRst.width(), readRst.height(),
                                   entries)
        calc.processCalculation()
        return sum_first_impact_pop

    # Sum second impact-population classes.
    def sum_second_impact_pop(self):
        # Get temporary directory.
        tempdir = self.gettempdir()

        # Get layers already created in temporary directory.
        first_area_list, second_area_list, first_pop_list, second_pop_list = self.pre_siam()
        entries = []

        for i in range(0, len(second_pop_list)):
            raster = second_pop_list[i]
            readRst = processing.getObject(str(raster))
            ras1 = QgsRasterCalculatorEntry()
            ras1.raster = readRst
            ras1.ref = "second_impact_pop_" + str(i + 1) + "@1"
            ras1.bandNumber = 1
            entries.append(ras1)

        reflist = " + ".join([ent.ref for ent in entries])
        formula = '(' + reflist + ')'

        readRst = QgsRasterLayer(second_pop_list[0])
        sum_second_impact_pop = tempdir + "/sum_second_impact_pop.tif"
        calc = QgsRasterCalculator(formula, sum_second_impact_pop, 'GTiff', readRst.extent(), readRst.width(), readRst.height(),
                                   entries)
        calc.processCalculation()
        return sum_second_impact_pop
    # Calculate siam.
    def siam(self):
        # Get temporary directory.
        tempdir = self.gettempdir()

        # Get weights
        area_weight = self.ui.lineEdit_1.text()
        pop_weight = self.ui.lineEdit_2.text()

        # Get layers
        sum_first_impact_area = self. sum_first_impact_area()
        sum_second_impact_area = self.sum_second_impact_area()
        sum_first_impact_pop = self.sum_first_impact_pop()
        sum_second_impact_pop = self.sum_second_impact_pop()

        # List layers for calculate.
        entries = []

        # First impact-area layer info.
        fileInfo_1a = QFileInfo(str(sum_first_impact_area))
        lyr_1a_Name = fileInfo_1a.baseName()
        impactlyr_1a = processing.getObject(str(sum_first_impact_area))
        first_ras_a = QgsRasterCalculatorEntry()
        first_ras_a.ref = lyr_1a_Name + '@1'
        first_ras_a.raster = impactlyr_1a
        first_ras_a.bandNumber = 1
        entries.append(first_ras_a)

        # Second impact-area layer info.
        fileInfo_2a = QFileInfo(str(sum_second_impact_area))
        lyr_2a_Name = fileInfo_2a.baseName()
        impactlyr_2a = processing.getObject(str(sum_second_impact_area))
        second_ras_a = QgsRasterCalculatorEntry()
        second_ras_a.ref = lyr_2a_Name + '@1'
        second_ras_a.raster = impactlyr_2a
        second_ras_a.bandNumber = 1
        entries.append(second_ras_a)

        # First impact-pop layer info.
        fileInfo_1p = QFileInfo(str(sum_first_impact_pop))
        lyr_1p_Name = fileInfo_1p.baseName()
        impactlyr_1p = processing.getObject(str(sum_first_impact_pop))
        first_ras_p = QgsRasterCalculatorEntry()
        first_ras_p.ref = lyr_1p_Name + '@1'
        first_ras_p.raster = impactlyr_1p
        first_ras_p.bandNumber = 1
        entries.append(first_ras_p)

        # Second impact-pop layer info.
        fileInfo_2p = QFileInfo(str(sum_second_impact_pop))
        lyr_2p_Name = fileInfo_2p.baseName()
        impactlyr_2p = processing.getObject(str(sum_second_impact_pop))
        second_ras_p = QgsRasterCalculatorEntry()
        second_ras_p.ref = lyr_2p_Name + '@1'
        second_ras_p.raster = impactlyr_2p
        second_ras_p.bandNumber = 1
        entries.append(second_ras_p)


        formula_siam = "(" + str(area_weight) + ' * ' + "("  + "\"" + first_ras_a.ref + "\"" + ' + ' + "\"" + second_ras_a.ref  + "\"" + "))" + ' + ' + "("+ str(pop_weight) + ' * ' + "("  + "\"" + first_ras_p.ref + "\"" + ' + ' + "\"" + second_ras_p.ref  + "\"" + "))"

        #print formula_siam

        output_siam = tempdir + "/SIAM.tif"
        calc_siam = QgsRasterCalculator(formula_siam, output_siam, 'GTiff', impactlyr_1a.extent(),
                                        impactlyr_1a.width(), impactlyr_1a.height(), entries)
        calc_siam.processCalculation()
        del entries

