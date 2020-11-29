from PyQt4 import QtGui
from PyQt4.QtGui import *
from qgis.core import *
from qgis.analysis import QgsRasterCalculator, QgsRasterCalculatorEntry
from qgis.utils import *
from ui.ui_Synthesis_dialog_base import Ui_SynClassDialogBase
import collections
from osgeo import gdal, osr
import numpy as np
import platform, tempfile
import struct

# create the dialog.
class SynthesisClassDialog(QtGui.QDialog):
    def __init__(self, iface):
        QtGui.QDialog.__init__(self)
        self.iface = iface

        # Set up the user interface from Designer.
        self.ui = Ui_SynClassDialogBase()
        self.ui.setupUi(self)
        self.ui.btnProcess.setDisabled(True)

        # Set up the signals.
        self.connectSignals()

    # Connections.
    def connectSignals(self):
        self.ui.cBox_Vulnerability.currentIndexChanged.connect(self.chek_fields)
        self.ui.cBox_Impact.currentIndexChanged.connect(self.chek_fields)
        self.ui.cBox_Footprint.currentIndexChanged.connect(self.chek_fields)
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
        self.ui.cBox_Vulnerability.clear()
        self.ui.cBox_Impact.clear()
        self.ui.cBox_Footprint.clear()

        # Clear fields.
        self.ui.linOutput.clear()

        # Add layers to combos.
        for layer in raster_list:
            self.ui.cBox_Vulnerability.addItem(layer.name(), layer)
            self.ui.cBox_Impact.addItem(layer.name(), layer)
            self.ui.cBox_Footprint.addItem(layer.name(), layer)

    # Check input fields.
    def chek_fields(self):
        # Defining layer input
        self.layer1 = self.ui.cBox_Vulnerability.itemData(self.ui.cBox_Vulnerability.currentIndex())
        self.layer2 = self.ui.cBox_Impact.itemData(self.ui.cBox_Impact.currentIndex())
        self.layer3 = self.ui.cBox_Footprint.itemData(self.ui.cBox_Footprint.currentIndex())

        # Make a list from raster layers in combo boxes.
        lst = []
        lst.extend([self.layer1, self.layer2, self.layer3])

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
        self.layer1 = self.ui.cBox_Vulnerability.currentText()
        self.layer2 = self.ui.cBox_Impact.currentText()
        self.layer3 = self.ui.cBox_Footprint.currentText()

        # Make a list from raster layers in combo boxes.
        lstLyr = []
        lstLyr.extend(
            [self.layer1, self.layer2, self.layer3])

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
                self.synthesis()    #############################################################################

    # Check for repeated item in a list.
    def all_same(self,items):
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

    # Called when "Process" button pressed.
    def btnProcessClicked(self):
        self.checkGeodata()

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

    # Create temporary folder for intermediate layers.
    def gettempdir(self):
        platform_type = platform.system()
        tempplatform = tempfile.gettempdir()
        tempdir = ''
        if platform_type == 'Linux':
            tempdir = str(tempplatform) + '/tmpSDSSRoad/synthesis'
        elif platform_type == 'Windows':
            tempdir = str(tempplatform) + '/tmpSDSSRoad/synthesis'
        elif platform_type == 'Darwin':
            tempdir = str(tempplatform) + '/tmpSDSSRoad/synthesis'
        return tempdir

    # Standarization All input raster layers.
    # Divide a raster layer to sum of its all pixel values (execute for each index - Yij maps).
    def normalization(self):
        tempdir = self.gettempdir()
        lst = self.chek_fields()
        lst_normalized = []

        for layer in lst[0:1]:
            lyrPath = layer.source()
            lyrName = layer.name()
            lyr2 = QgsRasterLayer(lyrPath, lyrName)

            # Get Min & Max values.
            extent = lyr2.extent()
            provider = lyr2.dataProvider()
            for band in range(1, layer.bandCount() + 1):
                stats = provider.bandStatistics(band, QgsRasterBandStats.All, extent, 0)
                minV = stats.minimumValue
                maxV = stats.maximumValue

            # Get raster info.
            entries = []
            ras1 = QgsRasterCalculatorEntry()
            ras1.ref = lyrName + '@1'
            ras1.raster = lyr2
            ras1.bandNumber = 1
            entries.append(ras1)

            formulaV = "(" + "\"" + ras1.ref + "\"" + ' - ' + str(minV) + ")" + ' / ' + "(" + str(maxV) + ' - ' + str(minV) + ")"
            output = tempdir + "/Stnd_V_%s.tif" % str(lyrName)
            lst_normalized.append(str(output))
            calc = QgsRasterCalculator(formulaV, output, 'GTiff', lyr2.extent(), lyr2.width(), lyr2.height(), entries)
            calc.processCalculation()

        for layer in lst[1:2]:
            lyrPath = layer.source()
            lyrName = layer.name()
            lyr2 = QgsRasterLayer(lyrPath, lyrName)

            # Get Min & Max values.
            extent = lyr2.extent()
            provider = lyr2.dataProvider()
            for band in range(1, layer.bandCount() + 1):
                stats = provider.bandStatistics(band, QgsRasterBandStats.All, extent, 0)
                minI = stats.minimumValue
                maxI = stats.maximumValue

            entries = []
            ras1 = QgsRasterCalculatorEntry()
            ras1.ref = lyrName + '@1'
            ras1.raster = lyr2
            ras1.bandNumber = 1
            entries.append(ras1)

            formulaI = "(" + "\"" + ras1.ref + "\"" + ' - ' + str(minI) + ")" + ' / ' + "(" + str(maxI) + ' - ' + str(minI) + ")"
            output = tempdir + "/Stnd_I_%s.tif" % str(lyrName)
            lst_normalized.append(str(output))
            calc = QgsRasterCalculator(formulaI, output, 'GTiff', lyr2.extent(), lyr2.width(), lyr2.height(), entries)
            calc.processCalculation()

        for layer in lst[2:]:
            lyrPath = layer.source()
            lyrName = layer.name()
            lyr2 = QgsRasterLayer(lyrPath, lyrName)

            # Get Min & Max values.
            extent = lyr2.extent()
            provider = lyr2.dataProvider()
            for band in range(1, layer.bandCount() + 1):
                stats = provider.bandStatistics(band, QgsRasterBandStats.All, extent, 0)
                minF = stats.minimumValue
                maxF = stats.maximumValue

            entries = []
            ras1 = QgsRasterCalculatorEntry()
            ras1.ref = lyrName + '@1'
            ras1.raster = lyr2
            ras1.bandNumber = 1
            entries.append(ras1)

            formulaF = "(" + "\"" + ras1.ref + "\"" + ' - ' + str(minF) + ")" + ' / ' + "(" + str(maxF) + ' - ' + str(minF) + ")"
            output = tempdir + "/Stnd_F_%s.tif" % str(lyrName)
            lst_normalized.append(str(output))
            calc = QgsRasterCalculator(formulaF, output, 'GTiff', lyr2.extent(), lyr2.width(), lyr2.height(), entries)
            calc.processCalculation()

        return lst_normalized

    # Classify normalized raster layer based on quantile classification.
    def quantile(self):
        tempdir = self.gettempdir()
        lst_norm = self.normalization()
        lst_quantiled = []

        for layer in lst_norm:
            dataset = gdal.Open(str(layer))
            lyrName = os.path.basename(layer)
            band = dataset.GetRasterBand(1)
            nodata = band.GetNoDataValue()
            array = dataset.ReadAsArray()

            new_array = array
            nan_array = array

            nan_array[array == nodata] = np.nan

            percentile_80 = np.nanpercentile(nan_array, 80)
            percentile_60 = np.nanpercentile(nan_array, 60)
            percentile_40 = np.nanpercentile(nan_array, 40)
            percentile_20 = np.nanpercentile(nan_array, 20)
            percentile_0 = np.nanpercentile(nan_array, 0)

            for i, v in enumerate(new_array):
                for j, element in enumerate(v):
                    if element <= percentile_0:
                        new_array[i, j] = 1
                    if element > percentile_0 and element <= percentile_20:
                        new_array[i, j] = 2
                    if element > percentile_20 and element <= percentile_40:
                        new_array[i, j] = 3
                    if element > percentile_40 and element <= percentile_60:
                        new_array[i, j] = 4
                    if element > percentile_60 and element <= percentile_80:
                        new_array[i, j] = 5
                    if element > percentile_80:
                        new_array[i, j] = 6

            new_array[new_array != new_array] = nodata

            geotransform = dataset.GetGeoTransform()
            wkt = dataset.GetProjection()

            # Create gtif file
            driver = gdal.GetDriverByName("GTiff")
            output_file = tempdir + "/quantile_%s.tif" % str(lyrName)
            lst_quantiled.append(str(output_file))
            #output_file = "/tmp/tmpSDSSRoad/synthesis/Buffer_quantile_2.tif"

            dst_ds = driver.Create(output_file,
                                   band.XSize,
                                   band.YSize,
                                   1,
                                   gdal.GDT_Int16)

            # writting output raster
            dst_ds.GetRasterBand(1).WriteArray(new_array)
            # setting nodata value
            dst_ds.GetRasterBand(1).SetNoDataValue(nodata)
            # setting extension of output raster
            # top left x, w-e pixel resolution, rotation, top left y, rotation, n-s pixel resolution
            dst_ds.SetGeoTransform(geotransform)
            # setting spatial reference of output raster
            srs = osr.SpatialReference()
            srs.ImportFromWkt(wkt)
            dst_ds.SetProjection(srs.ExportToWkt())
            # Close output raster dataset

            dataset = None
            dst_ds = None

        return lst_quantiled

    # Get three layers of quantile output then multiply first in 100, second in 10 and third in 1, then sum three new
    # mltiplied layers. Output layer showes three digits for each pixel. First digit showes Vulnerability class, second
    # digit showes Impact class and third digit showes Footprint class.
    def sum_quantiles(self):
        tempdir = self.gettempdir()
        lst_quantile = self.quantile()

        # List layers for calculate.
        entries = []

        for lyrPath in lst_quantile[0:1]:
            lyrName = os.path.basename(lyrPath)
            lyr1 = QgsRasterLayer(lyrPath, lyrName)

            # Get raster info.
            ras1 = QgsRasterCalculatorEntry()
            ras1.ref = lyrName + '@1'
            ras1.raster = lyr1
            ras1.bandNumber = 1
            entries.append(ras1)

        for lyrPath in lst_quantile[1:2]:
            lyrName = os.path.basename(lyrPath)
            lyr2 = QgsRasterLayer(lyrPath, lyrName)

            # Get raster info.
            ras2 = QgsRasterCalculatorEntry()
            ras2.ref = lyrName + '@1'
            ras2.raster = lyr2
            ras2.bandNumber = 1
            entries.append(ras2)

        for lyrPath in lst_quantile[2:3]:
            lyrName = os.path.basename(lyrPath)
            lyr3 = QgsRasterLayer(lyrPath, lyrName)

            # Get raster info.
            ras3 = QgsRasterCalculatorEntry()
            ras3.ref = lyrName + '@1'
            ras3.raster = lyr3
            ras3.bandNumber = 1
            entries.append(ras3)

        formula_sum = "((" + "\"" + ras1.ref + "\"" + '*' + '100' + ')' + ' + ' + '(' + "\"" + ras2.ref + "\"" + '*' + '10' + ")" + ' + ' + "(" + "\"" + ras3.ref + "\"" + '))'

        output = tempdir + "/VIF.tif"
        calc_sum = QgsRasterCalculator(formula_sum, output, 'GTiff', lyr1.extent(),
                                       lyr1.width(), lyr1.height(), entries)
        calc_sum.processCalculation()
        del entries
        return output

    # Sum the digits of a umber.
    def sum_digits(self, n):
            s = 0
            while n:
                s += n % 10
                n //= 10
            return s

    # Change pixel value to sum of the digits.
    def changeRasterValues(self, band):
        fmttypes = {'Byte': 'B', 'UInt16': 'H', 'Int16': 'h', 'UInt32': 'I', 'Int32': 'i', 'Float32': 'f',
                    'Float64': 'd'}

        data_type = band.DataType

        BandType = gdal.GetDataTypeName(band.DataType)

        raster = []

        for y in range(band.YSize):
            scanline = band.ReadRaster(0, y, band.XSize, 1, band.XSize, 1, data_type)
            values = struct.unpack(fmttypes[BandType] * band.XSize, scanline)
            raster.append(values)

        raster = [list(item) for item in raster]

        # changing raster values
        for i, item in enumerate(raster):
            for j, value in enumerate(item):
                if value >= 0:
                    sumdigits = self.sum_digits(value)
                    raster[i][j] = sumdigits
                elif value < 0:
                    raster[i][j] = -9999

        # transforming list in array
        raster = np.asarray(raster)

        return raster

    # Sum three digits of each pixel value in VIF layer.
    def sum_digits_vif(self):
        tempdir = self.gettempdir()
        layer = self.sum_quantiles()
        dataset = gdal.Open(str(layer))

        # Get projection
        prj = dataset.GetProjection()

        # setting band
        number_band = 1

        band = dataset.GetRasterBand(number_band)

        # Get raster metadata
        geotransform = dataset.GetGeoTransform()

        # Set name of output raster
        output_file = tempdir + "/sum_digits_VIF.tif"

        # Create gtif file with rows and columns from parent raster
        driver = gdal.GetDriverByName("GTiff")

        raster = self.changeRasterValues(band)

        dst_ds = driver.Create(output_file,
                               band.XSize,
                               band.YSize,
                               number_band,
                               band.DataType)

        # writting output raster
        dst_ds.GetRasterBand(number_band).WriteArray(raster)

        # setting extension of output raster
        # top left x, w-e pixel resolution, rotation, top left y, rotation, n-s pixel resolution
        dst_ds.SetGeoTransform(geotransform)

        # setting spatial reference of output raster
        srs = osr.SpatialReference(wkt=prj)
        dst_ds.SetProjection(srs.ExportToWkt())

        # Close output raster dataset
        dst_ds = None

        # Close main raster dataset
        dataset = None
        return output_file

    # Classify sum_digits_VIF.tif in 5 classes.
    def classify_result(self, band):
        fmttypes = {'Byte': 'B', 'UInt16': 'H', 'Int16': 'h', 'UInt32': 'I', 'Int32': 'i', 'Float32': 'f',
                    'Float64': 'd'}

        data_type = band.DataType

        BandType = gdal.GetDataTypeName(band.DataType)

        raster = []

        for y in range(band.YSize):
            scanline = band.ReadRaster(0, y, band.XSize, 1, band.XSize, 1, data_type)
            values = struct.unpack(fmttypes[BandType] * band.XSize, scanline)
            raster.append(values)

        raster = [list(item) for item in raster]

        # changing raster values
        for i, item in enumerate(raster):
            for j, value in enumerate(item):
                if 0 <= value <= 6:
                    raster[i][j] = 1
                elif 6 < value <= 9:
                    raster[i][j] = 2
                elif 9 < value <= 12:
                    raster[i][j] = 3
                elif 12 < value <= 15:
                    raster[i][j] = 4
                elif value > 15:
                    raster[i][j] = 5
                elif value < 0:
                    raster[i][j] = -9999

        # transforming list in array
        raster = np.asarray(raster)

        return raster

    # Classify and save final synthesis map.
    def synthesis(self):
        tempdir = self.gettempdir()
        layer = self.sum_digits_vif()
        dataset = gdal.Open(str(layer))

        # Get projection
        prj = dataset.GetProjection()

        # setting band
        number_band = 1

        band = dataset.GetRasterBand(number_band)

        # Get raster metadata
        geotransform = dataset.GetGeoTransform()

        # Set name of output raster
        output_file = tempdir + "/Synthesis.tif"

        # Create gtif file with rows and columns from parent raster
        driver = gdal.GetDriverByName("GTiff")

        raster = self.classify_result(band)

        dst_ds = driver.Create(output_file,
                               band.XSize,
                               band.YSize,
                               number_band,
                               band.DataType)

        # writting output raster
        dst_ds.GetRasterBand(number_band).WriteArray(raster)

        # setting extension of output raster
        # top left x, w-e pixel resolution, rotation, top left y, rotation, n-s pixel resolution
        dst_ds.SetGeoTransform(geotransform)

        # setting spatial reference of output raster
        srs = osr.SpatialReference(wkt=prj)
        dst_ds.SetProjection(srs.ExportToWkt())

        # Close output raster dataset
        dst_ds = None

        # Close main raster dataset
        dataset = None
