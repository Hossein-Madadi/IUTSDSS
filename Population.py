from PyQt4 import QtGui
from PyQt4.QtGui import *
from PyQt4.QtCore import *
from qgis.core import *
from qgis.utils import *
import numpy as np
from osgeo import gdal, ogr, gdalconst
import os, platform, tempfile, pickle, processing, itertools, osr
from ui.ui_Pop_Dialog_base import Ui_PopClassDialogBase
from qgis.analysis import QgsRasterCalculator, QgsRasterCalculatorEntry


# Create the dialog.
class PopClassDialog(QtGui.QDialog):
    def __init__(self, iface):
        QtGui.QDialog.__init__(self)
        self.iface = iface

        # Set up the user interface from Designer.
        self.ui = Ui_PopClassDialogBase()
        self.ui.setupUi(self)
        self.ui.btnProcess.setDisabled(True)

        # Set up the signals.
        self.connectSignals()

    # Connections.
    def connectSignals(self):
        self.ui.cBox_Pop_1.currentIndexChanged.connect(self.chek_fields)
        self.ui.cBox_Rst_1.currentIndexChanged.connect(self.chek_fields)
        self.ui.linOutput.textChanged.connect(self.chek_fields)
        self.ui.pButton_Output.clicked.connect(self.select_output_file)
        self.ui.btnClose.clicked.connect(self.close)
        self.ui.btnProcess.clicked.connect(self.btnProcessClicked)

    # Get point layer into a list.
    def getPoints(self):
        layers = QgsMapLayerRegistry.instance().mapLayers().values()
        point_list = []
        for layer in layers:
            if layer.type() == QgsMapLayer.VectorLayer and (
                            layer.wkbType() == QGis.WKBPoint or layer.wkbType() == QGis.WKBMultiPoint):
                point_list.append(layer)
        return point_list

    # Get raster layers into a list.
    def getRasters(self):
        layers = QgsMapLayerRegistry.instance().mapLayers().values()
        raster_list = []
        for layer in layers:
            if layer.type() == QgsMapLayer.RasterLayer:
                raster_list.append(layer)
        return raster_list

    # Add Population field to combo box.
    def loadPopField(self):
        layers = self.getPoints()
        self.ui.cBox_Pop_3.clear()
        selectedLayerIndex = self.ui.cBox_Pop_1.currentIndex()
        selectedLayer = layers[selectedLayerIndex]
        for field in selectedLayer.pendingFields():
            if field.type() == QVariant.Int or field.type() == QVariant.Double or \
                            field.type() == QVariant.LongLong or field.type() == QVariant.ULongLong:
                self.ui.cBox_Pop_3.addItem(field.name(), field)  # lists layer fields

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
        raster_list = self.getRasters()
        point_list = self.getPoints()

        # clear the combos.
        self.ui.cBox_Pop_1.clear()
        self.ui.cBox_Pop_3.clear()

        self.ui.cBox_Rst_1.clear()
        self.ui.lineEdit.clear()

        # Clear fields.
        self.ui.linOutput.clear()

        # Add layers to combos.
        for layer in point_list:
            self.ui.cBox_Pop_1.addItem(layer.name(), layer)
            self.pntLayer = self.ui.cBox_Pop_1.itemData(self.ui.cBox_Pop_1.currentIndex())
            if self.pntLayer is not None:
                self.loadPopField()

        for layer in raster_list:
            self.ui.cBox_Rst_1.addItem(layer.name(), layer)

    # Check input fields.
    def chek_fields(self):
        # Defining layer input
        self.layer1 = self.ui.cBox_Pop_1.itemData(self.ui.cBox_Pop_1.currentIndex())
        self.field = self.ui.cBox_Pop_3.itemData(self.ui.cBox_Pop_3.currentIndex())
        self.layer2 = self.ui.cBox_Rst_1.itemData(self.ui.cBox_Rst_1.currentIndex())
        self.radius = self.ui.lineEdit.text()
        self.outPutName = self.ui.linOutput.text()

        if not self.outPutName:
            if not (self.layer1 and self.layer2 and self.field and self.radius):
                self.ui.btnProcess.setEnabled(False)
            elif (self.layer1 and self.layer2 and self.field and self.radius) is not None:
                self.ui.btnProcess.setEnabled(False)
        elif self.outPutName is not None:
            if not (self.layer1 and self.layer2 and self.field and self.radius):
                self.ui.btnProcess.setEnabled(False)
            elif (self.layer1 and self.layer2 and self.field and self.radius) is not None:
                if not self.is_number(self.radius):
                    self.iface.messageBar().pushMessage("Error", "Radius value not numeric",
                                                        level=QgsMessageBar.CRITICAL,
                                                        duration=5)
                else:
                    self.ui.btnProcess.setEnabled(True)
    # Get Point Layer to action.
    def point_layer(self):
        for layer in self.getPoints():
            if layer.name() == self.ui.cBox_Pop_1.currentText():
                self.vlayer = layer
                return self.vlayer

    # Get Rater Layer to action.
    def grid_layer(self):
        for layer in self.getRasters():
            if layer.name() == self.ui.cBox_Rst_1.currentText():
                self.rlayer = layer
                return self.rlayer

    # Calculate Mean distance from each point to all points in search radius.
    def mean_distance(self):
        # Get active layer
        layer = self.point_layer()
        self.radius = float(self.ui.lineEdit.text())

        # Get features
        feats = [feat for feat in layer.getFeatures()]

        # Number of features
        n = len(feats)

        # Get range for features
        comb = range(n)

        # Introduce variables
        distances = [[] for i in range(n)]
        indexes = [[] for i in range(n)]

        # Calculate distance between points less than 2500 meter
        for i, j in itertools.combinations(comb, 2):

            dist = feats[i].geometry().distance(feats[j].geometry())

            if dist < self.radius:
                i_dist = distances[i].append(dist)
                i_index = indexes[i].append([i, j])
                j_dist = distances[j].append(dist)
                j_index = indexes[j].append([i, j])

        prov = layer.dataProvider()

        # Check if attribute is already there, return "-1" if not create it
        if layer.fieldNameIndex("MeanDist") == -1:
            prov.addAttributes([QgsField("MeanDist", QVariant.Double, "double", 10, 2)])
            layer.updateFields()
        else:
            pass

        # Calculate mean distance in search radius and add it to attribute
        layer.startEditing()
        for feature in layer.getFeatures():
            attrName = "MeanDist"
            for i, group in enumerate(distances):
                if i == feature.id():
                    # print feature.id(), np.mean(group)   # print for test results, it must remove.
                    layer.changeAttributeValue(feature.id(), layer.fieldNameIndex(attrName), float(np.mean(group)))

        layer.commitChanges()

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

    # Called when "Process" button pressed.
    def btnProcessClicked(self):
        self.checkGeodata()

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

    # Create Population map in raster.
    def population_mapping(self):
        tempdir = self.gettempdir()
        vlayer = self.point_layer()
        rlayer = self.grid_layer()
        radius = float(self.ui.lineEdit.text())
        pop_field = self.ui.cBox_Pop_3.currentText()

        def pixel2coord(x, y):
            xp = (pixelWidth * x) + originX + (pixelWidth / 2)
            yp = (pixelHeight * y) + originY + (pixelHeight / 2)
            return QgsPoint(xp, yp)

        # Get Point layer features
        feats = [feat for feat in vlayer.getFeatures()]

        # Get Pixel Size
        pixelWidth = rlayer.rasterUnitsPerPixelX()
        pixelHeight = rlayer.rasterUnitsPerPixelY()

        # extent of the Raster layer
        ext = rlayer.extent()

        originX, originY = (ext.xMinimum(), ext.yMinimum())
        src_cols = rlayer.width()
        src_rows = rlayer.height()

        filepath = rlayer.source()
        drive = gdal.GetDriverByName('GTiff')
        src_ds = gdal.Open(filepath)
        outBand = src_ds.GetRasterBand(1)

        pntRstList = []

        provider = rlayer.dataProvider()
        block = provider.block(1, ext, src_cols, src_rows)

        # Get X, Y pixel centers and then gey them values, if value is 1 then add X, Y
        # to a list.
        for i in range(0, src_cols):
            for j in range(0, src_rows):
                # if block.value(i,j) == 1:
                rspnt = pixel2coord(i, j)
                # rspnt = Pixel2world(geoTrans, i, j)
                pntRstList.append(rspnt)

        idxMeanDist = vlayer.fieldNameIndex("MeanDist")
        idxPopulaton = vlayer.fieldNameIndex(pop_field)

        # dictionary for store population data for each pixel based on points in its
        dictDist = {}
        dictPop = {}

        for rpoint in pntRstList:
            for ft in feats:
                vgeometry = ft.geometry()
                rgeometry = QgsGeometry.fromPoint(QgsPoint(rpoint[0], rpoint[1]))
                dist = vgeometry.distance(rgeometry)  # get distance between cell center to point
                # print ft.id(), rpoint, dist
                if dist < radius:
                    weightPop = (idxMeanDist ** 2 - dist ** 2) / (
                        idxMeanDist ** 2 + dist ** 2)  # abs((meanDist**2 - dist**2)/(meanDist**2 + dist**2))
                    pntPop = ft.attributes()[idxPopulaton]
                    keyR = rpoint
                    dictDist.setdefault(keyR, [])
                    dictDist[keyR].append(weightPop)
                    dictPop.setdefault(keyR, [])
                    dictPop[keyR].append(pntPop)
                else:
                    dictDist.setdefault(rpoint, [])
                    dictDist[rpoint].append(0)
                    dictPop.setdefault(rpoint, [])
                    dictPop[rpoint].append(0)

        # Sum all weight values for each population point
        bb = [sum(value) for value in zip(*dictDist.values())]
        # Divide weight for each pixel to sum of weights of each point population
        cc = {key: [0.0 if not y else ((x * 1.0) / y) for x, y in zip(value, bb)] for key, value in dictDist.items()}

        # Multiply cc to dictPop
        # Compare to dictionaris (cc and dictPop), because Keys in tow dictionaries
        # are same so values (in cc dictionary is scaled weights and in dictPop
        # dictionary is population)
        dd = {key: [0.0 if not y else x * y for x, y in zip(value, dictPop[key])] for key, value in cc.items() if
              key in dictPop}
        # self.userWarning("Print", str(dd))

        # Sum values for each key in dd, This final Value that shows population for each
        # pixel. Note: Sum population values for all pixels must be equal to all population.
        finalPop = {}
        for key, value in (dd.items()):
            finalPop.setdefault(key, [])
            finalPop[key].append(sum(value))

        ###############################################################################
        # Create memory layer based on cell centers for population calculated.
        mergeList = [[point, finalPop[point]] for point in finalPop]
        points = [QgsPoint(item[0][0], item[0][1]) for item in mergeList]
        values = [item[1] for item in mergeList]
        new_values = []

        for item in values:
            if type(item) == list:
                new_values.append(item[0])
            else:
                new_values.append(item)

        epsg = vlayer.crs().postgisSrid()
        uri = "Point?crs=epsg:" + str(epsg) + \
              "&field=id:integer&field=value:real""&index=yes"
        mem_layer = QgsVectorLayer(uri,
                                   'points',
                                   'memory')

        prov = mem_layer.dataProvider()
        feats = [QgsFeature() for i in range(len(points))]

        for i, feat in enumerate(feats):
            feat.setAttributes([i, new_values[i]])
            feat.setGeometry(QgsGeometry.fromPoint(points[i]))

        prov.addFeatures(feats)
        #######################################################################################
        # Save mem_layer
        pnt_temp = tempdir + "/Pop_pnt_temp.shp"
        crs = QgsCoordinateReferenceSystem("epsg:" + str(epsg))
        outLayer = QgsVectorFileWriter.writeAsVectorFormat(mem_layer,
                                                           pnt_temp,
                                                           "UTF-8",
                                                           crs,
                                                           "ESRI Shapefile")

        #######################################################################################
        # Create Final Population map in raster format.
        #### Define pixel_size and NoData value of new raster
        cellsize = pixelWidth
        NoData_value = -9999
        x_res = cellsize  # assuming these are the cell sizes
        y_res = cellsize  # change as appropriate

        # Chosen output raster extent from input model extent layer
        extLayer_extent = rlayer.extent()

        x_min = extLayer_extent.xMinimum()
        x_max = extLayer_extent.xMaximum()
        y_min = extLayer_extent.yMinimum()
        y_max = extLayer_extent.yMaximum()

        # Chosen Shape layer with GDAL
        pnt_ds = ogr.Open(pnt_temp)  #### This Line is very important for rasterize of point layer!!!!
        source_layer = pnt_ds.GetLayer()

        # Output file name
        rasterfileOutputName = self.ui.linOutput.text()  # tempdir + "/Population.tif"

        # Create Target - TIFF
        srs = osr.SpatialReference()
        srs.ImportFromWkt(rlayer.crs().toWkt())
        cols = int((x_max - x_min) / x_res)
        rows = int((y_max - y_min) / y_res)
        drive = gdal.GetDriverByName('GTiff')
        out_rst = drive.Create(rasterfileOutputName, cols, rows, 1, gdal.GDT_Float32)
        out_rst.SetGeoTransform((x_min, x_res, 0, y_max, 0, -y_res))
        out_rst.SetProjection(srs.ExportToWkt())
        band_rst = out_rst.GetRasterBand(1)
        band_rst.SetNoDataValue(NoData_value)
        band_rst.FlushCache()

        # Rasterize why is the burn value 0... isn't that the same as the background?
        gdal.RasterizeLayer(out_rst, [1], source_layer, options=["ATTRIBUTE=value"])

        # close bands
        band_rst = None

        # close rasters
        out_rst = None

        # Show new layer
        raster = QgsRasterLayer(rasterfileOutputName, 'Population')
        QgsMapLayerRegistry.instance().addMapLayer(raster)


