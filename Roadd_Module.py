# -*- coding: utf-8 -*-
"""
/***************************************************************************
 EviClass
                                 A QGIS plugin
 Environmental Vulnerability
                              -------------------
        begin                : 2017-08-23
        git sha              : $Format:%H$
        copyright            : (C) 2017 by Hossein Madadi
        email                : hosein.madadi@gmail.com
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
"""
from PyQt4.QtCore import *
from PyQt4.QtGui import *
from qgis.core import *
import resources
from SDSSRoad.dialogs.Evi_Module_dialog import EviClassDialog
import Vulnerability, Siam, Population, Footprint, Synthesis
import os.path, platform, sys, tempfile, shutil


class roadd_menu:
    def __init__(self, iface):
        self.iface = iface
        self.roadd_menu = None

        # Make folder for intermediate files
        tempplatform = tempfile.gettempdir()
        platform_type = platform.system()

        if platform_type == 'Linux':
            tempdir_v = str(tempplatform) + '/tmpSDSSRoad/vulnerability'
            if os.path.exists(tempdir_v):
                shutil.rmtree(tempdir_v)
            os.makedirs(tempdir_v)

            tempdir_p = str(tempplatform) + '/tmpSDSSRoad/population'
            if os.path.exists(tempdir_p):
                shutil.rmtree(tempdir_p)
            os.makedirs(tempdir_p)

            tempdir_s = str(tempplatform) + '/tmpSDSSRoad/siam'
            if os.path.exists(tempdir_s):
                shutil.rmtree(tempdir_s)
            os.makedirs(tempdir_s)

            tempdir_t = str(tempplatform) + '/tmpSDSSRoad/synthesis'
            if os.path.exists(tempdir_t):
                shutil.rmtree(tempdir_t)
            os.makedirs(tempdir_t)
        elif platform_type == 'Windows':
            tempdir_v = str(tempplatform) + '/tmpSDSSRoad/vulnerability'
            if os.path.exists(tempdir_v):
                shutil.rmtree(tempdir_v)
            os.makedirs(tempdir_v)

            tempdir_p = str(tempplatform) + '/tmpSDSSRoad/population'
            if os.path.exists(tempdir_p):
                shutil.rmtree(tempdir_p)
            os.makedirs(tempdir_p)

            tempdir_s = str(tempplatform) + '/tmpSDSSRoad/siam'
            if os.path.exists(tempdir_s):
                shutil.rmtree(tempdir_s)
            os.makedirs(tempdir_s)

            tempdir_t = str(tempplatform) + '/tmpSDSSRoad/synthesis'
            if os.path.exists(tempdir_t):
                shutil.rmtree(tempdir_t)
            os.makedirs(tempdir_t)
        elif platform_type == 'Darwin':
            tempdir_v = str(tempplatform) + '/tmpSDSSRoad/vulnerability'
            if os.path.exists(tempdir_v):
                shutil.rmtree(tempdir_v)
            os.makedirs(tempdir_v)

            tempdir_p = str(tempplatform) + '/tmpSDSSRoad/population'
            if os.path.exists(tempdir_p):
                shutil.rmtree(tempdir_p)
            os.makedirs(tempdir_p)

            tempdir_s = str(tempplatform) + '/tmpSDSSRoad/siam'
            if os.path.exists(tempdir_s):
                shutil.rmtree(tempdir_s)
            os.makedirs(tempdir_s)

            tempdir_t = str(tempplatform) + '/tmpSDSSRoad/synthesis'
            if os.path.exists(tempdir_t):
                shutil.rmtree(tempdir_t)
            os.makedirs(tempdir_t)
        else:
            iface.messageBar().pushInfo(u'There is not temporary directory in your platform system')

    def roadd_add_submenu(self, submenu):
        if self.roadd_menu != None:
            self.roadd_menu.addMenu(submenu)
        else:
            self.iface.addPluginToMenu("&roadd", submenu.menuAction())

    def initGui(self):
        # Uncomment the following two lines to have ROADD accessible from a top-level menu
        self.roadd_menu = QMenu(QCoreApplication.translate("roadd", "SDSSRoad"))
        self.iface.mainWindow().menuBar().insertMenu(self.iface.firstRightStandardMenu().menuAction(), self.roadd_menu)

        # Vulnerability Submenue
        self.vulnerability_menu = QMenu(QCoreApplication.translate("roadd", "&Vulnerability"))
        self.roadd_add_submenu(self.vulnerability_menu)

        icon = QIcon(os.path.dirname(__file__) + "/icons/mmqgis_animate_columns.png")
        self.vulnerability_action = QAction(icon, "Vulnerability", self.iface.mainWindow())
        QObject.connect(self.vulnerability_action, SIGNAL("triggered()"), self.vulnerability)
        self.vulnerability_menu.addAction(self.vulnerability_action)

        # SIAM Submenu
        self.siam_menu = QMenu(QCoreApplication.translate("roadd", "&SIAM"))
        self.roadd_add_submenu(self.siam_menu)

        icon = QIcon(os.path.dirname(__file__) + "/icons/mmqgis_attribute_join.png")
        self.population_action = QAction(icon, "Population", self.iface.mainWindow())
        QObject.connect(self.population_action, SIGNAL("triggered()"), self.population)
        self.siam_menu.addAction(self.population_action)

        icon = QIcon(os.path.dirname(__file__) + "/icons/mmqgis_merge.png")
        self.siam_action = QAction(icon, "SIAM", self.iface.mainWindow())
        QObject.connect(self.siam_action, SIGNAL("triggered()"), self.siam)
        self.siam_menu.addAction(self.siam_action)

        # Sustainability
        self.sustainability_menu = QMenu(QCoreApplication.translate("roadd", "&Sustainability"))
        self.roadd_add_submenu(self.sustainability_menu)

        icon = QIcon(os.path.dirname(__file__) + "/icons/mmqgis_buffers.png")
        self.sustainability_action = QAction(icon, "Ecological Footprint", self.iface.mainWindow())
        QObject.connect(self.sustainability_action, SIGNAL("triggered()"), self.sustainability)
        self.sustainability_menu.addAction(self.sustainability_action)

        # Synthesis
        self.synthesis_menu = QMenu(QCoreApplication.translate("roadd", "&Synthesis"))
        self.roadd_add_submenu(self.synthesis_menu)

        icon = QIcon(os.path.dirname(__file__) + "/icons/mmqgis_buffers.png")
        self.synthesis_action = QAction(icon, "Synthesis", self.iface.mainWindow())
        QObject.connect(self.synthesis_action, SIGNAL("triggered()"), self.synthesis)
        self.synthesis_menu.addAction(self.synthesis_action)

        # Help
        self.help_menu = QMenu(QCoreApplication.translate("roadd", "&Help"))
        self.roadd_add_submenu(self.help_menu)

        # About
        self.about_menu = QMenu(QCoreApplication.translate("roadd", "&About"))
        self.roadd_add_submenu(self.about_menu)

    def unload(self):
        if self.roadd_menu != None:
            self.iface.mainWindow().menuBar().removeAction(self.roadd_menu.menuAction())
        else:
            self.iface.removePluginMenu("&roadd", self.vulnerability_menu.menuAction())
            self.iface.removePluginMenu("&roadd", self.siam_menu.menuAction())
            self.iface.removePluginMenu("&roadd", self.sustainability_menu.menuAction())
            self.iface.removePluginMenu("&roadd", self.synthesis_menu.menuAction())
            self.iface.removePluginMenu("&roadd", self.help_menu.menuAction())
            self.iface.removePluginMenu("&roadd", self.about_menu.menuAction())

    def vulnerability(self):
        # show the dialog
        self.dlg = Vulnerability.EviClassDialog(self.iface)
        self.dlg.update_fields()
        self.dlg.show()

        # Run the dialog event loop
        result = self.dlg.exec_()

        # See if OK was pressed
        if result == 1:
            # do something useful (delete the line containing pass and
            # substitute with your code)
            Vulnerability.self.btnProcessClicked()

    def siam(self):
        # show the dialog
        self.dlg = Siam.SiamClassDialog(self.iface)
        self.dlg.update_fields()
        self.dlg.show()

        # Run the dialog event loop
        result = self.dlg.exec_()

        # See if OK was pressed
        if result == 1:
            # do something useful (delete the line containing pass and
            # substitute with your code)
            Siam.self.btnProcessClicked()

    def population(self):
        # show the dialog
        self.dlg = Population.PopClassDialog(self.iface)
        self.dlg.update_fields()
        self.dlg.show()

        # Run the dialog event loop
        result = self.dlg.exec_()

        # See if OK was pressed
        if result == 1:
            # do something useful (delete the line containing pass and
            # substitute with your code)
            Population.self.btnProcessClicked()

    def sustainability(self):

        # show the dialog
        self.dlg = Footprint.FootprintClassDialog(self.iface)
        self.dlg.update_fields()
        self.dlg.show()

        # Run the dialog event loop
        result = self.dlg.exec_()

        # See if OK was pressed
        if result == 1:
            # do something useful (delete the line containing pass and
            # substitute with your code)
            Footprint.self.btnProcessClicked()

    def synthesis(self):

        # show the dialog
        self.dlg = Synthesis.SynthesisClassDialog(self.iface)
        self.dlg.update_fields()
        self.dlg.show()

        # Run the dialog event loop
        result = self.dlg.exec_()

        # See if OK was pressed
        if result == 1:
            # do something useful (delete the line containing pass and
            # substitute with your code)
            Synthesis.self.btnProcessClicked()