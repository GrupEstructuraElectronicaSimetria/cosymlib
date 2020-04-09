#!/usr/bin/python3
import sys
from os.path import expanduser
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5 import QtOpenGL
from PyQt5.Qt3DExtras import Qt3DWindow, QFirstPersonCameraController
from PyQt5.Qt3DCore import QEntity
from cosymlib.shape import tools
from cosymlib import file_io, shape
import numpy as np

home = expanduser("~")


class Execute(QPushButton):
    def __init__(self, title, parent):
        super().__init__(title, parent)

    def mousePressEvent(self, e):

        super().mousePressEvent(e)

        if e.button() == Qt.LeftButton:

                reply = QMessageBox.question(self, 'Message',
                                             "Run?", QMessageBox.Yes |
                                             QMessageBox.No, QMessageBox.No)

                if reply == QMessageBox.Yes:
                    quit()


class Exit(QPushButton):
    def __init__(self, title, parent):
        super().__init__(title, parent)

    def mousePressEvent(self, e):

        super().mousePressEvent(e)

        if e.button() == Qt.LeftButton:

                reply = QMessageBox.question(self, 'Message',
                                             "Are you sure to quit?", QMessageBox.Yes |
                                             QMessageBox.No, QMessageBox.No)

                if reply == QMessageBox.Yes:
                    quit()


class MainWindow(QWidget):

    def __init__(self):
        super().__init__()

        self.filename = None
        self.directory = None
        self.central_atom = 0
        self.initUI()

    def initUI(self):

        x = 500
        y = 300
        self.resize(x, y)
        self.center()

        # Input file
        file = QPushButton('Input file', self)
        file.move(x-100, -4)
        file.clicked.connect(self.open)
        self.file_text = QLineEdit(self)
        self.file_text.move(2, 2)
        self.file_text.resize(x - 150, 20)

        self.save_button = QPushButton('Save file', self)
        self.save_button.move(x-100, 20)
        self.save_button.setEnabled(False)
        self.save_button.clicked.connect(self.save)

        self.savefile_text = QLineEdit(self)
        self.savefile_text.move(2, 25)
        self.savefile_text.resize(x - 150, 20)

        # Label selection
        self.central_atom_name = QCheckBox('Central atom', self)
        self.central_atom_name.move(30, 60)
        self.central_atom_name.toggle()
        self.central_atom_name.stateChanged.connect(self.changeVertices)
        self.central_atom_box = QLineEdit(self)
        self.central_atom_box.move(150, 60)
        self.central_atom_box.resize(30, 30)
        self.central_atom_box.textChanged.connect(self.set_central_atom)

        labels_name = QLabel('Label :', self)
        labels_name.move(200, 60)
        self.labels_box = QComboBox(self)
        self.labels_box.move(250, 60)
        self.label = str(self.labels_box.currentText())
        self.labels_box.activated[str].connect(self.labels_on)

        self.results = QPlainTextEdit(self)
        self.results.setReadOnly(True)
        self.results.setStyleSheet("""QPlainTextEdit {background-color: #333; color: #00FF00; font-family: Courier;}""")
        self.results.move(30, 100)
        self.results.resize(250, 175)

        execute = QPushButton("Run")
        exit = Exit("Exit", self)

        hbox = QHBoxLayout()
        hbox.addStretch(1)
        hbox.addWidget(execute)
        hbox.addWidget(exit)

        vbox = QVBoxLayout()
        vbox.addStretch(1)
        vbox.addLayout(hbox)

        self.setLayout(vbox)
        self.setWindowTitle('SHAPE')
        execute.clicked.connect(self.execute)

    def changeVertices(self, state):

        self.labels_box.clear()
        if state:
            self.central_atom_box.setReadOnly(False)
            if self.filename:
                self.labels_on(tools.get_structure_references(self.geometries[0].get_n_atoms() - 1))
        elif not state:
            # self.central_atom_box.clear()
            self.set_central_atom(0)
            self.central_atom_box.setReadOnly(True)
            if self.filename:
                self.labels_on(tools.get_structure_references(self.geometries[0].get_n_atoms()))

    def labels_on(self, vertices):
        for n in vertices:
            self.labels_box.addItem(n)
        self.label = str(self.labels_box.currentText())

    def open(self):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        self.filename, _ = QFileDialog.getOpenFileName(self, "Open", "", "All Files (*)", options=options)

        self.file_text.setText(str(self.filename))
        if self.filename:
            self.geometries = file_io.read_generic_structure_file(self.filename, read_multiple=True)
            self.changeVertices(self.central_atom_name.isChecked())
            # SecondWindow(self, local_geometry=[self.x(), self.y()]).show()

    def save(self):
        self.directory = QFileDialog.getExistingDirectory(self, "Select Directory")
        self.savefile_text.setText(str(self.directory))
        print(self.directory)
        # measures = [[shape.Shape(geometry).measure(self.label, self.central_atom)
        #              for geometry in self.geometries]]
        # names = [geometry.get_name() for geometry in self.geometries]
        # file_io.shape2file.write_shape_measure_data(measures, names, [self.label])

    def set_central_atom(self, text):
        try:
            self.central_atom = int(text)
        except ValueError:
            self.results.setPlainText('WARNING. Central atom must be an integer')
            self.results.repaint()
            # print('WARNING. Central atom must be an integer')

    def execute(self):
        if self.filename is None:
            self.results.setPlainText('WARNING. No input file selected')
            self.results.repaint()
            # print('WARNING. No input file selected, choose one before running')
        else:
            # if self.directory is None:
            text = "SHAPE's measure " + self.label + "\n"
            for geometry in self.geometries:
                text += '{:10} {:10.3f}'.format(geometry.name,
                                                shape.Shape(geometry).measure(self.label, self.central_atom)) + '\n'

            self.results.setPlainText(text)
            self.results.repaint()
            self.save_button.setEnabled(True)
            # else:
            #     measures = [[shape.Shape(geometry).measure(self.label, self.central_atom)
            #                 for geometry in self.geometries]]
            #     names = [geometry.get_name() for geometry in self.geometries]
            #     file_io.shape2file.write_shape_measure_data(measures, names, [self.label])

    def center(self):
        qr = self.frameGeometry()
        cp = QDesktopWidget().availableGeometry().center()
        qr.moveCenter(cp)
        self.move(qr.topLeft())


class SecondWindow(QMainWindow):
    def __init__(self, parent=None, local_geometry=[0, 0]):
        # super(SecondWindow, self).__init__(parent)
        QWidget.__init__(self, parent)
        self.setGeometry(local_geometry[0] + 500, local_geometry[1], 300, 400)
        self.setWindowTitle('Molecule View')
        self.pen = QPen(QColor(0,0,0))                      # set lineColor
        self.pen.setWidth(3)                                            # set lineWidth
        self.brush = QBrush(QColor(255,255,255,255))        # set fillColor
        self.polygon = self.createpoly(8, 150, 0)                         # polygon with n points, radius, angle of the first point


    def createpoly(self, n, r, s):
            polygon = QPolygonF()
            w = 360 / n  # angle per step
            for i in range(n):  # add the points of polygon
                t = w * i + s
                x = r * np.cos(np.radians(t))
                y = r * np.sin(np.radians(t))
                polygon.append(QPointF(self.width() / 2 + x, self.height() / 2 + y))

            return polygon

    def paintEvent(self, event):
        painter = QPainter(self)
        painter.setPen(self.pen)
        painter.setBrush(self.brush)
        painter.drawPolygon(self.polygon)

if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = MainWindow()
    ex.show()
    sys.exit(app.exec_())