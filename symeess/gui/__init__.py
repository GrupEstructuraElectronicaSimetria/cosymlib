#!/usr/bin/python3
import sys
from os.path import expanduser
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from symeess.shape import shape_tools
from symeess import file_io, shape2file, shape

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


class Example(QWidget):

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
        file.move(x-100, 0)
        file.clicked.connect(self.open)
        self.file_text = QLineEdit(self)
        self.file_text.move(2, 2)
        self.file_text.resize(x - 150, 20)

        # savefile = QPushButton('Save file', self)
        # savefile.move(x-100, 30)
        # savefile.clicked.connect(self.save)
        # self.savefile_text = QLineEdit(self)
        # self.savefile_text.move(2, 30)
        # self.savefile_text.resize(x - 150, 20)

        # Label selection
        vertices_name = QLabel('Vertices :', self)
        vertices_name.move(30, 60)
        vertices_box = QComboBox(self)
        vertices = ['2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '20', '24', '48', '60']
        for n in vertices:
            vertices_box.addItem(n)
        vertices_box.move(100, 60)
        vertices_box.activated[str].connect(self.vertices_on)

        labels_name = QLabel('Label :', self)
        labels_name.move(200, 60)
        self.labels_box = QComboBox(self)
        labels = shape_tools.get_structure_references(int(2))
        for n in labels:
            self.labels_box.addItem(n)
        self.labels_box.move(250, 60)
        self.label = str(self.labels_box.currentText())
        self.labels_box.activated[str].connect(self.labels_on)

        central_atom_name = QLabel('Central atom :', self)
        central_atom_name.move(350, 60)
        central_atom = QLineEdit(self)
        central_atom.move(450, 60)
        central_atom.resize(30, 30)
        central_atom.textChanged.connect(self.set_central_atom)

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

    def vertices_on(self, text):
        self.labels_box.clear()
        labels = shape_tools.get_structure_references(int(text))
        for n in labels:
            self.labels_box.addItem(n)
        self.labels_box.move(250, 60)
        self.labels_on()

    def labels_on(self):
        self.label = str(self.labels_box.currentText())

    def open(self):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        self.filename, _ = QFileDialog.getOpenFileName(self, "Open", "", "All Files (*)", options=options)

        self.file_text.setText(str(self.filename))
        if self.filename:
            self.geometries = file_io.read_input_file(self.filename)

    def save(self):
        self.directory = QFileDialog.getExistingDirectory(self, "Select Directory")
        self.savefile_text.setText(str(self.directory))

    def set_central_atom(self, text):
        try:
            self.central_atom = int(text)
        except ValueError:
            print('WARNING. Central atom must be an integer')

    def execute(self):
        if self.filename is None:
            print('WARNING. No input file selected, choose one before running')
        else:
            if self.directory is None:
                text = "SHAPE's measure " + self.label + "\n"
                for geometry in self.geometries:
                    text += '{:10} {:10.3f}'.format(geometry.get_name(),
                                                    shape.Shape(geometry).measure(self.label, self.central_atom)) + '\n'

                self.results.setPlainText(text)
                self.results.repaint()
            # else:
            #     measures = [shape.Shape(geometry).measure(self.label, self.central_atom)
            #                 for geometry in self.geometries]
            #     names = [geometry.get_name() for geometry in self.geometries]
            #     shape2file.write_shape_measure_data(measures, names, self.label)

    def center(self):
        qr = self.frameGeometry()
        cp = QDesktopWidget().availableGeometry().center()
        qr.moveCenter(cp)
        self.move(qr.topLeft())


if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = Example()
    ex.show()
    sys.exit(app.exec_())