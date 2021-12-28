import sys
from PyQt6.QtWidgets import QApplication, QWidget
from PyQt6.QtGui import QIcon


class Window(QWidget):
    def __init__(self):
        super().__init__()
        self.icon = QIcon('images/icon.png')
        self.setGeometry(200, 200, 700, 400)
        self.setWindowTitle("Induce Voltage GUI")
        self.setWindowIcon(QIcon(self.icon))
        self.setFixedSize(700, 400)
app = QApplication(sys.argv)
window = Window()
window.show()
sys.exit(app.exec())
