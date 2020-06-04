from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
#/############################################################################
#
# This module is based on an answer published in:
#
# http://stackoverflow.com/questions/11513132/embedding-ipython-qt-console-in-a-pyqt-application
#
# by Tim Rae
#
#############################################################################*/


import os
import sys


os.environ['QT_API'] = 'pyqt'
try:
    import sip
    sip.setapi("QString", 2)
    sip.setapi("QVariant", 2)
except:
    pass

from PyQt4 import QtGui 


import IPython
print( IPython.__version__)
if IPython.__version__.startswith("2"):
    QTCONSOLE = False
else:
    try:
        import qtconsole
        QTCONSOLE = True
    except ImportError:
        QTCONSOLE = False

if QTCONSOLE:
    try:
        from qtconsole.rich_ipython_widget import RichJupyterWidget as RichIPythonWidget
    except:
        from qtconsole.rich_ipython_widget import RichIPythonWidget
    from qtconsole.inprocess import QtInProcessKernelManager
else:
    from IPython.qt.console.rich_ipython_widget import RichIPythonWidget
    from IPython.qt.inprocess import QtInProcessKernelManager
from IPython.lib import guisupport

class QIPythonWidget(RichIPythonWidget):
    """ Convenience class for a live IPython console widget. We can replace the standard banner using the customBanner argument"""
    def __init__(self,customBanner=None,*args,**kwargs):
        super(QIPythonWidget, self).__init__(*args,**kwargs)
        if customBanner != None:
            self.banner = customBanner
        self.setWindowTitle(self.banner)
        self.kernel_manager = kernel_manager = QtInProcessKernelManager()
        kernel_manager.start_kernel()
        kernel_manager.kernel.gui = 'qt4'
        self.kernel_client = kernel_client = self._kernel_manager.client()
        kernel_client.start_channels()

        def stop():
            kernel_client.stop_channels()
            kernel_manager.shutdown_kernel()
            guisupport.get_app_qt4().exit()
        self.exit_requested.connect(stop)

    def updateNamespace(self,variableDict):
        """ Given a dictionary containing name / value pairs, push those variables to the IPython console widget """
        self.kernel_manager.kernel.shell.push(variableDict)
    def clearTerminal(self):
        """ Clears the terminal """
        self._control.clear()
    def showMessage(self,text):
        self._append_plain_text(text)
    def executeCommand(self,command):
        """ Execute a command in the frame of the console widget """
        self._execute(command,False)


class ExampleWidget(QtGui.QWidget):
    def __init__(self, parent=None):

        super(ExampleWidget, self).__init__(parent)
        layout = QtGui.QVBoxLayout(self)
        ipyConsole = QIPythonWidget(customBanner="Welcome to the embedded ipython console\n")
        layout.addWidget(ipyConsole)
        ## ipyConsole.showMessage(message)
        self.console = ipyConsole


if __name__ =="__main__":
    app  = QtGui.QApplication([])
    widget = ExampleWidget()
    
    widget.show()
    ## widget.console.printText(message)
    app.exec_()








