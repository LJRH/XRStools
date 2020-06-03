from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import os
import sys

from PyQt4 import Qt

class ProxyModel(Qt.QSortFilterProxyModel):

    def __init__(self, model, parent=None):
        Qt.QSortFilterProxyModel.__init__(self, parent)
        self.setSourceModel(model)
    
    def setRoots(self, roots):
        self.roots = roots
        self.base_root = os.path.commonprefix(roots)
        self.sourceModel().setRootPath(self.base_root)
        self.all_dirs = set()
        if 'posix' in os.name:
            for root in roots:
                last = root
                while last != '/':
                    self.all_dirs.add(last)
                    last = os.path.dirname(last)

    def get_root(self, path):
        for i, root in enumerate(self.roots):
            if root == path:
                return root, i
            if root in path:
                return root, -1
        return None, False

    def mapFromSource(self, sourceIndex):
#        return super(ProxyModel, self).mapFromSource(sourceIndex)
        model = sourceIndex.model()
        if model is None:
            return Qt.QModelIndex()
        path = model.filePath(sourceIndex)
        root, root_idx = self.get_root(path)
        if root_idx >= 0:
            return self.createIndex(root_idx, 0)
        return Qt.QModelIndex()

    def mapToSource(self, proxyIndex):
#        print 'maptosource', proxyIndex
        return super(ProxyModel, self).mapToSource(proxyIndex)
        
            
        
app = Qt.QApplication([])
view = Qt.QTreeView()

roots = ['/usr/local/include', '/usr/local/share']

model = Qt.QFileSystemModel()
proxy = ProxyModel(model, view)
proxy.setRoots(roots)
view.setModel(proxy)
view.show()

sys.exit(app.exec_())







