#include <QtGui/QtGui>
#include<stdio.h> 
QList<QStandardItem *> list;
 
class SortProxy : public QAbstractProxyModel
{
 Q_OBJECT
 
 public:
  SortProxy(QObject *parent = 0) : QAbstractProxyModel(parent), hideThem(false)
  {
    fixModel();
  }
  
  int rowCount(const QModelIndex &parent) const
  {
    QModelIndex sourceParent;
    if (parent.isValid())
      sourceParent = mapToSource(parent);
    int count = 0;
    QMapIterator<QPersistentModelIndex, QPersistentModelIndex> it(proxySourceParent);
    while (it.hasNext()) {
      it.next();
      if (it.value() == sourceParent)
	count++;
    }
    return count;
  }
  
  int columnCount(const QModelIndex &) const
  {
    return 1;
  }
  
  QModelIndex index(int row, int column, const QModelIndex &parent = QModelIndex()) const
  {
    QModelIndex sourceParent;
    if (parent.isValid())
      sourceParent = mapToSource(parent);
    QMapIterator<QPersistentModelIndex, QPersistentModelIndex> it(proxySourceParent);
    while (it.hasNext()) {
      it.next();
      if (it.value() == sourceParent && it.key().row() == row &&
	  it.key().column() == column)
	return it.key();
    }
    return QModelIndex();
  }
  
  QModelIndex parent(const QModelIndex &child) const
  {
    QModelIndex mi = proxySourceParent.value(child);
    if (mi.isValid())
      return mapFromSource(mi);
    return QModelIndex();
  }
  
  QModelIndex mapToSource(const QModelIndex &proxyIndex) const
  {
    if (!proxyIndex.isValid())
      return QModelIndex();
 return mapping.key(proxyIndex);
  }
  
  QModelIndex mapFromSource(const QModelIndex &sourceIndex) const
  {
    if (!sourceIndex.isValid())
      return QModelIndex();
    return mapping.value(sourceIndex);
  }
 
public slots:
  void hideEverythingButA1AndChildren()
  {
    hideThem = !hideThem;
    // Now we set up the proxy <-> source mappings
    emit layoutAboutToBeChanged();
    fixModel();
    emit layoutChanged();
  }
  
private:
  void fixModel()
  {
    mapping.clear();
    proxySourceParent.clear();
    for (int i=0;i<list.size();i++) {
      QStandardItem *si = list.at(i);
      if (hideThem) {
	if (!si->text().startsWith("A") || !si->parent())
	  continue;
	QModelIndex proxy = createIndex(si->row(), si->column(), si->index().internalPointer());
	mapping.insert(QPersistentModelIndex(si->index()), proxy);
	QModelIndex sourceParent;
	if (si->parent()->parent())
	  sourceParent = si->parent()->index();
	proxySourceParent.insert(proxy, sourceParent);
      } else {
	QModelIndex proxy = createIndex(si->row(), si->column(), si->index().internalPointer());
	mapping.insert(QPersistentModelIndex(si->index()), proxy);
	QModelIndex sourceParent;
	if (si->parent())
	  sourceParent = si->parent()->index();
	proxySourceParent.insert(proxy, sourceParent);
      }
    }
  }
  QMap<QPersistentModelIndex, QPersistentModelIndex> mapping;
  QMap<QPersistentModelIndex, QPersistentModelIndex> proxySourceParent;
  bool hideThem;
};

SortProxy *proxyModel = 0;

class Tree : public QTreeView
{
  Q_OBJECT
  
  public:
  Tree(QWidget *parent = 0) : QTreeView(parent)
  {
    QStandardItemModel *sourceModel = new QStandardItemModel(this);
    
    QStandardItem *parentA = sourceModel->invisibleRootItem();
    for (int i = 0; i < 2;i++) {
      itemA = new QStandardItem(QString("A %0").arg(i));
      parentA->appendRow(itemA);
      parentA = itemA;
      list.append(itemA);
    }
    itemA = new QStandardItem(QString("A 2"));
    parentA->appendRow(itemA);
    list.append(itemA);
    itemA3 = new QStandardItem(QString("A 3"));
    list.append(itemA3);
    parentA->appendRow(itemA3);
    itemA4 = new QStandardItem(QString("A 4"));
    list.append(itemA4);
    parentA->appendRow(itemA4);
    itemNonA = new QStandardItem(QString("Non A"));
    list.append(itemNonA);
    parentA->appendRow(itemNonA);
    
    QStandardItem *parentB = sourceModel->invisibleRootItem();
    for (int i = 0; i < 3;i++) {
      itemB = new QStandardItem(QString("B %0").arg(i));
      parentB->appendRow(itemB);
      parentB = itemB;
      list.append(itemB);
    }
    
    QStandardItem *parentC = sourceModel->invisibleRootItem();
    for (int i = 0; i < 3;i++) {
      itemC = new QStandardItem(QString("C %0").arg(i));
      parentC->appendRow(itemC);
      parentC = itemC;
      list.append(itemC);
      
    }
    
    proxyModel = new SortProxy(this);
    proxyModel->setSourceModel(sourceModel);
    setModel(proxyModel);
    expandAll();
  }
  QStandardItem *itemA;
  QStandardItem *itemA3;
  QStandardItem *itemA4;
  QStandardItem *itemNonA;
  QStandardItem *itemB;
  QStandardItem *itemC;
};


#include "main.moc"

int main(int argc, char **argv)
{
  QApplication app(argc, argv);
  QWidget widget;
  QPushButton *button = new QPushButton("Make only A1 A  children visible", &widget);
  printf(" Tree \n");
  Tree *tree = new Tree(&widget);
  printf(" tree OK \n");
  QVBoxLayout *lay = new QVBoxLayout(&widget);
  lay->addWidget(button);
  QObject::connect(button, SIGNAL (clicked()), proxyModel, SLOT (hideEverythingButA1AndChildren()));
  lay->addWidget(tree);
  printf(" SHOW \n");
  widget.show();
  return app.exec();
}
