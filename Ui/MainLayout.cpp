#include "MainLayout.h"
#define FIXED_RIGHT_ITEM_WIDGET 280
#define FIXED_BOTTOM_ITEM_HEIGHT 80
#define FIXED_SPACING 0

MainLayout::MainLayout(QWidget * parent, int margin, int spacing): QLayout(parent)
{
	setMargin(margin);
	setSpacing(spacing);
}

MainLayout::MainLayout(int spacing)
{
	setSpacing(spacing);
}

MainLayout::~MainLayout()
{
	QLayoutItem *l;
	while ((l = takeAt(0)))
		delete l;
}

void MainLayout::addItem(QLayoutItem *item)
{
	add(item, CENTER);
}

void MainLayout::addWidget(QWidget *widget, Position position)
{
	add(new QWidgetItem(widget), position);
}

Qt::Orientations MainLayout::expandingDirections() const
{
	return Qt::Horizontal | Qt::Vertical;
}

bool MainLayout::hasHeightForWidth() const
{
	return false;
}

int MainLayout::count() const
{
	return _list.size();
}

QLayoutItem* MainLayout::itemAt(int index) const
{
	ItemWrapper *wrapper = _list.value(index);
	if (wrapper)
		return wrapper->item;
	else
		return 0;
}

QSize MainLayout::minimumSize() const
{
	return calculateSize(MinimumSize);
}

void MainLayout::setGeometry(const QRect &rect)
{
	ItemWrapper *center = 0;
	int centerHeight = 0;
	int centerWidth = 0;
	int bottomWidth = 0;
	int bottomHeight = 0;	
	int rightWidth = 0;
	int i;

	QLayout::setGeometry(rect);

	for (int i = 0; i < _list.size(); ++i)
	{
		ItemWrapper *wrapper = _list.at(i);
		QLayoutItem *item = wrapper->item;
		Position position = wrapper->position;

		if (position == RIGHT)
		{
			item->setGeometry(QRect(item->geometry().x(), item->geometry().y(), FIXED_RIGHT_ITEM_WIDGET, rect.height()));
			rightWidth = FIXED_RIGHT_ITEM_WIDGET;
			item->setGeometry(QRect(rect.x() + rect.width() - rightWidth, rect.y(), FIXED_RIGHT_ITEM_WIDGET, rect.height()));
		}
	}

	centerWidth = rect.width() - rightWidth - FIXED_SPACING;
	bottomWidth = centerWidth;

	for (int i = 0; i < _list.size(); i++)
	{
		ItemWrapper *wrapper = _list.at(i);
		QLayoutItem *item = wrapper->item;
		Position position = wrapper->position;

		if (position == BOTTOM)
		{
			item->setGeometry(QRect(item->geometry().x(),
				item->geometry().y(), bottomWidth,
				FIXED_BOTTOM_ITEM_HEIGHT));

			bottomHeight = FIXED_BOTTOM_ITEM_HEIGHT;

			item->setGeometry(QRect(rect.x(),
				rect.y() + rect.height() - FIXED_BOTTOM_ITEM_HEIGHT,
				bottomWidth - FIXED_SPACING,
				FIXED_BOTTOM_ITEM_HEIGHT - FIXED_SPACING));
		}else if (position == CENTER)
		{
			center = wrapper;
		}

	}

	centerHeight = rect.height() - bottomHeight - FIXED_SPACING;

	if (center)
	{
		center->item->setGeometry(QRect(rect.x(), rect.y(),
			centerWidth,
			centerHeight));
	}
}

QSize MainLayout::sizeHint() const
{
	return calculateSize(SizeHint);
}

QLayoutItem* MainLayout::takeAt(int index)
{
	if (index >= 0 && index < _list.size()) {
		ItemWrapper *layoutStruct = _list.takeAt(index);
		return layoutStruct->item;
	}
	return 0;
}

void MainLayout::add(QLayoutItem * item, Position position)
{
	_list.append(new ItemWrapper(item, position));
}

void MainLayout::remove(Position p)
{
	for (int i = 0; i < _list.size(); i++)
	{
		ItemWrapper *wrapper = _list.at(i);
		QLayoutItem *item = wrapper->item;
		Position position = wrapper->position;
		if (position == p)
		{
			removeItem(item);
			_list.removeAt(i);
		}
	}
}

QSize MainLayout::calculateSize(SizeType sizeType) const
{
	QSize totalSize;
	int cw = 0;
	int ch = 0;
	int bw = 0;
	int bh = 0;
	int rw = 0;
	int rh = 0;


	for (int i = 0; i < _list.size(); ++i)
	{
		ItemWrapper *wrapper = _list.at(i);
		Position position = wrapper->position;
		QSize itemSize = wrapper->item->sizeHint();
		if (position == CENTER)
		{
			cw = itemSize.width();
			ch = itemSize.height();
		}
		else if (position == BOTTOM)
		{
			bw = itemSize.width();
			bh = itemSize.height();
		}
		else if (position == RIGHT)
		{
			rw = itemSize.width();
			rh = itemSize.height();
		}

	}

	totalSize.rwidth() = rw;
	if (cw > 0)totalSize.rwidth() += cw;
	else totalSize.rwidth() += bw;

	if (rh > 0)
		totalSize.rheight() = rh;
	else totalSize.rheight() = ch + bh;
	
	return totalSize;
}
