#pragma once
#ifndef DATA_SET_H
#define DATA_SET_H

#include <set>
#include <stdio.h>
#include "Commom/FileIO.h"

template <class T>
class DataSet
{
public:
	DataSet(const std::string & name) {
		_name = name;
	}

	DataSet(const DataSet & set)
	{
		_name = set._name;
		_data = set._data;
	}

	DataSet(const std::string & name, const std::set<T> & data)
	{
		_name = name;
		_data = data;
	}

	inline std::string getName() const
	{
		return _name;
	}

	inline int getNumElements() const
	{
		return _data.size();
	}

	inline void getElements(std::set<T> & data) const
	{
		data = _data;
	}

	inline bool isMember(T data) const
	{
		return (_data.find(data) != _data.end());
	}

	inline void insert(T data)
	{
		_data.insert(data);
	}

	inline void clear()
	{
		_data.clear();
	}

	bool readBinary(std::ifstream fin)
	{
		readStringAsBinary(fin, _name);
		readOneLevelSet(fin, _data);
	}

	bool writeBinary(std::ofstream fout)
	{
		writeStringAsBinary(fout, _name);
		writeOneLevelSet(fout, _data);
	}
protected:
	std::string _name;
	std::set<T> _data;
};
#endif