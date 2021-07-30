#pragma once
#ifndef FILE_IO_H
#define FILE_IO_H
#include <iostream>
#include <fstream>
#include <string>
#include <strstream>
#include <sstream>

static void getFilenameInfo(const std::string fullStr, std::string& dir, std::string& name, std::string& format)
{
	dir = fullStr.substr(0, fullStr.find_last_of("/") + 1);
	std::string temp = fullStr.substr(dir.length(), fullStr.length() - dir.length());
	name = temp.substr(0, temp.find_first_of("."));
	format = temp.substr(temp.find_first_of(".") + 1, temp.length() - name.length());
}

static std::string getPathDir(const std::string& filepath)
{
	if (filepath.find_last_of("/\\") != std::string::npos)
		return filepath.substr(0, filepath.find_last_of("/\\") + 1);
	return "";
}

static void writeStringAsBinary(std::ofstream& fout, std::string& s)
{
	int size = s.size();
	fout.write((char*)&size, sizeof(int));
	fout.write(s.c_str(), size);
}

static void readStringAsBinary(std::ifstream& fin, std::string& s)
{
	int size;
	fin.read((char*)&size, sizeof(int));
	char* buf = new char[size];
	fin.read(buf, size);
	s.append(buf, size);
}

#include <vector>
#include <set>

namespace StlBinaryIO
{	
	template <typename T>
	void writeOneLevelVector(std::ofstream& fout, std::vector<T>& vec)
	{
		int size = vec.size();
		fout.write((const char*)&size, sizeof(int));
		if (size != 0)
			fout.write((const char*)&vec[0], size * sizeof(T));
	}

	template <typename T>
	void readOneLevelVector(std::ifstream& fin, std::vector<T>& vec)
	{
		int size;
		fin.read((char*)&size, sizeof(int));
		if (size != 0)
		{
			vec.resize(size);
			fin.read((char*)&vec[0], size * sizeof(T));
		}
	}

	template <typename T>
	void writeTwoLevelVector(std::ofstream& fout, std::vector<std::vector<T>>& vec)
	{
		int outter_size = vec.size();
		if (outter_size != 0)
		{
			fout.write((const char*)&outter_size, sizeof(int));
			for (int i = 0; i < outter_size; i++)
			{
				int inner_size = vec[i].size();
				fout.write((const char*)&inner_size, sizeof(int));
				if (inner_size != 0)
					fout.write((const char*)&vec[i][0], inner_size * sizeof(T));
			}
		}

	}

	template <typename T>
	void readTwoLevelVector(std::ifstream& fin, std::vector<std::vector<T>>& vec)
	{
		int outter_size;
		fin.read((char*)&outter_size, sizeof(int));
		if (outter_size != 0)
		{
			vec.resize(outter_size);
			for (int i = 0; i < outter_size; i++)
			{
				int inner_size;
				fin.read((char*)&inner_size, sizeof(int));
				if (inner_size != 0)
				{
					vec[i].resize(inner_size);
					fin.read((char*)&vec[i][0], inner_size * sizeof(T));
				}
			}
		}
	}

	template <typename T>
	void writeOneLevelSet(std::ofstream& fout, std::set<T>& set)
	{
		int size = set.size();
		std::vector<T> vec;
		if (size != 0)
		{
			vec.resize(size);
			std::copy(set.begin(), set.end(), vec.begin());
		}
		writeOneLevelVector(fout, vec);
	}

	template <typename T>
	void readOneLevelSet(std::ifstream& fin, std::set<T>& set)
	{
		std::vector<T> vec;
		readOneLevelVector(fin, vec);
		if (vec.size() != 0)
			set = std::set<T>(vec.begin(), vec.end());
	}
}

namespace StlTextIO
{
	template <typename T>
	void writeOneLevelVector(std::ofstream& fout, std::vector<T>& vec)
	{
		fout << size << std::endl;
		for (int i = 0; i < size; i++)
			fout << vec[i] << std::endl;
	}

	template <typename T>
	void readOneLevelVector(std::ifstream& fin, std::vector<T>& vec)
	{
		int size;
		fin >> size;
		vec.resize(size);
		for (int i = 0; i < size; i++)
		{
			T data;
			fin >> data;
			vec[i] = data;
		}
	}

	template <typename T>
	void writeTwoLevelVector(std::ofstream& fout, std::vector<std::vector<T>>& vec)
	{
		int outter_size = vec.size();
		if (outter_size != 0)
		{
			fout << outter_size << std::endl;
			for (int i = 0; i < outter_size; i++)
			{
				int inner_size = vec[i].size();
				fout << inner_size;
				if (inner_size != 0)
				{
					for (int j = 0; j < inner_size, j++)
					{
						fout << " " << vec[i][j];
					}					
				} 
				fout << std::endl;
			}
		}
	}

	template <typename T>
	void readTwoLevelVector(std::ifstream& fin, std::vector<std::vector<T>>& vec)
	{
		int outter_size;
		fin >> outter_size;
		if (outter_size != 0)
		{
			vec.resize(outter_size);
			for (int i = 0; i < outter_size; i++)
			{
				int inner_size;
				fin >> inner_size;
				if (inner_size != 0)
				{
					vec[i].resize(inner_size);
					for (int j = 0; j < inner_size, j++)
					{
						int data;
						fin >> data;
						vec[i][j] = data;
					}
				}
			}
		}
	}

	template <typename T>
	void writeOneLevelSet(std::ofstream& fout, std::set<T>& set)
	{
		int size = set.size();
		std::vector<T> vec;
		if (size != 0)
		{
			vec.resize(size);
			std::copy(set.begin(), set.end(), vec.begin());
		}
		writeOneLevelVector(fout, vec);
	}

	template <typename T>
	void readOneLevelSet(std::ifstream& fin, std::set<T>& set)
	{
		std::vector<T> vec;
		readOneLevelVector(fin, vec);
		if (vec.size() != 0)
			set = std::set<T>(vec.begin(), vec.end());
	}
}

#endif