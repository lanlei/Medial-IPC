#pragma once
#ifndef BUFFER_SERIALIZATION_H
#define BUFFER_SERIALIZATION_H
#include <vector>

template<typename T>
struct OverallVectorBuffer
{
	int fragNum = 0;
	std::vector<T> buffer;
	std::vector<int> fragOffset;
	std::vector<int> fragSpan;
};

template<typename T>
struct FragmentVectorBufferPtr
{
	int fragId = 0;
	T* buffer = nullptr;
	int offset = 0;
	int span = 0;
};

template<typename T>
void registerToFragmentBuffer(std::vector<T>& source, FragmentVectorBufferPtr<T>& result)
{
	size_t size = source.size();
	result.offset = 0;
	result.span = size;
	if (result.buffer)
		free(result.buffer);
	result.buffer = (T*)malloc(sizeof(T) * size);
	std::copy(source.begin(), source.end(), result.buffer);
}

template<typename T>
void registerToOverallBuffer(std::vector<T>& source, OverallVectorBuffer<T>& target, FragmentVectorBufferPtr<T>& result)
{
	result.fragId = target.fragNum++;

	int targetSize = target.buffer.size();
	target.fragOffset.push_back(targetSize);

	int sourceSize = source.size();
	target.fragSpan.push_back(sourceSize);

	result.offset = target.fragOffset[result.fragId];
	result.span = target.fragSpan[result.fragId];

	target.buffer.resize(sourceSize + targetSize);

	std::copy(source.begin(), source.end(), target.buffer.data() + result.offset);
	result.buffer = target.buffer.data() + result.offset;
}

template<typename T>
void registerToOverallBuffer(std::vector<std::vector<T>>& source, OverallVectorBuffer<T>& target, std::vector<FragmentVectorBufferPtr<T>>& result)
{
	for (int i = 0; i < source.size(); i++)
	{
		FragmentVectorBufferPtr<T> result_;
		registerToOverallBuffer(source[i], target, result_);
		result.push_back(result_);
	}
	for (int i = 0; i < result.size(); i++)
	{
		result[i].buffer = target.buffer.data() + result[i].offset;
	}
}

template<typename T>
void updateOverallBufferLastFrag(std::vector<T>& source, OverallVectorBuffer<T>& target, FragmentVectorBufferPtr<T>& result)
{
	int fragId = result.fragId;
	int sourceSize = source.size();
	int targetSize = target.buffer.size();
	target.buffer.resize(targetSize - result.span + sourceSize);
	std::copy(source.begin(), source.end(), target.buffer.data() + result.offset);
	result.buffer = target.buffer.data() + result.offset;
	target.fragSpan[fragId] = sourceSize;
	result.span = target.fragSpan[fragId];
}

template<typename T>
void updateOverallBufferLastFrag(std::vector<std::vector<T>>& source, OverallVectorBuffer<T>& target, std::vector<FragmentVectorBufferPtr<T>>& result)
{
	result.resize(source.size());
	for (int i = 0; i < source.size(); i++)
	{
		FragmentVectorBufferPtr<T> result_;
		registerToOverallBuffer(source[i], target, result_);
		result[i] = result_;
	}
	for (int i = 0; i < result.size(); i++)
	{
		result[i].buffer = target.buffer.data() + result[i].offset;
	}
}

#endif