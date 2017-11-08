#include <iostream>
#include <vector>
#include <algorithm>
#include <string>
#include <fstream>
#include "simpson.hpp"

using std::vector;
using std::string;


// Function for reading CT scans as binary mode file
vector<unsigned char> readData(const string& fileName){	
	std::ifstream file(fileName, std::ios::in | std::ios::binary);
	vector<unsigned char> contents((std::istreambuf_iterator<char>(file)),
		std::istreambuf_iterator<char>());
	return contents;
}

// Second attempt at reading CT scans as binary file
vector<unsigned char> readData2(const string& fileName)
{
	vector<unsigned char> rawFileBuffer;
	std::ifstream in;	
	
	in.open(fileName, std::ios::binary);	
	in.seekg(0, std::ios::end);
	size_t fileSize = in.tellg();
	in.seekg(0, std::ios::beg);
	
	rawFileBuffer.resize(fileSize/sizeof(unsigned char));

	in.read((char*)rawFileBuffer.data(), fileSize);

	std::cout << fileSize << std::endl;
	return rawFileBuffer;
}

// Third attempt at reading CT scans as a binary file
static std::vector<char> ReadAllBytes(char const* filename)
{
    std::ifstream ifs(filename, std::ios::binary | std::ios::ate);
    std::ifstream::pos_type pos = ifs.tellg();

    std::vector<char>  result(pos);

    ifs.seekg(0, std::ios::beg);
    ifs.read(&result[0], pos);

    return result;
}

int main(){
	auto data = readData("head-8bit.raw");
	auto data2 = readData2("head-8bit.raw");
	auto data3 = ReadAllBytes("head-8bit.raw");
	for(auto c : data3){
		std::cout << c << " ";
	}
	return 0;
}
