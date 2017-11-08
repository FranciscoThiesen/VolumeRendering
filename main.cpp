#include <iostream>
#include <vector>
#include <algorithm>
#include <string>
#include <fstream>
#include "simpson.hpp"

using std::vector;
using std::string;


// Function for reading CT scans as binary mode file
vector<unsigned char> readData(const string& fileName)
{	
	std::ifstream file(fileName, std::ios::in | std::ios::binary);
	vector<unsigned char> contents((std::istreambuf_iterator<char>(file)),
		std::istreambuf_iterator<char>());
	return contents;
}

int main()
{
	vector<uint8_t> data = readData("head-8bit.raw");
	for(const int x : data){
		std::cout << x << std::endl;
	}
	return 0;
}
