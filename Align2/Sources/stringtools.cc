// stringtools.cc

#include <stringtools.h>
#include <cstdlib>
using namespace std;

void
strip(string &str)
{
	string temp = "";

	for (string::iterator pos = str.begin(); pos != str.end(); pos++)
		if (*pos != ' ' and *pos != '\n')
			temp += *pos;

	str = temp;
}


string
strip_and_return_string(string str)
{
	string temp = "";

	for (string::iterator pos = str.begin(); pos != str.end(); pos++)
		if (*pos != ' ' and *pos != '\n')
			temp += *pos;

	return temp;
}


int
strip_and_return_int(string str)
{
	string temp = "";

	for (string::iterator pos = str.begin(); pos != str.end(); pos++)
		if (*pos != ' ' and *pos != '\n')
			temp += *pos;

	int integ = atoi(temp.c_str());
	return integ;
}


unsigned int
strip_and_return_unsigned_int(string str)
{
	string temp = "";

	for (string::iterator pos = str.begin(); pos != str.end(); pos++)
		if (*pos != ' ' and *pos != '\n')
			temp += *pos;

	unsigned int integ = atoi(temp.c_str());
	return integ;
}


float
strip_and_return_float(string str)
{
	string temp = "";

	for (string::iterator pos = str.begin(); pos != str.end(); pos++)
	{
		if (*pos == ',')
			temp += '.';
		else
			if (*pos != ' ' and *pos != '\n')
				temp += *pos;
	}

	float doub = atof(temp.c_str());
	return doub;
}


vector<string>
split(string text_to_split, char delimiter)
{
	vector<string> vector_with_splited_text;
	string temp_text = "";

	for (string::iterator pos = text_to_split.begin(); pos != text_to_split.end(); pos++)
		if (*pos != delimiter and *pos !='\n')
			temp_text += *pos;
		else
		{
			vector_with_splited_text.push_back(temp_text);
			temp_text = "";
		}
	vector_with_splited_text.push_back(temp_text); // last piece

	return vector_with_splited_text;
}


int
seq_length(const string str)
{
	unsigned int size = 0;

	for (string::const_iterator pos = str.begin(); *pos != ' '; pos++)
		size++;

	return size;
}


/// Change each element of the string to upper case.
string
StringToUpper(string strToConvert)
{
	for (unsigned int i = 0; i < strToConvert.length(); i++)
		strToConvert[i] = toupper(strToConvert[i]);
	return strToConvert; // return the converted string
}


/// Change each element of the string to lower case.
string
StringToLower(string strToConvert)
{
	for (unsigned int i = 0; i < strToConvert.length(); i++)
		strToConvert[i] = tolower(strToConvert[i]);
	return strToConvert; // return the converted string
}
