#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
using namespace std;

#define SPJ_OK 0
#define SPJ_WA 1
#define SPJ_STD_IO_ERROR 2
#define SPJ_SRC_IO_ERROR 3
#define SPJ_SCR_IO_ERROR 4

struct SPJException{
	std::string message;
	int exitcode;
	SPJException() {}
	SPJException(const std::string &message, int exitcode) : message(message), exitcode(exitcode) {}
	void print() {
		std::cout << message << std::endl;
	}
};

long getFileSize(FILE *stream, int exitcode) {
	long curpos, ret;
	if ((curpos = ftell(stream)) < 0) {
		throw SPJException("Failed to read the file.", exitcode);
	}
	if (fseek(stream, 0L, SEEK_END) < 0) {
		throw SPJException("Failed to read the file.", exitcode);
	}
	if ((ret = ftell(stream)) < 0) {
		throw SPJException("Failed to read the file.", exitcode);
	}
	if (fseek(stream, 0L, curpos) < 0) {
		throw SPJException("Failed to read the file.", exitcode);
	}
	return ret;
}

double getScore(int argc, char **argv, int exitcode)  {
	if (argc == 6) {
		FILE *inscr = fopen(argv[4], "r");
		if (!inscr) {
			throw SPJException("Failed to read the score file.", exitcode);
		}
		double buffer;
		fscanf(inscr, "%lf", &buffer);
		fclose(inscr);
		return buffer;
	}
	return 0;
}

void putScore(int argc, char **argv, double score, int exitcode) {
	if (argc == 6) {
		FILE *ouscr = fopen(argv[5], "w");
		if (!ouscr) {
			throw SPJException("Failed to write the score file.", exitcode);
		}
		fprintf(ouscr, "%.20f", score);
		fclose(ouscr);
	}
}

int main(int argc, char *argv[]) {
	try{
		ifstream indat(argv[1]);
		ifstream instd(argv[2]);
		ifstream insrc(argv[3]);
		if (!indat) throw SPJException("Cannot open the input data.", SPJ_SRC_IO_ERROR);
		if (!instd) throw SPJException("Cannot open the stdfile.", SPJ_STD_IO_ERROR);
		if (!insrc) throw SPJException("Cannot open the srcfile.", SPJ_SRC_IO_ERROR);
		int size = 0; string ans;
		while (getline(insrc, ans)) 
        {
			string t;
            stringstream ss(ans);
            ss >> t >> t;
            if (t != "PASS.") 
            {
                instd.close();
                insrc.close();
                indat.close();
                putScore(argc, argv, 0.0, SPJ_SCR_IO_ERROR);
                return SPJ_WA;
            }
		}
		instd.close();
		insrc.close();
		indat.close();
		putScore(argc, argv, getScore(argc, argv, SPJ_SCR_IO_ERROR), SPJ_SCR_IO_ERROR);
		return SPJ_OK;
	} catch(SPJException e) {
		e.print();
		return e.exitcode;
	}
	return 0;
}
