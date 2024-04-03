#ifndef __QUALIFIER_H__
#define __QUALIFIER_H__
#include "commondefs.h"
#include "log.h"

using namespace std;
namespace NMSUPratt{

	class CQualifier{
	public:
		int npars;
		string qualname;
		vector<string> type;
		vector<string> parname;
		vector<string> value;
		char message[CLog::CHARLENGTH];
	};

	class CQualifiers{
	public:
		vector<CQualifier *> qualifier;
		int nqualifiers;
		void Read(string qfilename);
		void SetPars(CparameterMap *pmap,int iqualifier);
		void Print();
		char message[CLog::CHARLENGTH];
	};
}

#endif
