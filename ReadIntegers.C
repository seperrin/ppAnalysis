#if !defined(__CINT__) || defined(__MAKECINT__)

#include <set>
#include <Riostream.h>

#include "TROOT.h"
#include "TSystem.h"
#include "TString.h"
#include "TObjString.h"
#include "TObjArray.h"


#endif


//______________________________________________________________________________
void ReadIntegers(const char* filename,
									std::vector<int>& integers,
									Bool_t resetVector = kTRUE)
{
	/// Read integers from filename, where integers are either
	/// separated by "," or by return carriage
	
	if ( gSystem->AccessPathName(gSystem->ExpandPathName(filename))==kTRUE )
	{
		return;
	}
	std::ifstream in(gSystem->ExpandPathName(filename));
	int i;
	
	std::set<int> runset;
	
	if (!resetVector)
	{
		for ( std::vector<int>::size_type j = 0; j < integers.size(); ++ j )
		{
			runset.insert(integers[j]);
		}
	}
	
	char line[10000];
	
	in.getline(line,10000,'\n');
	
	TString sline(line);
	
	if (sline.Contains(","))
	{
		TObjArray* a = sline.Tokenize(",");
		TIter next(a);
		TObjString* s;
		while ( ( s = static_cast<TObjString*>(next()) ) )
		{
			runset.insert(s->String().Atoi());
		}
		delete a;
	}
	else
	{
		runset.insert(sline.Atoi());
		
		while ( in >> i )
		{
			runset.insert(i);
		}
	}
	
	integers.clear();
	
	for ( std::set<int>::const_iterator it = runset.begin(); it != runset.end(); ++it )
	{
		integers.push_back((*it));
	}
	
//	std::sort(integers.begin(),integers.end());
}
