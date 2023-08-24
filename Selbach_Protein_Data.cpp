// Selbach_Protein_Data.cpp : Defines the entry point for the console application.
//

#include <stdio.h>
#include <string.h>
#include <list>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <ctype.h>
#include<iterator>
#include "Protein_Information.h"
#include "Protein.h"
using namespace std;

int main()
{
	std::list<Protein> list;
	Protein p;
	p.Get_Protein_ID(list);
	
	/*std::list<Protein>::iterator it;

	for(it=list.begin(); it!=list.end(); it++)
	{
		
		cout<<it->protein_id<<" ";
		for(std::list<int>::iterator it2=it->peptide_id.begin(); it2!=it->peptide_id.end(); it2++)
		{
			cout<<(*it2)<<" ";

		}
		cout<<endl;
	}
	cout<<list.size()<<endl;*/

	cout<<"Done with Protein info"<<endl;
	Protein_Information pi;
	std::list<Protein_Information> listPi;
	pi.Get_Evidence_File_Info(listPi);

	/*std::multimap<string, Protein_Information> m;
	p.Get_Evidence_File_Info(m);*/

	cout<<"Done with Evidence content"<<endl;
	/*for(std::list<Protein_Information>::iterator it = listPi.begin(); it!=listPi.end(); it++)
	{
		cout<<it->Raw_File<<" "<<it->ratio<<" "<<it->Intensity_L<<" "<<it->Intensity_H<<" "<<it->unmod_pep_id<<" "<<it->mod_pep_id;
		
		cout<<endl;
	}*/
	cout<<"Evidence list size:"<<listPi.size()<<endl;

	
	cout<<"Start Ratio Calculation"<<endl;

	//pi.Calculate_Protein_Ratio(list,listPi);
	//pi.Calculate_Protein_Ratio_With_Evidence(list,listPi);
	
	pi.Calculate_Protein_Ratio_With_Evidence_Peptide_Ratio(list,listPi);

	//getchar();
	return 0;
}

