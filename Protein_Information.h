#include <string>
#include <list>
#include <map>
#include "Protein.h"

using namespace std;

#pragma once

class Protein_Information
{


public:
	string seq;
	string Raw_File;
	double ratio;
	long Intensity_L;
	long Intensity_H;
	//list<int> pep_id;
	int unmod_pep_id;
	int mod_pep_id;
	long scan_number;
	long evidence_number;
	double m_z;
	void Calculate_Protein_Ratio(list<Protein> &p, list<Protein_Information> &pi);
	void Calculate_Protein_Ratio_With_Evidence(list<Protein> &p, list<Protein_Information> &pi);
	void Calculate_Protein_Ratio_With_Evidence_Peptide_Ratio(list<Protein> &p, list<Protein_Information> &pi);
	void Cal_Median(list<double> &r1, list<double> &r2, list<double> &r3, double &r1_median, double &r2_median, double &r3_median);
	void Get_Evidence_File_Info(list<Protein_Information> & m);
	double Cal_Median_One(list<double> &r1);
	void Traverse_Map(std::multimap<string,double> &r_peptide_list, std::list<double> &r);
	Protein_Information(void);
	~Protein_Information(void);


};


