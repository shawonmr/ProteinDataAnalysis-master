#include <string>
#include <list>

using namespace std;

#pragma once

class Protein
{
public:
	string protein_id;
	list<string> protein_id_list;
	list<int> peptide_id;
	list<int> mod_peptide_id;
	list<int> evidence_id;
	void Get_Protein_ID(list<Protein> &protein_id);
	Protein(void);
	~Protein(void);
};

