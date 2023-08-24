#include "Protein.h"
#include <stdio.h>
#include <string.h>
#include <list>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <ctype.h>

Protein::Protein(void)
{
}


Protein::~Protein(void)
{
}

/*This function reads the input file and collects protein ids*/
void Protein::Get_Protein_ID(list<Protein> &protein_id)
{


	 //contains protein ids from proteinGroups
	int column_count, line_count = 0;
	string line, csvItem, protein;

	ifstream myfile;
	ofstream out,out2;

	myfile.open("proteinGroups01.csv");
	//out2.open("peptide_id_list.txt");

	if(!myfile)
	{
		cout<<"Can not open input csv file"<<endl;
		getchar();
		exit(-1);
	}


	line_count = 1;
    while (getline(myfile,line)) 
    {
		csvItem.clear();
		if(line_count >= 2)
		{
			column_count = 1;
			istringstream myline(line);
			Protein p;
			while(getline(myline, csvItem,','))
			{
					/*if(csvItem.size()>2)
					{*/
                       
						  if(column_count==1)
							{
									char *str = new char[csvItem.length()+1];

									char *pch;

									std::strcpy(str,csvItem.c_str());

									pch = strtok(str,":|;");
								    
									string cmp = pch;

									int count  = 1;
								    
									if(cmp == "IPI")
									{
										  while(pch != NULL) 
										  {
											//out<<pch<<" ";
											if(count == 2)
											{
												p.protein_id = pch;
												/*if(pch == "IPI00112645.1")
													out2<<pch<<endl;*/
												//cout<<pch<<"\n";
												break;
											}
											//cout<<pch<<"\n";
											pch = strtok (NULL,":|;");
											count ++;
										  }

									}

									std::strcpy(str,csvItem.c_str());

									pch = strtok(str,":|;");
								    
									cmp = pch;

									if(cmp == "IPI")
									{
										  while(pch != NULL) 
										  {
											  cmp =pch;

											  if(cmp[0]=='I' && cmp[1]=='P' && cmp[2]=='I' && cmp.length() > 3)
												  p.protein_id_list.push_back(cmp);
											  
											  pch = strtok (NULL,":|;");
											  
										  }

									}
									



									delete []str;
									//break;
								}
						    
						      if(column_count==32)
							  {

								    char *str = new char[csvItem.length()+1];

									char *pch;

									std::strcpy(str,csvItem.c_str());

									pch = strtok(str,";");									
									//out2<<"unmodified peptide id"<<endl;
									while(pch != NULL) 
								    {
										/*if(p.protein_id =="IPI00112645.1")
											out2<<atoi(pch)<<endl;*/

										p.peptide_id.push_back(atoi(pch));
											
											//cout<<pch<<"\n";
										pch = strtok (NULL,";");
											
								     }
									delete []str;
									//out2<<endl;
								}

								if(column_count==34)
							    {

								    char *str = new char[csvItem.length()+1];

									char *pch;

									std::strcpy(str,csvItem.c_str());

									pch = strtok(str,";");									
									
									//out2<<"modified peptide id"<<endl;

									while(pch != NULL) 
								    {
											
										p.mod_peptide_id.push_back(atoi(pch));
										/*if(p.protein_id =="IPI00112645.1")
											out2<<atoi(pch)<<endl;*/
											//cout<<pch<<"\n";
										pch = strtok (NULL,";");
											
								     }
									delete []str;
									
								}								  
								if(column_count==35)
								{
									 char *str = new char[csvItem.length()+1];

									char *pch;

									std::strcpy(str,csvItem.c_str());

									pch = strtok(str,";");									
									
									//out2<<"modified peptide id"<<endl;

									while(pch != NULL) 
								    {
											
										p.evidence_id.push_back(atoi(pch));
										/*if(p.protein_id =="IPI00112645.1")
											out2<<atoi(pch)<<endl;*/
											//cout<<pch<<"\n";
										pch = strtok (NULL,";");
											
								     }
									delete []str;

								}


					 
					/*}*/
					column_count++;
			
			}
			if(p.protein_id!="")
			  protein_id.push_back(p);		    
		
		}

		line_count++;
		
	}


	myfile.close();
	//out2.close();
	
}

