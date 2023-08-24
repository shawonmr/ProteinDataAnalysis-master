#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <list>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <ctype.h>
#include <algorithm>
#include <vector>
#include <math.h>
#include <iomanip>
#include <map>
#include "Protein_Information.h"


Protein_Information::Protein_Information(void)
{

}

Protein_Information::~Protein_Information(void)
{

}


/*This function trverses peptide maps to find out common H/L ratio belong to the peptide*/
void Protein_Information::Traverse_Map(std::multimap<string,double> &r_peptide_list, std::list<double> &r)
{

	std::list<double> new_r1;
	std::multimap<string, double>::iterator it;
	for(it=r_peptide_list.begin(); it!=r_peptide_list.end(); it = r_peptide_list.upper_bound(it->first))
	{
		std::pair<std::multimap<string, double>::iterator, std::multimap<string,double>::iterator> ret;
		
		ret = r_peptide_list.equal_range(it->first);
		
		for(std::multimap<string, double>::iterator it1 = ret.first; it1!=ret.second; it1++)
		{
			new_r1.push_back(it1->second);
		}
		r.push_back(Cal_Median_One(new_r1));
		//r_peptide_list.erase(r_peptide_list.find(it->first));
		new_r1.clear();

	}


}

/*This function calculates median to a given list*/
double Protein_Information :: Cal_Median_One(list<double> &r1)
{

	int N,i;
	double r1_median = 0;

  size_t size;

  N	= r1.size();

  vector<double> v1(N);
  
 
  if(r1.size() > 0)
  {
       i = 0;
	  for(list<double>::iterator t1 = r1.begin(); t1!=r1.end(); t1++)
       {
	      v1[i++] = (*t1); 
       }
  
	  size = v1.size();

	  std::sort(v1.begin(), v1.end());

	  if (size  % 2 == 0)
	  {
		  r1_median = (v1[size / 2 - 1] + v1[size / 2]) / 2.0;
	  }
	  else 
	  {
		  r1_median = v1[size / 2];
	  }
  }
  else
	  r1_median = 0;

  v1.clear();

  return r1_median;

}
void Protein_Information :: Cal_Median(list<double> &r1, list<double> &r2, list<double> &r3, double &r1_median, double &r2_median, double &r3_median)
{
  int N,i;
  size_t size;

  N	= r1.size();

  vector<double> v1(N);
  
 
  if(r1.size() > 0)
  {
       i = 0;
	  for(list<double>::iterator t1 = r1.begin(); t1!=r1.end(); t1++)
       {
	      v1[i++] = (*t1); 
       }
  
	  size = v1.size();

	  std::sort(v1.begin(), v1.end());

	  if (size  % 2 == 0)
	  {
		  r1_median = (v1[size / 2 - 1] + v1[size / 2]) / 2.0;
	  }
	  else 
	  {
		  r1_median = v1[size / 2];
	  }
  }
  else
	  r1_median = 0;

  v1.clear();

  N	= r2.size();

  vector<double> v2(N);
  
  if(r2.size() > 0)
  {
	  i = 0;
	  for(list<double>::iterator t1 = r2.begin(); t1!=r2.end(); t1++)
	  {
		  v2[i++] = (*t1); 
	  }
  
	  size = v2.size();

	  std::sort(v2.begin(), v2.end());

	  if (size  % 2 == 0)
	  {
		  r2_median = (v2[size / 2 - 1] + v2[size / 2]) / 2.0;
	  }
	  else 
	  {
		  r2_median = v2[size / 2];
	  }
  }
  else
	  r2_median = 0;

  v2.clear();
   
  N	= r3.size();

  vector<double> v3(N);
  
  if(r3.size() > 0)
  {
	  i = 0;
	  for(list<double>::iterator t1 = r3.begin(); t1!=r3.end(); t1++)
	  {
		  v3[i++] = (*t1); 
	  }
  
	  size = v3.size();

	  std::sort(v3.begin(), v3.end());

	  if (size  % 2 == 0)
	  {
		  r3_median = (v3[size / 2 - 1] + v3[size / 2]) / 2.0;
	  }
	  else 
	  {
		  r3_median = v3[size / 2];
	  }
  }

	  else
		  r3_median = 0;
  
  v3.clear();
  
}

/*This function calculates protein half-life while using lists of protien and protein information classes*/
void Protein_Information::Calculate_Protein_Ratio(list<Protein> &p, list<Protein_Information> &pi)
{
	ofstream out,out2;
	out.open("Protein_ratio_calculation.txt");
	out2.open("debug_info.txt");
	std::list<Protein>::iterator it;

	out<<"Protein ID"<<" "<<"R1(1.5h)"<<" "<<"R2(4.5h)"<<" "<<"R3(13.5)"<<" "<<"Half Life (h)"<<" "<<"R1 count"<<" "<<"R2 count"<<" "<<"R3 count"<<endl;
	for(it=p.begin(); it!=p.end(); it++)
	{
		
		//cout<<it->protein_id<<endl;
		
		if(it->protein_id=="IPI00109529.2")
		{
		

							std::list<double> r1;
							std::list<double> r2;
							std::list<double> r3;
							double r1_median = 0;
							double r2_median = 0;
							double r3_median = 0;
							int c1, c2, c3 = 0;

				   c1 = c2 = c3 = 0;

				for(std::list<Protein_Information>::iterator it3 = pi.begin(); it3!=pi.end(); it3++)
				{
                      
				         for(std::list<int>::iterator it2=it->peptide_id.begin(); it2!=it->peptide_id.end(); it2++)
				         {																			
							
								if((*it2)==it3->unmod_pep_id)
								{		
									out2<<it->protein_id<<" "<<it3->seq<<" "<<it3->Raw_File<<" "<<it3->scan_number<<" "<<it3->m_z<<" "<<it3->ratio<<" pep id "<<it3->unmod_pep_id<<endl;
										std::size_t found;
										found = it3->Raw_File.find("P1");
										if(found!=std::string::npos)
										{
											c1++;
											if(it3->ratio != 0 )
											{
												r1.push_back(it3->ratio);
												
											}
											

										}
										found = it3->Raw_File.find("P2");
										if(found!=std::string::npos)
										{
											c2++;
											if(it3->ratio != 0 )
											{
												r2.push_back(it3->ratio);
												
											}
											

										}
										found = it3->Raw_File.find("P3");
										if(found!=std::string::npos)
										{
											c3++;
											if(it3->ratio != 0 )
											{
												r3.push_back(it3->ratio);
												
												//out2<<it3->seq<<" "<<it3->Raw_File<<" "<<it3->scan_number<<" "<<it3->m_z<<" "<<it3->ratio<<" "<<it3->unmod_pep_id<<endl;
											}
										

										}

										//break;

								   }
								 								   

							  }

						 for(std::list<int>::iterator it2= it->mod_peptide_id.begin(); it2!=it->mod_peptide_id.end(); it2++)
				              {																			
							
                              
								
							  if((*it2)==it3->mod_pep_id)
								{															
								
									out2<<it->protein_id<<" "<<it3->seq<<" "<<it3->Raw_File<<" "<<it3->scan_number<<" "<<it3->m_z<<" "<<it3->ratio<<" mod pep id: "<<it3->mod_pep_id<<endl;
									    std::size_t found;
										found = it3->Raw_File.find("P1");
										if(found!=std::string::npos)
										{
											c1++;
											if(it3->ratio != 0 )
											{
												r1.push_back(it3->ratio);
												
											}
											

										}
										found = it3->Raw_File.find("P2");
										if(found!=std::string::npos)
										{
											c2++;
											if(it3->ratio != 0 )
											{
												r2.push_back(it3->ratio);
												
											}
											

										}
										found = it3->Raw_File.find("P3");
										if(found!=std::string::npos)
										{
											c3++;
											if(it3->ratio != 0 )
											{
												r3.push_back(it3->ratio);
												
												//out2<<it3->seq<<" "<<it3->Raw_File<<" "<<it3->scan_number<<" "<<it3->m_z<<" "<<it3->ratio<<" "<<it3->mod_pep_id<<endl;
											}
										

										}

										//break;

								   }
								 								   

							  }





								

						 }//Done ratio assigning


							if(it->protein_id=="IPI00109529.2")
							{
								out2<<it->protein_id<<endl;
								out2<<"r1"<<" ";
								for(list<double>::iterator i1 = r1.begin(); i1!=r1.end(); i1++)
								{
									out2<<(*i1)<<" ";
								}
								out2<<endl;

								out2<<"r2"<<" ";
								for(list<double>::iterator i1 = r2.begin(); i1!=r2.end(); i1++)
								{
									out2<<(*i1)<<" ";
								}
								out2<<endl;

								out2<<"r3"<<" ";
								for(list<double>::iterator i1 = r3.begin(); i1!=r3.end(); i1++)
								{
									out2<<(*i1)<<" ";
								}
								out2<<endl;



							}

							Protein_Information::Cal_Median(r1,r2,r3,r1_median,r2_median,r3_median);

							if(it->protein_id=="IPI00109529.2")
							{
								out2<<"r1 median "<<r1_median<<endl;
								out2<<"r2 median "<<r2_median<<endl;
								out2<<"r3 median "<<r3_median<<endl;

							}


							r1.clear();
							r2.clear();
							r3.clear();
							//now caluclate kdeg and half life

							double kdeg = (log(r1_median+1)*1.5+log(r2_median+1)*4.5+log(r3_median+1)*13.5)/204.75 - log(2.0)/27.5;
							double half_life = log(2.0)/kdeg;

							if(half_life > 0)
							{
							  out<<it->protein_id<<" "<<r1_median<<" "<<r2_median<<" "<<r3_median<<" "<<half_life<<" "<<c1<<" "<<c2<<" "<<c3<<endl;
		  
							}
		


							//cout<<endl;
	

		  }

    }

	out.close();
	//out2.close();
}

/*This function calculates protein half-life while using protein and protein information classes
with the evidence information from the input file*/
void Protein_Information::Calculate_Protein_Ratio_With_Evidence(list<Protein> &p, list<Protein_Information> &pi)
{

	ofstream out,out2;
	int three_ratio_counter = 0;
	out.open("Protein_ratio_calculation.txt");
	out2.open("debug_info.txt");
	std::list<Protein>::iterator it;

	out<<"Protein ID"<<"                                        "<<"R1(1.5h)"<<"     "<<"R2(4.5h)"<<"    "<<"R3(13.5)"<<"    "<<"Half Life (h)"<<"      "<<"R1 count"<<"  "<<"R2 count"<<"  "<<"R3 count"<<endl;
	for(it=p.begin(); it!=p.end(); it++)
	{
		
		cout<<it->protein_id<<endl;
		
		if(it->protein_id == "IPI00881242.1"||it->protein_id=="IPI00112641.2" || it->protein_id=="IPI00649471.1" || it->protein_id=="IPI00876083.1"||it->protein_id=="IPI00330591.1")
		{
		

							std::list<double> r1;
							std::list<double> r2;
							std::list<double> r3;
							double r1_median = 0;
							double r2_median = 0;
							double r3_median = 0;
							double t1 = 1.5;
							double t2 = 4.5;
							double t3 = 13.5;
							int c1, c2, c3 = 0;

				   c1 = c2 = c3 = 0;

				for(std::list<Protein_Information>::iterator it3 = pi.begin(); it3!=pi.end(); it3++)
				{
                      
				         for(std::list<int>::iterator it2=it->evidence_id.begin(); it2!=it->evidence_id.end(); it2++)
				         {																			
							
								if((*it2)==it3->evidence_number)
								{		
									out2<<it->protein_id<<" "<<it3->seq<<" "<<it3->Raw_File<<" "<<it3->scan_number<<" "<<it3->m_z<<" "<<it3->ratio<<" "<<it3->evidence_number<<endl;
										std::size_t found;
										found = it3->Raw_File.find("P1");
										if(found!=std::string::npos)
										{
											c1++;
											if(it3->ratio != 0 )
											{
												r1.push_back(it3->ratio);
												
											}
											

										}
										found = it3->Raw_File.find("P2");
										if(found!=std::string::npos)
										{
											c2++;
											if(it3->ratio != 0 )
											{
												r2.push_back(it3->ratio);
												
											}
											

										}
										found = it3->Raw_File.find("P3");
										if(found!=std::string::npos)
										{
											c3++;
											if(it3->ratio != 0 )
											{
												r3.push_back(it3->ratio);
												
												//out2<<it3->seq<<" "<<it3->Raw_File<<" "<<it3->scan_number<<" "<<it3->m_z<<" "<<it3->ratio<<" "<<it3->unmod_pep_id<<endl;
											}
										

										}

										//break;

								   }
								 								   

							  }

						




								

						 }//Done ratio assigning

                          
							if(it->protein_id == "IPI00881242.1"||it->protein_id=="IPI00112641.2" || it->protein_id=="IPI00649471.1" || it->protein_id=="IPI00876083.1"||it->protein_id=="IPI00330591.1")
							{
								out2<<it->protein_id<<endl;
								out2<<"r1"<<" ";
								for(list<double>::iterator i1 = r1.begin(); i1!=r1.end(); i1++)
								{
									out2<<(*i1)<<" ";
								}
								out2<<endl;

								out2<<"r2"<<" ";
								for(list<double>::iterator i1 = r2.begin(); i1!=r2.end(); i1++)
								{
									out2<<(*i1)<<" ";
								}
								out2<<endl;

								out2<<"r3"<<" ";
								for(list<double>::iterator i1 = r3.begin(); i1!=r3.end(); i1++)
								{
									out2<<(*i1)<<" ";
								}
								out2<<endl;



							}

							Protein_Information::Cal_Median(r1,r2,r3,r1_median,r2_median,r3_median);

							if(it->protein_id == "IPI00881242.1"||it->protein_id=="IPI00112641.2" || it->protein_id=="IPI00649471.1" || it->protein_id=="IPI00876083.1"||it->protein_id=="IPI00330591.1")
							{
								out2<<"r1 median "<<r1_median<<endl;
								out2<<"r2 median "<<r2_median<<endl;
								out2<<"r3 median "<<r3_median<<endl;

							}

							if(r1_median == 0)
								t1 = 0;
							if(r2_median == 0)
								t2 = 0;
							if(r3_median == 0)
								t3 = 0;


							if(r1_median != 0 && r2_median !=0 && r3_median !=0)
								three_ratio_counter ++;

							r1.clear();
							r2.clear();
							r3.clear();
							//now caluclate kdeg and half life
		
							

							double kdeg = (log(r1_median+1)*t1+log(r2_median+1)*t2+log(r3_median+1)*t3)/(t1*t1+t2*t2+t3*t3) - log(2.0)/27.5;
							double half_life = log(2.0)/kdeg;

							if(half_life > 0)
							{
							  
								
								list<string>::iterator it4;
								/*for(it4=it->protein_id_list.begin(); it4!=it->protein_id_list.end(); it4++)
								{

									out<<(*it4)<<" ";
								}
								out<<"\n";*/

								out<<it->protein_id<<" "<<r1_median<<" "<<r2_median<<" "<<r3_median<<" "<<half_life<<" "<<c1<<" "<<c2<<" "<<c3<<endl;
		  
							}
		


							//cout<<endl;
	

		  }

    }

	out.close();
	out2.close();

	cout<<"No of proteins with three ratios "<<three_ratio_counter<<endl;
}


/*Calculates protein half-life while using peptide ratio*/
void Protein_Information::Calculate_Protein_Ratio_With_Evidence_Peptide_Ratio(list<Protein> &p, list<Protein_Information> &pi)
{

	ofstream out,out2;
	int three_ratio_counter = 0;
	out.open("Protein_ratio_calculation.txt");
	out2.open("debug_info.txt");
	std::list<Protein>::iterator it;

	out<<"Protein ID"<<"   "<<"R1(1.5h)"<<"     "<<"R2(4.5h)"<<"    "<<"R3(13.5)"<<"    "<<"Half Life (h)"<<"      "<<"R1 count"<<"  "<<"R2 count"<<"  "<<"R3 count"<<endl;
	for(it=p.begin(); it!=p.end(); it++)
	{
		
		cout<<it->protein_id<<endl;
		
		if(it->protein_id=="IPI00112641.2" || it->protein_id=="IPI00881242.1" || it->protein_id=="IPI00649471.1" || it->protein_id=="IPI00876083.1" || it->protein_id=="IPI00330591.1")
		{
		

							std::list<double> r1;
							std::list<double> r2;
							std::list<double> r3;
							std::multimap<string,double> r1_peptide;
							std::multimap<string,double> r2_peptide;
							std::multimap<string,double> r3_peptide;
							double r1_median = 0;
							double r2_median = 0;
							double r3_median = 0;
							double t1 = 1.5;
							double t2 = 4.5;
							double t3 = 13.5;
							int c1, c2, c3 = 0;

				   c1 = c2 = c3 = 0;

				for(std::list<Protein_Information>::iterator it3 = pi.begin(); it3!=pi.end(); it3++)
				{
                      
				         for(std::list<int>::iterator it2=it->evidence_id.begin(); it2!=it->evidence_id.end(); it2++)
				         {																			
							
								if((*it2)==it3->evidence_number)
								{		
									out2<<it->protein_id<<" "<<it3->seq<<" "<<it3->Raw_File<<" "<<it3->scan_number<<" "<<it3->m_z<<" "<<it3->ratio<<" "<<it3->evidence_number<<endl;
										std::size_t found;
										found = it3->Raw_File.find("P1");
										if(found!=std::string::npos)
										{
											c1++;
											if(it3->ratio != 0 )
											{
												r1.push_back(it3->ratio);
												r1_peptide.insert(std::pair<string,double>(it3->seq,it3->ratio));
											}
											

										}
										found = it3->Raw_File.find("P2");
										if(found!=std::string::npos)
										{
											c2++;
											if(it3->ratio != 0 )
											{
												r2.push_back(it3->ratio);
												r2_peptide.insert(std::pair<string,double>(it3->seq,it3->ratio));
											}
											

										}
										found = it3->Raw_File.find("P3");
										if(found!=std::string::npos)
										{
											c3++;
											if(it3->ratio != 0 )
											{
												r3.push_back(it3->ratio);
												r3_peptide.insert(std::pair<string,double>(it3->seq,it3->ratio));
												//out2<<it3->seq<<" "<<it3->Raw_File<<" "<<it3->scan_number<<" "<<it3->m_z<<" "<<it3->ratio<<" "<<it3->unmod_pep_id<<endl;
											}
										

										}

										//break;

								   }
								 								   

							  }

						




								

						 }//Done ratio assigning

				 
				         std::list<double>new_r1;
						 std::list<double>new_r2;
						 std::list<double>new_r3;

						 Traverse_Map(r1_peptide, new_r1);
						 Traverse_Map(r2_peptide, new_r2);
						 Traverse_Map(r3_peptide, new_r3);
						  

                          
							
								out2<<it->protein_id<<endl;
								out2<<"r1"<<" ";
								for(list<double>::iterator i1 = new_r1.begin(); i1!=new_r1.end(); i1++)
								{
									out2<<(*i1)<<" ";
								}
								out2<<endl;

								out2<<"r2"<<" ";
								for(list<double>::iterator i1 = new_r2.begin(); i1!=new_r2.end(); i1++)
								{
									out2<<(*i1)<<" ";
								}
								out2<<endl;

								out2<<"r3"<<" ";
								for(list<double>::iterator i1 = new_r3.begin(); i1!=new_r3.end(); i1++)
								{
									out2<<(*i1)<<" ";
								}
								out2<<endl;



							

							Protein_Information::Cal_Median(new_r1,new_r2,new_r3,r1_median,r2_median,r3_median);

							
								out2<<"r1 median "<<r1_median<<endl;
								out2<<"r2 median "<<r2_median<<endl;
								out2<<"r3 median "<<r3_median<<endl;

							

							if(r1_median == 0)
								t1 = 0;
							if(r2_median == 0)
								t2 = 0;
							if(r3_median == 0)
								t3 = 0;


							if(r1_median != 0 && r2_median !=0 && r3_median !=0)
								three_ratio_counter ++;

							r1.clear();
							r2.clear();
							r3.clear();
							//now caluclate kdeg and half life
		
							

							double kdeg = (log(r1_median+1)*t1+log(r2_median+1)*t2+log(r3_median+1)*t3)/(t1*t1+t2*t2+t3*t3) - log(2.0)/27.5;
							double half_life = log(2.0)/kdeg;

							if(half_life > 0)
							{
							  
								
								list<string>::iterator it4;
								/*for(it4=it->protein_id_list.begin(); it4!=it->protein_id_list.end(); it4++)
								{

									out<<(*it4)<<" ";
								}
								out<<"\n";*/

								out<<it->protein_id<<" "<<r1_median<<" "<<r2_median<<" "<<r3_median<<" "<<half_life<<" "<<c1<<" "<<c2<<" "<<c3<<endl;
		  
							}
		


							//cout<<endl;
	

		  }

    }

	out.close();
	out2.close();

	cout<<"No of proteins with three ratios "<<three_ratio_counter<<endl;


}


/*This function retrieves evidence information from the input file for the proteins*/
void Protein_Information::Get_Evidence_File_Info(list<Protein_Information> &m)
{




	int column_count, line_count = 0;
	string line, csvItem, protein;

	ifstream myfile;
	ofstream out;

	myfile.open("evidence01.csv");

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
			string protein_name="";
			
			Protein_Information p;
			
			string raw_file_name="";
			list<string> name_list;
			double ratio;
			long intensity_low;
			long intensity_high;
			while(getline(myline, csvItem,','))
			{
					/*if(csvItem.size()>2)
					{*/
                          if(column_count == 1)
						  {
							  p.seq = csvItem.c_str();
						  }
						  if(column_count==12)
							{
									//char *str = new char[csvItem.length()+1];

									//char *pch;

									//strcpy(str,csvItem.c_str());

									//pch = strtok(str,":|;");
								 //   
									//string cmp = pch;

									//int count  = 1;
								 //   
									///*if(cmp == "IPI")
									//{*/
									//	  while(pch != NULL) 
									//	  {
									//		//out<<pch<<" ";
									//		  cmp = pch;
									//		  std::size_t found = cmp.find(".");
									//		  if(found!=std::string::npos)
									//		  {
									//			 protein_name=pch;
									//			 name_list.push_back(pch);
									//			//cout<<pch<<"\n";
									//			 //break;
									//		  }
									//		//cout<<pch<<"\n";
									//		pch = strtok (NULL,":|;");											
									//		count ++;
									//	  }

									///*}*/

									//delete []str;
									////break;
								}
						      
						       if(column_count==17)
							   {
								   p.Raw_File = csvItem.c_str();
								    
							   }
							   if(column_count==20)
							   {
								   p.m_z = atof(csvItem.c_str());
							   }
							   if(column_count==53)
							   {
								   p.ratio = atof(csvItem.c_str());
							   }
							   if(column_count==49)
							   {

								   p.scan_number = atol(csvItem.c_str());
							   }
							   if(column_count==57)
							   {
								   if(csvItem.size()>0)
								   {
								     if(isdigit(csvItem[0]))
								        p.Intensity_L = atol(csvItem.c_str());
									 else
                                        p.Intensity_L = 0;
								   }
								  else
                                   p.Intensity_L = 0;
							   
							   }
							   if(column_count==58)
							   {
								   if(csvItem.size()>0)
								   {
								      if(isdigit(csvItem[0]))
								         p.Intensity_H = atol(csvItem.c_str());
									  else
										  p.Intensity_H = 0;
								   }
								   else
									   p.Intensity_H = 0;
							   }
							   if(column_count==61)
							   {
								   p.evidence_number = atol(csvItem.c_str());

							   }
							   if(column_count==63)
							   {

								 //  char *str = new char[csvItem.length()+1];

									//char *pch;

									//std::strcpy(str,csvItem.c_str());

									//pch = strtok(str,";");									
									//
									//while(pch != NULL) 
								 //   {
									//		
									//	p.pep_id.push_back(atoi(pch));
									//		
									//		//cout<<pch<<"\n";
									//	pch = strtok (NULL,";");
									//		
								 //    }
									//delete []str;

								   p.unmod_pep_id = atoi(csvItem.c_str());

								}					
							     if(column_count==64)
							     {

								 //  char *str = new char[csvItem.length()+1];

									//char *pch;

									//std::strcpy(str,csvItem.c_str());

									//pch = strtok(str,";");									
									//
									//while(pch != NULL) 
								 //   {
									//		
									//	p.pep_id.push_back(atoi(pch));
									//		
									//		//cout<<pch<<"\n";
									//	pch = strtok (NULL,";");
									//		
								 //    }
									//delete []str;

									 p.mod_pep_id = atoi(csvItem.c_str());

								}				 

								 
							   
					 
					/*}*/
					column_count++;
			
			}
			
			//m.push_back(p);
			m.push_back(p);
		
		}

		line_count++;
		
	}


	myfile.close();






}
