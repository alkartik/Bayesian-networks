#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <unordered_map>
#include <unordered_set>
#include <cmath>
#include <ctime>
#include <iomanip>

// Format checker just assumes you have Alarm.bif and Solved_Alarm.bif (your file) in current directory
using namespace std;

// Our graph consists of a vector of nodes where each node is represented as follows:
class Graph_Node
{

private:
	string Node_Name;		// Variable name
	vector<int> Children;	// Children of a particular node - these are index of nodes in graph.
	vector<int> Parents; // Parents of a particular node- note these are names of parents
	int nvalues;			// Number of categories a variable represented by this node can take
	vector<string> values;	// Categories of possible values
	vector<float> CPT;		// conditional probability table as a 1-d array . Look for BIF format to understand its meaning

public:
	// Constructor- a node is initialised with its name and its categories
	Graph_Node(string name, int n, vector<string> vals)
	{
		Node_Name = name;

		nvalues = n;
		values = vals;
	}
	string get_name()
	{
		return Node_Name;
	}
	int size_of_CPT()
	{
		return CPT.size();
	}
	vector<int>* get_children()
	{
		return &Children;
	}
	vector<int>* get_Parents()
	{
		return &Parents;
	}
	float find_CP(int i)
	{
		return CPT[i];
	}
	vector<float> get_CPT()
	{
		return CPT;
	}
	int get_nvalues()
	{
		return nvalues;
	}
	vector<string>* get_values()
	{
		return &values;
	}
	void set_CPT(vector<float> new_CPT)
	{
		CPT.clear();
		CPT = new_CPT;
	}
	void set_Parents(vector<int> Parent_Nodes)
	{
		Parents.clear();
		Parents = Parent_Nodes;
	}
	int add_parent(int index){
		for (int i = 0; i < Parents.size(); i++)
		{
			if (Parents[i] == index)
				return 0;
		}
		Parents.push_back(index);
		return 1;
	}
	// add another node in a graph as a child of this node
	int add_child(int new_child_index)
	{
		for (int i = 0; i < Children.size(); i++)
		{
			if (Children[i] == new_child_index)
				return 0;
		}
		Children.push_back(new_child_index);
		return 1;
	}
};

// The whole network represted as a vector of nodes
class network
{

	vector<Graph_Node> Pres_Graph;

public:
	int totalline = 0;

	unordered_map<string, int> var_to_indx;
	vector<vector<string>> allvalues;
	int addNode(Graph_Node node)
	{
		Pres_Graph.push_back(node);
		return 0;
	}

	int netSize()
	{
		return Pres_Graph.size();
	}
	// get the index of node with a given name
	int get_index(string val_name)
	{
		vector<Graph_Node>::iterator vectorIt;
		int count = 0;
		for (vectorIt = Pres_Graph.begin(); vectorIt != Pres_Graph.end(); vectorIt++)
		{
			if (vectorIt->get_name().compare(val_name) == 0)
				return count;
			count++;
		}
		return -1;
	}
	// get the node at nth index
	vector<Graph_Node>::iterator get_nth_node(int n)
	{
		vector<Graph_Node>::iterator vectorIt=Pres_Graph.begin();
		vectorIt+=n;
		return vectorIt;
	}
	// get the iterator of a node with a given name
	vector<Graph_Node>::iterator search_node(string val_name)
	{
		vector<Graph_Node>::iterator vectorIt;
		for (vectorIt = Pres_Graph.begin(); vectorIt != Pres_Graph.end(); vectorIt++)
		{
			if (vectorIt->get_name().compare(val_name) == 0)
				return vectorIt;
		}
		// //cout << "node not found\n";
		return vectorIt;
	}
};

network read_network(string alarmfile)
{
	network Alarm;
	string line;
	int find = 0;
	int *totalline = &Alarm.totalline;
	ifstream myfile(alarmfile);
	string temp;
	string name;
	vector<string> values;
	if (myfile.is_open())
	{
		while (!myfile.eof())
		{
			stringstream ss;
			getline(myfile, line);
			(*totalline) += 1;

			ss.str(line);
			ss >> temp;

			if (temp.compare("variable") == 0)
			{

				ss >> name;
				getline(myfile, line);
				(*totalline) += 1;
				stringstream ss2;
				ss2.str(line);
				for (int i = 0; i < 4; i++)
				{
					ss2 >> temp;
				}
				values.clear();
				while (temp.compare("};") != 0)
				{
					values.push_back(temp);

					ss2 >> temp;
				}
				Graph_Node new_node(name, values.size(), values);
				Alarm.allvalues.push_back(values);
				Alarm.var_to_indx[name] = Alarm.netSize();
				int pos = Alarm.addNode(new_node);
			}
			else if (temp.compare("probability") == 0)
			{
				ss >> temp;
				ss >> temp;
				vector<Graph_Node>::iterator vectorIt;
				vector<Graph_Node>::iterator vectorIt1;
				vectorIt = Alarm.search_node(temp);
				int index = Alarm.get_index(temp);
				ss >> temp;
				values.clear();
				while (temp.compare(")") != 0)
				{
					vectorIt1 = Alarm.search_node(temp);
					vectorIt1->add_child(index);
					vectorIt->add_parent(Alarm.var_to_indx[temp]);
					ss >> temp;
				}\
				getline(myfile, line);
				stringstream ss2;

				ss2.str(line);
				ss2 >> temp;

				ss2 >> temp;

				vector<float> curr_CPT;
				string::size_type sz;
				while (temp.compare(";") != 0)
				{

					curr_CPT.push_back(atof(temp.c_str()));

					ss2 >> temp;
				}

				vectorIt->set_CPT(curr_CPT);
			}
			else
			{
			}
		}

		if (find == 1)

			myfile.close();
	}

	return Alarm;
}

class solver
{
public:
	string data_file; // name of data file
	string alarmfile;
	network *Alarm;
	unordered_map<int,vector<string>> Data;			 // stores patient data
	unordered_map<int, int> questionmark;	 // maps patient index to missing variable index
	unordered_map<int,vector<float>> prob_of_questmark; // stores probability of question mark being some value given parents.
	vector<vector<vector<string>>> var_to_allCPT;
	solver(network *Alarm, string alarmfile, string datafile)
	{
		this->alarmfile=alarmfile;
		this->data_file = datafile;
		this->Alarm = Alarm;
	}
	void CPT_init()
	{										// initializes CPT for each variable////////////////////////////////////////////////////////////
		long long Gsize = Alarm->netSize(); // size of graph
		for (long long i = 0; i < Gsize; i++)
		{
			vector<Graph_Node>::iterator Gnode = Alarm->get_nth_node(i);
			int CPTsize = Gnode->size_of_CPT();
			vector<float> CPT(CPTsize, 0);
			float nvalues = Gnode->get_nvalues();
			for (int i = 0; i < CPT.size(); i++)
			{
				CPT[i] = 1 / (nvalues);
			}
			Gnode->set_CPT(CPT);
		}
		return;
	}
	void Datafile_reader()
	{ ////////////////////////////////////////////////////
		ifstream datafile(this->data_file);
		string temp;
		string line;
		int index = 0;
		if (datafile.is_open())
		{

			while (!datafile.eof())
			{
				stringstream ss;
				getline(datafile, line);
				ss.str(line);
				vector<string> patient;
				for (int i = 0; i < 37; i++)
				{
					ss >> temp;
					if (temp == "\"?\"")
					{
						questionmark[index] = i;
					}
					patient.push_back(temp);
				}
				Data[index]=(patient);
				index++;
			}
			datafile.close();
		}
	}

	int CPT_index(vector<Graph_Node>::iterator &child, vector<string> values)
	{ ////////////////////////////////////////////
		unordered_map<string, int> *var_to_indx = &Alarm->var_to_indx;
		vector<vector<string>> *allvalues = &Alarm->allvalues;
		vector<int> *parent = child->get_Parents();
		int CPTindex = 0;
		int nvalues = child->get_nvalues();
		int CPTsize = child->size_of_CPT();
		int Cindex = (*var_to_indx)[child->get_name()];
		vector<string> *value = &(*allvalues)[Cindex];
		for (int j = 0; j < nvalues; j++)
		{
			if ((*value)[j] == values[0])
			{
				CPTindex += ((CPTsize / nvalues) * j);
				CPTsize /= nvalues;
				break;
			}
		}

		for (int i = 1; i < values.size(); i++)
		{
			int Pindex = (*parent)[i-1];
			value = &(*allvalues)[Pindex];
			for (int j = 0; j < (*value).size(); j++)
			{
				if ((*value)[j] == values[i])
				{
					CPTindex += ((CPTsize / (*value).size()) * j);
					CPTsize /= (*value).size();
					break;
				}
			}
		}
		return CPTindex;
	}
	float CP_finder(vector<Graph_Node>::iterator child, vector<string> &patient, pair<int, string> missing)
	{ //////////////////////////////////////////////
		vector<vector<string>> *allvalues = &Alarm->allvalues;
		unordered_map<string, int> *var_to_indx = &Alarm->var_to_indx;
		vector<int> *parent = child->get_Parents();
		vector<string> values;
		int Cindex = (*var_to_indx)[child->get_name()];
		if (Cindex == missing.first)
		{
			values.push_back(missing.second);
		}
		else
		{
			values.push_back(patient[Cindex]);
		}
		for (int i = 0; i < (*parent).size(); i++)
		{
			int Pindex = (*parent)[i];
			if (Pindex == missing.first)
			{
				values.push_back(missing.second);
			}
			else
			{
				values.push_back(patient[Pindex]);
			}
		}
		int CPTindex = CPT_index(child, values);
		return child->find_CP(CPTindex);
	}
	void CPT_to_dataset()
	{ /////////////////////////////////////////////////////
		unordered_map<string, int> *var_to_indx = &Alarm->var_to_indx;
		for (auto &data:Data)
		{
			int index = questionmark[data.first];
			vector<Graph_Node>::iterator datapoint = Alarm->get_nth_node(index);
			vector<string> *value = datapoint->get_values();
			int nvalues = datapoint->get_nvalues();
			vector<int> *children = datapoint->get_children();
			vector<float> prob;
			float total = 0;
			for (int v = 0; v < nvalues; v++)
			{
				float cp = CP_finder(datapoint, data.second, make_pair(index, (*value)[v]));
				// //cout<<cp<<" ";
				for (auto &child : *children)
				{
					vector<Graph_Node>::iterator Cnode = Alarm->get_nth_node(child);
					cp *= CP_finder(Cnode, data.second, make_pair(index, (*value)[v]));
					// //cout<<cp<<" ";
				}
				// //cout<<endl;
				total += cp;
				prob.push_back(cp);
			}
			for (auto &a : prob)
			{
				a /= (total);
			}
			prob_of_questmark[data.first]=(prob);
		}
	}
	bool patientChecker(vector<string> &patient, vector<string> &values, vector<Graph_Node>::iterator &child, bool num)
	{ ////////////////////////////////////////////////
		vector<int> *parent = child->get_Parents();
		unordered_map<string, int> *var_to_indx = &Alarm->var_to_indx;
		int Cindex = (*var_to_indx)[child->get_name()];
		if (num)
		{
			if (values[0] != patient[Cindex] && patient[Cindex] != "\"?\"")
			{
				// //cout<<"False"<<endl;
				return false;
			}
		}
		int i = 1;
		for (auto &Pindex : *parent)
		{
			if (values[i] != patient[Pindex] && patient[Pindex] != "\"?\"")
			{
				// //cout<<"False"<<endl;
				return false;
			}
			i++;
		}
		// //cout<<"True"<<endl;
		return true;
	}
	vector<vector<string>> values_Init(vector<Graph_Node>::iterator &child)

	{  ///////////////////////////////////////////////////////////////////////////////////
		unordered_map<string, int> *var_to_indx = &Alarm->var_to_indx;
		vector<vector<string>> *allvalues = &Alarm->allvalues;
		vector<int> *parent = child->get_Parents();
		int CPTsize = child->size_of_CPT();
		int nvalues = child->get_nvalues();
		vector<vector<string>> value_init(CPTsize, vector<string>((*parent).size() + 1, "khushal"));
		vector<string> *value = child->get_values();
		int index = 0;
		int blocksize = CPTsize;
		for (int i = 0; i < blocksize; i += blocksize / nvalues)  ////////////
		{
			string val = (*value)[index];
			if (index != (*value).size())
			{
				index++;
			}
			else
			{
				index = 0;
			}
			for (int j = 0; j < blocksize / nvalues; j++)
			{
				value_init[i + j][0] = val;
			}
		}
		blocksize /= nvalues;
		for (int p = 0; p < (*parent).size(); p++)
		{
			int start = 0;
			int Pindex = (*parent)[p];
			value = &(*allvalues)[Pindex];
			int size = (*value).size();
			while (start + blocksize <= CPTsize)
			{
				index = 0;
				for (int i = 0; i < blocksize; i += blocksize / size)
				{
					string val = (*value)[index];
					if (index != (*value).size())
					{
						index++;
					}
					else
					{
						index = 0;
					}
					for (int j = 0; j < blocksize / size; j++)
					{
						value_init[start + i + j][p + 1] = val;
					}
				}
				start += blocksize;
			}
			blocksize/=size;
		}
	return value_init;
	}
	void var_to_CPT(){
		int Gsize=Alarm->netSize();
		for (int i=0;i<Gsize;i++){
			vector<vector<string>> values;
			vector<Graph_Node>::iterator Gnode = Alarm->get_nth_node(i);
			values=values_Init(Gnode);
			var_to_allCPT.push_back(values);
		}
	}
	int value_index(int index, string val)
	{ ////////////////////////////////////////////////
		unordered_map<string, int> *var_to_indx = &Alarm->var_to_indx;
		vector<vector<string>> *allvalues = &Alarm->allvalues;
		vector<string> value = (*allvalues)[index];
		for (int k = 0; k < value.size(); k++)
		{
			if (value[k] == val)
			{
				return k;
			}
		}
		return -1;
	}
	float dataset_to_CP(vector<Graph_Node>::iterator &child, vector<string> &values)
	{ ///////////////////////////////////////
		vector<int> indexes;
		unordered_map<string, int> *var_to_indx = &Alarm->var_to_indx;
		vector<vector<string>> *allvalues = &Alarm->allvalues;
		int Cindex = (*var_to_indx)[child->get_name()];
		indexes.push_back(Cindex);
		vector<int> *parent = child->get_Parents();
		for (auto &Pindex : *parent)
		{
			indexes.push_back(Pindex);
		}
		float num = 0;
		float den = 0;
		for (auto &p:Data)
		{
			bool Ncheck = patientChecker(p.second, values, child, 1);
			bool Dcheck = patientChecker(p.second, values, child, 0);
			int Qindex = -1;
			for (int i = 0; i < values.size(); i++)
			{
				if (questionmark[p.first] == indexes[i])
				{
					Qindex = i;
				}
			}
			if (Ncheck)
			{
				if (Qindex != -1)
				{
					num += prob_of_questmark[p.first][value_index(indexes[Qindex], values[Qindex])];
				}
				else
				{
					num += 1;
				}
			}
			if (Dcheck)
			{
				if (Qindex != -1)
				{
					den += prob_of_questmark[p.first][value_index(indexes[Qindex], values[Qindex])];
				}
				else
				{
					den += 1;
				}
			}
		}
		// //cout<<"Done: "<<(num+0.0035)/(den+0.0035)<<endl;
		return (num) / (den+0.00001)+0.0001;
	}
	void datafile_to_CPT()
	{
		int Gsize = Alarm->netSize();
		for (int i = 0; i < Gsize; i++)
		{
			vector<Graph_Node>::iterator Gnode = Alarm->get_nth_node(i);
			// //cout<<"testing value_init!!"<<endl;
			vector<vector<string>> *value_init = &var_to_allCPT[i];
			// //cout<<"value_init worked!!"<<endl;
			//cout << i << " " << Gsize << endl;
			vector<float> CPT((*value_init).size());
			for (int j = 0; j < (*value_init).size(); j++)
			{
				// //cout<<"testing dataset_to_CP"<<endl;
				CPT[j] = dataset_to_CP(Gnode, (*value_init)[j]);
				if (CPT[j]>0.9999){
					CPT[j]=0.9999;
				}
				if (CPT[j]<=0.0001){
					CPT[j]=0.0001;
				}
				// //cout<<"tested dataset_to_CP"<<endl;
			}
			Gnode->set_CPT(CPT);
		}
		//cout << "datafile tested" << endl;
	}

	void CPTfilewriter(std::string InFile, std::string OutFile) ////////////////////////////////////////////////////
	{
		std::string line;
		std::ifstream infile(InFile);
		std::ofstream outfile(OutFile);
		std::string temp;
		std::string name;
		if (!infile.is_open())
			return;
		else
		{
			int line_no = 0;
			while (!infile.eof())
			{
				line_no++;
				std::stringstream ss;
				getline(infile, line);
				ss.str(line);
				ss >> temp;
				if (temp.compare("probability") == 0)
				{
					ss >> temp;
					ss >> temp;
					string name = temp;
					vector<float> CPT = (Alarm->get_nth_node((Alarm->var_to_indx)[name]))->get_CPT();
					outfile << '\n'<<line << std::endl;
					getline(infile, line);
					outfile << "    table ";
					for (int i = 0; i < CPT.size(); i++)
						outfile << std::fixed << std::setprecision(4) << CPT[i] << ' ';
					outfile << ";";
				}
				else
					if (line_no!=1)
					outfile << '\n'<<line;
					else
					outfile<<line;
			}
			infile.close();
		}
	}

	void converge()
	{
		time_t start, end;
		time(&start);
		//cout<<"hello1";
		CPT_init();
		//cout<<"hello1";

		var_to_CPT();
		//cout<<"hello1";

		Datafile_reader();
		//cout<<"Datafile reader working!!"<<endl;
		while (true)
		{
			double elapsed = time(&end) - start;
			//cout << "elapsed: " << elapsed << endl;
			if (elapsed >= 119)
			{
				break;
			}
			CPT_to_dataset();
			datafile_to_CPT();
			CPTfilewriter(alarmfile, "solved_alarm.bif");
		}
	}
};

int main(int argc, char* argv[]){
    string infile=string(argv[1]);
    string outfile=string(argv[2]);
	network Alarm;
		//cout<<"hello1";

	Alarm = read_network(infile); // working
		//cout<<"hello1";

	solver tester(&Alarm,infile, "records.dat");
	tester.converge();
	// Example: to do something
	// //cout << "Perfect! Hurrah! \n";
}