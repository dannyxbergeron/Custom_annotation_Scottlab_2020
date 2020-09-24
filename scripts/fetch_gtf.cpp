#include <iostream>
#include <fstream>
#include <string>
#include <unordered_map>

using namespace std;

void replace_biotype(string& line, const string& newBiotype)
{
    int found;
    string st, nd;
    string gene_biotype = "gene_biotype \"";
    if(newBiotype != "Unknown")
    {
        found = line.find(gene_biotype) + gene_biotype.size();
        st = line.substr(0, found);
        found = line.find('"', found);
        nd = line.substr(found);
        line = st + newBiotype + nd;
    }

    cout << line << endl;
}

int main(int argc, char** argv)
{
  ifstream infile(argv[1]);
  string line, id, biotype;
  int found;
  unordered_map<string, string> idsBiotypes;

  // Read the id file and create a map<pattern, biotype>
  while(getline(infile, line))
  {
      found = line.find('\t');
      id = line.substr(0, found);
      biotype = line.substr(found+1);
      idsBiotypes[id] = biotype;
  }

  // Close the input stream and open the gtf file to search for patterns
  infile.close();
  infile.clear();
  infile.open(argv[2]);

  // Look in each line to find the pattern and replace the gene_biotype
  while(getline(infile, line))
  {
    for(auto& m: idsBiotypes)
      if(line.find(m.first) != string::npos)
        replace_biotype(line, m.second);
  }

  return 0;
}
