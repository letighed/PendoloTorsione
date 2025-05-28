#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>

using namespace std;

struct risonanza {
    vector <double> tempo ;
    vector <double> forzante ;
    vector <double> pendolo ;
    vector <double> ampiezza ;
    vector <double> fase ;
};

int main(){

string file_input;
string file_output;

cout << "Nome file input" << endl;
cin >> file_input;

cout << "Nome file output" << endl;
cin >> file_output;

ifstream fin(file_input);
	if (!fin){
		cout << "error " << endl;
		return -1;
}


risonanza file;
/*
vector <double> tempo ;
vector <double> forzante ;
vector <double> pendolo ;
vector <double> ampiezza ;
vector <double> fase ;
*/
double temp, forz, pend, amp, fas ;


while (fin >> temp >> forz >> pend >> amp >> fas){
	if (forz != 0){
    file.tempo.push_back(temp);
	file.forzante.push_back(forz);
	file.pendolo.push_back(pend);
	file.ampiezza.push_back(amp);
	file.fase.push_back(fas);
}
}


// sopra il 80%
vector <int> indici_mass80 ;
double max =0;
for (auto c : file.pendolo){
if (c > max){
max = c;
}}

for (int i = 0; i < file.pendolo.size(); ++i) {
    if (file.pendolo[i] >= 0.8 * max) {
            indici_mass80.push_back(i);
        }
    }

// sotto l'80 %
vector <int> indici_min80 ;
double min = 0 ;
for (auto c : file.pendolo){
    if (c <  min ){
        min = c;
    }
}
for (int i = 0; i < file.pendolo.size(); ++i) {
    if (file.pendolo[i] <= 0.8 * min) {
            indici_min80.push_back(i);
    }
}


// output
ofstream fout(file_output);
	if (!fout){	
		cout << "error output " << endl;
		return -1;
}
    fout << "tempo\tforzante\tpendolo\tampiezza\tfase\n";
    for (int i : indici_mass80) {
        fout << file.tempo[i] << "\t"
             << file.forzante[i] << "\t"
             << file.pendolo[i] << "\t"
             << file.ampiezza[i] << "\t"
             << file.fase[i] << "\n";
    }
return 0;
}
